import re
from itertools import takewhile
from Bio import Align
from Bio.Seq import Seq
from django.db import transaction
from PinkStrawberry import models

# an attempt at imitating V-Quest's 'DNAPLOT' aligner without the imgt regions inserted
# alignment is local to cover some weird edgecases
# todo rewrite to account for IMGT introduced gaps, currently fails on low identity (< 60%)
class PairwiseAligner(Align.PairwiseAligner):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.match_score = 5
        self.mismatch_score = -4
        self.gap_score = -5000  # don't allow for gaps
        self.mode = 'local'


# global mode, gaps on right and left still allowed with penalty
class AltPairwiseAligner(Align.PairwiseAligner):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.match_score = 5
        self.mismatch_score = -4
        self.gap_score = -5000  # don't allow for gaps
        self.right_gap_score = -5
        self.left_gap_score = -5
        self.mode = 'global'


# align the given region and germline sequence using the PairwiseAligner class
# alignment is local to cover some weird edgecases, any coverage issues are compensated for by AltPairwiseAligner
def _align_germline(sequenceobj, startpos, regionseq, germlineseq, relative_coverage=1.0, d_region=False):
    alignment = PairwiseAligner().align(regionseq, germlineseq)[0]
    indent = calc_localindent(alignment)
    germstart = startpos + indent
    # check if the alignment has decent coverage this way, if not that realign with Alt
    germstop = indent + len(germlineseq)
    sizediff = len(germlineseq)/len(regionseq)
    coverage = (germstop-max(0, indent))/len(alignment.target)
    if coverage < sizediff*relative_coverage and coverage < 0.90:
        #print(f"{sequenceobj.id:<7} | {sequenceobj.Sequence_ID:<20} | {indent:<4}, {germstop:<4} | {coverage:.2f} | {sizediff:.2f} | {sizediff * relative_coverage:.2f} | {relative_coverage:.2f} | {regionseq}")
        # if d_region is True then get the D-region for the regionseq instead of the ND region, could just be mismatching N after all.
        if d_region:
            regionseq = sequenceobj.Nt_sequences.D_REGION
            startpos = sequenceobj.Nt_sequences.D_REGION_start

        #redo the alignment with Alt
        alignment = AltPairwiseAligner().align(regionseq, germlineseq)[0]
        germstart = startpos + calc_globalindent(alignment)
        sequenceobj.addcomment(f"redone alignment because of poor coverage {coverage:.2f}")
        sequenceobj.save()
    return (germstart-1, germlineseq)


def get_v_alignment(sequenceobj):
    return _align_germline(sequenceobj, sequenceobj.Nt_sequences.V_REGION_start,
                           sequenceobj.Nt_sequences.V_REGION,
                           sequenceobj.V_Germline.Nt_Sequence,
                           relative_coverage= 0.90)


def get_d_alignment(sequenceobj):
    if sequenceobj.islight():
        return None

    startpos = sequenceobj.Nt_sequences.V_REGION_end + 1
    return _align_germline(sequenceobj, startpos,
                           sequenceobj.Nt_sequences.V_D_J_REGION[startpos - sequenceobj.Nt_sequences.V_REGION_start:sequenceobj.Nt_sequences.J_REGION_start-sequenceobj.Nt_sequences.V_REGION_start],
                           sequenceobj.D_Germline.Nt_Sequence,
                           relative_coverage= 0.80,
                           d_region= True)


def get_j_alignment(sequenceobj):
    startpos = sequenceobj.Nt_sequences.V_REGION_end + 1
    if sequenceobj.isheavy():
        nxj_regionseq = sequenceobj.Nt_sequences.V_D_J_REGION[startpos-sequenceobj.Nt_sequences.V_REGION_start:]
    else:
        nxj_regionseq = sequenceobj.Nt_sequences.V_J_REGION[startpos-sequenceobj.Nt_sequences.V_REGION_start:]
    return _align_germline(sequenceobj, startpos,
                           nxj_regionseq,
                           sequenceobj.J_Germline.Nt_Sequence,
                           relative_coverage= 0.80)


# calculate indentation needed to make the start of the target and query match up. Negative means target indent, positive means query indent.
def calc_localindent(alignment, returnsequence=False):
    diff = alignment.coordinates[0][0] - alignment.coordinates[1][0]
    if not returnsequence:
        return diff

    if diff > 0:
        ans = [alignment.target, f'{diff*" "}{alignment.query}']
    elif diff < 0:
        ans = [f'{abs(diff)*" "}{alignment.target}', alignment.query]
    else:
        ans = [alignment.target, alignment.query]
    return ans

# calculate indentation needed to make the start of the target and query match up. Negative means target indent, positive means query indent.
def calc_globalindent(alignment):
    regionstart = sum(1 for i in takewhile(lambda x: x == -1, alignment.indices[0]))
    germstart = sum(1 for i in takewhile(lambda x: x == -1, alignment.indices[1]))
    diff = germstart - regionstart
    return diff


def align_germline(sequenceobj):
    v = get_v_alignment(sequenceobj)
    models.alignment.objects.get_or_create(name='V_REGION', position = v[0], sequence = v[1], linked_sequence=sequenceobj)

    j = get_j_alignment(sequenceobj)
    models.alignment.objects.get_or_create(name='J_REGION', position = j[0], sequence = j[1], linked_sequence=sequenceobj)

    if sequenceobj.isheavy():
        d = get_d_alignment(sequenceobj)
        models.alignment.objects.get_or_create(name='D_REGION', position=d[0], sequence=d[1], linked_sequence=sequenceobj)
    return


def show_alignment(sequencestr, alignmentsqueryset):
    print(''.join([f"{i:<{5}}" for i in range(0, len(sequencestr), 5)]))
    print(sequencestr)
    for alignmentobj in alignmentsqueryset:
        alignmentseq = alignmentobj.sequence if alignmentobj.position >= 0 else alignmentobj.sequence[abs(alignmentobj.position):]
        alignmentpos = alignmentobj.position if alignmentobj.position >= 0 else 0
        print(f"{alignmentobj.position * ' '}{alignmentseq} {alignmentpos}")
    return


def translate_alignment(sequenceobj, alignmentsqueryset):
    translated_alignments = []
    class tmpalign:
        def __init__(self, name, position, sequence):
            self.name = name
            self.position = position
            self.sequence = sequence
        def range(self):
            return (self.position, self.position + len(self.sequence))

    vframe = sequenceobj.Nt_sequences.V_REGION_reading_frame-1
    for alignmentobj in alignmentsqueryset:
        rfpos = alignmentobj.position - vframe
        rf = abs(rfpos%3 - 3) if rfpos%3 != 0 else 0
        pos = (rfpos+rf) // 3
        inframe = alignmentobj.sequence[rf:]
        if (_:=(len(alignmentobj.sequence)-rf)%3):
            inframe = inframe[:-_]
        translated = Seq(inframe).translate()
        translated_alignments.append(tmpalign(alignmentobj.name, pos, translated))
    return translated_alignments


def translate_range(start, stop, readingframe):
    newstart = -((start - readingframe) // -3)  # ceil, start at .3 or .6 implies only partial part, .0 full codon
    newstop  = ((stop - readingframe) // 3) - 1 + (((stop - readingframe) % 3) == 2)  # floor, but add one if modulo3 == 2 because .6 implies full codon
    return (newstart, newstop)


def get_alignmentdiff(sequencestr, alignmentobj):
    # catch for negative alignment positions
    alignmentseq = alignmentobj.sequence if alignmentobj.position >= 0 else alignmentobj.sequence[abs(alignmentobj.position):]
    alignmentpos = alignmentobj.position if alignmentobj.position >= 0 else 0

    # init list and start comparing
    diff = []
    for ind, a, b in zip(range(alignmentpos, alignmentpos+len(alignmentseq)),
                         sequencestr[alignmentpos:], alignmentseq):
        if a != b:
            diff.append({'original':b, 'changed':a, 'position':ind})
    return diff


# get all alignments, also D if heavychain
def get_germlinealignments(sequenceobj):
    vdj = [sequenceobj.alignment.get(name='V_REGION'),
              sequenceobj.alignment.get(name='J_REGION')]
    if sequenceobj.isheavy():
        vdj.insert(1, sequenceobj.alignment.get(name='D_REGION'))
    return vdj


def save_mutations(sequenceobj):
    # init storage to get the differences
    nt = dict()
    aa = dict()

    # get germline alignment objects & reading frame
    vdj = get_germlinealignments(sequenceobj)
    vframe = sequenceobj.Nt_sequences.V_REGION_reading_frame - 1

    #get the nt differences
    for i in vdj:
        nt[i.name[0]] = {'diff':get_alignmentdiff(sequenceobj.Nt_Sequence, i)}

    #get the aa differences
    for i in translate_alignment(sequenceobj, vdj):
        aa[i.name[0]] = {'diff':get_alignmentdiff(sequenceobj.AA_Sequence, i)}

    # get region range to make sure we don't store mutations/changes outside it
    nt['V']['range'] = (sequenceobj.Nt_sequences.V_REGION_start-1, sequenceobj.Nt_sequences.V_REGION_end)
    nt['J']['range'] = (sequenceobj.Nt_sequences.J_REGION_start - 1, sequenceobj.Nt_sequences.J_REGION_end)
    if sequenceobj.isheavy():
        nt['D']['range'] = ((sequenceobj.Nt_sequences.D_REGION_start-1, sequenceobj.Nt_sequences.D_REGION_end))

    for region, data in nt.items():
        ntrange = data['range']
        aa[region]['range'] = translate_range(ntrange[0], ntrange[1], vframe)

    # load aa into db first so linking with nucleotide mutations is easier
    with transaction.atomic():
        for region, data in aa.items():
            for i in data['diff']:
                # check for allowed ranges
                if not data['range'][0] <= i['position'] < data['range'][1]:
                    continue

                models.AA_change.objects.get_or_create(**i, sequence=sequenceobj)

    # load nt into db
    with transaction.atomic():
        for region, data in nt.items():
            for i in data['diff']:
                # check for allowed ranges
                if not data['range'][0] <= i['position'] < data['range'][1]:
                    continue

                # find position of related aa & change if any
                corrected = (i['position'] - vframe) // 3
                if (_ := sequenceobj.AA_changes.filter(position=corrected)):
                    i['aa_change'] = _.get()
                models.Nt_mutation.objects.get_or_create(**i, sequence=sequenceobj)

    return


def save_glycosites(sequenceobj):
    pattern = re.compile("N[^P][ST]")

    # check for matches in germline
    germline_starts = []
    germline_sequences = []
    vdj = translate_alignment(sequenceobj, get_germlinealignments(sequenceobj))
    for obj in vdj:
        for match in re.finditer(pattern, str(obj.sequence)):
            germline_starts.append(match.start()+obj.position)
            germline_sequences.append(match.group(0))


    # check for matches in sequence
    sequence_starts = []
    for match in re.finditer(pattern, sequenceobj.AA_Sequence):
        sequence_starts.append(match.start())
        models.Glycosite.objects.get_or_create(sequence=match.group(0), position=match.start(),
                                               germline=True if match.start() in germline_starts else False,
                                               Sequence_Identifier=sequenceobj)


    for pos in [x for x in germline_starts if x not in sequence_starts]:
        if sequenceobj.AA_changes.all().filter(position__range=(pos, pos+2)).exists():
            models.Glycosite.objects.get_or_create(sequence=germline_sequences[germline_starts.index(pos)],
                                                   position=pos,
                                                   germline=True,
                                                   insequence=False,
                                                   Sequence_Identifier=sequenceobj)
    return


def shortcut(sequenceobj):
    with transaction.atomic():
        align_germline(sequenceobj)
        save_mutations(sequenceobj)
        save_glycosites(sequenceobj)
    return
