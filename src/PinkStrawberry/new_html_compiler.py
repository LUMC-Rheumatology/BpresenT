from collections import defaultdict
from itertools import cycle
import re
from PinkStrawberry import Aligner

def __process_highlights(sequences, highlights):
    span_pattern = r"data-bs-html='true'>(.*)</span>"

    #toss to list
    for key, value in sequences.items():
        if not isinstance(value, list):
            sequences[key] = list(value)


    #start checking for highlights
    for key, sequence in sequences.items():
        if not highlights.get(key):
            continue

        for high in highlights.get(key, []):
            tooltiped = False

            if high['loc'] < 0:
                continue
            if sequence[high['loc']] == ' ':
                continue
            # location
            loc = high['loc']
            #colour
            colour = high['colour']
            #text
            text = high.get('text')
            #opacity
            opacity = high.get('opacity', '')

            #processing
            # check if there's already a span and work splice inside of it to prevent overlap later
            # todo rewrite this entire visualization to something dataclass / lambda based, this is WET enough for a slip 'n slide
            if sequence[loc].startswith("<span"):
                rematch = re.search(span_pattern, sequence[loc])
                # colouring
                sequence[loc] = sequence[loc][:rematch.start(1)] + \
                                f"<div style='background-color:{colour}{opacity};'>{sequence[loc][rematch.start(1):rematch.end(1)]}</div>" + \
                                sequence[loc][rematch.end(1):]
                if text:
                    sequence[loc] = sequence[loc][:82] + text + '<br>' + sequence[loc][82:]
            #colouring
            else:
                sequence[loc] = f"<div style='background-color:{colour}{opacity};'>{sequence[loc]}</div>"
                if text:
                    sequence[loc] = f"<span class='d-inline-block' tabindex='0' data-bs-toggle='tooltip'" \
                                    f" data-bs-title='{text}' data-bs-html='true'>{sequence[loc]}</span>"
    return sequences


def __rows_order(rows):
    rows_order = ('header', '>sequences', 'stepper')
    ordered_rows = dict()

    for rkey in rows_order:
        #dump sequences
        if rkey == '>sequences':
            for key, val in rows.items():
                if key.startswith('>'):
                    ordered_rows[key] = val
        #other vals
        elif (_ := rows.get(rkey)):
            ordered_rows[rkey] = _
    return ordered_rows


def visualizer_shortcut(sequenceobj, germlines=False, translated=False, linebr=80, predictions=None, massspec=None,
                        start=0, stop=None, region=None):
    sequenceid = '>'+sequenceobj.Sequence_ID
    header = construct_header(sequenceobj, translated=translated)

    # highlight mutations and glycosites
    seqhighlights = construct_mutation_highlights(sequenceobj, translated=translated) + \
                    construct_glycosite_highlights(sequenceobj, translated=translated)

    # get amino-acid / nucleotide specific data that can't be obtained with the translate argument
    if translated:
        vdjgerm = Aligner.translate_alignment(sequenceobj, Aligner.get_germlinealignments(sequenceobj))
        sequencestr = sequenceobj.AA_Sequence
    else:
        vdjgerm = Aligner.get_germlinealignments(sequenceobj)
        sequencestr = sequenceobj.Nt_Sequence

    sequences = {sequenceid: sequencestr}
    highlights = {sequenceid: seqhighlights}


    # set range to render and specify germline
    if region is not None:
        start, stop = [(start, stop) for reg, start, stop in sequenceobj.get_vdjrange(translated=translated) if reg == region][0]
        if germlines:
            vdjgerm = [aln for aln in vdjgerm if aln.name == f'{region}_REGION']
            assert len(vdjgerm) == 1
            germlinerange = vdjgerm[0].range()
            start = min(start, germlinerange[0])
            stop = max(stop, germlinerange[1])


    # add germlines if provided
    if germlines:
        vdjseq, vdjhigh = fuse_alignments(vdjgerm)
        sequences |= vdjseq
        highlights |= vdjhigh

    # add mass-spec petides if provided
    if massspec:
        massseq, masshigh = fuse_massspec(massspec)
        if not translated:
            massseq, masshigh = translated_transpose(sequenceobj, massseq, masshigh)
        sequences |= massseq
        highlights |= masshigh

    # add predicted peptides if provided
    if predictions:
        pepseq, pephigh = fuse_predictions(predictions)
        if not translated:
            pepseq, pephigh = translated_transpose(sequenceobj, pepseq, pephigh)
        sequences |= pepseq
        highlights |= pephigh

    # call the visualizer function
    formatted_sequence = visualize_to_html(sequences=sequences, highlights=highlights, header=header, linebr=linebr,
                                           start=start, stop=stop)

    return formatted_sequence


# transpose amino acid sequence and related highlight data to fit under nucleotide sequences
def translated_transpose(sequenceobj, seqdict, highdict=None):
    vframe = sequenceobj.Nt_sequences.V_REGION_reading_frame-1
    seqdict = {key:vframe*' ' + '  '.join(value) for key, value in seqdict.items()}
    if not highdict:
        return seqdict

    for key, value in highdict.items():
        for highargs in value:
            highargs.update({'loc':highargs['loc']*3+vframe})
    return seqdict, highdict


def construct_header(sequenceobj, translated=False):
    constructstr = ''
    indicator = '.'
    get_indicator = lambda: '-' if indicator == '.' else '.'

    ranges = [('V', sequenceobj.Nt_sequences.V_REGION_start-1, sequenceobj.Nt_sequences.V_REGION_end - 1),
              ('J', sequenceobj.Nt_sequences.J_REGION_start - 1, sequenceobj.Nt_sequences.J_REGION_end - 1)]
    if sequenceobj.isheavy():
        ranges.insert(1, ('D', sequenceobj.Nt_sequences.D_REGION_start-1, sequenceobj.Nt_sequences.D_REGION_end - 1))

    if translated == True:
        ranges = [(reg, *Aligner.translate_range(start, end, sequenceobj.Nt_sequences.V_REGION_reading_frame - 1)) for reg, start, end in ranges]

    for region, start, end in ranges:
        indicator = get_indicator()
        tmp = f"{constructstr[:start]:<{start}}{region}"
        constructstr = f"{tmp:{indicator}<{end}}|"

    return constructstr


def construct_mutation_highlights(sequenceobj, translated=False):
    highlights = []
    if translated:
        for obj in sequenceobj.AA_changes.all():
            related = ', '.join([f"#{_.position} {_.original} > {_.changed}" for _ in obj.mutations.all()])
            highlights.append(
                {"loc": obj.position, "colour": "#FFD4B8", "opacity":90,
                 "text": f"#{obj.position} {obj.original} > {obj.changed} ({related})"})
    else:
        for obj in sequenceobj.Nt_mutations.all():
            if (_:=obj.aa_change):
                related = f"(#{_.position} {_.original} > {_.changed})"
            else:
                related = ''
            highlights.append(
                {"loc":obj.position, "colour":"#FFD4B8", "opacity":90,
                 "text":f"#{obj.position} {obj.original} > {obj.changed} {related}"})

    return highlights


def construct_glycosite_highlights(sequenceobj, translated=False):
    highlights = []
    if translated:
        for obj in sequenceobj.Glycosite.all():
            colour = ("#97E5D7" if obj.insequence else "#7bd6ed") if obj.germline else "#A8FFEF"
            tooltip = f"#{obj.position} {obj.sequence}{(' also germline' if obj.insequence else ' only germline') if obj.germline else ''}"
            for i in range(len(obj.sequence)):
                highlights.append(
                    {"loc": obj.position+i, "colour": colour, "text": tooltip})
    else:
        vframe = sequenceobj.Nt_sequences.V_REGION_reading_frame - 1
        for obj in sequenceobj.Glycosite.all():
            colour = ("#97E5D7" if obj.insequence else "#7bd6ed") if obj.germline else "#A8FFEF"
            tooltip = f"#{obj.position} {obj.sequence}{(' also germline' if obj.insequence else ' only germline') if obj.germline else ''}"
            for i in range(len(obj.sequence)*3):
                highlights.append(
                    {"loc": obj.position*3+i+vframe, "colour": colour, "text": tooltip})

    return highlights


# get multiple querysets on one line
def fuse_predictions(queryset):
    sequences = defaultdict(str)
    highlights = defaultdict(list)
    pos_tracker = defaultdict(set)
    # loop over predictions of peptides in queryset
    for obj in queryset:
        # get related germline prediction

        # check if we can add it to the current line or need a new one
        pepkey = ''
        positions = set(obj.offset+i for i in range(0,len(obj.peptide)))
        while pos_tracker[pepkey] & positions:
            if not pepkey:
                pepkey = 1
            pepkey += 1
        pos_tracker[pepkey] |= positions
        # produce the highlights
        for loc in positions:
            highlights[f">prediction {pepkey}"].append(
                {"loc": loc, "colour": 'None',
                 "text": f"Allele: {obj.allele}<br>"
                         f"Affinity: {obj.affinity}<br>"
                         f"Percentile rank: {obj.percentile_rank}"}
            )
        # insert the peptide into the sequence, build if needed
        tmpstr = sequences[f">prediction {pepkey}"]
        sequences[f">prediction {pepkey}"] = f"{tmpstr[:obj.offset]:<{obj.offset}}{obj.peptide}{tmpstr[obj.offset+len(obj.peptide):]}"

    # sort highlights to add colours that help differentiate neighbours
    highlights = {key:sorted(value, key=lambda x: x["loc"]) for key, value in highlights.items()}
    #start the colour loop
    colourgen = cycle(['#C4DFDF', '#E3F4F4', '#D2E9E9', '#F8F6F4'])
    for highlightlists in highlights.values():
        lasttext, colour = None, None
        for highlight in highlightlists:
            if highlight["text"] != lasttext:
                lasttext = highlight["text"]
                colour = next(colourgen)
            highlight["colour"] = colour
    return sequences, highlights



# basically the same as fuse_predictions, just altered for mass-spec peptides.
# Should really turn this whole 3x copything into a  function with lambdas to make highlights and colouring more flexible
def fuse_massspec(queryset):
    sequences = defaultdict(str)
    highlights = defaultdict(list)
    pos_tracker = defaultdict(set)
    # loop over predictions of peptides in queryset
    for obj in queryset:
        # get related germline prediction

        # check if we can add it to the current line or need a new one
        pepkey = ''
        positions = set(obj.sequence_start+i for i in range(0,len(obj.peptide)))
        while pos_tracker[pepkey] & positions:
            if not pepkey:
                pepkey = 1
            pepkey += 1
        pos_tracker[pepkey] |= positions
        # produce the highlights
        for loc in positions:
            highlights[f">mass-spec {pepkey}"].append(
                {"loc": loc, "colour": 'None',
                 "text": f"Allele: {obj.allele}<br>"
                         f"Confidence: {obj.confidence}<br>"
                         f"Quality q-val: {obj.quality_qval}"}
            )
        # insert the peptide into the sequence, build if needed
        tmpstr = sequences[f">mass-spec {pepkey}"]
        sequences[f">mass-spec {pepkey}"] = f"{tmpstr[:obj.sequence_start]:<{obj.sequence_start}}{obj.peptide}{tmpstr[obj.sequence_start+len(obj.peptide):]}"

    # sort highlights to add colours that help differentiate neighbours
    highlights = {key:sorted(value, key=lambda x: x["loc"]) for key, value in highlights.items()}
    #start the colour loop
    colourgen = cycle(['#DFCCFB', '#BEADFA', '#D0BFFF'])
    for highlightlists in highlights.values():
        lasttext, colour = None, None
        for highlight in highlightlists:
            if highlight["text"] != lasttext:
                lasttext = highlight["text"]
                colour = next(colourgen)
            highlight["colour"] = colour
    return sequences, highlights


# slightly modified fuse_predictions, should rewrite it to a more general function with lambdas if possible
def fuse_alignments(alignmentset, colours = None):
    sequences = defaultdict(str)
    highlights = defaultdict(list)
    pos_tracker = defaultdict(set)
    # loop over alignments
    for obj in alignmentset:
        # check if we can add it to the current line or need a new one
        pepkey = ''
        positions = set(obj.position+i for i in range(0,len(obj.sequence)))
        while pos_tracker[pepkey] & positions:
            if not pepkey:
                pepkey = 1
            pepkey += 1
        pos_tracker[pepkey] |= positions
        # produce the highlights
        for loc in positions:
            highlights[f">germline {pepkey}"].append(
                {"loc": loc, "colour": 'None',
                 "text": f"{obj.name}"}
            )
        # insert the alignment into the sequence, build if needed
        tmpstr = sequences[f">germline {pepkey}"]
        sequencestr = obj.sequence if obj.position >= 0 else obj.sequence[abs(obj.position):]
        sequencepos = obj.position if obj.position >= 0 else 0
        sequences[f">germline {pepkey}"] = f"{tmpstr[:sequencepos]:<{sequencepos}}{sequencestr}{tmpstr[sequencepos+len(sequencestr):]}"

    # sort highlights to add colours that help differentiate neighbours
    highlights = {key:sorted(value, key=lambda x: x["loc"]) for key, value in highlights.items()}
    #start the colour loop
    if colours is None:
        colours = ['None']
    colourgen = cycle(colours)
    for highlightlists in highlights.values():
        lasttext, colour = None, None
        for highlight in highlightlists:
            if highlight["text"] != lasttext:
                if highlight["text"][0] in ('V', 'D', 'J'):
                    colour = {'V': '#C1EFE7', 'D': '#FDF7EB', 'J': '#E4F3E8'}[highlight["text"][0]]
                else:
                    lasttext = highlight["text"]
                    colour = next(colourgen)
            highlight["colour"] = colour
    return sequences, highlights


# construct a 'html' version of the given sequences and highlights
def visualize_to_html(sequences, highlights=None, header=None, linebr=50, stepper=5, start=0, stop=None):
    rows = dict()

    #load sequences into dict
    for ind, sequence in sequences.items():
        ind = ind if ind[0] == '>' else f">{ind}"
        rows[f"{ind}"] = sequence

    #create numbering stepper
    maxlen = max(map(len, sequences.values()))
    maxseqidlen = max(map(len, sequences.keys())) + 2
    if stepper > 0:
        rows['stepper'] = ''.join([f"{i:<{stepper}}" for i in range(0, maxlen, stepper)])

    #add header if provided
    if header:
        rows['header'] = header

    #process highlights
    if highlights:
        rows = __process_highlights(rows, highlights)

    #get the stop position if it hasn't been defined
    if stop is None:
        stop = maxlen

    formatted = ''
    rows = __rows_order(rows)
    for i in range(start, stop, linebr):
        for key, row in rows.items():
            # slice from original value
            limit = stop if i+linebr > stop else i+linebr
            split = ''.join(row[i: limit])
            # skip if it's an empty line
            if set(split) == {' '} or split == '':
                continue
            # name row for sequences
            if key.startswith('>'):
                split = f"{key:<{maxseqidlen}}{split}"
            else:
                split = f"{maxseqidlen*' '}{split}"
            #finish line format
            formatted += f"{split}\n"
        formatted += "\n"
    return formatted


if __name__ == '__main__':
    vis = visualize_to_html(
        sequences= {">seq1":"--------------gcctccaccaagggcccatcggtcttccccctggcaccctcctccaagagcacctctggggcacagcggccctgggctgcctggtcaaggactgcctccaccaagggcccatcggtcttccccctggcaccctcctccaagagcacctctggggcacagcggccctgggctgcctggtcaaggact----------------------",
                    ">seq2":"cgaactgtggctgcaccatctgtcttcatcttcccgccatctgatgagcagttgaaatctggaactgcctctgttgtgtgcctgctgaataacttctatcccagagaggccaaagtacagtggaaggtggataacgccctccaatcgggtaactcccaggagagtgtcacagagcaggacagcaaggacagcacctacagcctcagcagcaccctgaccgct"},
        highlights= {">seq1":[{"loc":0, "colour":'yellow', "text":'n'},{"loc":155, "colour":'cyan', "text":'n'}]}
    )
    print(vis)