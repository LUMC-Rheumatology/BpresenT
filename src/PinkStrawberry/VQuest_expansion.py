import itertools
import math
import re
from collections import Counter
from django.db import transaction
from django.core.exceptions import ObjectDoesNotExist
from Bio import Align
from PinkStrawberry import models, Aligner


def link_germline(Sequence_obj):
    fetchcode = lambda x: x.split(" ")[1]
    latest_db = models.IMGT_Gene_DB.objects.latest()
    Summary = Sequence_obj.Summary
    Sequence_obj.V_Germline = models.Germline.objects.get(IMGT_Gene_DB=latest_db, Gene_and_allele=fetchcode(Sequence_obj.V_GENE_and_allele))
    Sequence_obj.J_Germline = models.Germline.objects.get(IMGT_Gene_DB=latest_db, Gene_and_allele=fetchcode(Summary.J_GENE_and_allele))
    if Sequence_obj.isheavy():
        Sequence_obj.D_Germline = models.Germline.objects.get(IMGT_Gene_DB=latest_db, Gene_and_allele=fetchcode(Summary.D_GENE_and_allele))
    Sequence_obj.save()
    return

# shows region but misses germline similarity, use Aligner.save_glycosites instead.
def findglycosite(Sequence_obj):
    pattern = re.compile("N[^P][ST]")
    # get AA sequence of V(D)J-Region
    sequence = _ if (_ := Sequence_obj.AA_sequences.V_D_J_REGION) != None else Sequence_obj.AA_sequences.V_J_REGION
    regions = Sequence_obj.get_regions(aa_mode=True)
    # check for matches
    matches = re.finditer(pattern, sequence)
    for match in matches:
        regionconstruct = set()
        for region, flags in regions.items():
            if flags['start'] <= match.start() < flags['stop']:
                regionconstruct.add(region)
            if flags['start'] <= match.end() < flags['stop']:
                regionconstruct.add(region)
        newobj = models.Glycosite(Sequence=match.group(0), Start=match.start(),
                                  Stop=match.end(), Sequence_Identifier=Sequence_obj,
                                  Region=', '.join(regionconstruct))
        newobj.save()
    return


def determine_heavytype(sequenceobj):
    if sequenceobj.Nt_Constant is None:
        return None

    if len(sequenceobj.Nt_Constant) < 3:
        return None

    aligner = Align.PairwiseAligner(scoring='blastn')
    aligner.right_gap_score = 0
    sequencedict = {"IgG1":"gcctccaccaagggcccatc",
                    "IgG3":"gcttccaccaagggcccatc",
                    "IgA0":"gcatccccgaccagccccaa",
                    "IgM0":"gggagtgcatccgccccaac",
                    }
    scoring = {subtype : aligner.align(sequenceobj.Nt_Constant.upper(), sequence.upper())[0] for subtype, sequence in sequencedict.items()}
    heavytype = max(scoring.items(), key=lambda x: x[1].score)[0][-2]
    return heavytype


def determine_subtype(sequenceobj):
    if sequenceobj.Nt_Constant is None:
        return None

    if len(sequenceobj.Nt_Constant) < 3:
        return None

    aligner = Align.PairwiseAligner(scoring='blastn')
    aligner.right_gap_score = 0
    sequencedict = {"IgG1":"gcctccaccaagggcccatcggtcttccccctggcaccctcctccaagagcacctctgggggcacagcggccctgggctg",
                    "IgG2":"gcctccaccaagggcccatcggtcttccccctggcgccctgctccaggagcacctccgagagcacagccgccctgggctg",
                    "IgG3":"gcttccaccaagggcccatcggtcttccccctggcgccctgctccaggagcacctctgggggcacagcggccctgggctg",
                    "IgG4":"gcttccaccaagggcccatccgtcttccccctggcgccctgctccaggagcacctccgagagcacagccgccctgggctg",
                    }

    scoring = {subtype : aligner.align(sequenceobj.Nt_Constant.upper(), sequence.upper())[0] for subtype, sequence in sequencedict.items()}
    highscore = max(scoring.values(), key = lambda x: x.score).score
    isotypes = [subtype for subtype, aln in scoring.items() if aln.score == highscore]
    return isotypes


# Depricated, full sequences in use now
def _start_aligner():
    global ALIGNER
    ALIGNER = Align.PairwiseAligner()
    ALIGNER.mode = 'local'
    ALIGNER.match_score = 1
    ALIGNER.mismatch_score = -0.2
    ALIGNER.open_gap_score = -1
    ALIGNER.extend_gap_score = -0.5
    ALIGNER.target_end_gap_score = 0.0
    ALIGNER.query_end_gap_score = 0.0
    print("started aligner!")
    return

#continuation with _start_aligner
def depricated_determine_subtype(sequence, type_):
#check if sequence is empty, if so then return None
    if sequence == None:
        return None

    #start aligner
    if "ALIGNER" not in globals():
        _start_aligner()

    predicted = Counter()
    matchrounds = []  # declare to make sure to skip if there's no patterns available
    if type_ == 'G': #could make the if a dict with lookup to be more flexible but for now it's fine
        matchrounds = (({"start": 35, "stop": 52, "pattern": "tcctccaag", "subtype": ["IgG1", "IgG3"]},
                        {"start": 35, "stop": 52, "pattern": "tgctccagg", "subtype": ["IgG2", "IgG4"]}),
                       ({"start": 54, "stop": 67, "pattern": "gggggc", "subtype": ["IgG1"]},
                        {"start": 54, "stop": 67, "pattern": "gagagc", "subtype": ["IgG2", "IgG3", "IgG4"]}),
                       )

    for round_ in matchrounds:
        scores = []
        for params in round_:
            target = sequence[params["start"]:params["stop"]]
            if target == '': #triggers if the sequence isn't long enough
                continue
            alignments = ALIGNER.align(target, params["pattern"])
            alignment = alignments[0]
            #print(f"Score: {alignment.score}\tAlignments: {len(alignments)}", alignment, sep="\n")
            if alignment.score < math.floor(len(params["pattern"])/2):
                continue #triggers if the alignment score was poor (half or less than half of a match
            scores.append((alignment.score, params["subtype"]))

        scores_length = len(scores)
        if scores_length == 1:
            predicted.update(scores[0][1])
        elif scores_length > 1:   # doesn't work for ranking if PATTERN has differing len, could div by PATTERN len for better accuracy.
            scores.sort(key=lambda _: -_[0])
            #Honestly this line is garbage, should be rewritten with python unpacking behavior in loops to be more readable. It's fast tho!
            predicted.update(itertools.chain.from_iterable(map(lambda _: _[1], scores[:sum((1 for score in map(lambda _: _[0], scores) if score == scores[0][0]))])))

    predicted_length = len(predicted)
    if predicted_length == 1:
        return (subtype for subtype, _ in predicted.items())
    elif predicted_length > 1:
        predicted = predicted.most_common()
        return (subtype for subtype, _ in predicted[:sum((1 for _, score in predicted if score == predicted[0][1]))])
    return None

#c4>t, L2; L2 ctg 4-6 [agct 2-5]>L ttg|c59>a,S20>Y(+ - -); S20 tcc 58-60>Y tac|
#index refers to that of the IMGT gapped sequence matched to FR and CDR regions, this func translates it back to the original ungapped index
#DEVNOTE: if this ends up supplementing aligner.py, the positions might need to be transposed to the entire sequence
def parse_mutations(sequence_obj):
    #stuff to correct for IMGT gap magic for nt
    nt_gapcorrector = dict()
    gapcount = 0
    gappedseq = sequence_obj.IMGT_gapped_nt_sequences.V_D_J_REGION \
                if sequence_obj.isheavy() else\
                sequence_obj.IMGT_gapped_nt_sequences.V_J_REGION
    for i, nt in enumerate(gappedseq):
        nt_gapcorrector[i] = gapcount
        if nt == '.':
            gapcount += 1

    #stuff to correct for IMGT gap magic for AA
    aa_gapcorrector = dict()
    gapcount = 0
    gappedseq = sequence_obj.IMGT_gapped_AA_sequences.V_D_J_REGION \
                if sequence_obj.isheavy() else\
                sequence_obj.IMGT_gapped_AA_sequences.V_J_REGION
    for i, nt in enumerate(gappedseq):
        aa_gapcorrector[i] = gapcount
        if nt == '.':
            gapcount += 1

    for region in ("FR1_IMGT", "CDR1_IMGT", "FR2_IMGT", "CDR2_IMGT", "FR3_IMGT", "CDR3_IMGT"):
        regionname = region[:-5]
        mutstr = getattr(sequence_obj.V_REGION_mutation_and_AA_change_table, region)
        if mutstr is None: continue

        for change in mutstr[:-1].split('|'):
            change = [_.split('>') for _ in change.split(";")[0].split(',')]
            aa = None
            #triggers if mutation caused an AA change:
            if len(change) == 2:
                aa_data = change[1]
                #triggers if no change happened
                if len(aa_data) == 1:
                    aa_data[0] = aa_data[0].lstrip(' ')
                    aa = models.IMGT_AA_change.objects.get_or_create(
                        aa_original = aa_data[0][0],
                        aa_pos = int(aa_data[0][1:])-1-aa_gapcorrector[int(aa_data[0][1:])],
                        aa_changed = aa_data[0][0],
                        hydrophobicity = True,
                        volume = True,
                        chemical_characteristics = True,
                        similarity = 4,
                        region = regionname,
                        sequence = sequence_obj,
                    )[0]
                else:
                    aa = models.IMGT_AA_change.objects.get_or_create(
                        aa_original = aa_data[0][0],
                        aa_pos = int(aa_data[0][1:])-1-aa_gapcorrector[int(aa_data[0][1:])],
                        aa_changed = aa_data[1][0],
                        hydrophobicity = (hydrophobicity := True if aa_data[1][2] == "+" else False),
                        volume = (volume := True if aa_data[1][3] == "+" else False),
                        chemical_characteristics = (chemical_characteristics := True if aa_data[1][4] == "+" else False),
                        similarity = sum((hydrophobicity, volume, chemical_characteristics)),
                        region = regionname,
                        sequence = sequence_obj,
                    )[0]
            nt_data = change[0]
            models.IMGT_Nt_mutation.objects.get_or_create(
                nt_original = nt_data[0][0],
                nt_pos = int(nt_data[0][1:])-1-nt_gapcorrector[int(nt_data[0][1:])],
                nt_mutated = nt_data[1],
                sequence = sequence_obj,
                aa_change = aa,
            )


# main calling func
def generate_metadata(VQuest_Run):
    late_removal = []
    query = VQuest_Run.Sequence.all()

    with transaction.atomic():
        for obj, properties in zip(query, query.values(
                'id', 'Nt_sequences__V_REGION_reading_frame', 'Summary__Sequence',
                'AA_sequences__V_D_J_REGION', 'AA_sequences__V_J_REGION', 'Nt_sequences__V_J_REGION',
                'Nt_sequences__V_D_J_REGION',
                'Nt_sequences__V_D_J_REGION_start', 'Nt_sequences__V_J_REGION_start',
                'Nt_sequences__V_D_J_REGION_end', 'Nt_sequences__V_J_REGION_end', 'Sequence_ID')):

            # The sequence is always displayed in the "sense" orientation whatever the orientation at the submission.
            # The results provided in the other spreadsheets of the Excel file correspond to the IMGT/V-QUEST analysis
            # of the sequence in the sense orientation.
            # from: https://www.imgt.org/IMGT_vquest/user_guide

            # split info in sequence ID (Possile Fluff / Patient ID / BCR pair / Type) and make relations
            sequence_info = properties.get("Sequence_ID").replace('__', '_').split("_")
            if len(sequence_info) not in (3, 4):
                VQuest_Run.addcomment(f"Invalid entry: {obj.Sequence_ID}, faulty sequence ID")
                late_removal.append(obj)
                continue

            # bind to BCR and Patient, and create objects for them if they don't already exist
            patientmatch = re.search("PT[0-9]+$", sequence_info[-3], flags=re.IGNORECASE)
            patientvalues = {"Identifier": patientmatch[0]}

            if not models.Patient.objects.filter(**patientvalues).exists():
                newpatient = models.Patient(**patientvalues)
                newpatient.save()
            obj.Patient = models.Patient.objects.get(**patientvalues)

            bcrvalues = {"Identifier": sequence_info[-2], "Run": VQuest_Run}
            if not models.BCR.objects.filter(**bcrvalues).exists():
                newbcr = models.BCR(**bcrvalues)
                newbcr.save()
            obj.BCR = models.BCR.objects.get(**bcrvalues)

            # remove insertions from summary sequence
            Nt_Sequence = ''.join((i if i.lower() == i else '' for i in properties['Summary__Sequence']))

            # translate to amino acids
            AA_Sequence = codontranslate(Nt_Sequence, properties['Nt_sequences__V_REGION_reading_frame'])

            # load translated and cleaned nt sequence into obj
            obj.Nt_Sequence = Nt_Sequence
            obj.AA_Sequence = AA_Sequence

            # check if we're dealing with a heavy or light chain and specify if we need VJ or VDJ region
            region = "V_D_J_REGION" if properties['Nt_sequences__V_D_J_REGION_start'] != None else "V_J_REGION"

            # splice and translate where needed
            reading_frame = properties.get('Nt_sequences__V_REGION_reading_frame')
            VDJ_start = properties.get(
                f'Nt_sequences__{region}_start') - 1  # -1 to compensate for index starting at 0 in python
            VDJ_stop = VDJ_start + (len(properties.get(f'AA_sequences__{region}')) * 3)

            Nt_Leading = Nt_Sequence[:VDJ_start]
            Nt_Constant = Nt_Sequence[VDJ_stop:]
            AA_Leading = codontranslate(Nt_Leading, reading_frame)
            AA_Constant = codontranslate(Nt_Constant)
            obj.AA_Leading = AA_Leading
            obj.AA_Constant = AA_Constant
            obj.Nt_Leading = Nt_Leading
            obj.Nt_Constant = Nt_Constant

            # determine chain type, check it, and load into obj
            rawtype = sequence_info[-1]
            if rawtype.lower() in ('lambda', 'kappa', 'λ', 'κ'):
                cleantype = {'lambda': 'L', 'kappa': 'K', 'λ':'L', 'κ':'K'}[rawtype]
            elif len(rawtype) >= 3: #could write this as regex to cover the edgecases?
                if (subtype := models.Subtype.objects.filter(Name__exact=rawtype)).exists():
                    obj.Subtype.add(subtype[0])
                cleantype = rawtype[2]
            else:
                cleantype = rawtype

            obj.Type = cleantype.upper()

            # check for mismatch in IMGT V_GENE_and_allele
            if 'H' if obj.isheavy() else obj.Type != obj.imgttype():
                if obj.imgttype() in ('L', 'K'):  # True if there's a mismatch and the IMGT chain is light
                    obj.addcomment(f"Chain type mismatch, {obj.Type} -> {obj.imgttype()}")
                    obj.Type = obj.imgttype()
                elif obj.islight():  # True if there's a mismatch in chain type and the input chain was light
                    heavytype = determine_heavytype(obj)
                    if heavytype is None:
                        obj.addcomment(
                            f"Chain type mismatch, heavy chain type could not be determined, {obj.Type} -> X"
                        )
                        obj.Type = 'X'
                    else:
                        obj.addcomment(f"Chain type mismatch, {obj.Type} -> {heavytype}")
                        obj.Type = heavytype

            # check for mismatch in constant domain for Ig(A/G/M)
            if obj.isheavy():
                heavytype = determine_heavytype(obj)
                if obj.Type != heavytype and heavytype is not None:
                    obj.addcomment(f"Chain type mismatch, {obj.Type} -> {heavytype}")
                    obj.Type = heavytype

                # catch the unidentified D-Region germline here because we have to establish proper heavytype first
                if obj.Summary.D_GENE_and_allele is None:
                    VQuest_Run.addcomment(f"Invalid entry: {obj.Sequence_ID}, missing D-Region germline. ")
                    late_removal.append(obj)
                    continue


            #for IgG, check if we can determine it's subtype(s)
            if obj.get_subtype is None and obj.Type == 'G':
                subtypes = determine_subtype(obj)
                if subtypes is not None:
                    for subtype in subtypes:
                        obj.Subtype.add(models.Subtype.objects.get(Name=subtype))

            #link the germline and allele cells to the gene-db
            try:
                link_germline(obj)
            except ObjectDoesNotExist as err:
                VQuest_Run.addcomment(f"Outdated entry: {obj.Sequence_ID}, no matching germline ID found in local Gene-DB. ")
                late_removal.append(obj)
                continue

            # Parse mutations in v-quest files
            parse_mutations(obj)

            # BpresenT align, find mutations and glycosites
            Aligner.shortcut(obj)

            # Save all generated data!
            obj.save()

    #remove all objects that weren't able to be properly sanitized
    for obj in late_removal:
        obj.delete()
    VQuest_Run.save()
    return

def update_germline(vquest_obj):
    for sequence_obj in vquest_obj.Sequence.all():
        link_germline(sequence_obj)

# realistically could be replaced by the Bio.Seq.Seq(foo).translate()
def codontranslate(seq, frameshift = 1, reverse = False):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X', #stopcodon at 3 and 4
        'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W', #stopcodon at 3
    }

    table2 = {
        'CU': 'L', 'GU': 'V', 'UC': 'S', 'CC': 'P',
        'AC': 'T', 'GC': 'A', 'CG': 'R', 'GG': 'G'
    }

    if reverse:
        seq = seq[::-1]
    seq = seq[frameshift-1::].upper()
    tmp = ""
    if len(seq) >= 3:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if len(codon) == 3:
                if codon.__contains__("N"):
                    print("N hit", codon)
                    tmp += table2.get(codon[:2], "X")
                    print(tmp)
                else:
                    tmp += table[codon]
    return tmp