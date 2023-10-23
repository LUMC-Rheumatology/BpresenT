import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import chain, combinations, cycle
from Bio import Align
import numpy as np
import base64
from io import BytesIO
from django.db.models import Count
from django.db.models.functions import Substr, Length
from PinkStrawberry import models

plt.switch_backend('AGG') # switch to prevent screen output

def getgraph():
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches="tight", dpi=300, transparent=True)
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png)
    graph = graph.decode('utf-8')
    buffer.close()
    plt.close()
    return graph


def placeholder():
    return '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'


def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_format


def peptidedistribution(patientobj):
    #fetch related peptides
    predicted_found = models.Predicted_Found.objects.filter(sequence__Patient=patientobj)
    predicted = models.Prediction.objects.filter(sequence__Patient=patientobj, germline=False).exclude(id__in=predicted_found.values_list('predicted', flat=True))
    massspec = models.Ms_peptide.objects.filter(sequence__Patient=patientobj).exclude(id__in=predicted_found.values_list('found', flat=True))
    # set up the data for plotting
    labels = 'both', 'predicted', 'found',
    data = [len(predicted_found), len(predicted), len(massspec)]
    # check if we actually have enough data
    if sum(data) == 0:
        return placeholder()
    #generate plot and turn into datastream
    plt.figure()
    plt.pie(data, labels=labels, textprops={'weight':'bold'},
       wedgeprops={"linewidth": 0.5, "edgecolor": "white"},  autopct=autopct_format(data))
    #plt.tight_layout()
    graph = getgraph()
    return graph


def isotypedistribution(patientobj):
    #fetch related peptides
    queryout = (models.Sequence.objects.filter(Patient=patientobj)
              .values('Type')
              .annotate(tcount=Count('Type'))
              .order_by()
              )
    # set up the data for plotting
    labels = [_ if (_:=i['Type']) in ('K', 'L') else f"Ig{_}" for i in queryout]
    data = [i['tcount'] for i in queryout]
    #generate plot and turn into datastream
    plt.figure()
    plt.pie(data, labels=labels, textprops={'weight':'bold'},
       wedgeprops={"linewidth": 0.5, "edgecolor": "white"},  autopct=autopct_format(data))
    graph = getgraph()
    return graph


def germlinedistribution(patientobj, isotype):
    if isotype not in ('H', 'L', 'K'):
        raise Exception("given isotype unknown")

    graphlist = []
    if isotype == 'H':
        querystart = models.Sequence.objects.filter(Patient=patientobj, D_Germline__isnull=False)
    elif isotype in ('L', 'K'):
        querystart = models.Sequence.objects.filter(Patient=patientobj, D_Germline__isnull=True)\
            .exclude(id__in=[i.id for i in models.Sequence.objects.filter(Patient=patientobj, D_Germline__isnull=True) if i.imgttype() != isotype])

    for region in ('V', 'D', 'J') if isotype == 'H' else ('V', 'J'):
        #fetch related data
        queryout = (querystart
                  .values(f'{region}_Germline__Gene_and_allele')
                  .annotate(gene=Substr(f'{region}_Germline__Gene_and_allele', 1, Length(f'{region}_Germline__Gene_and_allele')-3))
                  .annotate(gcount=Count(f'{region}_Germline__Gene_and_allele'))
                  .order_by(f'{region}_Germline__Gene_and_allele')
                  )

        fig, ax = plt.subplots()

        keys = [i[f'{region}_Germline__Gene_and_allele'] for i in queryout]
        data = [i['gcount'] for i in queryout]

        ax.grid(zorder=0, axis='y')

        fontsizer = min(270/max((len(queryout),1)), 13)

        plt.xticks(rotation=50, ha='right', fontsize=fontsizer)
        ax.bar(keys, data, color=mcolors.TABLEAU_COLORS, width=0.8, align='center', zorder=3)
        ax.set_ylabel('occurances')

        graphlist.append((getgraph(), region, sum(queryout.values_list('gcount', flat=True))))
    return graphlist


def mutation_glycosilationdistribution(patientobj):
    mutationcount, glycosilationcount = [], []
    for i in patientobj.Sequence.all():
        mutationcount.append(i.Nt_mutations.count())
        glycosilationcount.append(i.Glycosite.count())

    plt.figure(figsize=(8, 1.5))
    plt.scatter(mutationcount, glycosilationcount)
    plt.xlabel('#Mutations')
    plt.ylabel('#Glycosylation')
    plt.yticks(np.arange(0, max(glycosilationcount)+1, 1))

    graph = getgraph()
    return graph


###### pairwise comparison that was introduced too late and never saw completion ######
def patient_similarity():
    aligner = Align.PairwiseAligner(scoring='blastp')
    aligner.mode = 'global'

    counter = 0
    patients = models.Patient.objects.all()[1:4]
    for pt1, pt2 in combinations(patients, 2):
        for pt1seqdata in pt1.Sequence.all():
            pt1seq = (_ if (_ := pt1seqdata.AA_sequences.V_D_J_REGION) else pt1seqdata.AA_sequences.V_J_REGION).upper()
            print(pt1.Identifier, pt1seqdata.Sequence_ID)
            for pt2seqdata in pt2.Sequence.all():
                #check H/V/K same
                if pt2seqdata.imgttype() != pt1seqdata.imgttype():
                    continue
                pt2seq = (_ if (_ := pt2seqdata.AA_sequences.V_D_J_REGION) else pt2seqdata.AA_sequences.V_J_REGION).upper()
                alignment = aligner.align(pt1seq, pt2seq)[0]
                scoreadjust = alignment.score / ((len(pt1seq)+len(pt2seq)/2))
                print(f"{counter} | {alignment.score:>{6}} | {scoreadjust:>5.2f} | {matched(alignment):.2f}% | {pt1seqdata.Sequence_ID} & {pt2seqdata.Sequence_ID}")
                #print(alignment)
                counter += 1
            print(end='\n')


def matched(alignobj):
    target = build_alinged(alignobj.target, alignobj.indices[0])
    query  = build_alinged(alignobj.query, alignobj.indices[1])
    miss = 0
    hit = 0
    for t, q in zip(target, query):
        if '-' in (t, q):
            continue
        if t == q:
            hit += 1
        else:
            miss += 1
    return hit / (hit+miss) * 100


def build_alinged(seq, ind):
    string = ''
    for i in ind:
        string += seq[i] if i != -1 else '-'
    return string