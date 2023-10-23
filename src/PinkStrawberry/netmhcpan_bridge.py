import os
from datetime import datetime
from mhctools import NetMHCpan
from PinkStrawberry import models
from django.db import transaction
import pandas as pd
import numpy as np

# call everything needed to go straight to database
def fullride(sequenceset, mhcalleles, lengths):
    results = apply_predictor(sequenceset, mhcalleles, lengths)
    results = get_best_binders(results)
    prediction_to_db(results)
    return


def find_optimal_mhc(peptideobj, mhcalleles):
    sequence = peptideobj.peptide
    results = predict(mhcalleles, {"peptide":sequence.center(8, 'X')}, max(len(sequence), 8))
    return results.iloc[0]['allele']

# change all known changes back into their original amino acids
def build_unmutated_query(sequence_queryset):
    germline_query = dict()
    for sequence_obj in sequence_queryset:
        user_seq = list(sequence_obj.AA_Sequence)
        for change in sequence_obj.AA_changes.all():
            user_seq[change.position] = change.original
        #get the VXJ range to make sure we don't predict anything outside that
        ranges = sequence_obj.get_vdjrange(translated=True)
        user_seq = user_seq[ranges[0][1]:ranges[-1][2]]
        germline_query[f"GL{sequence_obj.id}"] = ''.join(user_seq).replace('*', 'X')

    return germline_query


def get_best_binders(results, head=30):
    # prepare df to load stuff into
    returndf = pd.DataFrame()

    # get all sequences from the dataframe
    seqlist = np.unique(np.char.lstrip(results['source_sequence_name'].unique().astype(str), "GL"))

    # make dataframe for each sequence
    for seqid in seqlist:
        seqdf = results[(results['source_sequence_name'].str.endswith(seqid))]
        # make dataframe for each unique peptide length
        for length in seqdf.length.unique():
            lendf = seqdf[(seqdf['length'] == length)]
            # get 30 best %rank non-germline peptides
            bestdf = lendf[(~lendf['source_sequence_name'].str.startswith('GL'))].head(head)
            # create a new df that includes the associated germline sequences
            totalseqdf = bestdf.copy()
            for pos, entry in bestdf.iterrows():
                totalseqdf = pd.concat([totalseqdf, lendf[
                                   (lendf['source_sequence_name'].str.startswith('GL')) &
                                   (lendf['offset'] == entry['offset']) &
                                   (lendf['length'] == entry['length']) &
                                   (lendf['allele'] == entry['allele'])
                               ]])
            # add the coupeling of best performing non-germline peptides and their germline versions to output
            returndf = pd.concat([returndf, totalseqdf])
    return returndf.sort_values(by=['offset', 'length', 'peptide'])


def apply_predictor(sequence_queryset, alleles=None, peptide_length=None):
    user_query = {
        sequence_obj.id : (sequence_obj.AA_sequences.V_D_J_REGION if sequence_obj.isheavy() else sequence_obj.AA_sequences.V_J_REGION).replace('*', 'X').center(8, 'X')
        for sequence_obj in sequence_queryset}

    germline_query = build_unmutated_query(sequence_queryset)

    query = user_query | germline_query

    print(f"predicting {len(query)} sequence(s)...")
    results = predict(alleles, query, peptide_length)
    print("done predicting!")
    return results


"""
For a list of allowed alleles, check https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/MHC_allele_names.txt
alleles structured as list: ["HLA-A02:18", HLA-A03:293, ...]
protein_sequences structured as dict: {"1L2Y" : "NLYIQWLKDGGPSSGRPPPS",}
"""
def predict(alleles, protein_sequences, peptide_length=None):
    if peptide_length is None:
        peptide_length = [9, 10]
    if alleles is None:
        alleles = ["A*02:01"]

    predictor = NetMHCpan(alleles=alleles)
    predictor.allow_X_in_peptides = True
    df = predictor.predict_subsequences(protein_sequences, peptide_lengths=peptide_length)\
            .to_dataframe().sort_values(by=['percentile_rank'], ascending=True).astype({'source_sequence_name': 'string'})
    df['allele'] = df['allele'].str.replace('*', '', regex=True)
    return df


def prediction_to_db(df):
    # show timing in terminal
    print("starting transaction")
    total, cur, start = len(df), 0, datetime.now()
    ##
    target_germline_matcher = dict()
    constructkey = lambda pred: str(pred.source_sequence_name).lstrip('GL')+\
                                str(len(pred.peptide))+\
                                str(pred.allele)+\
                                str(pred.offset)
    with transaction.atomic():
        for index, binding_prediction in df.iterrows():
            # show timing in terminal
            cur += 1
            runtime = datetime.now()-start
            if cur % 25 == 0: print(f"{runtime.seconds // 3600:>0{2}}:{(runtime.seconds // 60) % 60:>0{2}}:{runtime.seconds % 60:>0{2}}\t{(cur/total)*100:.2f}%\tpeptides waiting: {len(target_germline_matcher)}")
            # load the prediction into the database, sequence is associated on save

            prediction_obj = models.Prediction.objects.get_or_create(
                **{'peptide': binding_prediction.peptide, 'allele': models.MHC_allele.objects.get(name=binding_prediction.allele),
                   'score': binding_prediction.score, 'percentile_rank': binding_prediction.percentile_rank,
                   'affinity': binding_prediction.affinity, 'offset': 0,
                   'mhcpansequence': binding_prediction.source_sequence_name}
                )[0]

            prediction_obj.offset = prediction_obj.sequence.get_vdjrange(translated=True)[0][1] + binding_prediction.offset
            prediction_obj.save()

            match_identifier = constructkey(binding_prediction)

            #create dict if we haven't seen this one yet
            if match_identifier not in target_germline_matcher:
                target_germline_matcher[match_identifier] = dict()

            if str(binding_prediction.source_sequence_name).startswith('GL'):
                target_germline_matcher[match_identifier] = target_germline_matcher[match_identifier] | {"Germline_peptide":prediction_obj}
            else:
                target_germline_matcher[match_identifier] = target_germline_matcher[match_identifier] | {"Target_peptide":prediction_obj}

            if len(target_germline_matcher[match_identifier].keys()) == 2:
                sequence_obj = models.Sequence.objects.get(id=target_germline_matcher[match_identifier]["Target_peptide"].mhcpansequence)
                comparison = models.Peptide_comparison.objects.get_or_create(
                            **target_germline_matcher[match_identifier],
                            sequence=sequence_obj,
                            )[0]
                comparison.update_diff()
                comparison.save()
                del target_germline_matcher[match_identifier]
        # show timing in terminal
        print(f"{(cur / total) * 100:.2f}%\tpeptides waiting: {len(target_germline_matcher)}")
        ###
    return


if __name__ == '__main__':
    df = predict(["HLA-A02:01", "HLA-A01:01"],
                {"1_IgG4":
                 "VQLVQSAAEVKKPGESLKISCKGSGYNFTNNWIGWVRQMPGNGLEWMGIIDPRGSNTKYSPSLEGQVTFSADKSINTAYLQWSSLKASDTAIYFCARRGGKDNVWGDW",})
