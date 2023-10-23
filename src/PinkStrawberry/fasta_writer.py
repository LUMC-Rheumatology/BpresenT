from io import StringIO

def write(sequencequeryset, translated=False):
    fastafile = StringIO(newline='\n')

    for seq in sequencequeryset:
        id = f"{seq.Sequence_ID}"
        cell = f"cell={seq.BCR}"
        patient = f"patient={seq.Patient}"
        germlines = f"{' '.join([reg+'-region='+obj.Gene_and_allele for reg, obj in seq.get_germlines().items() if obj is not None])}"
        germline_ranges = ' '.join([f"{regrange[0]}-range={regrange[1]}-{regrange[2]}" for regrange in seq.get_vdjrange(translated=translated)])

        if translated:
            massspec = f"ms_peptides={[(pep.sequence_start if pep.sequence_start is not None else 'NF', pep.peptide, pep.quality_qval) for pep in seq.Ms_peptide.all()]}"  if seq.Ms_peptide.all().exists() else ''
            sequence = seq.AA_Sequence

        else:
            vrf = seq.Nt_sequences.V_REGION_reading_frame - 1
            massspec = f"ms_peptides={[(pep.sequence_start*3+vrf if pep.sequence_start is not None else 'NF', pep.peptide, pep.quality_qval) for pep in seq.Ms_peptide.all()]}" if seq.Ms_peptide.all().exists() else ''
            sequence = seq.Nt_Sequence

        entry = f">{id} {cell} {patient} {germlines} {germline_ranges} {massspec}\n" \
                f"{sequence}\n"
        fastafile.write(entry)

    fastafile.seek(0)
    response = fastafile.getvalue()
    fastafile.close()
    print(response)
    return response