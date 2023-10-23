from urllib import request
from PinkStrawberry import models
from django.db import transaction

"""
The IMGT/GENE-DB FASTA header contains 15 fields separated by '|':

1. IMGT/LIGM-DB accession number(s)
2. IMGT gene and allele name
3. species (may be followed by an "_" and the name of the strain, breed or isolate, if defined)
4. IMGT gene and allele functionality
5. exon(s), region name(s), or extracted label(s)
6. start and end positions in the IMGT/LIGM-DB accession number(s)
7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
8. codon start, or 'NR' (not relevant) for non coding labels
9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
12. number of amino acids (AA): this field indicates that the sequence is in amino acids
13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
14. partial (if it is)
15. reverse complementary (if it is)
"""
def fetch():
    with request.urlopen('https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP') as response_aa_db, \
         request.urlopen("https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP") as response_nt_db, \
         request.urlopen('https://www.imgt.org/download/GENE-DB/RELEASE') as response_release:
        aa_db = response_aa_db.read().decode("utf-8")
        nt_db = response_nt_db.read().decode("utf-8")
        release = response_release.read().decode("utf-8")

        #read through nt_db first to get the nucleotide sequences
        nt_db = parse_genedb_fasta(nt_db)

        #read through aa_db to make the proper objects from that
        main_db = parse_genedb_fasta(aa_db)

        #loop through aa_db to match against
        for head in main_db.keys():
            main_db[head]['AA_Sequence'] = main_db[head]['Sequence']
            main_db[head]['Nt_Sequence'] = nt_db[head]['Sequence']
            del main_db[head]['Sequence']

        return release, main_db.values()


#this one was written for the AA db but for a proper overview it's also used on the nt database
def parse_genedb_fasta(flat_fasta):
    generator = iter(flat_fasta.split("\n"))
    entries = dict()
    sequence, data = '', dict()
    header = generator.__next__().split("|")  # make sure first header is already in there

    for entry in generator:
        if entry.startswith(">") or entry == "":
            # conclude last sequence is finished and add it to the list
            if header[2] == "Homo sapiens":
                for key, value in zip(("Accesion_number", "Gene_and_allele", "Species",
                                       "Functionality", "Labels", "Accesion_start_end",
                                       "Accesion_nucleotide_nb", "Codon_start",
                                       "Fiveprime_changed", "Threeprime_changed",
                                       "Sequencing_changed", "AA_nb", "Length",
                                       "Partial", "Reverse_complementary"), header):
                    data[key] = value
                data['Sequence'] = ''.join((i if i != '*' else 'X' for i in sequence))
                entries['|'.join(header[:8])] = data

                data = dict()
                sequence = ''

            # header of new sequence
            header = entry[1:].split("|")
        else:
            sequence += entry

    return entries


def to_db(release, entries):
    # create new IMGT_Gene_DB object if there's been a new release to make sure old links don't get lost
    releaseobj = models.IMGT_Gene_DB.objects.get_or_create(release=release)[0]

    with transaction.atomic():
        for entry in entries:
            models.Germline.objects.get_or_create(**entry, IMGT_Gene_DB=releaseobj)
    return


def update():
    a = fetch()
    to_db(*a)
    return
