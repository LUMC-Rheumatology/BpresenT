import json
import os
import re
from django.db import transaction

from Bafstu import settings
from PinkStrawberry import models


def _wrangle_tsv(file, filename):
    _value_with_value_between_brackets = re.compile('^\d+ \(\d+\)')
    _value_between_brackets = re.compile('\(\d+\)$')
    entries, invalid_list = [], []
    fileheaders = _make_headers_compliant(file.readline().rstrip("\t\n").split("\t"))
    for line in file.readlines():
        row = dict()
        for header, value in zip(fileheaders, line.strip("\n").split("\t")):
            if filename in [
                "V_REGION_nt_mutation_statistics",
                "V_REGION_AA_change_statistics",
            ]:
                if re.match(_value_with_value_between_brackets, value):
                    value = re.sub(_value_between_brackets, '', value)
                if len(value) == 0:
                    value = '0'
            if filename == "Summary":
                if value == "X":
                    # only occurs in super poor recognition by V-Quest
                    print("found invalid val:")
                    print(line)
                    print(header, value)
                    invalid_list.append(row["Sequence_ID"])
            if value == '':
                value = None


            row[header] = value
        entries.append(row)
    return entries, invalid_list


def _make_headers_compliant(headers):
    compliant_headers = []
    for header in headers:
        if header.endswith('-') or header.endswith('+'):  # check for those pesky ++-
            if sum(True for _ in header[-3:] if _ in ('+', '-')) == 3:  # make sure it's ONLY those pesky -++
                header = header[:-3] + ''.join(
                    [{'+': 'P', '-': 'N'}[_] for _ in header[-3:] if _ in ('+', '-')])
        for target, change in (('-', '_'),
                               (' ', '_'),
                               ('/', ''),
                               ('(', ''),
                               (')', ''),
                               ('>', '_to_'),
                               ('%', 'percent'),
                               ('5prime', 'fiveprime_'),
                               ('3prime', 'threeprime_'),
                               ("5'", 'fiveprime_'),
                               ("3'", 'threeprime_'),
                               ("__", "_"),):
            header = header.replace(target, change)
        compliant_headers.append(header)
    return compliant_headers


def _remove_entries(VQuest_output_dirty, delete_list):
    print("remlist", sorted(delete_list))
    VQuest_output_clean = dict()
    for file, data in VQuest_output_dirty.items():
        VQuest_output_clean[file] = data.copy()
        print("parsing for deletion", file)
        for row in data:
            if row["Sequence_ID"] in delete_list:
                #print(file, row["Sequence_ID"], '< DEL')
                VQuest_output_clean[file].remove(row)
            else:
                pass
                #print(file, row["Sequence_ID"])
    return VQuest_output_clean

#for django db obj
def database_loader(VQuest_obj):
    VQuest_output = VQuest_object_parser(VQuest_obj)
    with transaction.atomic():
        # load generic sequence data
        for entry in VQuest_output["Summary"]:
            obj = models.Sequence(Sequence_ID=entry["Sequence_ID"], Sequence_number=entry["Sequence_number"],
                                  V_DOMAIN_Functionality=entry["V_DOMAIN_Functionality"],
                                  V_GENE_and_allele=entry["V_GENE_and_allele"], Run=VQuest_obj)
            obj.save()

        # load file data
        for file, data in VQuest_output.items():
            print(f"loading {file} into db")
            try:
                for entry in data:
                    entry["Sequence_Identifier"] = models.Sequence.objects.get(
                        Sequence_number=entry["Sequence_number"], Run=VQuest_obj)
                    # clear match sequences
                    for key in ('Sequence_number', 'Sequence_ID', 'V_DOMAIN_Functionality', 'V_GENE_and_allele'):
                        entry.pop(key)
                    model = getattr(models, file)(**entry)  # unpacking is important here, django save does not handle a dict like the forms do!
                    model.save()
            except ValueError as error:
                print(entry["Sequence_Identifier"].Sequence_ID)
                for key, value in entry.items():
                    print(f"{key:<30} | {value}")
                raise error
    return


#for django db obj
def VQuest_object_parser(VQuest_obj):
    print(f"parsing {VQuest_obj.name}")
    search_files = ("Summary", "IMGT_gapped_nt_sequences", "Nt_sequences", "IMGT_gapped_AA_sequences",
                    "AA_sequences", "Junction", "V_REGION_mutation_and_AA_change_table",
                    "V_REGION_nt_mutation_statistics", "V_REGION_AA_change_statistics", "V_REGION_mutation_hotspots")
    VQuest_output = dict()
    delete_set = set()
    for filename in search_files:
        with getattr(VQuest_obj, filename).open('r') as file:
            VQuest_output[filename], invalid_ids = _wrangle_tsv(file, filename)
        delete_set.update(invalid_ids)

    VQuest_output = _remove_entries(VQuest_output, delete_set)
    for id in delete_set:
        VQuest_obj.addcomment(f"Invalid entry: {id}, missing IMGT region")
    VQuest_obj.save()
    return VQuest_output