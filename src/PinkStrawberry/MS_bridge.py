import re
from PinkStrawberry import models
from django.db import transaction
from openpyxl import load_workbook

def associate_peptides():
    pass

def read_msexcel(file, worksheets=None):
    #worksheets to read
    if not worksheets:
        worksheets = ['DATA_received_CPM', ]

    wb = load_workbook(filename = file, data_only = True)
    main_dict = dict()
    for ws in wb.worksheets:
        if ws.title not in worksheets:
            continue
        main_dict[ws.title] = dict()
        subtable_title = False #declare to make sure we toss everything before the first subtable
        generator_func = ws.values
        for row in generator_func:
            #check for a title value
            (_ := set(row)).discard(None)
            if len(_) == 1:
                subtable_title = _.pop()
                headers = next(generator_func)
                main_dict[ws.title][subtable_title] = []
            if len(_) > 1 and subtable_title:
                main_dict[ws.title][subtable_title].append(dict(zip(headers, row)))
    return main_dict


def to_db(excel):
    with transaction.atomic():
        for _, worksheet in excel.items():
            for title, subtable in worksheet.items():
                for entry in subtable:
                    dbkwargs = {'peptide': entry.get('Sequence'),
                                'quality_qval':entry.get('Qvality q-value'),
                                'Comment':entry.get('Modifications')}
                    dbkwargs['confidence'] = {'High':2, 'Medium':1, 'Low':0}[entry.get('Confidence')]
                    for masterprotein in entry.get('Positions in Master Proteins').split(';'):
                        masterprotein = masterprotein.strip()
                        search = re.search(r'\[.+]', masterprotein)
                        pos = search.group(0)[1:-1].split('-')
                        dbkwargs['start'] = pos[0]
                        dbkwargs['stop'] = pos[1]
                        dbkwargs['sheet_sequence'] = masterprotein[:search.start()]
                        models.Ms_peptide.objects.get_or_create(**dbkwargs)
    return
