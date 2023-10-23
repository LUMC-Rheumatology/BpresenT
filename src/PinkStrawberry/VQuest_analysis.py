import pathlib
import os
import re
from collections import defaultdict, Counter


def parse_tsv(path):
    with open(path) as file:
        entries = []
        fileheader = file.readline().rstrip("\t\n").split("\t")
        for line in file.readlines():
            data = dict()
            for cell, col in zip(line.strip("\n").split("\t"), fileheader):
                data[col] = cell
            entries.append(data)
    return entries


def IMGT_output_parser(path):
    pattern = re.compile("^(?!11)[0-9]([0-9]|_)_?")  # regex explainer: https://regex101.com/r/YLVliU/1
    IMGT_output = dict()
    for link in os.listdir(path):
        if re.match(pattern, link):
            filename = re.sub(pattern, '', link)
            IMGT_output[filename] = parse_tsv(pathlib.Path.joinpath(path, link))

    return IMGT_output


def visualize_IMGT_output_parser(IMGT_output: dict):
    for file, data in IMGT_output.items():
        print(f"{file}")
        for key, value in data[0].items():
            print(f"{key}: {value}")
        print("\n\n")


def check_line_synchronisation(IMGT_output):
    Summary = IMGT_output.pop("Summary.txt")
    length = len(Summary)

    for line in range(0, length):
        check = [Summary[line].get(key).replace(" (see comment)", "") for key in ("Sequence number", "Sequence ID", "V-DOMAIN Functionality", "V-GENE and allele")]

        for file, data in IMGT_output.items():
            compare = [data[line].get(key) for key in ("Sequence number", "Sequence ID", "V-DOMAIN Functionality", "V-GENE and allele")]

            if check != compare:
                print(f"{file} has a mismatch on {compare} with check {check}")


def headercount(IMGT_output):
    counter = defaultdict(lambda: [0, []])
    for file, data in IMGT_output.items():
        for key in data[0].keys():
            counter[key][0] += 1
            counter[key][1].append(file)

    counter = [(key, value[0], value[1]) for key, value in counter.items()]
    counter.sort(key=lambda x: x[1])

    for header in counter:
        print(header)


def headercount_withvalues(IMGT_output):
    counter = defaultdict(lambda: [0, [], []])
    for file, data in IMGT_output.items():
        for key, value in data[0].items():
            counter[key][0] += 1
            counter[key][1].append(file)
            counter[key][2].append(value)

    counter = [(key, value) for key, value in counter.items()]
    counter.sort(key=lambda x: x[1][0])

    for header, value in counter:
        formatlen = max((len(_) for _ in value[1]))+2

        print(f"{value[0]}, {header}")
        for file, cell in zip(value[1], value[2]):
            print(f"{file:<{formatlen}}: {cell}")
        print(end="\n\n")


def modelmaker(IMGT_output, show_reasoning=False, sanitize_data=True):
    floatpattern = re.compile("^\d+\.\d+$")
    datapredictions = dict()

    #store headers for easier access and proper commenting after sanitizing
    headers = dict()
    for key, data in IMGT_output.items():
        headers[key] = list(data[0].keys())

    if sanitize_data:
        IMGT_output = sanitize_output(IMGT_output)

        newheaders = dict()
        for key, value in headers.items(): #make sure the titles still match up
            newkey = key[:-4].replace('-', '_').replace(' ', '_')
            newheaders[newkey] = value
        headers = newheaders

    #start actual prediction
    for file, data in IMGT_output.items():
        datapredictions[file] = dict()
        for key in data[0].keys():
            occurances = {"str": {"counter":0, "example":''}, "int": {"counter":0, "example":''},
                          "float": {"counter":0, "example":''}, "bool": {"counter":0, "example":''},
                          "none": {"counter":0, "example":''}}
            occurances["str"]["length"] = 0 #make sure we catch str length too

            for row in data:
                cell = row[key]
                if cell == '' or cell == None:
                    occurances["none"]["counter"] += 1
                elif cell.isdecimal():
                    if cell in ('0', '1'):  # <!> It's never a proper bool if it's a decimal, manually checked
                        occurances["bool"]["counter"] += 1
                        occurances["bool"]["example"] = cell
                    else:
                        occurances["int"]["counter"] += 1
                        occurances["int"]["example"] = cell
                elif re.match(floatpattern, cell):
                    occurances["float"]["counter"] += 1
                    occurances["float"]["example"] = cell
                else:
                    occurances["str"]["counter"] += 1
                    if len(cell) > occurances["str"]["length"]: #catch longest string length for CharField
                        occurances["str"]["length"] = len(cell)
                        occurances["str"]["example"] = cell

            #conclude probable datatype
            predictedtype = None

            for candidate, field in {"float":"FloatField()", "int":"IntegerField()", "bool":"BooleanField()", "str":"TextField()"}.items():
                if occurances[candidate]["counter"] > 0:
                    predictedtype = field
                    break

            #translate to django model datatype
            if predictedtype == "TextField()":
                if occurances["str"]["length"] == 1:
                    predictedtype = "CharField(1)"

                ##making this dynamic gave me anxiety so it's just all TextField now.
                #if maxlen <= 50:
                #    predictedtype = f"CharField(max_length={maxlen+10 if maxlen != 1 else maxlen})"

            if predictedtype == None:
                predictedtype = "none"
            else:
                predictedtype = 'models.'+predictedtype

            #store prediction for display purposes
            datapredictions[file][key] = (occurances, predictedtype)

    for file, data in datapredictions.items():
        #start showing file
        print(f"READIN --> {file}")

        #strip first 4 cols as they're all sequence identifiers and are replaced by the sequence table ID ref
        data = list(data.items())[4:] #side effect is that we lose key lookup but we don't use it anyway, could call dict() around it if we end up needing it
        header_iter = 4 #start keeping up with header dict position

        #show generated model
        print("\ngenerated model")
        print(f"# model for file {file}")
        print(f"class {file}(models.Model):")
        print(f"\tSequence_Identifier = models.ForeignKey(Sequence, on_delete=models.CASCADE, related_name='{file}')")
        for key, value in data:
            field = value[1]
            topad = f"\t{key} = {field}"
            original_column = " # "+headers[file][header_iter] if key != headers[file][header_iter] else ''
            print(f"{topad:<55}{original_column}")
            header_iter += 1 #trust me it's better than eumerate at this point
        print(f"\n\tclass Meta:\n"
              f"\t\tget_latest_by = 'id'")

        if show_reasoning:
            input()
            print("\nreasoning")
            # split into rows of 2 columns
            x1 = []
            x2 = []
            i = 0
            for entry in data: #this is dumb at this stage but the commitment has been made and I'm not rewriting it
                if i == 2:
                    x1.append(x2)  # .copy() probably isn't needed but python is weird with list refs
                    x2 = []
                    i = 0
                x2.append(entry)
                i += 1
            x1.append(x2)
            #outline format of reasoning
            for splits in x1:
                transformedlists = ['' for _ in range(7)]  # force 6 dimensions because of header + 5 possible datatypes
                for tups in splits:
                    transformedlists[0] += f"{tups[0]:<60}"
                    transformedlists[1] += f"Conclusion: {tups[1][1]:<48}"
                    for index, key in enumerate(("str", "int", "float", "bool", "none"), start=2):
                        transformedlists[index] += f"{key:<5} | {tups[1][0][key]['counter']:<3} | {tups[1][0][key]['example'][:45]:<46}"

                #print reasoning
                for line in transformedlists:
                    print(line)
                print(end="\n")
        print(end="\n\n\n")
        input() #wait to show next file


def sanitize_output(IMGT_output):
    clean_IMGT_output = dict()
    value_with_value_between_brackets = re.compile('^\d+ \(\d+\)')
    value_between_brackets = re.compile('\(\d+\)$')
    delete_list = []

    for file, data in IMGT_output.items():
        file = file[:-4].replace('-', '_').replace(' ', '_')
        clean_IMGT_output[file] = []
        for entry in data:
            newentry = dict()
            for header, value in entry.items():
                #clean headers for model matching
                if header.endswith('-') or header.endswith('+'):  # check for those pesky ++-
                    if sum([True for _ in header[-3:] if _ in ('+', '-')]) == 3:  # make sure it's ONLY those pesky -++
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

                #clean values for model loadin
                if file == "V_REGION_nt_mutation_statistics" or file == "V_REGION_AA_change_statistics": # extra value between brackets is only a thing for this file
                    if re.match(value_with_value_between_brackets, value):
                        value = re.sub(value_between_brackets, '', value)
                    if len(value) == 0:
                        value = '0'

                if value == '':
                    value = None
                elif value == "X":
                    sequenceid = newentry["Sequence_ID"]
                    delete_list.append(sequenceid)
                    value = '0'

                newentry[header] = value
            clean_IMGT_output[file].append(newentry)

    for file, data in clean_IMGT_output.items():
        for entry in data:
            if entry["Sequence_ID"] in delete_list:
                print("removed", entry["Sequence_ID"])
                clean_IMGT_output[file].remove(entry)

    return clean_IMGT_output, delete_list


if __name__ == "__main__":
    static = pathlib.Path.joinpath(pathlib.Path("/mnt/c"), "Users", "Mila", "PycharmProjects", "Bafstu", "PinkStrawberry", "static")
    a = IMGT_output_parser(pathlib.Path.joinpath(static, 'IMGT_output'))
    #print(a["Summary.txt"][0])
    #visualize_IMGT_output_parser(a)
    #check_line_synchronisation(a)
    #headercount_withvalues(a)
    modelmaker(a, show_reasoning=True, sanitize_data=False)
    b = sanitize_output(a)
