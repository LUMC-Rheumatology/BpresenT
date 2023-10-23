import numpy as np
from django.shortcuts import render, redirect
from django.contrib import messages
from django.conf import settings
from django.http import HttpResponse
from rest_framework import generics
from .serializers import *
import re
from pathlib import Path
from PinkStrawberry import models, forms, VQuest_bridge, VQuest_expansion, netmhcpan_bridge, new_html_compiler, \
    MS_bridge, stats, fasta_writer


# Create your views here.
def landing(request):
    return render(request, "landing.html", {})


# "Frontend" basic view
def base(request): # just for extending, user shouldn't go here but it's nice to keep around for testing
    print(request.resolver_match)
    return render(request, "frontend/base.html", {})


def projects(request):
    if request.method == 'POST':
        # submit new project to be created
        modal_show = True
        form = forms.project_form(request.POST)
        if form.is_valid():
            if models.Project.objects.filter(name=request.POST.get('name')).exists():
                form.add_error('name', 'Already exists!')
            else:
                form.save()
                form = forms.project_form
                modal_show = False
        projectdata = models.Project.objects.all().order_by('-created_at')
        return render(request, "frontend/projects.html",
                      {"projects": projectdata, 'form': form, 'modal_show':modal_show, 'search_form':forms.search_form})

    elif request.method == 'GET':
        form = forms.project_form
        projectsdata = models.Project.objects.all().order_by('-created_at')
        return render(request, "frontend/projects.html",
                {"projects":projectsdata, 'form': form, 'modal_show':False, 'search_form':forms.search_form})


def project_detail(request, pk):
    projectdata = models.Project.objects.get(id=pk)
    if request.method == 'POST':
        if request.POST.get('name'):
            if (form := forms.vquest_form(request.POST)).is_valid():
                numstartpattern = re.compile("^[0-9]([0-9]|_)_?")
                kwargs = dict()
                requiredfiles = ['1_Summary.txt', '2_IMGT-gapped-nt-sequences.txt', '3_Nt-sequences.txt',
                                 '4_IMGT-gapped-AA-sequences.txt', '5_AA-sequences.txt', '6_Junction.txt',
                                 '7_V-REGION-mutation-and-AA-change-table.txt', '8_V-REGION-nt-mutation-statistics.txt',
                                 '9_V-REGION-AA-change-statistics.txt', '10_V-REGION-mutation-hotspots.txt',
                                 '11_Parameters.txt']
                for file in request.FILES.getlist('file_field'):
                    if (filename := file.__str__()) in requiredfiles:
                        requiredfiles.pop(requiredfiles.index(filename))
                        filename = re.sub(numstartpattern, '', filename)[:-4].replace('-', '_').replace(' ', '_')
                        kwargs[filename] = file

                if requiredfiles:
                    form.add_error("file_field", [f"Missing files: {', '.join(requiredfiles)}"])
                    return render(request, "frontend/project_detail.html",
                                  {"project": projectdata, 'vquest_form': form, 'modal_show': True})

                kwargs['Project'] = projectdata
                form.cleaned_data.pop('file_field') #clear files as those are in kwargs
                obj = models.VQuest_Run(**kwargs, **form.cleaned_data)
                obj.save()
                try:
                    VQuest_bridge.database_loader(obj)
                    VQuest_expansion.generate_metadata(obj)
                except Exception as err:
                    obj.delete()
                    raise err
            else:
                print(form)

        return render(request, "frontend/project_detail.html",
                      {"project": projectdata, 'vquest_form': forms.vquest_form, 'modal_show': False,
                       "vquestruns":projectdata.VQuest_Run.all()})

    elif request.method == 'GET':
        return render(request, "frontend/project_detail.html",
                {"project":projectdata, 'vquest_form':forms.vquest_form, 'modal_show':False,
                       "vquestruns":projectdata.VQuest_Run.all()})


def sequences(request):
    if request.method == 'POST':
        search_form = forms.search_form(request.POST)
        print(request.POST)
        if search_form.is_valid():
            cleankwargs = {key : value for key, value in search_form.cleaned_data.items() if (value not in ('', None, [], ['']))}
            print("search with arguments:", cleankwargs, sep="\n", end="\n\n")
            query = models.Sequence.objects.all()
            #handle project selection
            if cleankwargs.get('Project', False):
                project = cleankwargs.pop('Project')
                query = project.get_sequences()
                search_form.fields['Patient'].queryset = models.Patient.objects.filter(id__in=query.values_list('Patient', flat = True))
            #check for foreign key associations and trickle down models for the ModelChoiceField
            fkslen = len(fks := ('Patient', 'BCR'))
            for index, fk in enumerate(fks, start=1):
                if cleankwargs.get(fk, False):
                    query = query.filter(**{fk:cleankwargs.pop(fk)})
                    if index < fkslen:
                        # find next fk field in fks,
                        # check what objects are associated with the sequences still in the query for that fk field,
                        # give those as options for the dropdown
                        search_form.fields[nfk].queryset = getattr(models, (nfk := fks[index])).objects.filter(id__in=query.values_list(nfk, flat = True))

            #parse for MutipleChoiceFields
            mcfs = ('Type',)
            for mcf in mcfs:
                if cleankwargs.get(mcf, False):
                    print(temp := cleankwargs.pop(mcf))
                    query = query.filter(**{f"{mcf}__in" : temp})

            #compile LIKE %val% query for remaining kwargs
            querykwargs = {f"{key}__contains" : value for key, value in cleankwargs.items() if (value not in ('', None))}
            query = query.filter(**querykwargs)

            #store the queryset so we can download it next step if we need to

        else:
            #triggers if search form was entered wrong somehow
            print("Faulty search request!", request.POST, search_form, sep="\n", end="\n\n")
            query = models.Sequence.objects.all()

        # check if any of the specialized buttons were pressed
        if 'getfastant' in request.POST:
            fasta = fasta_writer.write(query)
            return HttpResponse(fasta, content_type='text/plain; charset=utf-8')
        if 'getfastaaa' in request.POST:
            fasta = fasta_writer.write(query, translated=True)
            return HttpResponse(fasta, content_type='text/plain; charset=utf-8')

        return render(request, "frontend/sequences.html",
                      {"sequences": query.order_by('id'), 'search_form': search_form})

    if request.method == 'GET':
        search_form = forms.search_form
        query = models.Sequence.objects.all()
        return render(request, "frontend/sequences.html",
                      {"sequences": query.order_by('id'), 'search_form': search_form})


def sequence_detail(request, id):
    if request.method == 'GET':
        sequenceobj = models.Sequence.objects.get(id=id)
        subtype = sequenceobj.get_subtype()
        data = {"Project": sequenceobj.Run.Project,
                "V-Quest": sequenceobj.Run,
                "Sequence ID": sequenceobj.Sequence_ID,
                "Isotype": (f"Ig{sequenceobj.Type}" if sequenceobj.isheavy() else sequenceobj.Type)+
                           (f": {subtype}" if subtype is not None else ''),
                "Germline": sequenceobj.get_germlines().values(),
                "Functionality": sequenceobj.V_DOMAIN_Functionality,
                "Comments": sequenceobj.Comment,
                "Glycosites": sequenceobj.Glycosite.all().count(),
                "Mutations": sequenceobj.IMGT_Nt_mutation.all().count(),
                }
        peptides = [peptide_obj for peptide_obj in sorted([_.Target_peptide for _ in sequenceobj.Peptide_comparison.all()],  key=lambda _: _.percentile_rank)]
        header = ">"+sequenceobj.Sequence_ID
        nt_sequence = sequenceobj.Nt_Sequence
        associated_glycosites = sequenceobj.Glycosite.all()

        # get visualized alignments of all regions separately
        visualizer_dict = dict()
        regions = ('V', 'D', 'J') if sequenceobj.isheavy() else ('V', 'J')
        for region in regions:
            for notation in ('nt', 'aa'):
                visualizer_dict[f'germline_alignment_{region}_{notation}'] = \
                    new_html_compiler.visualizer_shortcut(sequenceobj, germlines=True, region=region, linebr=50,
                                                          translated=True if notation=='aa' else False)

        return render(request, "frontend/sequence_detail.html",
                      {"sequence_obj":sequenceobj, "data":data, "header":header, "nt_sequence":nt_sequence,
                       "associated_glycosites":associated_glycosites, "peptides":peptides,
                       **visualizer_dict,
                       },
                      )


def patient_detail(request, id):
    patientobj = models.Patient.objects.get(Identifier=id)
    # fetch data
    patientseqs = patientobj.Sequence.all()
    mutationcount = [len(i.Nt_mutations.all()) for i in patientseqs]
    glycocount = [len(i.Glycosite.all()) for i in patientseqs]
    data = {'ID':patientobj.Identifier,
            '#Sequences': len(patientseqs),
            '#Mutations mean':f"{np.mean(mutationcount):.2f}",
            '#Mutations median':f"{np.median(mutationcount):.2f}",
            '#Glycosite mean': f"{np.mean(glycocount):.2f}",
            '#Glycosite median': f"{np.median(glycocount):.2f}",
            }

    # start charts
    peptidepiechart = stats.peptidedistribution(patientobj)
    isotypepiechart = stats.isotypedistribution(patientobj)
    h_germlinebarcharts = stats.germlinedistribution(patientobj, "H")
    k_germlinebarcharts = stats.germlinedistribution(patientobj, "K")
    l_germlinebarcharts = stats.germlinedistribution(patientobj, "L")
    mut_glychart = stats.mutation_glycosilationdistribution(patientobj)
    return render(request, 'frontend/patient_detail.html',
                  {'data': data,
                   'peptidepiechart': peptidepiechart,
                   'isotypepiechart': isotypepiechart,
                   'h_germlinebarcharts': h_germlinebarcharts,
                   'k_germlinebarcharts': k_germlinebarcharts,
                   'l_germlinebarcharts': l_germlinebarcharts,
                   'mutglychart': mut_glychart,
                   'object': patientobj})


def germline_detail(request, gene_and_allele):
    if request.method == 'GET':
        labels = f"{gene_and_allele[3]}-REGION"
        try:
            germline = models.Germline.objects.get(Gene_and_allele=gene_and_allele, Labels=labels,
                                                   IMGT_Gene_DB=models.IMGT_Gene_DB.objects.latest())
        except models.Germline.DoesNotExist as error:
            messages.error(request, "No related germline object found")
            return redirect("sequences")

        germline_data = {field: getattr(germline, field) for field in
                         ['Accesion_number', 'Gene_and_allele', 'Species', 'Functionality', 'Labels',
                          'Accesion_start_end', 'Accesion_nucleotide_nb', 'Codon_start', 'Fiveprime_changed',
                          'Threeprime_changed', 'Sequencing_changed', 'AA_nb', 'Length', 'Partial',
                          'Reverse_complementary', 'Nt_Sequence']}
        germline_sequence = germline_data.pop('Nt_Sequence')+'\n'
        germline_header = '>'+'|'.join((val for val in germline_data.values()))+'|\n'
        associated_sequences = getattr(germline, f"Sequence{gene_and_allele[3]}").all()
        return render(request, "frontend/germline_detail.html",
                      {"germline":germline, "germline_data":germline_data,
                       "germline_sequence":germline_sequence, "germline_header":germline_header,
                       "associated_sequences":associated_sequences},)


def visualizer(request, id):
    # check if there's one already loaded this session
    if id == 0:
        if (_ := request.session.get("lastvisualized")):
            return redirect("visualizer/id", id=_)
        else:
            return redirect('sequences')

    sequenceobj = models.Sequence.objects.get(id=id)
    predictions = sequenceobj.Prediction.filter(germline=False).order_by('percentile_rank')
    massspec = sequenceobj.Ms_peptide.filter(sequence_start__isnull=False).order_by('quality_qval')
    request.session['lastvisualized'] = id
    formatted_sequence_nt = new_html_compiler.visualizer_shortcut(sequenceobj, germlines=True, translated=False)
    formatted_sequence_aa = new_html_compiler.visualizer_shortcut(sequenceobj, germlines=True, translated=True)

    if request.method == 'GET':
        form = forms.visualizer_form()
        form.fields['predicted_peptides'].queryset = predictions
        form.fields['mass_spec_peptides'].queryset = massspec
        return render(request, "frontend/visualizer.html",
                      {"form": form, "sequenceobj":sequenceobj,
                       "formatted_sequence_nt":formatted_sequence_nt,
                       "formatted_sequence_aa":formatted_sequence_aa})

    if request.method == 'POST':
        filledform = forms.visualizer_form(request.POST)
        filledform.fields['predicted_peptides'].queryset = predictions
        filledform.fields['mass_spec_peptides'].queryset = massspec

        if filledform.is_valid():
            germbool = filledform.cleaned_data.get('germline')
            predictedqueryset = filledform.cleaned_data.get('predicted_peptides')
            massqueryset = filledform.cleaned_data.get('mass_spec_peptides')
            linebreak = filledform.cleaned_data.get('linebreak')
            formatted_sequence_nt = new_html_compiler.visualizer_shortcut(sequenceobj,
                                        germlines=germbool, translated=False, predictions=predictedqueryset,
                                        linebr=linebreak, massspec=massqueryset)
            formatted_sequence_aa = new_html_compiler.visualizer_shortcut(sequenceobj,
                                        germlines=germbool, translated=True,  predictions=predictedqueryset,
                                        linebr=linebreak, massspec=massqueryset)
        else:
            pass
        form = filledform
        return render(request, "frontend/visualizer.html",
                      {"form": form, "sequenceobj":sequenceobj,
                       "formatted_sequence_nt":formatted_sequence_nt,
                       "formatted_sequence_aa":formatted_sequence_aa})


def netmhcpan(request):
    form = forms.netmhcpan_form
    if request.method == 'GET':
        return render(request, "frontend/netmhcpan.html",
                      {"form": form})

    if request.method == 'POST':
        filledform = form(request.POST)
        if filledform.is_valid():
            mhc_alleles = filledform.data['allele'].replace(' ', '').split(';')
            lengths = list(map(int, filledform.data['length'].replace(' ', '').split(';')))
            runset = filledform.cleaned_data.get('run')
            sequenceset = filledform.cleaned_data.get('sequences')

            # fuse choices from runs into one
            if runset:
                generator = (run.Sequence.all() for run in runset)
                mainset = next(generator)
                for otherset in generator:
                    mainset = mainset | otherset
                if sequenceset:
                    mainset = mainset | sequenceset
            else:
                mainset = sequenceset
            netmhcpan_bridge.fullride(mainset, mhc_alleles, lengths)
            messages.success(request, "Peptides done predicting!")
            return redirect("sequence/id", id=mainset[0].id)

        else:
            pass
        form = filledform
        return render(request, "frontend/netmhcpan.html",
                      {"form":form})

#return .txt of the netmhcpan MHC allele list
def allelelist(request):
    textfile = open(Path.joinpath(settings.BASE_DIR, "PinkStrawberry/static/netmhcpan/services.healthtech.dtu.dk_services_NetMHCpan-4.1_MHC_allele_names.txt")).read()
    return HttpResponse(textfile, content_type='text/plain')


def ms_overview(request):
    form = forms.ms_excel_form
    queryset = models.Ms_peptide.objects.filter(sequence__isnull=True)
    if request.method == 'POST':
        if request.FILES:
            try:
                exceldata = request.FILES['excelsheet']
                excel = MS_bridge.read_msexcel(exceldata)
                MS_bridge.to_db(excel)
            except Exception as error:
                messages.error(request, f"Encountered error: {error}")
        elif (_ := request.POST.get('ms_peptide')):
            models.Ms_peptide.objects.get(id=_).delete()
    return render(request, "frontend/ms_overview.html", {'form':form, 'queryset':queryset})


def ms_form(request, id):
    obj = models.Ms_peptide.objects.get(id=id)
    peptide_form = forms.ms_peptide_form(instance=obj)
    form = forms.ms_addpeptide_form

    if request.method == 'GET':
        if (_ := request.session.get('lastpredicted')):
            print(_)
            form = form(initial={'sequence': models.Sequence.objects.get(id=_[0]),
                                 'test_allele': _[1]}
                        )
        return render(request, "frontend/ms_form.html", {"peptide_form":peptide_form, "form":form, 'object':obj})
    if request.method == 'POST':
        data = forms.ms_addpeptide_form(request.POST)
        if data.is_valid():
            # load last submission into session for easy repeat
            request.session['lastpredicted'] = (data.data['sequence'], data.data['test_allele'])
            _ = request.session.get('lastpredicted')
            form = form(initial={'sequence': models.Sequence.objects.get(id=_[0]),'test_allele': _[1]})
            # handle allele and sequence input into model
            test_allele = data.data['test_allele'].replace(' ', '').split(';')
            obj.possible_allele.set(models.MHC_allele.objects.filter(name__in=test_allele))
            # test for optimal allele
            optimal_allele = netmhcpan_bridge.find_optimal_mhc(obj, test_allele)
            obj.allele = models.MHC_allele.objects.get(name=optimal_allele)
            obj.mhcpansequence = models.Sequence.objects.get(id=data.data['sequence'])
            obj.save()
            messages.success(request, 'Peptide submission successful')
            peptide_form = forms.ms_peptide_form(instance=obj)
        else:
            form = data

        return render(request, "frontend/ms_form.html", {"peptide_form": peptide_form, "form":form, 'object':obj})


def peptides_overview(request):
    form = forms.peptidesearch_form
    if request.method == 'GET':
        queryset = models.Prediction.objects.filter(germline=False).order_by('percentile_rank')
        pairedset = list(zip(map(lambda x: models.Sequence.objects.get(id=x.mhcpansequence), queryset), queryset))
        return render(request, "frontend/peptides_overview.html",
                      {"queryset":pairedset, "peptidesearch_form":form},)

    elif request.method == 'POST':
        data = form(request.POST)
        if data.is_valid():
            sequenceobj = data.cleaned_data.get('sequence')
        queryset = models.Prediction.objects.filter(germline=False).order_by('percentile_rank')
        pairedset = list(zip(map(lambda x: models.Sequence.objects.get(id=x.mhcpansequence), queryset), queryset))
        return render(request, "frontend/peptides_overview.html",
                      {"queryset":pairedset, "peptidesearch_form":form},)


def settings_regular(request):
    return render(request, "frontend/base.html",
                  {})

def settings_advanced(request):
    return render(request, "frontend/base.html",
                  {})






# REST frame
class sequence_list(generics.ListCreateAPIView):
    queryset = models.Sequence.objects.all()
    serializer_class = SequenceSerializer

class vquest_list(ExampleViewSet):
    queryset = models.VQuest_Run.objects.all()
    serializer_class = VquestSerializer





"""
class run_list(generics.ListCreateAPIView):
    queryset = RunIdentifiers.objects.all()
    serializer_class = RunSerializer


class sequence_details(generics.RetrieveUpdateDestroyAPIView):
    lookup_field = 'id'
    serializer_class = SequenceSerializer

    def get_queryset(self):
        return Sequences.objects.filter(Run_Identifier_id=self.kwargs["runid"], id=self.kwargs["id"])


class run(generics.ListCreateAPIView):
    lookup_field = "Run_Identifier_id"
    serializer_class = SequenceSerializer

    def get_queryset(self):
        return Sequences.objects.filter(Run_Identifier_id=self.kwargs["Run_Identifier_id"])


class FileUploadView(APIView):
    parser_classes = (FileUploadParser,)

    def put(self, request, filename, format=None):
        file_obj = request.FILES['file']
        print(file_obj)
        return Response(status=204)

#path discovery assister
@api_view(['GET'])
def api_root(request, format=None):
    return Response({
        'run': reverse('run', request=request, format=format),
        #'run': reverse('run_detail', request=request, format=format),
    })
"""