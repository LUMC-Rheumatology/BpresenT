from rest_framework import serializers
from .models import *

class SequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sequence
        fields = '__all__'


class VquestSerializer(serializers.ModelSerializer):
    class Meta:
        model = VQuest_Run
        fields = '__all__'

# vquest run overview
from django.http import FileResponse
from rest_framework import viewsets, renderers
from rest_framework.decorators import action

class PassthroughRenderer(renderers.BaseRenderer):
    """
        Return data as-is. View should supply a Response.
    """
    media_type = ''
    format = ''
    def render(self, data, accepted_media_type=None, renderer_context=None):
        return data


class ExampleViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = VQuest_Run.objects.all()
    serializer_class = VquestSerializer

    @action(methods=['get'], detail=True, renderer_classes=(PassthroughRenderer,))
    def download(self, *args, **kwargs):
        instance = self.get_object()

        # get an open file handle (I'm just using a file attached to the model for this example):
        file_handle = instance.file.open()

        # send file
        response = FileResponse(file_handle, content_type='whatever')
        response['Content-Length'] = instance.file.size
        response['Content-Disposition'] = 'attachment; filename="%s"' % instance.file.name

        return response

"""
class SequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sequences
        fields = ['id', 'Sequence_Number', 'Sequence_ID', 'V_DOMAIN_Functionality', 'V_GENE_and_allele',
                  'submitted_sequence', 'submitted_sequence_translated', 'IMGT_sequence', 'Run_Identifier']


class HyperSequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sequences
        fields = ['url', 'id', 'Sequence_Number', 'Sequence_ID', 'V_DOMAIN_Functionality', 'V_GENE_and_allele',
                  'submitted_sequence', 'submitted_sequence_translated', 'IMGT_sequence', 'Run_Identifier']


class RunSerializer(serializers.ModelSerializer):
    class Meta:
        model = RunIdentifiers
        fields = ['id', 'run_name', 'description']


class HyperRunSerializer(serializers.HyperlinkedModelSerializer):
    sequences = serializers.HyperlinkedRelatedField(many=True, view_name='run', read_only=True)

    class Meta:
        model = RunIdentifiers
        fields = ['url', 'id', 'run_name', 'description', 'sequences']
        extra_kwargs = {
            'url': {'view_name': 'run_detail', 'lookup_field': '<int:Run_Identifier_id>'},
        }


class PredictionSerializer(serializers.ModelSerializer):

    class Meta:
        model = Predictions
        fields = ['id', 'peptide', 'allele', 'score', 'percentile_rank', 'affinity', 'offset', 'sequence']


class MutationSerializer(serializers.ModelSerializer):

    class Meta:
        model = Mutations
        fields = ['id', 'aa_original', 'aa_pos', 'aa_mutated', 'hydrophobicity', 'volume', 'chemical_characteristics', 'aa_changes', 'region', 'sequence']
"""