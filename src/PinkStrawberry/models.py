import math
import re

from django.db import models
from itertools import chain
from Bio.Seq import Seq
from PinkStrawberry.Aligner import translate_range

# Create your models here.

class Project(models.Model):
    name = models.TextField()
    description = models.TextField(null=True)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.name

    def get_sequences(self):
        generator = (run.Sequence.all() for run in self.VQuest_Run.all())
        try:
            mainset = next(generator)
        except StopIteration:
            mainset = set()
        else:
            for otherset in generator:
                mainset = mainset | otherset
        return mainset

    class Meta:
        get_latest_by = 'id'

# function to create the directory structure for storing VQuest runs, upload_to= automatically targets the project's
# media folder and this func takes it from there. Could probably be done inside the class definition.
def _vquest_path(instance, filename):
    vector = instance.created_at
    return f"vquest/{vector.year}-{vector.month}-{vector.day}/{filename}"


class VQuest_Run(models.Model):
    name = \
        models.TextField()
    description = \
        models.TextField(null=True)
    created_at = \
        models.DateTimeField(auto_now_add=True)
    Summary = \
        models.FileField(upload_to=_vquest_path)
    IMGT_gapped_nt_sequences = \
        models.FileField(upload_to=_vquest_path)
    Nt_sequences = \
        models.FileField(upload_to=_vquest_path)
    IMGT_gapped_AA_sequences = \
        models.FileField(upload_to=_vquest_path)
    AA_sequences = \
        models.FileField(upload_to=_vquest_path)
    Junction = \
        models.FileField(upload_to=_vquest_path)
    V_REGION_mutation_and_AA_change_table = \
        models.FileField(upload_to=_vquest_path)
    V_REGION_nt_mutation_statistics = \
        models.FileField(upload_to=_vquest_path)
    V_REGION_AA_change_statistics = \
        models.FileField(upload_to=_vquest_path)
    V_REGION_mutation_hotspots = \
        models.FileField(upload_to=_vquest_path)
    Parameters = \
        models.FileField(upload_to=_vquest_path)
    Project = models.ForeignKey(Project, on_delete=models.CASCADE,
                                related_name='VQuest_Run')
    Comment = models.TextField(null=True)

    def __str__(self):
        return f"{self.name}"

    #func used to add comments
    def addcomment(self, comment):
        if self.Comment is None:
            self.Comment = ''
        self.Comment = self.Comment + (comment+'\n')
        return

    class Meta:
        get_latest_by = 'id'


# models for the imgt Gene_DB comparison
class IMGT_Gene_DB(models.Model):
    release = models.TextField()
    downloaded_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"IMGT-Gene-DB {self.release}"

    class Meta:
        get_latest_by = 'id'


class Germline(models.Model):
    IMGT_Gene_DB = models.ForeignKey(IMGT_Gene_DB, on_delete=models.CASCADE,
                                            related_name='Germline')
    Accesion_number = models.TextField(null=True)
    Gene_and_allele = models.TextField(null=True)
    Species = models.TextField(null=True)
    Functionality = models.TextField(null=True)
    Labels = models.TextField(null=True)
    Accesion_start_end = models.TextField(null=True)
    Accesion_nucleotide_nb = models.TextField(null=True)
    Codon_start = models.TextField(null=True)
    Fiveprime_changed = models.TextField(null=True)
    Threeprime_changed = models.TextField(null=True)
    Sequencing_changed = models.TextField(null=True)
    AA_nb = models.TextField(null=True)
    Length = models.TextField(null=True)
    Partial = models.TextField(null=True)
    Reverse_complementary = models.TextField(null=True)
    AA_Sequence = models.TextField(null=True)
    Nt_Sequence = models.TextField(null=True)

    def __str__(self):
        return self.Gene_and_allele

    def translate(self, frame:int=1):
        if (frame := abs(frame)) != frame:
            return Seq(self.Nt_Sequence[::-1][frame - 1:]).translate()
        else:
            return Seq(self.Nt_Sequence[frame - 1:]).translate()


#coupeling stuff
class BCR(models.Model):
    Identifier = models.TextField()
    Mass = models.FloatField(null=True) # not implemented due to lack of data
    Run = models.ForeignKey(VQuest_Run, on_delete=models.CASCADE,
                            related_name='BCR')

    def __str__(self):
        return f"{self.Identifier}"


class Patient(models.Model):
    Identifier = models.TextField(unique=True)
    Comment = models.TextField(null=True)

    def __str__(self):
        return f"{self.Identifier}"

    class Meta:
        get_latest_by = 'id'


class Subtype(models.Model):
    Name = models.TextField(unique=True)

    def __str__(self):
        return self.Name


class Sequence(models.Model):
    Run = models.ForeignKey(VQuest_Run, on_delete=models.CASCADE,
                            related_name='Sequence')
    BCR = models.ForeignKey(BCR, on_delete=models.CASCADE,
                            related_name='Sequence', null=True)
    Patient = models.ForeignKey(Patient, on_delete=models.CASCADE,
                            related_name='Sequence', null=True)
    Sequence_number = models.IntegerField()
    Sequence_ID = models.CharField(max_length=25)
    V_DOMAIN_Functionality = models.CharField(max_length=20)
    V_GENE_and_allele = models.TextField()
    Type = models.CharField(max_length=1, null=True)
    Subtype = models.ManyToManyField(Subtype)  # MtM because subtype can be imprecise
    AA_Sequence = models.TextField(null=True)
    AA_Leading = models.TextField(null=True)
    AA_Constant = models.TextField(null=True)
    Nt_Sequence = models.TextField(null=True)
    Nt_Leading = models.TextField(null=True)
    Nt_Constant = models.TextField(null=True)
    Comment = models.TextField(null=True)
    V_Germline = models.ForeignKey(Germline, on_delete=models.SET_NULL,
                                   related_name='SequenceV', null=True)
    D_Germline = models.ForeignKey(Germline, on_delete=models.SET_NULL,
                                   related_name='SequenceD', null=True)
    J_Germline = models.ForeignKey(Germline, on_delete=models.SET_NULL,
                                   related_name='SequenceJ', null=True)

    #X is for not defined, should be avoided with determine_heavychain func from IMGT_analysis
    def chaintype(self):
        translation = {'X': 'heavy', 'G': 'heavy', 'M': 'heavy', 'A': 'heavy', 'D': 'heavy', 'E': 'heavy', 'K': 'light', 'L': 'light'}
        return translation[self.Type]

    def isheavy(self):
        return True if self.chaintype() == 'heavy' else False

    def islight(self):
        return True if self.chaintype() == 'light' else False

    #func used to isolate heavy/lambda/kappa from v gene i.e. (Homsap IG -> H <- V3-11*05 F)
    def imgttype(self):
        return self.V_GENE_and_allele[9]

    # get_vdjrange() is preferred
    # compensates for index difference in python, adjusts to the beginning of FR1 (V-Region start)
    def get_regions(self, aa_mode=False):
        region_locations = dict()
        table = self.Nt_sequences.__dict__
        startval = table["FR1_IMGT_start"] # complies with python index starting at 0 if subtracted
        region_names = ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4')
        for region in region_names:
            region_locations[region] = {'start': table.get(f"{region}_IMGT_start")-startval,
                                        'stop': table.get(f"{region}_IMGT_end")-startval+1} # +1 for python splicing ease
        if aa_mode:
            region_locations = {region: {flag: int(loc/3)
                                for flag, loc in locs.items()}
                                for region, locs in region_locations.items()}
        return region_locations


    def get_vdjrange(self, translated=False):
        ranges = [('V', self.Nt_sequences.V_REGION_start - 1, self.Nt_sequences.V_REGION_end),
                  ('J', self.Nt_sequences.J_REGION_start - 1, self.Nt_sequences.J_REGION_end)]
        if self.isheavy():
            ranges.insert(1,
                          ('D', self.Nt_sequences.D_REGION_start - 1, self.Nt_sequences.D_REGION_end))

        if not translated:
            return ranges

        ntrange = ranges.copy()
        ranges = []
        for key, start, stop in ntrange:
            ranges.append((key, *translate_range(start, stop, self.Nt_sequences.V_REGION_reading_frame)))
        return ranges


    def get_germlines(self):
        return {"V":self.V_Germline, "D":self.D_Germline, "J":self.J_Germline}


    def get_subtype(self):
        try:
            val = ', '.join(_.Name for _ in self.Subtype.all())
        except AttributeError:
            val = None
        return val if len(val) != 0 else None

    #func used to add comments
    def addcomment(self, comment):
        if self.Comment is None:
            self.Comment = ''
        self.Comment = self.Comment + (comment+'\n')
        return

    class Meta:
        get_latest_by = 'id'


class alignment(models.Model):
    name = models.TextField()
    position = models.IntegerField()
    sequence = models.TextField()
    linked_sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE, related_name="alignment")

    def range(self):
        return (self.position, self.position + len(self.sequence))


class IMGT_AA_change(models.Model):
    _similarities = (
        (0, "Very dissimilar"),
        (1, "Dissimilar"),
        (2, "Similar"),
        (3, "Very similar"),
        (4, "Identical")
    )
    aa_original = models.CharField(max_length=1)
    aa_pos = models.IntegerField()
    aa_changed = models.CharField(max_length=1)
    hydrophobicity = models.BooleanField()
    volume = models.BooleanField()
    chemical_characteristics = models.BooleanField()
    similarity = models.IntegerField(choices=_similarities)
    region = models.CharField(max_length=10)
    sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE, related_name="IMGT_AA_change")

    class Meta:
        get_latest_by = 'id'


class IMGT_Nt_mutation(models.Model):
    nt_original = models.CharField(max_length=1)
    nt_pos = models.IntegerField()
    nt_mutated = models.CharField(max_length=1)
    sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE, related_name="IMGT_Nt_mutation")
    aa_change = models.ForeignKey(IMGT_AA_change, on_delete=models.DO_NOTHING, related_name="IMGT_Nt_mutation", null=True)


class Glycosite(models.Model):
    Sequence_Identifier = models.ForeignKey(Sequence, on_delete=models.CASCADE,
                                            related_name='Glycosite')
    sequence = models.CharField(max_length=4)
    position = models.IntegerField()
    germline = models.BooleanField()
    insequence = models.BooleanField(default=1)

    class Meta:
        get_latest_by = 'id'


class AA_change(models.Model):
    original = models.CharField(max_length=1)
    position = models.IntegerField()
    changed = models.CharField(max_length=1)
    sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE, related_name="AA_changes")

    def get_region(self):
        for region, start, stop in self.sequence.get_vdjrange(translated=True):
            if start <= self.position < stop:
                return region


class Nt_mutation(models.Model):
    original = models.CharField(max_length=1)
    position = models.IntegerField()
    changed = models.CharField(max_length=1)
    aa_change = models.ForeignKey(AA_change, on_delete=models.DO_NOTHING, related_name="mutations", null=True)
    sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE, related_name="Nt_mutations")

    def get_region(self):
        for region, start, stop in self.sequence.get_vdjrange():
            if start <= self.position < stop:
                return region



class MHC_allele(models.Model): # MHC alleles for easy coupling
    name = models.TextField(primary_key=True)

    def __str__(self):
        return self.name

    def __repr__(self):
        return f"MHC|{self.name}"



class Ms_peptide(models.Model): # Peptides found through Mass-spec (user input)
    _confidences = (
        (0, "Low"),
        (1, "Medium"),
        (2, "High"),
    )
    peptide = models.TextField()
    start = models.IntegerField()
    stop = models.IntegerField()
    confidence = models.IntegerField(choices=_confidences)
    allele = models.ForeignKey(MHC_allele, on_delete=models.DO_NOTHING, related_name="Ms_peptides", null=True)
    possible_allele = models.ManyToManyField(MHC_allele, related_name="Ms_peptides_possible")
    quality_qval = models.FloatField()
    Comment = models.TextField(null=True) # modifications go here
    sequence_start = models.IntegerField(null=True)
    sequence = models.ForeignKey(Sequence, on_delete=models.DO_NOTHING,
                                related_name='Ms_peptide', null=True)
    sheet_sequence = models.TextField()

    def save(self, *args, **kwargs):
        if self.sequence:
            for predobj in Prediction.objects.filter(peptide__contains=self.peptide, sequence=self.sequence.id):
                Predicted_Found.objects.get_or_create(predicted=predobj, found=self)
            research = re.search(self.peptide, self.sequence.AA_Sequence)
            if research:
                self.sequence_start = research.start()

        super(Ms_peptide, self).save(*args, **kwargs)


class Prediction(models.Model):  # NetMHCpan4.1 predictions
    _binders = (
        (2, 'Strong binder'),
        (1, 'Weak binder'),
        (0, 'No binder'))
    peptide = models.TextField()
    allele = models.ForeignKey(MHC_allele, on_delete=models.DO_NOTHING, related_name="Pred_peptides", null=True)
    score = models.FloatField()
    percentile_rank = models.FloatField()
    affinity = models.FloatField()
    offset = models.IntegerField()
    mhcpansequence = models.TextField()
    germline = models.BooleanField(null=True)
    binder = models.IntegerField(choices=_binders, null=True)
    sequence = models.ForeignKey(Sequence, on_delete=models.DO_NOTHING,
                                related_name='Prediction', null=True)

    def save(self, *args, **kwargs):
        if self.percentile_rank < 0.50:
            self.binder = 2
        elif self.percentile_rank < 2:
            self.binder = 1
        else:
            self.binder = 0
        if self.mhcpansequence[:2] == 'GL':
            self.germline = True
            self.sequence = Sequence.objects.get(id=self.mhcpansequence[2:])
        else:
            self.germline = False
            self.sequence = Sequence.objects.get(id=self.mhcpansequence)
        super(Prediction, self).save(*args, **kwargs)


class Peptide_comparison(models.Model):
    peptide = models.TextField(null=True)
    allele = models.ForeignKey(MHC_allele, on_delete=models.DO_NOTHING, related_name="Pred_peptides_comparison", null=True)
    score_diff = models.FloatField(null=True)
    percentile_rank_diff = models.FloatField(null=True)
    affinity_diff = models.FloatField(null=True)
    Target_peptide = models.ForeignKey(Prediction, on_delete=models.CASCADE,
                                related_name='comparison', null=True)
    Germline_peptide = models.ForeignKey(Prediction, on_delete=models.CASCADE,
                                related_name='GLcomparison', null=True)
    sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE,
                                related_name='Peptide_comparison', null=True)

    def update_diff(self):
        if not (self.Target_peptide is None or self.Germline_peptide is None):
            self.peptide = self.Target_peptide.peptide
            self.allele = self.Target_peptide.allele
            self.score_diff = self.Target_peptide.score - self.Germline_peptide.score
            self.percentile_rank_diff = self.Target_peptide.percentile_rank - self.Germline_peptide.percentile_rank
            self.affinity_diff = self.Target_peptide.affinity - self.Germline_peptide.affinity



class Predicted_Found(models.Model):
    _qualities = (
        (0, "Low"),
        (1, "Medium"),
        (2, "High"),
    )
    predicted = models.ForeignKey(Prediction, on_delete=models.CASCADE,
                                   related_name='pairing')
    found = models.ForeignKey(Ms_peptide, on_delete=models.CASCADE,
                                   related_name='pairing')
    quality = models.IntegerField(choices=_qualities, null=True)
    sequence = models.ForeignKey(Sequence, on_delete=models.CASCADE,
                                related_name='predicted_found', null=True)

    def save(self, *args, **kwargs):
        if self.predicted.sequence == self.found.sequence:
            self.sequence = self.found.sequence
        super(Predicted_Found, self).save(*args, **kwargs)



# FILE MODELS https://www.imgt.org/IMGT_vquest/user_guide
# model for file V-REGION-mutation-hotspots.txt
class V_REGION_mutation_hotspots(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE,
                                            related_name='V_REGION_mutation_hotspots')
    ata = models.TextField()  # (a/t)a
    tat = models.TextField()  # t(a/t)
    aggctat = models.TextField()  # (a/g)g(c/t)(a/t)
    atagcct = models.TextField()  # (a/t)(a/g)c(c/t)

    class Meta:
        get_latest_by = 'id'


# model for file Summary.txt
class Summary(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE, related_name='Summary')
    V_REGION_score = models.IntegerField()  # V-REGION score
    V_REGION_identity_percent = models.FloatField()  # V-REGION identity %
    V_REGION_identity_nt = models.TextField()  # V-REGION identity nt
    V_REGION_identity_percent_with_insdel_events = models.FloatField(null=True)  # V-REGION identity % (with ins/del events)
    V_REGION_identity_nt_with_insdel_events = models.TextField(null=True)  # V-REGION identity nt (with ins/del events)
    J_GENE_and_allele = models.TextField()  # J-GENE and allele
    J_REGION_score = models.IntegerField()  # J-REGION score
    J_REGION_identity_percent = models.FloatField()  # J-REGION identity %
    J_REGION_identity_nt = models.TextField()  # J-REGION identity nt
    D_GENE_and_allele = models.TextField(null=True)  # D-GENE and allele
    D_REGION_reading_frame = models.IntegerField(null=True)  # D-REGION reading frame
    CDR1_IMGT_length = models.IntegerField()  # CDR1-IMGT length
    CDR2_IMGT_length = models.IntegerField()  # CDR2-IMGT length
    CDR3_IMGT_length = models.IntegerField()  # CDR3-IMGT length
    CDR_IMGT_lengths = models.TextField()  # CDR-IMGT lengths
    FR_IMGT_lengths = models.TextField()  # FR-IMGT lengths
    AA_JUNCTION = models.TextField()  # AA JUNCTION
    JUNCTION_frame = models.TextField()  # JUNCTION frame
    Orientation = models.CharField(max_length=1)
    V_DOMAIN_Functionality_comment = models.TextField(null=True)  # V-DOMAIN Functionality comment
    V_REGION_potential_insdel = models.TextField(null=True)  # V-REGION potential ins/del
    J_GENE_and_allele_comment = models.TextField(null=True)  # J-GENE and allele comment
    V_REGION_insertions = models.TextField(null=True)  # V-REGION insertions
    V_REGION_deletions = models.TextField(null=True)  # V-REGION deletions
    Sequence = models.TextField()
    fiveprime_trimmed_n_nb = models.IntegerField()  # 5prime trimmed-n nb
    threeprime_trimmed_n_nb = models.IntegerField()  # 3prime trimmed-n nb
    Analysed_sequence_length = models.IntegerField()  # Analysed sequence length
    Sequence_analysis_category = models.TextField()  # Sequence analysis category

    class Meta:
        get_latest_by = 'id'


# model for file IMGT-gapped-nt-sequences.txt
class IMGT_gapped_nt_sequences(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE, related_name='IMGT_gapped_nt_sequences')
    J_GENE_and_allele = models.TextField()  # J-GENE and allele
    D_GENE_and_allele = models.TextField(null=True)  # D-GENE and allele
    V_D_J_REGION = models.TextField(null=True)  # V-D-J-REGION
    V_J_REGION = models.TextField(null=True)  # V-J-REGION
    V_REGION = models.TextField()  # V-REGION
    FR1_IMGT = models.TextField()  # FR1-IMGT
    CDR1_IMGT = models.TextField()  # CDR1-IMGT
    FR2_IMGT = models.TextField()  # FR2-IMGT
    CDR2_IMGT = models.TextField()  # CDR2-IMGT
    FR3_IMGT = models.TextField()  # FR3-IMGT
    CDR3_IMGT = models.TextField()  # CDR3-IMGT
    JUNCTION = models.TextField()
    J_REGION = models.TextField()  # J-REGION
    FR4_IMGT = models.TextField()  # FR4-IMGT

    class Meta:
        get_latest_by = 'id'


# model for file Nt-sequences.txt
class Nt_sequences(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE, related_name='Nt_sequences')
    J_GENE_and_allele = models.TextField()  # J-GENE and allele
    D_GENE_and_allele = models.TextField(null=True)  # D-GENE and allele
    V_D_J_REGION = models.TextField(null=True)  # V-D-J-REGION
    V_J_REGION = models.TextField(null=True)  # V-J-REGION
    V_REGION = models.TextField()  # V-REGION
    FR1_IMGT = models.TextField(null=True)  # FR1-IMGT
    CDR1_IMGT = models.TextField(null=True)  # CDR1-IMGT
    FR2_IMGT = models.TextField(null=True)  # FR2-IMGT
    CDR2_IMGT = models.TextField(null=True)  # CDR2-IMGT
    FR3_IMGT = models.TextField(null=True)  # FR3-IMGT
    CDR3_IMGT = models.TextField(null=True)  # CDR3-IMGT
    JUNCTION = models.TextField(null=True)
    threeprime_V_REGION = models.TextField(null=True)  # 3'V-REGION
    N_D_J_REGION = models.TextField(null=True)  # (N-D)-J-REGION
    N_D_REGION = models.TextField(null=True)  # (N-D)-REGION
    Pthreeprime_V = models.TextField(null=True)  # P3'V
    N_REGION = models.TextField(null=True)  # N-REGION
    N1_REGION = models.TextField(null=True)  # N1-REGION
    Pfiveprime_D = models.TextField(null=True)  # P5'D
    D_REGION = models.TextField(null=True)  # D-REGION
    Pthreeprime_D = models.TextField(null=True)  # P3'D
    Pfiveprime_D1 = models.TextField(null=True)  # P5'D1
    D1_REGION = models.TextField(null=True)  # D1-REGION
    Pthreeprime_D1 = models.TextField(null=True)  # P3'D1
    N2_REGION = models.TextField(null=True)  # N2-REGION
    Pfiveprime_D2 = models.TextField(null=True)  # P5'D2
    D2_REGION = models.TextField(null=True)  # D2-REGION
    Pthreeprime_D2 = models.TextField(null=True)  # P3'D2
    N3_REGION = models.TextField(null=True)  # N3-REGION
    Pfiveprime_D3 = models.TextField(null=True)  # P5'D3
    D3_REGION = models.TextField(null=True)  # D3-REGION
    Pthreeprime_D3 = models.TextField(null=True)  # P3'D3
    N4_REGION = models.TextField(null=True)  # N4-REGION
    Pfiveprime_J = models.TextField(null=True)  # P5'J
    fiveprime_J_REGION = models.TextField(null=True)  # 5'J-REGION
    D_J_REGION = models.TextField(null=True)  # D-J-REGION
    J_REGION = models.TextField()  # J-REGION
    FR4_IMGT = models.TextField(null=True)  # FR4-IMGT
    V_D_J_REGION_start = models.IntegerField(null=True)  # V-D-J-REGION start
    V_D_J_REGION_end = models.IntegerField(null=True)  # V-D-J-REGION end
    V_J_REGION_start = models.IntegerField(null=True)  # V-J-REGION start
    V_J_REGION_end = models.IntegerField(null=True)  # V-J-REGION end
    V_REGION_start = models.IntegerField()  # V-REGION start
    V_REGION_end = models.IntegerField()  # V-REGION end
    FR1_IMGT_start = models.IntegerField(null=True)  # FR1-IMGT start
    FR1_IMGT_end = models.IntegerField(null=True)  # FR1-IMGT end
    CDR1_IMGT_start = models.IntegerField(null=True)  # CDR1-IMGT start
    CDR1_IMGT_end = models.IntegerField(null=True)  # CDR1-IMGT end
    FR2_IMGT_start = models.IntegerField(null=True)  # FR2-IMGT start
    FR2_IMGT_end = models.IntegerField(null=True)  # FR2-IMGT end
    CDR2_IMGT_start = models.IntegerField(null=True)  # CDR2-IMGT start
    CDR2_IMGT_end = models.IntegerField(null=True)  # CDR2-IMGT end
    FR3_IMGT_start = models.IntegerField(null=True)  # FR3-IMGT start
    FR3_IMGT_end = models.IntegerField(null=True)  # FR3-IMGT end
    CDR3_IMGT_start = models.IntegerField(null=True)  # CDR3-IMGT start
    CDR3_IMGT_end = models.IntegerField(null=True)  # CDR3-IMGT end
    JUNCTION_start = models.IntegerField(null=True)  # JUNCTION start
    JUNCTION_end = models.IntegerField(null=True)  # JUNCTION end
    threeprime_V_REGION_start = models.IntegerField(null=True)  # 3'V-REGION start
    threeprime_V_REGION_end = models.IntegerField(null=True)  # 3'V-REGION end
    N_D_J_REGION_start = models.IntegerField(null=True)  # (N-D)-J-REGION start
    N_D_J_REGION_end = models.IntegerField(null=True)  # (N-D)-J-REGION end
    N_D_REGION_start = models.IntegerField(null=True)  # (N-D)-REGION start
    N_D_REGION_end = models.IntegerField(null=True)  # (N-D)-REGION end
    Pthreeprime_V_start = models.IntegerField(null=True)  # P3'V start
    Pthreeprime_V_end = models.IntegerField(null=True)  # P3'V end
    N_REGION_start = models.IntegerField(null=True)  # N-REGION start
    N_REGION_end = models.IntegerField(null=True)  # N-REGION end
    N1_REGION_start = models.IntegerField(null=True)  # N1-REGION start
    N1_REGION_end = models.IntegerField(null=True)  # N1-REGION end
    Pfiveprime_D_start = models.IntegerField(null=True)  # P5'D start
    Pfiveprime_D_end = models.IntegerField(null=True)  # P5'D end
    D_REGION_start = models.IntegerField(null=True)  # D-REGION start
    D_REGION_end = models.IntegerField(null=True)  # D-REGION end
    Pthreeprime_D_start = models.IntegerField(null=True)  # P3'D start
    Pthreeprime_D_end = models.IntegerField(null=True)  # P3'D end
    Pfiveprime_D1_start = models.IntegerField(null=True)  # P5'D1 start
    Pfiveprime_D1_end = models.IntegerField(null=True)  # P5'D1 end
    D1_REGION_start = models.IntegerField(null=True)  # D1-REGION start
    D1_REGION_end = models.IntegerField(null=True)  # D1-REGION end
    Pthreeprime_D1_start = models.IntegerField(null=True)  # P3'D1 start
    Pthreeprime_D1_end = models.IntegerField(null=True)  # P3'D1 end
    N2_REGION_start = models.IntegerField(null=True)  # N2-REGION start
    N2_REGION_end = models.IntegerField(null=True)  # N2-REGION end
    Pfiveprime_D2_start = models.IntegerField(null=True)  # P5'D2 start
    Pfiveprime_D2_end = models.IntegerField(null=True)  # P5'D2 end
    D2_REGION_start = models.IntegerField(null=True)  # D2-REGION start
    D2_REGION_end = models.IntegerField(null=True)  # D2-REGION end
    Pthreeprime_D2_start = models.IntegerField(null=True)  # P3'D2 start
    Pthreeprime_D2_end = models.IntegerField(null=True)  # P3'D2 end
    N3_REGION_start = models.IntegerField(null=True)  # N3-REGION start
    N3_REGION_end = models.IntegerField(null=True)  # N3-REGION end
    Pfiveprime_D3_start = models.IntegerField(null=True)  # P5'D3 start
    Pfiveprime_D3_end = models.IntegerField(null=True)  # P5'D3 end
    D3_REGION_start = models.IntegerField(null=True)  # D3-REGION start
    D3_REGION_end = models.IntegerField(null=True)  # D3-REGION end
    Pthreeprime_D3_start = models.IntegerField(null=True)  # P3'D3 start
    Pthreeprime_D3_end = models.IntegerField(null=True)  # P3'D3 end
    N4_REGION_start = models.IntegerField(null=True)  # N4-REGION start
    N4_REGION_end = models.IntegerField(null=True)  # N4-REGION end
    Pfiveprime_J_start = models.IntegerField(null=True)  # P5'J start
    Pfiveprime_J_end = models.IntegerField(null=True)  # P5'J end
    fiveprime_J_REGION_start = models.IntegerField(null=True)  # 5'J-REGION start
    fiveprime_J_REGION_end = models.IntegerField(null=True)  # 5'J-REGION end
    D_J_REGION_start = models.IntegerField(null=True)  # D-J-REGION start
    D_J_REGION_end = models.IntegerField(null=True)  # D-J-REGION end
    J_REGION_start = models.IntegerField()  # J-REGION start
    J_REGION_end = models.IntegerField()  # J-REGION end
    FR4_IMGT_start = models.IntegerField(null=True)  # FR4-IMGT start
    FR4_IMGT_end = models.IntegerField(null=True)  # FR4-IMGT end
    V_REGION_reading_frame = models.IntegerField()  # V-REGION reading frame
    V_REGION_partial_fiveprime_missing_nt_nb = models.IntegerField(null=True)  # V-REGION partial 5prime missing nt nb
    V_REGION_uncertain_nt_nb = models.IntegerField(null=True)  # V-REGION uncertain nt nb
    J_REGION_partial_threeprime_missing_nt_nb = models.IntegerField(null=True)  # J-REGION partial 3prime missing nt nb

    class Meta:
        get_latest_by = 'id'


# model for file IMGT-gapped-AA-sequences.txt
class IMGT_gapped_AA_sequences(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE, related_name='IMGT_gapped_AA_sequences')
    J_GENE_and_allele = models.TextField()  # J-GENE and allele
    D_GENE_and_allele = models.TextField(null=True)  # D-GENE and allele
    V_D_J_REGION = models.TextField(null=True)  # V-D-J-REGION
    V_J_REGION = models.TextField(null=True)  # V-J-REGION
    V_REGION = models.TextField()  # V-REGION
    FR1_IMGT = models.TextField()  # FR1-IMGT
    CDR1_IMGT = models.TextField()  # CDR1-IMGT
    FR2_IMGT = models.TextField()  # FR2-IMGT
    CDR2_IMGT = models.TextField()  # CDR2-IMGT
    FR3_IMGT = models.TextField()  # FR3-IMGT
    CDR3_IMGT = models.TextField()  # CDR3-IMGT
    JUNCTION = models.TextField()
    J_REGION = models.TextField()  # J-REGION
    FR4_IMGT = models.TextField()  # FR4-IMGT

    class Meta:
        get_latest_by = 'id'


# model for file AA-sequences.txt
class AA_sequences(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE, related_name='AA_sequences')
    J_GENE_and_allele = models.TextField()  # J-GENE and allele
    D_GENE_and_allele = models.TextField(null=True)  # D-GENE and allele
    V_D_J_REGION = models.TextField(null=True)  # V-D-J-REGION
    V_J_REGION = models.TextField(null=True)  # V-J-REGION
    V_REGION = models.TextField()  # V-REGION
    FR1_IMGT = models.TextField(null=True)  # FR1-IMGT
    CDR1_IMGT = models.TextField(null=True)  # CDR1-IMGT
    FR2_IMGT = models.TextField(null=True)  # FR2-IMGT
    CDR2_IMGT = models.TextField(null=True)  # CDR2-IMGT
    FR3_IMGT = models.TextField(null=True)  # FR3-IMGT
    CDR3_IMGT = models.TextField(null=True)  # CDR3-IMGT
    JUNCTION = models.TextField()
    J_REGION = models.TextField()  # J-REGION
    FR4_IMGT = models.TextField(null=True)  # FR4-IMGT

    class Meta:
        get_latest_by = 'id'


# model for file Junction.txt
class Junction(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE, related_name='Junction')
    J_GENE_and_allele = models.TextField()  # J-GENE and allele
    D_GENE_and_allele = models.TextField(null=True)  # D-GENE and allele
    JUNCTION_frame = models.TextField(null=True)  # JUNCTION frame
    JUNCTION = models.TextField(null=True)
    JUNCTION_with_frameshift = models.TextField(null=True)  # JUNCTION (with frameshift)
    threeprime_V_REGION = models.TextField(null=True)  # 3'V-REGION
    Pthreeprime_V = models.TextField(null=True)  # P3'V
    N_REGION = models.TextField(null=True)  # N-REGION
    N1_REGION = models.TextField(null=True)  # N1-REGION
    Pfiveprime_D = models.TextField(null=True)  # P5'D
    D_REGION = models.TextField(null=True)  # D-REGION
    Pthreeprime_D = models.TextField(null=True)  # P3'D
    Pfiveprime_D1 = models.TextField(null=True)  # P5'D1
    D1_REGION = models.TextField(null=True)  # D1-REGION
    Pthreeprime_D1 = models.TextField(null=True)  # P3'D1
    N2_REGION = models.TextField(null=True)  # N2-REGION
    Pfiveprime_D2 = models.TextField(null=True)  # P5'D2
    D2_REGION = models.TextField(null=True)  # D2-REGION
    Pthreeprime_D2 = models.TextField(null=True)  # P3'D2
    N3_REGION = models.TextField(null=True)  # N3-REGION
    Pfiveprime_D3 = models.TextField(null=True)  # P5'D3
    D3_REGION = models.TextField(null=True)  # D3-REGION
    Pthreeprime_D3 = models.TextField(null=True)  # P3'D3
    N4_REGION = models.TextField(null=True)  # N4-REGION
    Pfiveprime_J = models.TextField(null=True)  # P5'J
    fiveprime_J_REGION = models.TextField(null=True)  # 5'J-REGION
    JUNCTION_nt_nb = models.IntegerField(null=True)  # JUNCTION-nt nb
    threeprime_V_REGION_nt_nb = models.IntegerField(null=True)  # 3'V-REGION-nt nb
    Pthreeprime_V_nt_nb = models.IntegerField(null=True)  # P3'V-nt nb
    N_REGION_nt_nb = models.IntegerField(null=True)  # N-REGION-nt nb
    N1_REGION_nt_nb = models.IntegerField(null=True)  # N1-REGION-nt nb
    Pfiveprime_D_nt_nb = models.IntegerField(null=True)  # P5'D-nt nb
    D_REGION_nt_nb = models.IntegerField(null=True)  # D-REGION-nt nb
    Pthreeprime_D_nt_nb = models.IntegerField(null=True)  # P3'D-nt nb
    Pfiveprime_D1_nt_nb = models.IntegerField(null=True)  # P5'D1-nt nb
    D1_REGION_nt_nb = models.IntegerField(null=True)  # D1-REGION-nt nb
    Pthreeprime_D1_nt_nb = models.IntegerField(null=True)  # P3'D1-nt nb
    N2_REGION_nt_nb = models.IntegerField(null=True)  # N2-REGION-nt nb
    Pfiveprime_D2_nt_nb = models.IntegerField(null=True)  # P5'D2-nt nb
    D2_REGION_nt_nb = models.IntegerField(null=True)  # D2-REGION-nt nb
    Pthreeprime_D2_nt_nb = models.IntegerField(null=True)  # P3'D2-nt nb
    N3_REGION_nt_nb = models.IntegerField(null=True)  # N3-REGION-nt nb
    Pfiveprime_D3_nt_nb = models.IntegerField(null=True)  # P5'D3-nt nb
    D3_REGION_nt_nb = models.IntegerField(null=True)  # D3-REGION-nt nb
    Pthreeprime_D3_nt_nb = models.IntegerField(null=True)  # P3'D3-nt nb
    N4_REGION_nt_nb = models.IntegerField(null=True)  # N4-REGION-nt nb
    Pfiveprime_J_nt_nb = models.IntegerField(null=True)  # P5'J-nt nb
    fiveprime_J_REGION_nt_nb = models.IntegerField(null=True)  # 5'J-REGION-nt nb
    threeprime_V_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 3'V-REGION trimmed-nt nb
    fiveprime_D_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 5'D-REGION trimmed-nt nb
    threeprime_D_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 3'D-REGION trimmed-nt nb
    fiveprime_D1_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 5'D1-REGION trimmed-nt nb
    threeprime_D1_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 3'D1-REGION trimmed-nt nb
    fiveprime_D2_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 5'D2-REGION trimmed-nt nb
    threeprime_D2_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 3'D2-REGION trimmed-nt nb
    fiveprime_D3_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 5'D3-REGION trimmed-nt nb
    threeprime_D3_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 3'D3-REGION trimmed-nt nb
    fiveprime_J_REGION_trimmed_nt_nb = models.IntegerField(null=True)  # 5'J-REGION trimmed-nt nb
    threeprime_V_REGION_mut_nt_nb = models.IntegerField(null=True)  # 3'V-REGION mut-nt nb
    D_REGION_mut_nt_nb = models.IntegerField(null=True)  # D-REGION mut-nt nb
    D1_REGION_mut_nt_nb = models.IntegerField(null=True)  # D1-REGION mut-nt nb
    D2_REGION_mut_nt_nb = models.IntegerField(null=True)  # D2-REGION mut-nt nb
    D3_REGION_mut_nt_nb = models.IntegerField(null=True)  # D3-REGION mut-nt nb
    fiveprime_J_REGION_mut_nt_nb = models.IntegerField(null=True)  # 5'J-REGION mut-nt nb
    D_REGION_reading_frame = models.IntegerField(null=True)  # D-REGION reading frame
    Ngc = models.TextField(null=True)
    CDR3_IMGT_length = models.IntegerField(null=True)  # CDR3-IMGT length
    Molecular_mass = models.FloatField(null=True)  # Molecular mass
    pI = models.FloatField(null=True)
    threeprime_V_REGION_accepted_mut_nb = models.IntegerField(null=True)  # 3'V-REGION accepted mut nb
    D_REGION_accepted_mut_nb = models.IntegerField(null=True)  # D-REGION accepted mut nb
    fiveprime_J_REGION_accepted_mut_nb = models.IntegerField(null=True)  # 5'J-REGION accepted mut nb
    Accepted_D_GENE_nb = models.IntegerField(null=True)  # Accepted D-GENE nb
    CDR3_IMGT = models.TextField(null=True)  # CDR3-IMGT
    CDR3_IMGT_nt_nb = models.IntegerField(null=True)  # CDR3-IMGT-nt nb
    CDR3_IMGT_with_frameshift = models.TextField(null=True)  # CDR3-IMGT (with frameshift)
    CDR3_IMGT_AA = models.TextField(null=True)  # CDR3-IMGT (AA)
    CDR3_IMGT_AA_with_frameshift = models.TextField(null=True)  # CDR3-IMGT (AA) (with frameshift)
    JUNCTION_AA = models.TextField(null=True)  # JUNCTION (AA)
    JUNCTION_AA_with_frameshift = models.TextField(null=True)  # JUNCTION (AA) (with frameshift)
    JUNCTION_decryption = models.TextField(null=True)  # JUNCTION decryption

    class Meta:
        get_latest_by = 'id'


# model for file V-REGION-mutation-and-AA-change-table.txt
class V_REGION_mutation_and_AA_change_table(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE,
                                            related_name='V_REGION_mutation_and_AA_change_table')
    V_REGION = models.TextField(null=True)  # V-REGION
    FR1_IMGT = models.TextField(null=True)  # FR1-IMGT
    CDR1_IMGT = models.TextField(null=True)  # CDR1-IMGT
    FR2_IMGT = models.TextField(null=True)  # FR2-IMGT
    CDR2_IMGT = models.TextField(null=True)  # CDR2-IMGT
    FR3_IMGT = models.TextField(null=True)  # FR3-IMGT
    CDR3_IMGT = models.TextField(null=True)  # CDR3-IMGT

    class Meta:
        get_latest_by = 'id'


# model for file V-REGION-nt-mutation-statistics.txt
class V_REGION_nt_mutation_statistics(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE,
                                            related_name='V_REGION_nt_mutation_statistics')
    V_REGION_Nb_of_positions_including_IMGT_gaps_nt = models.IntegerField()  # V-REGION Nb of positions including IMGT gaps (nt)
    V_REGION_Nb_of_nucleotides = models.IntegerField()  # V-REGION Nb of nucleotides
    V_REGION_Nb_of_identical_nucleotides = models.IntegerField()  # V-REGION Nb of identical nucleotides
    V_REGION_Nb_of_mutations = models.IntegerField()  # V-REGION Nb of mutations
    V_REGION_Nb_of_silent_mutations = models.IntegerField()  # V-REGION Nb of silent mutations
    V_REGION_Nb_of_nonsilent_mutations = models.IntegerField()  # V-REGION Nb of nonsilent mutations
    V_REGION_a_to_g = models.IntegerField()  # V-REGION a>g
    V_REGION_g_to_a = models.IntegerField()  # V-REGION g>a
    V_REGION_c_to_t = models.IntegerField()  # V-REGION c>t
    V_REGION_t_to_c = models.IntegerField()  # V-REGION t>c
    V_REGION_a_to_c = models.IntegerField()  # V-REGION a>c
    V_REGION_c_to_a = models.IntegerField()  # V-REGION c>a
    V_REGION_a_to_t = models.IntegerField()  # V-REGION a>t
    V_REGION_t_to_a = models.IntegerField()  # V-REGION t>a
    V_REGION_g_to_c = models.IntegerField()  # V-REGION g>c
    V_REGION_c_to_g = models.IntegerField()  # V-REGION c>g
    V_REGION_g_to_t = models.IntegerField()  # V-REGION g>t
    V_REGION_t_to_g = models.IntegerField()  # V-REGION t>g
    FR1_IMGT_Nb_of_positions_including_IMGT_gaps_nt = models.IntegerField()  # FR1-IMGT Nb of positions including IMGT gaps (nt)
    FR1_IMGT_Nb_of_nucleotides = models.IntegerField()  # FR1-IMGT Nb of nucleotides
    FR1_IMGT_Nb_of_identical_nucleotides = models.IntegerField()  # FR1-IMGT Nb of identical nucleotides
    FR1_IMGT_Nb_of_mutations = models.IntegerField()  # FR1-IMGT Nb of mutations
    FR1_IMGT_Nb_of_silent_mutations = models.IntegerField()  # FR1-IMGT Nb of silent mutations
    FR1_IMGT_Nb_of_nonsilent_mutations = models.IntegerField()  # FR1-IMGT Nb of nonsilent mutations
    FR1_IMGT_a_to_g = models.IntegerField()  # FR1-IMGT a>g
    FR1_IMGT_g_to_a = models.IntegerField()  # FR1-IMGT g>a
    FR1_IMGT_c_to_t = models.IntegerField()  # FR1-IMGT c>t
    FR1_IMGT_t_to_c = models.IntegerField()  # FR1-IMGT t>c
    FR1_IMGT_a_to_c = models.IntegerField()  # FR1-IMGT a>c
    FR1_IMGT_c_to_a = models.IntegerField()  # FR1-IMGT c>a
    FR1_IMGT_a_to_t = models.IntegerField()  # FR1-IMGT a>t
    FR1_IMGT_t_to_a = models.IntegerField()  # FR1-IMGT t>a
    FR1_IMGT_g_to_c = models.IntegerField()  # FR1-IMGT g>c
    FR1_IMGT_c_to_g = models.IntegerField()  # FR1-IMGT c>g
    FR1_IMGT_g_to_t = models.IntegerField()  # FR1-IMGT g>t
    FR1_IMGT_t_to_g = models.IntegerField()  # FR1-IMGT t>g
    CDR1_IMGT_Nb_of_positions_including_IMGT_gaps_nt = models.IntegerField()  # CDR1-IMGT Nb of positions including IMGT gaps (nt)
    CDR1_IMGT_Nb_of_nucleotides = models.IntegerField()  # CDR1-IMGT Nb of nucleotides
    CDR1_IMGT_Nb_of_identical_nucleotides = models.IntegerField()  # CDR1-IMGT Nb of identical nucleotides
    CDR1_IMGT_Nb_of_mutations = models.IntegerField()  # CDR1-IMGT Nb of mutations
    CDR1_IMGT_Nb_of_silent_mutations = models.IntegerField()  # CDR1-IMGT Nb of silent mutations
    CDR1_IMGT_Nb_of_nonsilent_mutations = models.IntegerField()  # CDR1-IMGT Nb of nonsilent mutations
    CDR1_IMGT_a_to_g = models.IntegerField()  # CDR1-IMGT a>g
    CDR1_IMGT_g_to_a = models.IntegerField()  # CDR1-IMGT g>a
    CDR1_IMGT_c_to_t = models.IntegerField()  # CDR1-IMGT c>t
    CDR1_IMGT_t_to_c = models.IntegerField()  # CDR1-IMGT t>c
    CDR1_IMGT_a_to_c = models.IntegerField()  # CDR1-IMGT a>c
    CDR1_IMGT_c_to_a = models.IntegerField()  # CDR1-IMGT c>a
    CDR1_IMGT_a_to_t = models.IntegerField()  # CDR1-IMGT a>t
    CDR1_IMGT_t_to_a = models.IntegerField()  # CDR1-IMGT t>a
    CDR1_IMGT_g_to_c = models.IntegerField()  # CDR1-IMGT g>c
    CDR1_IMGT_c_to_g = models.IntegerField()  # CDR1-IMGT c>g
    CDR1_IMGT_g_to_t = models.IntegerField()  # CDR1-IMGT g>t
    CDR1_IMGT_t_to_g = models.IntegerField()  # CDR1-IMGT t>g
    FR2_IMGT_Nb_of_positions_including_IMGT_gaps_nt = models.IntegerField()  # FR2-IMGT Nb of positions including IMGT gaps (nt)
    FR2_IMGT_Nb_of_nucleotides = models.IntegerField()  # FR2-IMGT Nb of nucleotides
    FR2_IMGT_Nb_of_identical_nucleotides = models.IntegerField()  # FR2-IMGT Nb of identical nucleotides
    FR2_IMGT_Nb_of_mutations = models.IntegerField()  # FR2-IMGT Nb of mutations
    FR2_IMGT_Nb_of_silent_mutations = models.IntegerField()  # FR2-IMGT Nb of silent mutations
    FR2_IMGT_Nb_of_nonsilent_mutations = models.IntegerField()  # FR2-IMGT Nb of nonsilent mutations
    FR2_IMGT_a_to_g = models.IntegerField()  # FR2-IMGT a>g
    FR2_IMGT_g_to_a = models.IntegerField()  # FR2-IMGT g>a
    FR2_IMGT_c_to_t = models.IntegerField()  # FR2-IMGT c>t
    FR2_IMGT_t_to_c = models.IntegerField()  # FR2-IMGT t>c
    FR2_IMGT_a_to_c = models.IntegerField()  # FR2-IMGT a>c
    FR2_IMGT_c_to_a = models.IntegerField()  # FR2-IMGT c>a
    FR2_IMGT_a_to_t = models.IntegerField()  # FR2-IMGT a>t
    FR2_IMGT_t_to_a = models.IntegerField()  # FR2-IMGT t>a
    FR2_IMGT_g_to_c = models.IntegerField()  # FR2-IMGT g>c
    FR2_IMGT_c_to_g = models.IntegerField()  # FR2-IMGT c>g
    FR2_IMGT_g_to_t = models.IntegerField()  # FR2-IMGT g>t
    FR2_IMGT_t_to_g = models.IntegerField()  # FR2-IMGT t>g
    CDR2_IMGT_Nb_of_positions_including_IMGT_gaps_nt = models.IntegerField()  # CDR2-IMGT Nb of positions including IMGT gaps (nt)
    CDR2_IMGT_Nb_of_nucleotides = models.IntegerField()  # CDR2-IMGT Nb of nucleotides
    CDR2_IMGT_Nb_of_identical_nucleotides = models.IntegerField()  # CDR2-IMGT Nb of identical nucleotides
    CDR2_IMGT_Nb_of_mutations = models.IntegerField()  # CDR2-IMGT Nb of mutations
    CDR2_IMGT_Nb_of_silent_mutations = models.IntegerField()  # CDR2-IMGT Nb of silent mutations
    CDR2_IMGT_Nb_of_nonsilent_mutations = models.IntegerField()  # CDR2-IMGT Nb of nonsilent mutations
    CDR2_IMGT_a_to_g = models.IntegerField()  # CDR2-IMGT a>g
    CDR2_IMGT_g_to_a = models.IntegerField()  # CDR2-IMGT g>a
    CDR2_IMGT_c_to_t = models.IntegerField()  # CDR2-IMGT c>t
    CDR2_IMGT_t_to_c = models.IntegerField()  # CDR2-IMGT t>c
    CDR2_IMGT_a_to_c = models.IntegerField()  # CDR2-IMGT a>c
    CDR2_IMGT_c_to_a = models.IntegerField()  # CDR2-IMGT c>a
    CDR2_IMGT_a_to_t = models.IntegerField()  # CDR2-IMGT a>t
    CDR2_IMGT_t_to_a = models.IntegerField()  # CDR2-IMGT t>a
    CDR2_IMGT_g_to_c = models.IntegerField()  # CDR2-IMGT g>c
    CDR2_IMGT_c_to_g = models.IntegerField()  # CDR2-IMGT c>g
    CDR2_IMGT_g_to_t = models.IntegerField()  # CDR2-IMGT g>t
    CDR2_IMGT_t_to_g = models.IntegerField()  # CDR2-IMGT t>g
    FR3_IMGT_Nb_of_positions_including_IMGT_gaps_nt = models.IntegerField()  # FR3-IMGT Nb of positions including IMGT gaps (nt)
    FR3_IMGT_Nb_of_nucleotides = models.IntegerField()  # FR3-IMGT Nb of nucleotides
    FR3_IMGT_Nb_of_identical_nucleotides = models.IntegerField()  # FR3-IMGT Nb of identical nucleotides
    FR3_IMGT_Nb_of_mutations = models.IntegerField()  # FR3-IMGT Nb of mutations
    FR3_IMGT_Nb_of_silent_mutations = models.IntegerField()  # FR3-IMGT Nb of silent mutations
    FR3_IMGT_Nb_of_nonsilent_mutations = models.IntegerField()  # FR3-IMGT Nb of nonsilent mutations
    FR3_IMGT_a_to_g = models.IntegerField()  # FR3-IMGT a>g
    FR3_IMGT_g_to_a = models.IntegerField()  # FR3-IMGT g>a
    FR3_IMGT_c_to_t = models.IntegerField()  # FR3-IMGT c>t
    FR3_IMGT_t_to_c = models.IntegerField()  # FR3-IMGT t>c
    FR3_IMGT_a_to_c = models.IntegerField()  # FR3-IMGT a>c
    FR3_IMGT_c_to_a = models.IntegerField()  # FR3-IMGT c>a
    FR3_IMGT_a_to_t = models.IntegerField()  # FR3-IMGT a>t
    FR3_IMGT_t_to_a = models.IntegerField()  # FR3-IMGT t>a
    FR3_IMGT_g_to_c = models.IntegerField()  # FR3-IMGT g>c
    FR3_IMGT_c_to_g = models.IntegerField()  # FR3-IMGT c>g
    FR3_IMGT_g_to_t = models.IntegerField()  # FR3-IMGT g>t
    FR3_IMGT_t_to_g = models.IntegerField()  # FR3-IMGT t>g
    CDR3_IMGT_Nb_of_positions_including_IMGT_gaps_nt = models.IntegerField()  # CDR3-IMGT Nb of positions including IMGT gaps (nt)
    CDR3_IMGT_Nb_of_nucleotides = models.IntegerField()  # CDR3-IMGT Nb of nucleotides
    CDR3_IMGT_Nb_of_identical_nucleotides = models.IntegerField()  # CDR3-IMGT Nb of identical nucleotides
    CDR3_IMGT_Nb_of_mutations = models.IntegerField()  # CDR3-IMGT Nb of mutations
    CDR3_IMGT_Nb_of_silent_mutations = models.IntegerField(null=True)  # CDR3-IMGT Nb of silent mutations
    CDR3_IMGT_Nb_of_nonsilent_mutations = models.IntegerField(null=True)  # CDR3-IMGT Nb of nonsilent mutations
    CDR3_IMGT_a_to_g = models.IntegerField(null=True)  # CDR3-IMGT a>g
    CDR3_IMGT_g_to_a = models.IntegerField(null=True)  # CDR3-IMGT g>a
    CDR3_IMGT_c_to_t = models.IntegerField(null=True)  # CDR3-IMGT c>t
    CDR3_IMGT_t_to_c = models.IntegerField(null=True)  # CDR3-IMGT t>c
    CDR3_IMGT_a_to_c = models.IntegerField(null=True)  # CDR3-IMGT a>c
    CDR3_IMGT_c_to_a = models.IntegerField(null=True)  # CDR3-IMGT c>a
    CDR3_IMGT_a_to_t = models.IntegerField(null=True)  # CDR3-IMGT a>t
    CDR3_IMGT_t_to_a = models.IntegerField(null=True)  # CDR3-IMGT t>a
    CDR3_IMGT_g_to_c = models.IntegerField(null=True)  # CDR3-IMGT g>c
    CDR3_IMGT_c_to_g = models.IntegerField(null=True)  # CDR3-IMGT c>g
    CDR3_IMGT_g_to_t = models.IntegerField(null=True)  # CDR3-IMGT g>t
    CDR3_IMGT_t_to_g = models.IntegerField(null=True)  # CDR3-IMGT t>g

    class Meta:
        get_latest_by = 'id'


# model for file V-REGION-AA-change-statistics.txt
class V_REGION_AA_change_statistics(models.Model):
    Sequence_Identifier = models.OneToOneField(Sequence, on_delete=models.CASCADE,
                                            related_name='V_REGION_AA_change_statistics')
    V_REGION_Nb_of_positions_including_IMGT_gaps_AA = models.IntegerField()  # V-REGION Nb of positions including IMGT gaps (AA)
    V_REGION_Nb_of_AA = models.IntegerField()  # V-REGION Nb of AA
    V_REGION_Nb_of_identical_AA = models.IntegerField()  # V-REGION Nb of identical AA
    V_REGION_Nb_of_AA_changes = models.IntegerField()  # V-REGION Nb of AA changes
    V_REGION_PPP = models.IntegerField()  # V-REGION +++
    V_REGION_PPN = models.IntegerField()  # V-REGION ++-
    V_REGION_PNP = models.IntegerField()  # V-REGION +-+
    V_REGION_PNN = models.IntegerField()  # V-REGION +--
    V_REGION_NPN = models.IntegerField()  # V-REGION -+-
    V_REGION_NNP = models.IntegerField()  # V-REGION --+
    V_REGION_NNN = models.IntegerField()  # V-REGION ---
    V_REGION_Very_similar = models.IntegerField()  # V-REGION Very similar
    V_REGION_Similar = models.IntegerField()  # V-REGION Similar
    V_REGION_Dissimilar = models.IntegerField()  # V-REGION Dissimilar
    V_REGION_Very_dissimilar = models.IntegerField()  # V-REGION Very dissimilar
    FR1_IMGT_Nb_of_positions_including_IMGT_gaps_AA = models.IntegerField()  # FR1-IMGT Nb of positions including IMGT gaps (AA)
    FR1_IMGT_Nb_of_AA = models.IntegerField()  # FR1-IMGT Nb of AA
    FR1_IMGT_Nb_of_identical_AA = models.IntegerField()  # FR1-IMGT Nb of identical AA
    FR1_IMGT_Nb_of_AA_changes = models.IntegerField()  # FR1-IMGT Nb of AA changes
    FR1_IMGT_PPP = models.IntegerField()  # FR1-IMGT +++
    FR1_IMGT_PPN = models.IntegerField()  # FR1-IMGT ++-
    FR1_IMGT_PNP = models.IntegerField()  # FR1-IMGT +-+
    FR1_IMGT_PNN = models.IntegerField()  # FR1-IMGT +--
    FR1_IMGT_NPN = models.IntegerField()  # FR1-IMGT -+-
    FR1_IMGT_NNP = models.IntegerField()  # FR1-IMGT --+
    FR1_IMGT_NNN = models.IntegerField()  # FR1-IMGT ---
    FR1_IMGT_Very_similar = models.IntegerField()  # FR1-IMGT Very similar
    FR1_IMGT_Similar = models.IntegerField()  # FR1-IMGT Similar
    FR1_IMGT_Dissimilar = models.IntegerField()  # FR1-IMGT Dissimilar
    FR1_IMGT_Very_dissimilar = models.IntegerField()  # FR1-IMGT Very dissimilar
    CDR1_IMGT_Nb_of_positions_including_IMGT_gaps_AA = models.IntegerField()  # CDR1-IMGT Nb of positions including IMGT gaps (AA)
    CDR1_IMGT_Nb_of_AA = models.IntegerField()  # CDR1-IMGT Nb of AA
    CDR1_IMGT_Nb_of_identical_AA = models.IntegerField()  # CDR1-IMGT Nb of identical AA
    CDR1_IMGT_Nb_of_AA_changes = models.IntegerField()  # CDR1-IMGT Nb of AA changes
    CDR1_IMGT_PPP = models.IntegerField()  # CDR1-IMGT +++
    CDR1_IMGT_PPN = models.IntegerField()  # CDR1-IMGT ++-
    CDR1_IMGT_PNP = models.IntegerField()  # CDR1-IMGT +-+
    CDR1_IMGT_PNN = models.IntegerField()  # CDR1-IMGT +--
    CDR1_IMGT_NPN = models.IntegerField()  # CDR1-IMGT -+-
    CDR1_IMGT_NNP = models.IntegerField()  # CDR1-IMGT --+
    CDR1_IMGT_NNN = models.IntegerField()  # CDR1-IMGT ---
    CDR1_IMGT_Very_similar = models.IntegerField()  # CDR1-IMGT Very similar
    CDR1_IMGT_Similar = models.IntegerField()  # CDR1-IMGT Similar
    CDR1_IMGT_Dissimilar = models.IntegerField()  # CDR1-IMGT Dissimilar
    CDR1_IMGT_Very_dissimilar = models.IntegerField()  # CDR1-IMGT Very dissimilar
    FR2_IMGT_Nb_of_positions_including_IMGT_gaps_AA = models.IntegerField()  # FR2-IMGT Nb of positions including IMGT gaps (AA)
    FR2_IMGT_Nb_of_AA = models.IntegerField()  # FR2-IMGT Nb of AA
    FR2_IMGT_Nb_of_identical_AA = models.IntegerField()  # FR2-IMGT Nb of identical AA
    FR2_IMGT_Nb_of_AA_changes = models.IntegerField()  # FR2-IMGT Nb of AA changes
    FR2_IMGT_PPP = models.IntegerField()  # FR2-IMGT +++
    FR2_IMGT_PPN = models.IntegerField()  # FR2-IMGT ++-
    FR2_IMGT_PNP = models.IntegerField()  # FR2-IMGT +-+
    FR2_IMGT_PNN = models.IntegerField()  # FR2-IMGT +--
    FR2_IMGT_NPN = models.IntegerField()  # FR2-IMGT -+-
    FR2_IMGT_NNP = models.IntegerField()  # FR2-IMGT --+
    FR2_IMGT_NNN = models.IntegerField()  # FR2-IMGT ---
    FR2_IMGT_Very_similar = models.IntegerField()  # FR2-IMGT Very similar
    FR2_IMGT_Similar = models.IntegerField()  # FR2-IMGT Similar
    FR2_IMGT_Dissimilar = models.IntegerField()  # FR2-IMGT Dissimilar
    FR2_IMGT_Very_dissimilar = models.IntegerField()  # FR2-IMGT Very dissimilar
    CDR2_IMGT_Nb_of_positions_including_IMGT_gaps_AA = models.IntegerField()  # CDR2-IMGT Nb of positions including IMGT gaps (AA)
    CDR2_IMGT_Nb_of_AA = models.IntegerField()  # CDR2-IMGT Nb of AA
    CDR2_IMGT_Nb_of_identical_AA = models.IntegerField()  # CDR2-IMGT Nb of identical AA
    CDR2_IMGT_Nb_of_AA_changes = models.IntegerField()  # CDR2-IMGT Nb of AA changes
    CDR2_IMGT_PPP = models.IntegerField()  # CDR2-IMGT +++
    CDR2_IMGT_PPN = models.IntegerField()  # CDR2-IMGT ++-
    CDR2_IMGT_PNP = models.IntegerField()  # CDR2-IMGT +-+
    CDR2_IMGT_PNN = models.IntegerField()  # CDR2-IMGT +--
    CDR2_IMGT_NPN = models.IntegerField()  # CDR2-IMGT -+-
    CDR2_IMGT_NNP = models.IntegerField()  # CDR2-IMGT --+
    CDR2_IMGT_NNN = models.IntegerField()  # CDR2-IMGT ---
    CDR2_IMGT_Very_similar = models.IntegerField()  # CDR2-IMGT Very similar
    CDR2_IMGT_Similar = models.IntegerField()  # CDR2-IMGT Similar
    CDR2_IMGT_Dissimilar = models.IntegerField()  # CDR2-IMGT Dissimilar
    CDR2_IMGT_Very_dissimilar = models.IntegerField()  # CDR2-IMGT Very dissimilar
    FR3_IMGT_Nb_of_positions_including_IMGT_gaps_AA = models.IntegerField()  # FR3-IMGT Nb of positions including IMGT gaps (AA)
    FR3_IMGT_Nb_of_AA = models.IntegerField()  # FR3-IMGT Nb of AA
    FR3_IMGT_Nb_of_identical_AA = models.IntegerField()  # FR3-IMGT Nb of identical AA
    FR3_IMGT_Nb_of_AA_changes = models.IntegerField()  # FR3-IMGT Nb of AA changes
    FR3_IMGT_PPP = models.IntegerField()  # FR3-IMGT +++
    FR3_IMGT_PPN = models.IntegerField()  # FR3-IMGT ++-
    FR3_IMGT_PNP = models.IntegerField()  # FR3-IMGT +-+
    FR3_IMGT_PNN = models.IntegerField()  # FR3-IMGT +--
    FR3_IMGT_NPN = models.IntegerField()  # FR3-IMGT -+-
    FR3_IMGT_NNP = models.IntegerField()  # FR3-IMGT --+
    FR3_IMGT_NNN = models.IntegerField()  # FR3-IMGT ---
    FR3_IMGT_Very_similar = models.IntegerField()  # FR3-IMGT Very similar
    FR3_IMGT_Similar = models.IntegerField()  # FR3-IMGT Similar
    FR3_IMGT_Dissimilar = models.IntegerField()  # FR3-IMGT Dissimilar
    FR3_IMGT_Very_dissimilar = models.IntegerField()  # FR3-IMGT Very dissimilar
    CDR3_IMGT_Nb_of_positions_including_IMGT_gaps_AA = models.IntegerField()  # CDR3-IMGT Nb of positions including IMGT gaps (AA)
    CDR3_IMGT_Nb_of_AA = models.IntegerField()  # CDR3-IMGT Nb of AA
    CDR3_IMGT_Nb_of_identical_AA = models.IntegerField()  # CDR3-IMGT Nb of identical AA
    CDR3_IMGT_Nb_of_AA_changes = models.IntegerField()  # CDR3-IMGT Nb of AA changes
    CDR3_IMGT_PPP = models.IntegerField()  # CDR3-IMGT +++
    CDR3_IMGT_PPN = models.IntegerField()  # CDR3-IMGT ++-
    CDR3_IMGT_PNP = models.IntegerField()  # CDR3-IMGT +-+
    CDR3_IMGT_PNN = models.IntegerField()  # CDR3-IMGT +--
    CDR3_IMGT_NPN = models.IntegerField()  # CDR3-IMGT -+-
    CDR3_IMGT_NNP = models.IntegerField()  # CDR3-IMGT --+
    CDR3_IMGT_NNN = models.IntegerField()  # CDR3-IMGT ---
    CDR3_IMGT_Very_similar = models.IntegerField()  # CDR3-IMGT Very similar
    CDR3_IMGT_Similar = models.IntegerField()  # CDR3-IMGT Similar
    CDR3_IMGT_Dissimilar = models.IntegerField()  # CDR3-IMGT Dissimilar
    CDR3_IMGT_Very_dissimilar = models.IntegerField()  # CDR3-IMGT Very dissimilar

    class Meta:
        get_latest_by = 'id'
