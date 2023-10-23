# Generated by Django 4.1.3 on 2023-02-14 12:48

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0010_alter_mutations_sequence_alter_mutations_bases_aa_and_more'),
    ]

    operations = [
        migrations.CreateModel(
            name='Run_Identifier',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('run_name', models.TextField()),
                ('description', models.TextField()),
                ('created_at', models.DateTimeField(auto_now_add=True)),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='Sequence',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('Sequence_number', models.IntegerField()),
                ('Sequence_ID', models.CharField(max_length=25)),
                ('V_DOMAIN_Functionality', models.CharField(max_length=20)),
                ('V_GENE_and_allele', models.TextField()),
                ('Run', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='Sequence', to='PinkStrawberry.run_identifier')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='V_REGION_nt_mutation_statistics',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('V_REGION_Nb_of_positions_including_IMGT_gaps_nt', models.IntegerField()),
                ('V_REGION_Nb_of_nucleotides', models.IntegerField()),
                ('V_REGION_Nb_of_identical_nucleotides', models.IntegerField()),
                ('V_REGION_Nb_of_mutations', models.IntegerField()),
                ('V_REGION_Nb_of_silent_mutations', models.IntegerField()),
                ('V_REGION_Nb_of_nonsilent_mutations', models.IntegerField()),
                ('V_REGION_a_to_g', models.IntegerField()),
                ('V_REGION_g_to_a', models.IntegerField()),
                ('V_REGION_c_to_t', models.IntegerField()),
                ('V_REGION_t_to_c', models.IntegerField()),
                ('V_REGION_a_to_c', models.IntegerField()),
                ('V_REGION_c_to_a', models.IntegerField()),
                ('V_REGION_a_to_t', models.IntegerField()),
                ('V_REGION_t_to_a', models.IntegerField()),
                ('V_REGION_g_to_c', models.IntegerField()),
                ('V_REGION_c_to_g', models.IntegerField()),
                ('V_REGION_g_to_t', models.IntegerField()),
                ('V_REGION_t_to_g', models.IntegerField()),
                ('FR1_IMGT_Nb_of_positions_including_IMGT_gaps_nt', models.IntegerField()),
                ('FR1_IMGT_Nb_of_nucleotides', models.IntegerField()),
                ('FR1_IMGT_Nb_of_identical_nucleotides', models.IntegerField()),
                ('FR1_IMGT_Nb_of_mutations', models.IntegerField()),
                ('FR1_IMGT_Nb_of_silent_mutations', models.IntegerField()),
                ('FR1_IMGT_Nb_of_nonsilent_mutations', models.IntegerField()),
                ('FR1_IMGT_a_to_g', models.IntegerField()),
                ('FR1_IMGT_g_to_a', models.IntegerField()),
                ('FR1_IMGT_c_to_t', models.IntegerField()),
                ('FR1_IMGT_t_to_c', models.IntegerField()),
                ('FR1_IMGT_a_to_c', models.IntegerField()),
                ('FR1_IMGT_c_to_a', models.IntegerField()),
                ('FR1_IMGT_a_to_t', models.IntegerField()),
                ('FR1_IMGT_t_to_a', models.IntegerField()),
                ('FR1_IMGT_g_to_c', models.IntegerField()),
                ('FR1_IMGT_c_to_g', models.IntegerField()),
                ('FR1_IMGT_g_to_t', models.IntegerField()),
                ('FR1_IMGT_t_to_g', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_positions_including_IMGT_gaps_nt', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_nucleotides', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_identical_nucleotides', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_mutations', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_silent_mutations', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_nonsilent_mutations', models.IntegerField()),
                ('CDR1_IMGT_a_to_g', models.IntegerField()),
                ('CDR1_IMGT_g_to_a', models.IntegerField()),
                ('CDR1_IMGT_c_to_t', models.IntegerField()),
                ('CDR1_IMGT_t_to_c', models.IntegerField()),
                ('CDR1_IMGT_a_to_c', models.IntegerField()),
                ('CDR1_IMGT_c_to_a', models.IntegerField()),
                ('CDR1_IMGT_a_to_t', models.IntegerField()),
                ('CDR1_IMGT_t_to_a', models.IntegerField()),
                ('CDR1_IMGT_g_to_c', models.IntegerField()),
                ('CDR1_IMGT_c_to_g', models.IntegerField()),
                ('CDR1_IMGT_g_to_t', models.IntegerField()),
                ('CDR1_IMGT_t_to_g', models.IntegerField()),
                ('FR2_IMGT_Nb_of_positions_including_IMGT_gaps_nt', models.IntegerField()),
                ('FR2_IMGT_Nb_of_nucleotides', models.IntegerField()),
                ('FR2_IMGT_Nb_of_identical_nucleotides', models.IntegerField()),
                ('FR2_IMGT_Nb_of_mutations', models.IntegerField()),
                ('FR2_IMGT_Nb_of_silent_mutations', models.IntegerField()),
                ('FR2_IMGT_Nb_of_nonsilent_mutations', models.IntegerField()),
                ('FR2_IMGT_a_to_g', models.IntegerField()),
                ('FR2_IMGT_g_to_a', models.IntegerField()),
                ('FR2_IMGT_c_to_t', models.IntegerField()),
                ('FR2_IMGT_t_to_c', models.IntegerField()),
                ('FR2_IMGT_a_to_c', models.IntegerField()),
                ('FR2_IMGT_c_to_a', models.IntegerField()),
                ('FR2_IMGT_a_to_t', models.IntegerField()),
                ('FR2_IMGT_t_to_a', models.IntegerField()),
                ('FR2_IMGT_g_to_c', models.IntegerField()),
                ('FR2_IMGT_c_to_g', models.IntegerField()),
                ('FR2_IMGT_g_to_t', models.IntegerField()),
                ('FR2_IMGT_t_to_g', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_positions_including_IMGT_gaps_nt', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_nucleotides', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_identical_nucleotides', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_mutations', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_silent_mutations', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_nonsilent_mutations', models.IntegerField()),
                ('CDR2_IMGT_a_to_g', models.IntegerField()),
                ('CDR2_IMGT_g_to_a', models.IntegerField()),
                ('CDR2_IMGT_c_to_t', models.IntegerField()),
                ('CDR2_IMGT_t_to_c', models.IntegerField()),
                ('CDR2_IMGT_a_to_c', models.IntegerField()),
                ('CDR2_IMGT_c_to_a', models.IntegerField()),
                ('CDR2_IMGT_a_to_t', models.IntegerField()),
                ('CDR2_IMGT_t_to_a', models.IntegerField()),
                ('CDR2_IMGT_g_to_c', models.IntegerField()),
                ('CDR2_IMGT_c_to_g', models.IntegerField()),
                ('CDR2_IMGT_g_to_t', models.IntegerField()),
                ('CDR2_IMGT_t_to_g', models.IntegerField()),
                ('FR3_IMGT_Nb_of_positions_including_IMGT_gaps_nt', models.IntegerField()),
                ('FR3_IMGT_Nb_of_nucleotides', models.IntegerField()),
                ('FR3_IMGT_Nb_of_identical_nucleotides', models.IntegerField()),
                ('FR3_IMGT_Nb_of_mutations', models.IntegerField()),
                ('FR3_IMGT_Nb_of_silent_mutations', models.IntegerField()),
                ('FR3_IMGT_Nb_of_nonsilent_mutations', models.IntegerField()),
                ('FR3_IMGT_a_to_g', models.IntegerField()),
                ('FR3_IMGT_g_to_a', models.IntegerField()),
                ('FR3_IMGT_c_to_t', models.IntegerField()),
                ('FR3_IMGT_t_to_c', models.IntegerField()),
                ('FR3_IMGT_a_to_c', models.IntegerField()),
                ('FR3_IMGT_c_to_a', models.IntegerField()),
                ('FR3_IMGT_a_to_t', models.IntegerField()),
                ('FR3_IMGT_t_to_a', models.IntegerField()),
                ('FR3_IMGT_g_to_c', models.IntegerField()),
                ('FR3_IMGT_c_to_g', models.IntegerField()),
                ('FR3_IMGT_g_to_t', models.IntegerField()),
                ('FR3_IMGT_t_to_g', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_positions_including_IMGT_gaps_nt', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_nucleotides', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_identical_nucleotides', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_mutations', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_silent_mutations', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_nonsilent_mutations', models.IntegerField()),
                ('CDR3_IMGT_a_to_g', models.IntegerField()),
                ('CDR3_IMGT_g_to_a', models.IntegerField()),
                ('CDR3_IMGT_c_to_t', models.IntegerField()),
                ('CDR3_IMGT_t_to_c', models.IntegerField()),
                ('CDR3_IMGT_a_to_c', models.IntegerField()),
                ('CDR3_IMGT_c_to_a', models.IntegerField()),
                ('CDR3_IMGT_a_to_t', models.IntegerField()),
                ('CDR3_IMGT_t_to_a', models.IntegerField()),
                ('CDR3_IMGT_g_to_c', models.IntegerField()),
                ('CDR3_IMGT_c_to_g', models.IntegerField()),
                ('CDR3_IMGT_g_to_t', models.IntegerField()),
                ('CDR3_IMGT_t_to_g', models.IntegerField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='V_REGION_nt_mutation_statistics', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='V_REGION_mutation_hotspots',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('ata', models.TextField()),
                ('tat', models.TextField()),
                ('aggctat', models.TextField()),
                ('atagcct', models.TextField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='V_REGION_mutation_hotspots', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='V_REGION_mutation_and_AA_change_table',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('V_REGION', models.TextField()),
                ('FR1_IMGT', models.TextField()),
                ('CDR1_IMGT', models.TextField()),
                ('FR2_IMGT', models.TextField()),
                ('CDR2_IMGT', models.TextField()),
                ('FR3_IMGT', models.TextField()),
                ('CDR3_IMGT', models.TextField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='V_REGION_mutation_and_AA_change_table', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='V_REGION_AA_change_statistics',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('V_REGION_Nb_of_positions_including_IMGT_gaps_AA', models.IntegerField()),
                ('V_REGION_Nb_of_AA', models.IntegerField()),
                ('V_REGION_Nb_of_identical_AA', models.IntegerField()),
                ('V_REGION_Nb_of_AA_changes', models.IntegerField()),
                ('V_REGION_PPP', models.IntegerField()),
                ('V_REGION_PPN', models.IntegerField()),
                ('V_REGION_PNP', models.IntegerField()),
                ('V_REGION_PNN', models.IntegerField()),
                ('V_REGION_NPN', models.IntegerField()),
                ('V_REGION_NNP', models.IntegerField()),
                ('V_REGION_NNN', models.IntegerField()),
                ('V_REGION_Very_similar', models.IntegerField()),
                ('V_REGION_Similar', models.IntegerField()),
                ('V_REGION_Dissimilar', models.IntegerField()),
                ('V_REGION_Very_dissimilar', models.IntegerField()),
                ('FR1_IMGT_Nb_of_positions_including_IMGT_gaps_AA', models.IntegerField()),
                ('FR1_IMGT_Nb_of_AA', models.IntegerField()),
                ('FR1_IMGT_Nb_of_identical_AA', models.IntegerField()),
                ('FR1_IMGT_Nb_of_AA_changes', models.IntegerField()),
                ('FR1_IMGT_PPP', models.IntegerField()),
                ('FR1_IMGT_PPN', models.IntegerField()),
                ('FR1_IMGT_PNP', models.IntegerField()),
                ('FR1_IMGT_PNN', models.IntegerField()),
                ('FR1_IMGT_NPN', models.IntegerField()),
                ('FR1_IMGT_NNP', models.IntegerField()),
                ('FR1_IMGT_NNN', models.IntegerField()),
                ('FR1_IMGT_Very_similar', models.IntegerField()),
                ('FR1_IMGT_Similar', models.IntegerField()),
                ('FR1_IMGT_Dissimilar', models.IntegerField()),
                ('FR1_IMGT_Very_dissimilar', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_positions_including_IMGT_gaps_AA', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_AA', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_identical_AA', models.IntegerField()),
                ('CDR1_IMGT_Nb_of_AA_changes', models.IntegerField()),
                ('CDR1_IMGT_PPP', models.IntegerField()),
                ('CDR1_IMGT_PPN', models.IntegerField()),
                ('CDR1_IMGT_PNP', models.IntegerField()),
                ('CDR1_IMGT_PNN', models.IntegerField()),
                ('CDR1_IMGT_NPN', models.IntegerField()),
                ('CDR1_IMGT_NNP', models.IntegerField()),
                ('CDR1_IMGT_NNN', models.IntegerField()),
                ('CDR1_IMGT_Very_similar', models.IntegerField()),
                ('CDR1_IMGT_Similar', models.IntegerField()),
                ('CDR1_IMGT_Dissimilar', models.IntegerField()),
                ('CDR1_IMGT_Very_dissimilar', models.IntegerField()),
                ('FR2_IMGT_Nb_of_positions_including_IMGT_gaps_AA', models.IntegerField()),
                ('FR2_IMGT_Nb_of_AA', models.IntegerField()),
                ('FR2_IMGT_Nb_of_identical_AA', models.IntegerField()),
                ('FR2_IMGT_Nb_of_AA_changes', models.IntegerField()),
                ('FR2_IMGT_PPP', models.IntegerField()),
                ('FR2_IMGT_PPN', models.IntegerField()),
                ('FR2_IMGT_PNP', models.IntegerField()),
                ('FR2_IMGT_PNN', models.IntegerField()),
                ('FR2_IMGT_NPN', models.IntegerField()),
                ('FR2_IMGT_NNP', models.IntegerField()),
                ('FR2_IMGT_NNN', models.IntegerField()),
                ('FR2_IMGT_Very_similar', models.IntegerField()),
                ('FR2_IMGT_Similar', models.IntegerField()),
                ('FR2_IMGT_Dissimilar', models.IntegerField()),
                ('FR2_IMGT_Very_dissimilar', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_positions_including_IMGT_gaps_AA', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_AA', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_identical_AA', models.IntegerField()),
                ('CDR2_IMGT_Nb_of_AA_changes', models.IntegerField()),
                ('CDR2_IMGT_PPP', models.IntegerField()),
                ('CDR2_IMGT_PPN', models.IntegerField()),
                ('CDR2_IMGT_PNP', models.IntegerField()),
                ('CDR2_IMGT_PNN', models.IntegerField()),
                ('CDR2_IMGT_NPN', models.IntegerField()),
                ('CDR2_IMGT_NNP', models.IntegerField()),
                ('CDR2_IMGT_NNN', models.IntegerField()),
                ('CDR2_IMGT_Very_similar', models.IntegerField()),
                ('CDR2_IMGT_Similar', models.IntegerField()),
                ('CDR2_IMGT_Dissimilar', models.IntegerField()),
                ('CDR2_IMGT_Very_dissimilar', models.IntegerField()),
                ('FR3_IMGT_Nb_of_positions_including_IMGT_gaps_AA', models.IntegerField()),
                ('FR3_IMGT_Nb_of_AA', models.IntegerField()),
                ('FR3_IMGT_Nb_of_identical_AA', models.IntegerField()),
                ('FR3_IMGT_Nb_of_AA_changes', models.IntegerField()),
                ('FR3_IMGT_PPP', models.IntegerField()),
                ('FR3_IMGT_PPN', models.IntegerField()),
                ('FR3_IMGT_PNP', models.IntegerField()),
                ('FR3_IMGT_PNN', models.IntegerField()),
                ('FR3_IMGT_NPN', models.IntegerField()),
                ('FR3_IMGT_NNP', models.IntegerField()),
                ('FR3_IMGT_NNN', models.IntegerField()),
                ('FR3_IMGT_Very_similar', models.IntegerField()),
                ('FR3_IMGT_Similar', models.IntegerField()),
                ('FR3_IMGT_Dissimilar', models.IntegerField()),
                ('FR3_IMGT_Very_dissimilar', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_positions_including_IMGT_gaps_AA', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_AA', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_identical_AA', models.IntegerField()),
                ('CDR3_IMGT_Nb_of_AA_changes', models.IntegerField()),
                ('CDR3_IMGT_PPP', models.IntegerField()),
                ('CDR3_IMGT_PPN', models.IntegerField()),
                ('CDR3_IMGT_PNP', models.IntegerField()),
                ('CDR3_IMGT_PNN', models.IntegerField()),
                ('CDR3_IMGT_NPN', models.IntegerField()),
                ('CDR3_IMGT_NNP', models.IntegerField()),
                ('CDR3_IMGT_NNN', models.IntegerField()),
                ('CDR3_IMGT_Very_similar', models.IntegerField()),
                ('CDR3_IMGT_Similar', models.IntegerField()),
                ('CDR3_IMGT_Dissimilar', models.IntegerField()),
                ('CDR3_IMGT_Very_dissimilar', models.IntegerField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='V_REGION_AA_change_statistics', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='Summary',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('V_REGION_score', models.IntegerField()),
                ('V_REGION_identity_percent', models.FloatField()),
                ('V_REGION_identity_nt', models.TextField()),
                ('V_REGION_identity_percent_with_insdel_events', models.FloatField()),
                ('V_REGION_identity_nt_with_insdel_events', models.TextField()),
                ('J_GENE_and_allele', models.TextField()),
                ('J_REGION_score', models.IntegerField()),
                ('J_REGION_identity_percent', models.FloatField()),
                ('J_REGION_identity_nt', models.TextField()),
                ('D_GENE_and_allele', models.TextField()),
                ('D_REGION_reading_frame', models.IntegerField()),
                ('CDR1_IMGT_length', models.IntegerField()),
                ('CDR2_IMGT_length', models.IntegerField()),
                ('CDR3_IMGT_length', models.IntegerField()),
                ('CDR_IMGT_lengths', models.TextField()),
                ('FR_IMGT_lengths', models.TextField()),
                ('AA_JUNCTION', models.TextField()),
                ('JUNCTION_frame', models.TextField()),
                ('Orientation', models.CharField(max_length=1)),
                ('V_DOMAIN_Functionality_comment', models.TextField()),
                ('V_REGION_potential_insdel', models.TextField()),
                ('J_GENE_and_allele_comment', models.TextField()),
                ('V_REGION_insertions', models.TextField()),
                ('V_REGION_deletions', models.TextField()),
                ('Sequence', models.TextField()),
                ('fiveprime_trimmed_n_nb', models.IntegerField()),
                ('threeprime_trimmed_n_nb', models.IntegerField()),
                ('Analysed_sequence_length', models.IntegerField()),
                ('Sequence_analysis_category', models.TextField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='Summary', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='Nt_sequences',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('J_GENE_and_allele', models.TextField()),
                ('D_GENE_and_allele', models.TextField()),
                ('V_D_J_REGION', models.TextField()),
                ('V_J_REGION', models.TextField()),
                ('V_REGION', models.TextField()),
                ('FR1_IMGT', models.TextField()),
                ('CDR1_IMGT', models.TextField()),
                ('FR2_IMGT', models.TextField()),
                ('CDR2_IMGT', models.TextField()),
                ('FR3_IMGT', models.TextField()),
                ('CDR3_IMGT', models.TextField()),
                ('JUNCTION', models.TextField()),
                ('threeprime_V_REGION', models.TextField()),
                ('N_D_J_REGION', models.TextField()),
                ('N_D_REGION', models.TextField()),
                ('Pthreeprime_V', models.TextField()),
                ('N_REGION', models.TextField()),
                ('N1_REGION', models.TextField()),
                ('Pfiveprime_D', models.TextField()),
                ('D_REGION', models.TextField()),
                ('Pthreeprime_D', models.TextField()),
                ('Pfiveprime_D1', models.TextField()),
                ('D1_REGION', models.TextField()),
                ('Pthreeprime_D1', models.TextField()),
                ('N2_REGION', models.TextField()),
                ('Pfiveprime_D2', models.TextField()),
                ('D2_REGION', models.TextField()),
                ('Pthreeprime_D2', models.TextField()),
                ('N3_REGION', models.TextField()),
                ('Pfiveprime_D3', models.TextField()),
                ('D3_REGION', models.TextField()),
                ('Pthreeprime_D3', models.TextField()),
                ('N4_REGION', models.TextField()),
                ('Pfiveprime_J', models.TextField()),
                ('fiveprime_J_REGION', models.TextField()),
                ('D_J_REGION', models.TextField()),
                ('J_REGION', models.TextField()),
                ('FR4_IMGT', models.TextField()),
                ('V_D_J_REGION_start', models.IntegerField()),
                ('V_D_J_REGION_end', models.IntegerField()),
                ('V_J_REGION_start', models.IntegerField()),
                ('V_J_REGION_end', models.IntegerField()),
                ('V_REGION_start', models.IntegerField()),
                ('V_REGION_end', models.IntegerField()),
                ('FR1_IMGT_start', models.IntegerField()),
                ('FR1_IMGT_end', models.IntegerField()),
                ('CDR1_IMGT_start', models.IntegerField()),
                ('CDR1_IMGT_end', models.IntegerField()),
                ('FR2_IMGT_start', models.IntegerField()),
                ('FR2_IMGT_end', models.IntegerField()),
                ('CDR2_IMGT_start', models.IntegerField()),
                ('CDR2_IMGT_end', models.IntegerField()),
                ('FR3_IMGT_start', models.IntegerField()),
                ('FR3_IMGT_end', models.IntegerField()),
                ('CDR3_IMGT_start', models.IntegerField()),
                ('CDR3_IMGT_end', models.IntegerField()),
                ('JUNCTION_start', models.IntegerField()),
                ('JUNCTION_end', models.IntegerField()),
                ('threeprime_V_REGION_start', models.IntegerField()),
                ('threeprime_V_REGION_end', models.IntegerField()),
                ('N_D_J_REGION_start', models.IntegerField()),
                ('N_D_J_REGION_end', models.IntegerField()),
                ('N_D_REGION_start', models.IntegerField()),
                ('N_D_REGION_end', models.IntegerField()),
                ('Pthreeprime_V_start', models.IntegerField()),
                ('Pthreeprime_V_end', models.IntegerField()),
                ('N_REGION_start', models.IntegerField()),
                ('N_REGION_end', models.IntegerField()),
                ('N1_REGION_start', models.IntegerField()),
                ('N1_REGION_end', models.IntegerField()),
                ('Pfiveprime_D_start', models.IntegerField()),
                ('Pfiveprime_D_end', models.IntegerField()),
                ('D_REGION_start', models.IntegerField()),
                ('D_REGION_end', models.IntegerField()),
                ('Pthreeprime_D_start', models.IntegerField()),
                ('Pthreeprime_D_end', models.IntegerField()),
                ('Pfiveprime_D1_start', models.IntegerField()),
                ('Pfiveprime_D1_end', models.IntegerField()),
                ('D1_REGION_start', models.IntegerField()),
                ('D1_REGION_end', models.IntegerField()),
                ('Pthreeprime_D1_start', models.IntegerField()),
                ('Pthreeprime_D1_end', models.IntegerField()),
                ('N2_REGION_start', models.IntegerField()),
                ('N2_REGION_end', models.IntegerField()),
                ('Pfiveprime_D2_start', models.IntegerField()),
                ('Pfiveprime_D2_end', models.IntegerField()),
                ('D2_REGION_start', models.IntegerField()),
                ('D2_REGION_end', models.IntegerField()),
                ('Pthreeprime_D2_start', models.IntegerField()),
                ('Pthreeprime_D2_end', models.IntegerField()),
                ('N3_REGION_start', models.IntegerField()),
                ('N3_REGION_end', models.IntegerField()),
                ('Pfiveprime_D3_start', models.IntegerField()),
                ('Pfiveprime_D3_end', models.IntegerField()),
                ('D3_REGION_start', models.IntegerField()),
                ('D3_REGION_end', models.IntegerField()),
                ('Pthreeprime_D3_start', models.IntegerField()),
                ('Pthreeprime_D3_end', models.IntegerField()),
                ('N4_REGION_start', models.IntegerField()),
                ('N4_REGION_end', models.IntegerField()),
                ('Pfiveprime_J_start', models.IntegerField()),
                ('Pfiveprime_J_end', models.IntegerField()),
                ('fiveprime_J_REGION_start', models.IntegerField()),
                ('fiveprime_J_REGION_end', models.IntegerField()),
                ('D_J_REGION_start', models.IntegerField()),
                ('D_J_REGION_end', models.IntegerField()),
                ('J_REGION_start', models.IntegerField()),
                ('J_REGION_end', models.IntegerField()),
                ('FR4_IMGT_start', models.IntegerField()),
                ('FR4_IMGT_end', models.IntegerField()),
                ('V_REGION_reading_frame', models.IntegerField()),
                ('V_REGION_partial_fiveprime_missing_nt_nb', models.IntegerField()),
                ('V_REGION_uncertain_nt_nb', models.IntegerField()),
                ('J_REGION_partial_threeprime_missing_nt_nb', models.IntegerField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='Nt_sequences', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='Junction',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('J_GENE_and_allele', models.TextField()),
                ('D_GENE_and_allele', models.TextField()),
                ('JUNCTION_frame', models.TextField()),
                ('JUNCTION', models.TextField()),
                ('JUNCTION_with_frameshift', models.TextField()),
                ('threeprime_V_REGION', models.TextField()),
                ('Pthreeprime_V', models.TextField()),
                ('N_REGION', models.TextField()),
                ('N1_REGION', models.TextField()),
                ('Pfiveprime_D', models.TextField()),
                ('D_REGION', models.TextField()),
                ('Pthreeprime_D', models.TextField()),
                ('Pfiveprime_D1', models.TextField()),
                ('D1_REGION', models.TextField()),
                ('Pthreeprime_D1', models.TextField()),
                ('N2_REGION', models.TextField()),
                ('Pfiveprime_D2', models.TextField()),
                ('D2_REGION', models.TextField()),
                ('Pthreeprime_D2', models.TextField()),
                ('N3_REGION', models.TextField()),
                ('Pfiveprime_D3', models.TextField()),
                ('D3_REGION', models.TextField()),
                ('Pthreeprime_D3', models.TextField()),
                ('N4_REGION', models.TextField()),
                ('Pfiveprime_J', models.TextField()),
                ('fiveprime_J_REGION', models.TextField()),
                ('JUNCTION_nt_nb', models.IntegerField()),
                ('threeprime_V_REGION_nt_nb', models.IntegerField()),
                ('Pthreeprime_V_nt_nb', models.IntegerField()),
                ('N_REGION_nt_nb', models.IntegerField()),
                ('N1_REGION_nt_nb', models.IntegerField()),
                ('Pfiveprime_D_nt_nb', models.IntegerField()),
                ('D_REGION_nt_nb', models.IntegerField()),
                ('Pthreeprime_D_nt_nb', models.IntegerField()),
                ('Pfiveprime_D1_nt_nb', models.IntegerField()),
                ('D1_REGION_nt_nb', models.IntegerField()),
                ('Pthreeprime_D1_nt_nb', models.IntegerField()),
                ('N2_REGION_nt_nb', models.IntegerField()),
                ('Pfiveprime_D2_nt_nb', models.IntegerField()),
                ('D2_REGION_nt_nb', models.IntegerField()),
                ('Pthreeprime_D2_nt_nb', models.IntegerField()),
                ('N3_REGION_nt_nb', models.IntegerField()),
                ('Pfiveprime_D3_nt_nb', models.IntegerField()),
                ('D3_REGION_nt_nb', models.IntegerField()),
                ('Pthreeprime_D3_nt_nb', models.IntegerField()),
                ('N4_REGION_nt_nb', models.IntegerField()),
                ('Pfiveprime_J_nt_nb', models.IntegerField()),
                ('fiveprime_J_REGION_nt_nb', models.IntegerField()),
                ('threeprime_V_REGION_trimmed_nt_nb', models.IntegerField()),
                ('fiveprime_D_REGION_trimmed_nt_nb', models.IntegerField()),
                ('threeprime_D_REGION_trimmed_nt_nb', models.IntegerField()),
                ('fiveprime_D1_REGION_trimmed_nt_nb', models.IntegerField()),
                ('threeprime_D1_REGION_trimmed_nt_nb', models.IntegerField()),
                ('fiveprime_D2_REGION_trimmed_nt_nb', models.IntegerField()),
                ('threeprime_D2_REGION_trimmed_nt_nb', models.IntegerField()),
                ('fiveprime_D3_REGION_trimmed_nt_nb', models.IntegerField()),
                ('threeprime_D3_REGION_trimmed_nt_nb', models.IntegerField()),
                ('fiveprime_J_REGION_trimmed_nt_nb', models.IntegerField()),
                ('threeprime_V_REGION_mut_nt_nb', models.IntegerField()),
                ('D_REGION_mut_nt_nb', models.IntegerField()),
                ('D1_REGION_mut_nt_nb', models.IntegerField()),
                ('D2_REGION_mut_nt_nb', models.IntegerField()),
                ('D3_REGION_mut_nt_nb', models.IntegerField()),
                ('fiveprime_J_REGION_mut_nt_nb', models.IntegerField()),
                ('D_REGION_reading_frame', models.IntegerField()),
                ('Ngc', models.TextField()),
                ('CDR3_IMGT_length', models.IntegerField()),
                ('Molecular_mass', models.FloatField()),
                ('pI', models.FloatField()),
                ('threeprime_V_REGION_accepted_mut_nb', models.IntegerField()),
                ('D_REGION_accepted_mut_nb', models.IntegerField()),
                ('fiveprime_J_REGION_accepted_mut_nb', models.IntegerField()),
                ('Accepted_D_GENE_nb', models.IntegerField()),
                ('CDR3_IMGT', models.TextField()),
                ('CDR3_IMGT_nt_nb', models.IntegerField()),
                ('CDR3_IMGT_with_frameshift', models.TextField()),
                ('CDR3_IMGT_AA', models.TextField()),
                ('CDR3_IMGT_AA_with_frameshift', models.TextField()),
                ('JUNCTION_AA', models.TextField()),
                ('JUNCTION_AA_with_frameshift', models.TextField()),
                ('JUNCTION_decryption', models.TextField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='Junction', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='IMGT_gapped_nt_sequences',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('J_GENE_and_allele', models.TextField()),
                ('D_GENE_and_allele', models.TextField()),
                ('V_D_J_REGION', models.TextField()),
                ('V_J_REGION', models.TextField()),
                ('V_REGION', models.TextField()),
                ('FR1_IMGT', models.TextField()),
                ('CDR1_IMGT', models.TextField()),
                ('FR2_IMGT', models.TextField()),
                ('CDR2_IMGT', models.TextField()),
                ('FR3_IMGT', models.TextField()),
                ('CDR3_IMGT', models.TextField()),
                ('JUNCTION', models.TextField()),
                ('J_REGION', models.TextField()),
                ('FR4_IMGT', models.TextField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='IMGT_gapped_nt_sequences', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='IMGT_gapped_AA_sequences',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('J_GENE_and_allele', models.TextField()),
                ('D_GENE_and_allele', models.TextField()),
                ('V_D_J_REGION', models.TextField()),
                ('V_J_REGION', models.TextField()),
                ('V_REGION', models.TextField()),
                ('FR1_IMGT', models.TextField()),
                ('CDR1_IMGT', models.TextField()),
                ('FR2_IMGT', models.TextField()),
                ('CDR2_IMGT', models.TextField()),
                ('FR3_IMGT', models.TextField()),
                ('CDR3_IMGT', models.TextField()),
                ('JUNCTION', models.TextField()),
                ('J_REGION', models.TextField()),
                ('FR4_IMGT', models.TextField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='IMGT_gapped_AA_sequences', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.CreateModel(
            name='AA_sequences',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('J_GENE_and_allele', models.TextField()),
                ('D_GENE_and_allele', models.TextField()),
                ('V_D_J_REGION', models.TextField()),
                ('V_J_REGION', models.TextField()),
                ('V_REGION', models.TextField()),
                ('FR1_IMGT', models.TextField()),
                ('CDR1_IMGT', models.TextField()),
                ('FR2_IMGT', models.TextField()),
                ('CDR2_IMGT', models.TextField()),
                ('FR3_IMGT', models.TextField()),
                ('CDR3_IMGT', models.TextField()),
                ('JUNCTION', models.TextField()),
                ('J_REGION', models.TextField()),
                ('FR4_IMGT', models.TextField()),
                ('Sequence_Identifier', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='AA_sequences', to='PinkStrawberry.sequence')),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
    ]
