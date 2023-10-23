# Generated by Django 4.1.3 on 2023-05-01 21:33

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0048_prediction'),
    ]

    operations = [
        migrations.CreateModel(
            name='Peptide_comparison',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('peptide', models.TextField(null=True)),
                ('allele', models.CharField(max_length=20, null=True)),
                ('score_diff', models.FloatField(null=True)),
                ('percentile_rank_diff', models.FloatField(null=True)),
                ('affinity_diff', models.FloatField(null=True)),
                ('Germline_peptide', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='GLcomparison', to='PinkStrawberry.prediction')),
                ('Target_peptide', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='comparison', to='PinkStrawberry.prediction')),
                ('sequence', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='Peptide_comparison', to='PinkStrawberry.sequence')),
            ],
        ),
    ]
