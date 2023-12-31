# Generated by Django 4.1.3 on 2023-05-01 15:09

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0047_rename_mutated_aa_change_changed'),
    ]

    operations = [
        migrations.CreateModel(
            name='Prediction',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('peptide', models.TextField()),
                ('allele', models.CharField(max_length=20)),
                ('score', models.FloatField()),
                ('percentile_rank', models.FloatField()),
                ('affinity', models.FloatField()),
                ('offset', models.IntegerField()),
                ('sequence', models.TextField()),
            ],
        ),
    ]
