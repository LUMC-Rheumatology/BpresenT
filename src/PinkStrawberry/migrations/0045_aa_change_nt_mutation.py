# Generated by Django 4.1.3 on 2023-05-01 13:50

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0044_alter_imgt_aa_change_sequence_and_more'),
    ]

    operations = [
        migrations.CreateModel(
            name='AA_change',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('nt_original', models.CharField(max_length=1)),
                ('nt_pos', models.IntegerField()),
                ('nt_mutated', models.CharField(max_length=1)),
                ('sequence', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='AA_changes', to='PinkStrawberry.sequence')),
            ],
        ),
        migrations.CreateModel(
            name='Nt_mutation',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('original', models.CharField(max_length=1)),
                ('pos', models.IntegerField()),
                ('mutated', models.CharField(max_length=1)),
                ('aa_change', models.ForeignKey(null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='mutations', to='PinkStrawberry.aa_change')),
                ('sequence', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='Nt_mutations', to='PinkStrawberry.sequence')),
            ],
        ),
    ]
