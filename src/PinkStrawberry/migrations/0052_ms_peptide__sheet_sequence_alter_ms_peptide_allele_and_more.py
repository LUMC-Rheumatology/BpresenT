# Generated by Django 4.1.3 on 2023-08-04 14:32

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0051_alter_ms_peptide_confidence'),
    ]

    operations = [
        migrations.AddField(
            model_name='ms_peptide',
            name='_sheet_sequence',
            field=models.TextField(default='TEST'),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='ms_peptide',
            name='allele',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='ms_peptide',
            name='confidence',
            field=models.IntegerField(choices=[(0, 'Low'), (1, 'Medium'), (2, 'High')]),
        ),
        migrations.AlterField(
            model_name='ms_peptide',
            name='sequence',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='Ms_peptide', to='PinkStrawberry.sequence'),
        ),
    ]
