# Generated by Django 4.1.3 on 2023-08-01 14:37

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0050_ms_peptide'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ms_peptide',
            name='confidence',
            field=models.IntegerField(choices=[(0, 'LOW'), (1, 'MEDIUM'), (2, 'HIGH')]),
        ),
    ]