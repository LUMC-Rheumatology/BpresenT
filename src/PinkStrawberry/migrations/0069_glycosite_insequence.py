# Generated by Django 4.1.3 on 2023-08-24 14:14

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0068_ms_peptide_sequence_start'),
    ]

    operations = [
        migrations.AddField(
            model_name='glycosite',
            name='insequence',
            field=models.BooleanField(default=1),
        ),
    ]
