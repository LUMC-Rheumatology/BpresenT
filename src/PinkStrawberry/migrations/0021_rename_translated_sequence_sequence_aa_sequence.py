# Generated by Django 4.1.3 on 2023-02-27 10:10

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0020_sequence_patient'),
    ]

    operations = [
        migrations.RenameField(
            model_name='sequence',
            old_name='Translated_Sequence',
            new_name='AA_Sequence',
        ),
    ]
