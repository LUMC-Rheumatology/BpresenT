# Generated by Django 4.1.3 on 2023-05-01 14:43

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0046_rename_nt_mutated_aa_change_mutated_and_more'),
    ]

    operations = [
        migrations.RenameField(
            model_name='aa_change',
            old_name='mutated',
            new_name='changed',
        ),
    ]
