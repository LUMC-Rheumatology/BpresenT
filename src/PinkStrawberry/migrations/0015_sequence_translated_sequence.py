# Generated by Django 4.1.3 on 2023-02-20 11:04

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0014_alter_nt_sequences_pfiveprime_j_end_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='sequence',
            name='Translated_Sequence',
            field=models.TextField(null=True),
        ),
    ]