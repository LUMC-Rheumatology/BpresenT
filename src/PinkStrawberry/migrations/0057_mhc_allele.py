# Generated by Django 4.1.3 on 2023-08-08 13:56

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0056_ms_peptide_possible_allele'),
    ]

    operations = [
        migrations.CreateModel(
            name='MHC_allele',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.TextField()),
            ],
        ),
    ]