# Generated by Django 4.1.3 on 2023-01-31 10:21

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0008_sequences_submitted_sequence_translated'),
    ]

    operations = [
        migrations.AlterField(
            model_name='sequences',
            name='IMGT_sequence',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='sequences',
            name='submitted_sequence',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='sequences',
            name='submitted_sequence_translated',
            field=models.TextField(blank=True, null=True),
        ),
    ]