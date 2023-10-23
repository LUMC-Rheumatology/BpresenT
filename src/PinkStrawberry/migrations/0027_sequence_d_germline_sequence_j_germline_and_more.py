# Generated by Django 4.1.3 on 2023-03-11 13:25

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0026_alter_patient_options_remove_patient_created_at'),
    ]

    operations = [
        migrations.AddField(
            model_name='sequence',
            name='D_Germline',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='SequenceD', to='PinkStrawberry.germline'),
        ),
        migrations.AddField(
            model_name='sequence',
            name='J_Germline',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='SequenceJ', to='PinkStrawberry.germline'),
        ),
        migrations.AddField(
            model_name='sequence',
            name='V_Germline',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='SequenceV', to='PinkStrawberry.germline'),
        ),
    ]