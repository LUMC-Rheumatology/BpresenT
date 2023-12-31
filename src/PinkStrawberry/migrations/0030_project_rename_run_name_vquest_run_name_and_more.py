# Generated by Django 4.1.3 on 2023-03-17 22:13

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0029_alter_aa_sequences_sequence_identifier_and_more'),
    ]

    operations = [
        migrations.CreateModel(
            name='Project',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.TextField()),
                ('description', models.TextField(null=True)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
            ],
            options={
                'get_latest_by': 'id',
            },
        ),
        migrations.RenameField(
            model_name='vquest_run',
            old_name='run_name',
            new_name='name',
        ),
        migrations.AlterField(
            model_name='bcr',
            name='Run',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='BCR', to='PinkStrawberry.vquest_run'),
        ),
        migrations.AlterField(
            model_name='sequence',
            name='Run',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='Sequence', to='PinkStrawberry.vquest_run'),
        ),
        migrations.DeleteModel(
            name='Run_Identifier',
        ),
        migrations.AddField(
            model_name='vquest_run',
            name='Project',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='VQuest_Run', to='PinkStrawberry.project'),
        ),
    ]
