# Generated by Django 4.1.3 on 2023-03-17 22:48

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('PinkStrawberry', '0031_vquest_run_comment'),
    ]

    operations = [
        migrations.AlterField(
            model_name='vquest_run',
            name='Project',
            field=models.ForeignKey(default=None, on_delete=django.db.models.deletion.CASCADE, related_name='VQuest_Run', to='PinkStrawberry.project'),
            preserve_default=False,
        ),
    ]
