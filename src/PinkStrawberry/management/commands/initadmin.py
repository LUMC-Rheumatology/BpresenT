from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from django.db import IntegrityError

class Command(BaseCommand):
    help = "Create an account, admin:admin, for the admin panel"

    def handle(self, *args, **options):
        try:
            superuser = User.objects.create_superuser(
                username="admin",
                email="noreply@lumc.nl",
                password="admin")
            superuser.save()
            self.stdout.write(self.style.SUCCESS("succesfully created admin:admin account"))
        except IntegrityError:
            self.stdout.write(self.style.ERROR("superuser admin:admin already exists"))
