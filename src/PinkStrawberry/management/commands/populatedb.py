from django.core.management.base import BaseCommand
from django.db import IntegrityError, transaction
from django.conf import settings
from pathlib import Path
from PinkStrawberry.models import MHC_allele

class Command(BaseCommand):
    help = "Populate BpresenT database with required information"

    def handle(self, *args, **options):
        textvec = open(Path.joinpath(
            settings.BASE_DIR, "PinkStrawberry/static/netmhcpan/services.healthtech.dtu.dk_services_NetMHCpan-4.1_MHC_allele_names.txt"))\
            .readlines()
        with transaction.atomic():
            for line in textvec:
                MHC_allele.objects.get_or_create(name=line.strip())

        self.stdout.write(self.style.SUCCESS("succesfully populated database"))


