from django.core.management.base import BaseCommand
from PinkStrawberry import IMGT_GENE_DB_Loader as loader

class Command(BaseCommand):
    help = "Update the BpresenT database with the newest version of IMGT/Gene-DB"

    def handle(self, *args, **options):
        loader.update()
        self.stdout.write(self.style.SUCCESS("succesfully updated IMGT/Gene-DB"))
