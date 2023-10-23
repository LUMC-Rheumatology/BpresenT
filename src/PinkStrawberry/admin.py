from django.contrib import admin
from django.apps import apps

# Customized admin site
class LocalAdminSite(admin.AdminSite):
    site_header = 'BpresenT Database Manager'
    final_catch_all_view = False
    login_template = "templates/admin/login.html"



# Register your models here.
class ShowAllAttr(admin.ModelAdmin):

    def __init__(self, model, admin_site):
        self.list_display = [field.name for field in model._meta.fields if field.name != "id"]
        super(ShowAllAttr, self).__init__(model, admin_site)

admin_site = LocalAdminSite()
for model in list(apps.get_app_config("PinkStrawberry").get_models()):
    admin.site.register(model, ShowAllAttr)

