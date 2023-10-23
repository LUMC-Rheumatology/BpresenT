"""Bafstu URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from PinkStrawberry import urls
from django.conf.urls.static import static
from django.conf import settings
import os

#Show call logo
if True:
    _cyan, _end = '\033[96m', '\033[0m'
    print(f"starting server from {settings.BASE_DIR}")
    print(f"""
                 )\._.,--....,'``.
   .{_cyan}B{_end}--.        /;   _{_cyan}presenT{_end}_\  (`._ ,.
  `=,-,-'~~~   `----(,_..'--(,_..'`-.;.'
    """)
    del _cyan, _end


urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('PinkStrawberry.urls')),
] + static(settings.STATIC_URL)
