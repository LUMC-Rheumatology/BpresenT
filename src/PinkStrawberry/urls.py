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
from django.urls import path, re_path
from PinkStrawberry.admin import admin_site
from rest_framework.urlpatterns import format_suffix_patterns
from .views import *


#temp rest patterns for demo
restpatterns = [
    path('api', sequence_list.as_view()),
    path('vquest', ExampleViewSet.as_view({'get': 'list'}))
]


frontpatterns = [
    # "frontend' stuff
    path('base', base, name="base"),
    path('', projects),
    path('projects', projects, name="projects"),
    path('project/<int:pk>', project_detail, name="project/pk"),
    path('sequences', sequences, name="sequences"),
    path('sequence/<int:id>', sequence_detail, name="sequence/id"),
    path('patient/<str:id>', patient_detail, name="patient/id"),
    path('germline/<str:gene_and_allele>', germline_detail, name="germline/gene_and_allele"),
    path('visualizer/<int:id>', visualizer, name="visualizer/id"),
    path('peptides/netmhcpan', netmhcpan, name="peptides/netmhcpan"),
    path('peptides/netmhcpan/allelelist.txt', allelelist, name="peptides/netmhcpan/allelelist"),
    path('peptides/netmhcpan', netmhcpan, name="peptides/netmhcpan"),
    path('peptides/msoverview', ms_overview, name='peptides/msoverview'),
    path('peptides/msform/<int:id>', ms_form, name="peptides/msform/id"),
    path('peptides/overview', peptides_overview, name="peptides/overview"),
    path('settings/regular', settings_regular, name="settings/regular"),
]

urlpatterns = frontpatterns + format_suffix_patterns(restpatterns)
