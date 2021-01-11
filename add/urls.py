"""add URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.2/topics/http/urls/
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

    path('',include('travel.urls'))
"""
from datetime import datetime
from django.urls import path
from django.contrib import admin
from django.contrib.auth.views import LoginView, LogoutView
from app import forms, views



# does this gets uploaded?

urlpatterns = [
    path('', views.home, name='home'),
    path('mandq/', views.mandq, name='mandq'),
    path('randp/', views.randp, name='randp'),
    path('survival/', views.survival, name='survival'),
    path('phfortests/', views.phfortests, name='phfortests'),
    path('contact/', views.contact, name='contact'),
    path('mandq/add1/', views.add1, name='add1'),
    path('randp/add2/', views.add2, name='add2'),
    path('survival/add3/', views.add3, name='add3'),
    path('mandq/add1/anova2/', views.anova2, name = 'twowayanova'),
    path('mandq/add1/ttest/', views.ttest, name = 'ttest'),
    path('mandq/add1/anova/', views.anova, name = 'anova'),
    path('randp/add2/chisq/', views.chisq, name = 'fishers'),
    path('randp/add2/multicat/', views.chisq, name = 'multicat'),
    path('randp/add2/chisqwithsex/', views.chisq, name = 'chisqwithsex'),
    path('survival/add3/logrank/', views.logrank, name = 'logrank'),
    path('survival/add3/logrankandsex/', views.logrank, name = 'logrankandsex'),
    path('login/',
         LoginView.as_view
         (
             template_name='app/login.html',
             authentication_form=forms.BootstrapAuthenticationForm,
             extra_context=
             {
                 'title': 'Log in',
                 'year' : datetime.now().year,
             }
         ),
         name='login'),
    path('logout/', LogoutView.as_view(next_page='/'), name='logout'),
    path('admin/', admin.site.urls),
]
