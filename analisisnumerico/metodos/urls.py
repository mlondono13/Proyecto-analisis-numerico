# mi_aplicacion/urls.py

from django.urls import path
from . import views

urlpatterns = [
    path('', views.inicio, name='inicio'),  # Página de inicio
    path('capitulo_1/', views.calcular_cap1, name='capitulo_1'),
    path('capitulo_2/', views.calcular_cap2, name='capitulo_2'),
    path('capitulo_3/', views.calcular_cap3, name='capitulo_3'),
]
