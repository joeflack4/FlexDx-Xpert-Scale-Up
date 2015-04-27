from django.conf.urls import url
from django.conf import settings
from django.conf.urls.static import static

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

from async import views as async_views

urlpatterns = [
    # Examples:
    # url(r'^$', 'xpert_ru.views.home', name='home'),
    # url(r'^xpert_ru/', include('flexdxtb.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
    url(r'^$', async_views.home_page),
    url(r'^index.html$', async_views.home_page),
    url(r'^about/$', async_views.about_page),
    url(r'^model/$', async_views.model_page),
    url(r'^help/$', async_views.help_page),
    url(r'^homebrew/model/$', async_views.homebrew_model_page),
    url(r'^bargraph.png$', async_views.bargraph),
    url(r'^bargraph_uncert.png$', async_views.bargraph_uncert),
    url(r'^dgraph1.png$', async_views.dgraph1),
    url(r'^dgraph2.png$', async_views.dgraph2),
    url(r'^dgraph3.png$', async_views.dgraph3),
    url(r'^dgraph4.png$', async_views.dgraph4),
    url(r'^excel/$', async_views.excel_sheet),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
  #Remove for production
