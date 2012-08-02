"""
WSGI config for gtc project.

This module contains the WSGI application used by Django's development server
and any production WSGI deployments. It should expose a module-level variable
named ``application``. Django's ``runserver`` and ``runfcgi`` commands discover
this application via the ``WSGI_APPLICATION`` setting.

Usually you will have the standard Django WSGI application here, but it also
might make sense to replace the whole Django WSGI application with a custom one
that later delegates to the Django one. For example, you could introduce WSGI
middleware here, or combine a Django application with an application of another
framework.

"""
import django
import os

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "gtc.settings")

# This application object is used by any WSGI server configured to use this
# file. This includes Django's development server, if the WSGI_APPLICATION
# setting points here.
#from django.core.wsgi import get_wsgi_application
#application = get_wsgi_application()
import django.core.handlers.wsgi
_application = django.core.handlers.wsgi.WSGIHandler()

def application(environ, start_response):
    if 'HTTP_X_FORWARDED_PROTOCOL' in environ:
        environ['wsgi.url_scheme'] = environ.get('HTTP_X_FORWARDED_PROTOCOL', 'http')
    if 'HTTP_X_FORWARDED_HOST' in environ:
        environ['HTTP_X_FORWARDED_HOST'] = environ['HTTP_X_FORWARDED_HOST'].split(',')[0].strip()
    return _application(environ, start_response)
# Apply WSGI middleware here.
# from helloworld.wsgi import HelloWorldApplication
# application = HelloWorldApplication(application)
