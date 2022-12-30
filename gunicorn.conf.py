workers = 2
threads = 2
bind = '0.0.0.0:5000'
accesslog = '/var/log/gunicorn_access.log'
errorlog = '/var/log/gunicorn_error.log'
loglevel = 'warning'
