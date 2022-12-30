from flask import Flask
import os

app = Flask('cnv_ana')
app.config.from_pyfile('settings.py')
app.config['WTF_I18N_ENABLED'] = False
app.config['UPLOAD_PATH'] = os.path.join(app.root_path, 'uploads')

app.config['RAW_PATH'] = os.path.join('/data/', 'rawdata')
app.config['OUT_PATH'] = os.path.join('/data/', 'output')
app.config['DB_PATH'] = os.path.join('/data/', 'db','hg19')

app.config['BIN_PATH'] = os.path.join('/cnv_ana/', 'bin', 'CNVana_batch.py')
app.config['BASE_PATH'] = os.path.join('/cnv_ana/', 'bin','baseline')



from cnv_ana import views