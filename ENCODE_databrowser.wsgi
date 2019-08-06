#!/var/www/encode-ui-sandbox/encode-ui-sandbox/bin/python
activate_this = '/var/www/encode-ui-sandbox/encode-ui-sandbox/bin/activate_this.py'
with open(activate_this) as file_:
    exec(file_.read(), dict(__file__=activate_this))
import sys
print(sys.path)
sys.path.insert(0,"/var/www/encode-ui-sandbox")
import os
cwd = os.getcwd()
print(cwd)
from app import server as application
