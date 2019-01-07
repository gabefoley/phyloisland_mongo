from flask import Flask
from flask_mongoengine import MongoEngine
from flask_uploads import UploadSet, configure_uploads, ALL
from flask_user import login_required, UserManager
from flask_bootstrap import Bootstrap
from mongoengine import connect
import configs.mongoconfig
import users

app = Flask(__name__)
Bootstrap(app)
app.config.from_pyfile('configs/mongoconfig.py')

# Connect the database
db = MongoEngine(app)
connect(configs.mongoconfig.MONGODB_DB)



# Setup the Login Manager
# lm = LoginManager()
# lm.init_app(app)
# lm.login_view = 'login'

# Add the Upload directory
allfiles = UploadSet('all', ALL)
configure_uploads(app, allfiles)

# Import views down here because we need to have already initialised
from views import *

if __name__ == "__main__":
    app.run(debug=True)
