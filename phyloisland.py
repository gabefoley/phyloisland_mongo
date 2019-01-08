from flask import Flask
from flask_mongoengine import MongoEngine
from flask_uploads import UploadSet, configure_uploads, ALL
from flask_bootstrap import Bootstrap
from mongoengine import connect
import configs.mongoconfig

app = Flask(__name__)
Bootstrap(app)
app.config.from_pyfile('configs/mongoconfig.py')

# Connect the database
db = MongoEngine(app)
connect(configs.mongoconfig.MONGODB_DB)

# Add the Upload directory
allfiles = UploadSet('all', ALL)
configure_uploads(app, allfiles)

# Import views down here because we need to have already initialised
from views import *

if __name__ == "__main__":
    app.run(debug=True)
