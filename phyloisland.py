from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_uploads import UploadSet, configure_uploads, ALL
from flask_bootstrap import Bootstrap

app = Flask(__name__)
Bootstrap(app)
app.config.from_pyfile('config.py')

# Connect the database
db = SQLAlchemy(app)

# Add the Upload directory
allfiles = UploadSet('all', ALL)
configure_uploads(app, allfiles)

# Import views down here because we need to have already initialised
from views import *

if __name__ == "__main__":
    app.run(debug=True)
