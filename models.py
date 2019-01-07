from phyloisland import app, db
from flask_user import UserMixin, UserManager

from flask_mongoengine import *
from flask_mongoengine.wtf import model_form


# Create models

class User(db.Document, UserMixin):
    active = db.BooleanField(default=True)

    # User authentication information
    username = db.StringField(default='')
    email = db.StringField(max_length=30)
    password = db.StringField()
    email_confirmed_at = db.DateTimeField()

    # User information
    first_name = db.StringField(default='')
    last_name = db.StringField(default='')

    # Customisable information
    page_size = db.IntField()


    # Relationships
    roles = db.ListField(db.StringField(), default=[])

# Setup Flask-User and specify the User data-model
user_manager = UserManager(app, db, User)