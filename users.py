from phyloisland import app, db
from flask_user import UserMixin, UserManager


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

    # Relationships
    roles = db.ListField(db.StringField(), default=[])


# Setup Flask-User and specify the User data-model
user_manager = UserManager(app, db, User)