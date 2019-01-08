from phyloisland import app, db
from flask_user import UserMixin, UserManager


class User(db.DynamicDocument, UserMixin):
    """
    User model
    """
    active = db.BooleanField(default=True)

    # User authentication information
    username = db.StringField(default='', unique=True)
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


class SequenceRecords(db.DynamicDocument):
    """
    Class for storing Sequence records
    """
    name = db.StringField()
    species = db.StringField()
    description = db.StringField()
    sequence = db.StringField()
    references = db.ListField(db.StringField(), default=list)


class GenomeRecords(db.DynamicDocument):
    """
    Class for storing Genome records
    """
    name = db.StringField()
    species = db.StringField()
    strain = db.StringField()
    description = db.StringField()
    sequence = db.StringField()
    present = db.DictField()
    references = db.ListField(db.StringField(), default=list)
