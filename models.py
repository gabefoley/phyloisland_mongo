from flask_mongoengine import *
from flask_mongoengine.wtf import model_form


# Create models

class User(DynamicDocument):
    username = mongoengine.StringField()
    password = mongoengine.StringField()

