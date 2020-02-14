# Names to use for our parameters
# MONGO_DBNAME = "phyloisland_mongo_2019"
# MONGO_URI = 'mongodb://pi:abc123@localhost/' + MONGO_DBNAME
#
# MONGODB_DB = 'phyloisland_mongo_2019_dream_again'
# MONGODB_DB = 'phyloisland_mongo_2019_dream_again_no_delete'
# MONGODB_DB = 'check'

MONGODB_DB = 'raptor'

# MONGODB_HOST= 'localhost'
# MONGODB_USERNAME = 'pi'
# MONGODB_PASSWORD = 'abc123'


# User / login configuration
USER_APP_NAME = "Phylo Island" # Shown in email templates and page footers
CSRF_ENABLED = False
USER_ENABLE_EMAIL = False # Require email authentication?
USER_ENABLE_USERNAME = True  # Enable username authentication
USER_REQUIRE_RETYPE_PASSWORD = True  # Require retyping of password
USER_ENABLE_CHANGE_PASSWORD = True
USER_AFTER_REGISTER_ENDPOINT = 'user.login'
USER_EMAIL_SENDER_EMAIL = 'phyloisland@gmail.com'
USER_EMAIL_SENDER_NAME = 'Phylo Island'
UPLOADS_ALL_DEST = 'static/uploads'
UPLOADED_ALL_DEST = 'static/uploads'

# Mail server configuration
MAIL_SERVER = 'smtp.gmail.com'
MAIL_PORT = 587
MAIL_USE_TLS = True
MAIL_USE_SSL = True
MAIL_USERNAME = 'phyloisland@gmail.com'
MAIL_PASSWORD = 'phyloislandpass'


SECRET_KEY = 'developmentkey2020askyourselfimportantquestions'

