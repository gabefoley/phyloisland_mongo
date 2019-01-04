# Names to use for our parameters
MONGO_DBNAME = "phyloisland_mongo_2019"

DEBUG = True
# MONGO_URI = 'mysql://pi:@localhost/' + bio_server_name

MONGO_URI = 'mongodb://pi:abc123@localhost/' + MONGO_DBNAME
UPLOADS_ALL_DEST = 'static/uploads'
UPLOADED_ALL_DEST = 'static/uploads'

SECRET_KEY = 'developmentkey'

