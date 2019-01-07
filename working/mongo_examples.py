#
#
# @app.route('/add')
# def add():
#     user = mongo.db.users
#     user.insert({'name' : 'Gabe', 'mood' : 'happy'})
#     user.insert({'name' : 'Amanda', 'mood' : 'sad'})
#     user.insert({'name' : 'Fred', 'mood' : 'ecstatic'})
#     return 'Updated user'
#
# @app.route('/find')
# def find():
#     user = mongo.db.users
#     result = user.find_one({'name' : 'Amanda'})
#     return 'You found ' + result['name'] + 'and she is ' +  result['mood']
#
# @app.route('/update')
# def update():
#     user = mongo.db.users
#     result = user.find_one({'name' : 'Amanda'})
#     result['mood'] = 'crazy'
#     user.save(result)
#
#     return 'Updated'
