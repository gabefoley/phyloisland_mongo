# # class SplodgeView(ModelView):
# #     column_list = ('name', 'splodgeword')
# #
# #     form = forms.SplodgeForm
#
# global passwordtochange
# global hopscotch
# passwordtochange = ['name', 'password']
# @app.route('/namechange/')
# def namechange():
#     passwordtochange = ['password']
#     hopscotch = False
#     return 'Changed?'
#
#
# @app.route('/namechangeback/')
# def namechangeback():
#     passwordtochange = ['name, password']
#     hopscotch = True
#     print
#
#     return 'Changed back?'

# class BlockedView(BaseView):
#     @expose("/", methods=('GET', 'POST'))
#
#     def is_accessible(self):
#         return not current_user.is_authenticated
#
#     def inaccessible_callback(self, name, **kwargs):
#         return render_template('login.html', title='login', form=form)
#
# hopscotch = True
#
# # class UserView(ModelView):
# #     @property
# #     def column_list(self):
# #         pass
# #     @property
# #     def _list_columns(self):
# #         return self.get_list_columns()
# #
# #     @_list_columns.setter
# #     def _list_columns(self, value):
# #         pass
# #
# #     def user(self):
# #         print('And now the password is ' + passwordtochange)
# #
# #
# #         return self.render('user.html')
# #
# #     # _handle_view called every request
# #     def _handle_view(self, name, **kwargs):
# #         print ('in the view')
# #         hopscotch = random.choice([True, False])
# #         print ('hopscotch is ' + str(hopscotch))
# #         if not hopscotch:
# #             print ('got to the false claim')
# #             # self._list_columns = ['name', 'password']
# #
# #             self._refresh_cache()
# #             # return super(UserView, self)._handle_view(name, **kwargs)
# #
# #         # re-scaffold views every request
# #         self._refresh_cache()
# #
# #         return super(UserView, self)._handle_view(name, **kwargs)
# #
# #     # _refresh_cache called once when view is added to admin interface
# #     def _refresh_cache(self):
# #         # do not _refresh_cache outside of a request context
# #         # if not hopscotch:
# #         #     # init members with empty tuples to avoid instantiation error
# #         #     self._list_columns = ()
# #         #     return
# #         self._list_columns = ['name', 'password']
# #
# #         super(UserView, self)._refresh_cache()
