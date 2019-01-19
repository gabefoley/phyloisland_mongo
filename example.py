class MyModelView(ModelView):

    @login_required
    @expose("/", methods=('GET', 'POST'))
    def index_view(self):
        self.form_columns = ('name', 'profile')

        return super(ModelView, self).index_view()


    def is_accessible(self):
        if (not current_user.is_active or not
        current_user.is_authenticated):
            return False
        return True

    @property
    def _list_columns(self):
        return self.get_list_columns()

    @_list_columns.setter
    def _list_columns(self, value):
        pass

    @action('generate_profile', 'Generate a profile from these sequences')
    def action_generate_profile(self, ids):
        pass


    @property
    def _action_form_class(self):
        return self.get_action_form()

    @_action_form_class.setter
    def _action_form_class(self, value):
        pass

    form_columns = ('name')