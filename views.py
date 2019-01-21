from phyloisland import app, allfiles
import models
import forms
import utilities
import getGenomes
import custom_filters
from flask import render_template, flash, request, session, send_file, has_app_context
from flask_login import login_required, current_user
from flask_admin import Admin, expose, BaseView
from flask_admin.actions import action
from flask_admin.contrib.mongoengine import ModelView, filters
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from markupsafe import Markup
import gettext
import timeit
import time
import os
from urllib.error import HTTPError
import random
import io

from flask_admin.model.template import EndpointLinkRowAction, LinkRowAction

from wtforms import SelectField

ref_names = ['A1', 'A2', 'Chitinase', 'TcB', 'TcC']


class UploadView(BaseView):

    """
    View for uploading files
    """

    @login_required
    @expose("/", methods=('GET', 'POST'))
    def upload(self):
        form = forms.UploadForm()
        print ('is it validated - ', form.validate())

        if request.method == "POST" and not form.validate():
            print ('not validated')
            print (form.genome_type)
            print (request.get_json())


        if request.method == 'POST' and form.validate():

            if request.method == 'POST':
                names = request.get_json()
                print ('names coming in')
                print (request)
                print (request.data)
                print (names)
                # for name in names:
                #     print (name)

            start_time = timeit.default_timer()

            # Get the information from the upload form
            filename = allfiles.save(request.files['file'])
            seq_type = form.type.data
            add_seq = form.add_sequence.data
            add_genome = form.add_genome.data
            single = form.single_genome.data
            genome_type=form.genome_type.data
            representative = form.representative.data
            assembly = form.assembly.data
            genbank = form.genbank.data
            failed_genomes = []

            print ('seq type is ', seq_type)

            try:

                if seq_type == "protein" or seq_type == "nucleotide":
                    print('here')
                    seq_records = utilities.read_fasta("static/uploads/" + filename)

                    if not seq_records:
                        print("Couldn't find any sequences in the uploaded file.")
                        flash("Couldn't find any sequences in the uploaded file.", category='error')

                    else:
                        if add_seq:
                            utilities.addSequence(seq_records)

                        if add_genome:
                            # TODO: Correct this call
                            getGenomes.add_genome(seq_records, seq_type)

                elif seq_type == "species":

                    species_names = utilities.readLinesFromFile("static/uploads/" + filename)
                    for species_name in species_names:

                        print("Species name is ", species_name)

                        destinations = [genome_type]

                        #TODO: Once the checkbox selection is dynamic we can just add freely here
                        if representative and genome_type not in ["representative genome", "assembly", "genbank"]:
                            destinations.append("representative genome")

                        if assembly and genome_type not in ["assembly", "genbank"]:
                            destinations.append("assembly")

                        print("Destinations is ", destinations)

                        genome_results = getGenomes.add_genome(species_name, destinations, single=single)

                        if genome_results and genome_results != "Fail":
                            utilities.add_genome(genome_results)
                        else:

                            # Couldn't find it in RefSeq, let's try genbank
                            if genbank:
                                destinations = ["genbank"]
                                genome_results = getGenomes.add_genome(species_name, destinations, single=single)

                                if genome_results and genome_results != "Fail":

                                    utilities.add_genome(genome_results)
                                else:
                                    failed_genomes.append(species_name)

                            else:
                                failed_genomes.append(species_name)

                elif seq_type == "genome":
                    pass
                    # genome_names = readLinesFromFile("static/uploads/" + filename)
                    #
                    # for name in genome_names:
                    #     genome_ids = {}
                    #     genome_results = {}
                    #     genome_ids = readLinesFromFile("static/uploads/" + filename)
                    #     genome_results = mapToGenome.get_full_genome([name])
                    #     if genome_results and genome_results != 'in_database':
                    #         addGenome(genome_results)
                    #     elif genome_results != 'in_database' and search_shotgun:
                    #         print("\nWe didn't identify any genome records for %s. Attempting to search for shotgun "
                    #               "sequenced genomes \n" % (
                    #                 name))
                    #
                    #         shotgun_id_dict = mapToGenome.get_shotgun_id_dict_from_id(name)
                    #         if shotgun_id_dict:
                    #             # print("Skipping getting shotgun")
                    #             # shotgun_results = None
                    #             genome_results = mapToGenome.get_shotgun_genome(shotgun_id_dict)
                    #             if genome_results:
                    #                 addGenome(genome_results)
                    #
                    #         else:
                    #             print("Couldn't find a shotgun sequence for %s \n" % (name))
                    #
                    #
                    #     elif genome_results != 'in_database' and not search_shotgun:
                    #         print(
                    #             "\nWe didn't identify any genome records for %s. And we are not attempting to search
                    # for shotgun sequenced genomes \n" % (
                    #                 name))

                elif seq_type == "profile":
                    print ('it is a profile')

                    hmm_path = "static/uploads/" + filename
                    while not os.path.exists(hmm_path):
                        time.sleep(1)
                    if os.path.isfile(hmm_path):
                        file = open(hmm_path, 'rb')

                        utilities.save_profile(file)

            except HTTPError as error:
                flash("There was a HTTP error. Please try again", category='error')

            seconds = timeit.default_timer() - start_time
            minutes, seconds = divmod(seconds, 60)
            hours, minutes = divmod(minutes, 60)

            periods = [('hours', hours), ('minutes', minutes), ('seconds', seconds)]
            time_string = ', '.join('{} {}'.format(value, name)
                                    for name, value in periods
                                    if value)

            print('\nFINISHED GETTING RECORDS: Time taken was {} \n'.format(time_string))
            if failed_genomes:
                flash ("The following genomes failed to map - " + " ".join(x for x in failed_genomes), category='error')
                print ("The following genomes failed to map - " + " ".join(x for x in failed_genomes))


        elif request.method == 'POST':
            print ('not validated')
            print (form.genome_type)


        return self.render("upload_admin.html", form=form)

        def is_accessible(self):
            if (not current_user.is_active or not
            current_user.is_authenticated):
                return False
            return True


class SetupView(BaseView):

    @login_required
    @expose("/", methods=('GET', 'POST'))
    def setup(self):
        form = forms.SetupForm()
        del_form = forms.DeleteForm()
        current = models.User.objects().get(username=str(current_user.username))
        prefs = {'page_size' : current.page_size, 'references' : current_user.references}

        if request.method == "POST":
            if form.submit.data:
                try:
                    models.User.objects().get(username=str(current_user.username)).update(page_size=form.page_size.data)
                    print ('data')
                    models.User.objects().get(username=str(current_user.username)).update(page_size=form.page_size.data)

                    for reference in request.form.getlist('references'):
                        models.User.objects().update(add_to_set__references=reference)

                    flash("User preferences updated", category='success')
                    return self.render('setup.html', form=form, del_form=del_form, prefs=prefs)

                except Exception as e:
                    print(e)
                    flash("Something went wrong", category='error')
                    return self.render('setup.html', form=form, del_form=del_form, prefs=prefs)

            else:
                return self.render('setup.html', form=form, del_form=del_form, prefs=prefs)

        elif request.method == "GET":
            print ('here')
            return self.render('setup.html', form=form, del_form=del_form, prefs=prefs)


class SequenceRecordsView(ModelView):


    edit_modal = True
    can_create = False
    can_view_details = True
    form_edit_rules = ('name', 'species', 'description', 'sequence')

    # list_template = 'list.html'
    # create_template = 'create.html'
    # edit_template = 'edit.html'

    column_formatters = {
        'sequence': custom_filters.seqdescription_formatter,
    }

    @action('generate_profile', 'Generate a profile from these sequences')
    def action_generate_profile(self, ids):
        try:
            query = models.SequenceRecords.objects(id in ids)
            print (query)
            align_list = []
            for record in query.all():
                align_record = SeqRecord(Seq(record.sequence, generic_protein), id=str(record.name) + "_" + "seq")
                align_list.append(align_record)

            utilities.createProfile(align_list)


        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext(ex))




class GenomeRecordsView(ModelView):
    form_edit_rules = ('name', 'species', 'strain', 'description')

    column_list = (
        'name', 'species', 'strain', 'description', 'sequence', 'Genome Overview')

    @login_required
    @expose("/", methods=('GET', 'POST'))
    def index_view(self):
        print (session.get('page_size'))
        current = models.User.objects().get(username=str(current_user.username))
        print (current.page_size)
        self.edit_modal = True
        self.can_create = False
        self.can_view_details = True

        self.page_size = current.page_size

        self.form_edit_rules = ('name', 'species', 'strain', 'description')

        self.column_list = (
            'name', 'species', 'strain', 'description', 'sequence', 'hits' 'Genome Overview')

        return super(ModelView, self).index_view()


    def is_accessible(self):
        if (not current_user.is_active or not
        current_user.is_authenticated):
            return False
        return True

    def _download_formatter(self, context, model, name):

        return Markup(
            "<a href='{url}' target='_blank'>View Genome</a>".format(url=self.get_url('download_genome_overview', id=model.name)))


    column_formatters = {
        'sequence': custom_filters.seqdescription_formatter,
        'Genome Overview': _download_formatter,

    }

    for name in ref_names:
        exec('@action("' + name + '_check_with_profile", "Check for ' + name + ' region with a profile")\ndef ' + name +
        '_check_with_profile(self, profile_ids):\n utilities.check_with_profile(profile_ids, "' + name + '")')

#
# class ProfileView(ModelView):
#     """
#     View of the Profile database for storing HMM Profiles """
#     column_list = (
#         'name', 'a1_profile_ref', 'a2_profile_ref', 'pore_profile_ref', 'chitinase_profile_ref',
# 'region1_profile_ref',
#         'region2_profile_ref',
#         'region3_profile_ref', 'region4_profile_ref', 'download')
#     form_columns = ('name', 'profile')
#
#     form_extra_fields = {'profile': models.BlobUploadField(
#         label='File',
#         allowed_extensions=['pdf', 'doc', 'docx', 'xls', 'xlsx', 'png', 'jpg', 'jpeg', 'gif', 'hmm', 'fa', 'fasta'],
#         size_field='size',
#         filename_field='filename',
#         mimetype_field='mimetype'
#     )}
#
#     def _download_formatter(self, context, model, name):
#         return Markup(
#             "<a href='{url}' target='_blank'>Download</a>".format(url=self.get_url('download_blob', id=model.uid)))
#
#     column_formatters = {
#         'download': _download_formatter,
#     }


class GenomeOverviewView(BaseView):

    @login_required
    @expose("/", methods=('GET', 'POST'))
    def genomeoverview(self):
        form = forms.GenomeOverviewSelectForm()
        # form.genome.choices = [(genome.name, genome.species) for genome in models.GenomeRecords.query.all()]
        return self.render('genomeoverview.html', form=form)


class GenomeDetailView(BaseView):

    @login_required
    @expose("/", methods=('GET', 'POST'))
    def genomedetail(self):
        return self.render('genomedetail.html')

class UserView(ModelView):

    @login_required
    @expose("/", methods=('GET', 'POST'))
    def index_view(self):
        self.edit_modal = True
        self.can_view_details = True

        self.can_create = False
        self.column_list = ['username']

        choice = bool(random.getrandbits(1))
        print ("Choice is ", choice)
        if choice:
            self.column_list = ['username']
        else:
            self.column_list = ['username', 'password']
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



class UserFormView(BaseView):

    @login_required
    @expose("/")
    def userform(self):
        form = forms.UserForm()
        return self.render('user.html', form=form)

class ProfileView(ModelView):
    """
    View of the Profile database for storing HMM Profiles """

    can_create = False


    @login_required
    @expose("/", methods=('GET', 'POST'))
    def index_view(self):

        def form_rules():
            action_rules = []

            if not has_app_context():
                ref_names = []
            else:
                current = models.User.objects().get(username=str(current_user.username))
                ref_names = current.references

            print('ref names is ', ref_names)

            for name in ref_names:

                action_rules.append(
                    '@action("' + name + '_set_reference", "Set this profile as the ' + name + ' reference")\ndef ' + name +
                    '_set_reference(self, ids):\n print("the name here is" + ids)')

            return action_rules

        self.column_list = (
            'name', 'references')
        # 'name', 'references', 'download')


        self.form_columns = ('name', 'profile')
        self.form_extra_fields = {'profile': models.BlobUploadField(
            label='File',
            allowed_extensions=['pdf', 'doc', 'docx', 'xls', 'xlsx', 'png', 'jpg', 'jpeg', 'gif', 'hmm', 'fa', 'fasta'],
            size_field='size',
            filename_field='filename',
            mimetype_field='mimetype'
        )}

        print ('here now in it')

        @action('doggydog', 'Generate a doggydog')
        def kittycat(self, ids):
            print('doggydog')

        self.action_form()
        # self.list_template = 'custom_list.html'
        #
        # ar = form_rules()
        # #
        # ref_names = ['region1', 'poopcity']
        # #
        # for name in ref_names:
        #     exec(
        #         '@action("' + name + '_set_reference", "Set this profile as the ' + name + ' reference")\ndef ' + name +
        #         '_set_reference(self, ids):\n print("the name here is" + ids)')
        #
        # print ('here is ar ')
        # print (ar)
        # for rule in ar:
        #     exec(rule)

        def on_model_change(self, form, model, is_created):
            print('change')

        def _download_formatter(self, context, model, name):

            return Markup(
                "<a href='{url}' target='_blank'>Download</a>".format(url=self.get_url('download_blob', id=model.name)))

        self.column_formatters = {
            'download': _download_formatter,
        }

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


    @property
    def _action_form_class(self):
        return self.get_action_form()

    @_action_form_class.setter
    def _action_form_class(self, value):
        pass


    for name in ref_names:
        exec('@action("' + name + '_set_reference", "Set this profile as the ' + name + ' reference")\ndef ' + name +
        '_set_reference(self, profile_ids):\n utilities.set_profile_as_reference(profile_ids, "' + name + '")')



class DocumentationView(BaseView):

    @expose("/")
    def documentation(self):
        return self.render('documentation.html')



@app.route("/overview_<id>", methods=['GET'])
def download_genome_overview(id):
    """
    Route for downloading profiles from the Profile view
    :param id: Profile to download
    :return:
    """
    print (id)
    record = models.GenomeRecords.objects().get(name=id)
    return send_file(
        record.genome_overview,
        attachment_filename=id + '.png')

@app.route("/<id>", methods=['GET'])
def download_blob(id):
    """
    Route for downloading profiles from the Profile view
    :param id: Profile to download
    :return:
    """
    print (id)
    profile = models.Profile.objects().get(name=id)
    return send_file(
        io.BytesIO(profile.profile),
        attachment_filename=id + '.hmm')

@app.errorhandler(404)
def page_not_found(e):
        return render_template("404.html")


@app.errorhandler(500)
def error_encountered(e):
    return render_template("500.html", error=e)


class MyModelView(ModelView):

    @property
    def _list_columns(self):
        return self.get_list_columns()

    @_list_columns.setter
    def _list_columns(self, value):
        pass

    @property
    def _action_form_class(self):
        return self.get_action_form()

    @_action_form_class.setter
    def _action_form_class(self, value):
        pass

    @property
    def column_list(self):
        if not has_app_context():
            @action('kittycat', 'Generate a kittycat')
            def kittycat(self, ids):
                print('kitty cat')

            return ['team', 'project_name', 'approve']
        else:

            @action('doggydog', 'Generate a doggydog')
            def kittycat(self, ids):
                print('doggydog')
            return ['team', 'project_name']

admin = Admin(app, 'Phylo Island', base_template='layout.html', url='/', template_mode='bootstrap3')

admin.add_view(UserView(model=models.User, endpoint='user'))
admin.add_view(SetupView(name='Setup', endpoint='setup'))
admin.add_view(UploadView(name='Upload', endpoint='upload_admin'))
admin.add_view(SequenceRecordsView(model=models.SequenceRecords, endpoint="sequence_records"))
admin.add_view(GenomeRecordsView(model=models.GenomeRecords, endpoint="genome_records"))
admin.add_view(ProfileView(model=models.Profile, name='Profiles', endpoint='profiles'))
# admin.add_view(MyModelView(model=models.Profile, name='Profiles', endpoint='models'))

admin.add_view(GenomeOverviewView(name='Genome Overview', endpoint='genomeoverview'))
admin.add_view(GenomeDetailView(name='Genome Detail', endpoint='genomedetail'))
admin.add_view(DocumentationView(name='Documentation & FAQ', endpoint='documentation'))
