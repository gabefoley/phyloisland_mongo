from phyloisland import app, allfiles
import models
import forms
import utilities
import getGenomes
import custom_filters
from flask import render_template, flash, request, session
from flask_login import login_required, current_user
from flask_admin import Admin, expose, AdminIndexView, BaseView
from flask_admin.contrib.mongoengine import ModelView, filters
import timeit
import time
import os
from urllib.error import HTTPError


class UploadView(BaseView):

    """
    View for uploading files
    """

    @login_required
    @expose("/", methods=('GET', 'POST'))
    def upload(self):
        form = forms.UploadForm()

        if request.method == 'POST' and form.validate():

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

                    hmm_path = "static/uploads/" + filename
                    while not os.path.exists(hmm_path):
                        time.sleep(1)
                    if os.path.isfile(hmm_path):
                        file = open(hmm_path, 'rb')

                        utilities.saveProfile(file)

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
                print("List of failed genomes - ", failed_genomes)
        return self.render("upload_admin.html", form=form)


class SetupView(BaseView):

    @login_required
    @expose("/", methods=('GET', 'POST'))
    def setup(self):
        form = forms.SetupForm()
        current = models.User.objects().get(username=str(current_user.username))
        if request.method == "POST":
            if form.submit.data:
                try:
                    models.User.objects().get(username=str(current_user.username)).update(page_size=form.page_size.data)
                    print ('data')
                    models.User.objects().get(username=str(current_user.username)).update(page_size=form.page_size.data)

                    for reference in request.form.getlist('references'):
                        models.User.objects().update(add_to_set__references=reference)

                    flash("User preferences updated", category='success')
                    return self.render('setup.html', form=form)

                except Exception as e:
                    print(e)
                    flash("Something went wrong", category='error')
                    return self.render('setup.html', form=form)

            else:
                return self.render('setup.html', form=form)

        elif request.method == "GET":
            return self.render('setup.html', form=form)


class SequenceRecordsView(ModelView):


    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True
    # list_template = 'list.html'
    # create_template = 'create.html'
    # edit_template = 'edit.html'

    column_formatters = {
        'sequence': custom_filters.seqdescription_formatter,
    }


class GenomeRecordsView(ModelView):

    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True
    page_size = 2

    column_formatters = {
        'sequence': custom_filters.seqdescription_formatter,
    }


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


class MyAdminIndexView(AdminIndexView):
    def is_accessible(self):
        return current_user.is_authenticated


class MyHomeView(AdminIndexView):
    @expose("/")
    def index(self):
        if 'username' in session:
            return 'You are logged in as ' + session['username']
        else:
            form = forms.LoginForm()
            return self.render('admin/index.html', form=form)


class UserView(ModelView):
    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True
    # column_list = ['name', 'password', 'thomas']


class UserFormView(BaseView):

    @login_required
    @expose("/")
    def userform(self):
        form = forms.UserForm()
        return self.render('user.html', form=form)


@app.errorhandler(404)
def page_not_found(e):
        return render_template("404.html")


@app.errorhandler(500)
def error_encountered(e):
    return render_template("500.html", error=e)

admin = Admin(app, 'Phylo Island', base_template='layout.html', url='/', template_mode='bootstrap3')

admin.add_view(UserView(model=models.User, endpoint='user'))
admin.add_view(SetupView(name='Setup', endpoint='setup'))
admin.add_view(UploadView(name='Upload', endpoint='upload_admin'))
admin.add_view(SequenceRecordsView(model=models.SequenceRecords, endpoint="sequence_records"))
admin.add_view(GenomeRecordsView(model=models.GenomeRecords, endpoint="genome_records"))
# admin.add_view(ProfileView(model=models.Profile, session=db.session, name='Profiles'))
admin.add_view(GenomeOverviewView(name='Genome Overview', endpoint='genomeoverview'))
admin.add_view(GenomeDetailView(name='Genome Detail', endpoint='genomedetail'))
