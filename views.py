from phyloisland import db, app, allfiles
import models
import forms
import utilities
import getGenomes
import custom_filters
from flask import render_template, flash, request
from flask_admin import Admin, expose, AdminIndexView, BaseView
from flask_admin.actions import action
from flask_admin.contrib.sqla import ModelView, filters
from markupsafe import Markup
import timeit
import time
import os
from urllib.error import HTTPError

# Setup a default value for page size
selected_page_size = 50


@app.route('/dashboard/')
def dashboard():
    flash ("flash test")
    return render_template("dashboard.html")


@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html")


@app.errorhandler(500)
def error_encountered(e):
    return render_template("500.html", error=e)


class MyHomeView(AdminIndexView):
    @expose("/")
    def index(self):
        return self.render('index.html')

class UploadView(BaseView):

    """
    View for uploading files
    """

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
            # search_shotgun = form.search_shotgun.data
            single = form.single_genome.data
            genome_type=form.genome_type.data
            representative = form.representative.data
            assembly = form.assembly.data
            genbank = form.genbank.data
            failed_genomes = []

            # # Create the initial seqRecords
            # phyloisland.seqDict = {}
            # phyloisland.unmappable = []

            try:

                if seq_type == "protein" or seq_type == "nucleotide":
                    seq_records = utilities.read_fasta("static/uploads/" + filename)

                    if not seq_records:
                        print("Couldn't find any sequences in the uploaded file.")
                        flash("Couldn't find any sequences in the uploaded file.")

                    else:
                        if add_seq:
                            utilities.addSequence(seq_records)

                        if add_genome:
                            #TODO: Correct this call
                            getGenomes.add_genome(seq_records, seq_type)




                elif seq_type == "species":

                    species_names = utilities.readLinesFromFile("static/uploads/" + filename)

                    for species_name in species_names:

                        print ("Species name is ", species_name)

                        destinations = [genome_type]

                        #TODO: Once the checkbox selection is dynamic we can just add freely here
                        if representative and genome_type not in ["representative genome", "assembly", "genbank"]:
                            destinations.append("representative genome")

                        if assembly and genome_type not in ["assembly", "genbank"]:
                            destinations.append("assembly")

                        print ("Destinations is ", destinations)

                        genome_results = getGenomes.add_genome2(species_name, destinations, single=single)

                        if genome_results and genome_results != "Fail":
                            print ("GENOME RESULTS IS ")
                            print (genome_results)
                            utilities.addGenome(genome_results)
                        else:

                            # Couldn't find it in RefSeq, let's try genbank
                            if genbank:
                                destinations = ["genbank"]
                                genome_results = getGenomes.add_genome2(species_name, destinations, single=single)

                                if genome_results and genome_results != "Fail":
                                    print("GENOME RESULTS IS ")
                                    print(genome_results)
                                    utilities.addGenome(genome_results)
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
                    #             "\nWe didn't identify any genome records for %s. And we are not attempting to search for shotgun sequenced genomes \n" % (
                    #                 name))

                elif seq_type == "profile":

                    hmm_path = "static/uploads/" + filename
                    while not os.path.exists(hmm_path):
                        time.sleep(1)
                    if os.path.isfile(hmm_path):
                        file = open(hmm_path, 'rb')

                        utilities.saveProfile(file)

            except HTTPError as error:
                flash("There was a HTTP error. Please try again")

            seconds = timeit.default_timer() - start_time
            minutes, seconds = divmod(seconds, 60)
            hours, minutes = divmod(minutes, 60)

            periods = [('hours', hours), ('minutes', minutes), ('seconds', seconds)]
            time_string = ', '.join('{} {}'.format(value, name)
                                    for name, value in periods
                                    if value)

            print('\nFINISHED GETTING RECORDS: Time taken was {} \n'.format( time_string))
            if failed_genomes:
                print ("List of failed genomes - ", failed_genomes)
        return self.render("upload_admin.html", form=form)

class SetupView(BaseView):
    @expose("/", methods=('GET', 'POST'))

    def setup(self):
        form = forms.SetupForm()
        if request.method == "POST":
            if form.submit.data:
                try:
                    print (form.name.data)
                    print (form.page_size.data)
                    print (form.region_1_name.data)

                    user = models.User("gabe")
                    setup_submission = models.User(form.name.data, int(form.page_size.data), form.region_1_name.data)

                    print ('got here')
                    db.session.add(setup_submission)
                    print ('added')
                    db.session.commit()
                    print ('commit')
                    flash ("User preferences updated")
                    return self.render('setup.html', form=form)


                except Exception as e:
                    print (e)
                    flash ("Something went wrong")
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

    def _seqdescription_formatter(view, context, model, name):

        if model.sequence:
            return model.sequence[:15] + "..."
        else:
            return model.sequence

    column_formatters = {
        'sequence': _seqdescription_formatter,
    }


class GenomeRecordsView(ModelView):
    column_list = (
        'name', 'species', 'strain', 'description', 'a1_ref', 'a2_ref', 'pore_ref', 'sequence', 'a1', 'a1_length',
        'a1_loc',
        'a2',
        'a2_length', 'a2_loc', 'overlap', 'distance', 'pore', 'pore_length', 'pore_loc', 'pore_within_a2',
        'chitinase', 'chitinase_length', 'chitinase_loc', 'chitinase_distance_from_a2', 'region1_ref', 'region2_ref',
        'region3_ref', 'region4_ref', 'region1', 'region1_length', 'region1_loc', 'region2',
        'region2_length', 'region2_loc', 'region3', 'region3_length', 'region3_loc', 'region4',
        'region4_length', 'region4_loc')
    column_searchable_list = ['name', 'species', 'a1', 'a2', 'overlap']
    create_modal = True
    edit_modal = True
    can_create = False
    can_view_details = True
    page_size = selected_page_size

    def _a1description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.a1:
            return model.a1[:15] + "..."
        else:
            return model.a1

    def _a2description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.a2:
            return model.a2[:15] + "..."
        else:
            return model.a2

    def _chitinase_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.chitinase:
            return model.chitinase[:15] + "..."
        else:
            return model.chitinase

    def _pore_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.pore:
            return model.pore[:15] + "..."
        else:
            return model.pore

    def _region1description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.region1:
            return model.region1[:15] + "..."
        else:
            return model.region1

    def _region2description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.region2:
            return model.region2[:15] + "..."
        else:
            return model.region2

    def _region3description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.region3:
            return model.region3[:15] + "..."
        else:
            return model.region3

    def _region4description_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window

        if model.region1:
            return model.region4[:15] + "..."
        else:
            return model.region4

    def _seqdescription_formatter(view, context, model, name):
        # Format your string here e.g show first 20 characters
        # can return any valid HTML e.g. a link to another view to show the detail or a popup window


        if model.sequence:
            return model.sequence[:15] + "..."
        else:
            return model.sequence

    column_formatters = {
        'a1': _a1description_formatter,
        'a2': _a2description_formatter,
        'pore': _pore_formatter,
        'chitinase': _chitinase_formatter,
        'region1': _region1description_formatter,
        'region2': _region2description_formatter,
        'region3': _region3description_formatter,
        'region4': _region4description_formatter,

        'sequence': _seqdescription_formatter,
    }

    column_filters = ('name',
                      'species',
                      'strain',
                      'description',
                      'a1',
                      'a2',
                      filters.FilterLike(models.GenomeRecords.overlap, 'Overlap',
                                         options=(('True', 'True'), ('False', 'False'))),
                      'sequence',
                      custom_filters.GetUniqueSpecies(
                          models.GenomeRecords.name, 'Get unique species'),
                      filters.IntGreaterFilter(models.GenomeRecords.distance, 'Distance greater than'),
                      filters.IntSmallerFilter(models.GenomeRecords.distance, 'Distance smaller than'),
                      filters.IntGreaterFilter(models.GenomeRecords.chitinase_distance_from_a2,
                                               'Chitinase to A2 distance greater than'),
                      filters.IntSmallerFilter(models.GenomeRecords.chitinase_distance_from_a2,
                                               'Chitinase to A2 distance smaller than')

                      )


class ProfileView(ModelView):
    """
    View of the Profile database for storing HMM Profiles """
    column_list = (
        'name', 'a1_profile_ref', 'a2_profile_ref', 'pore_profile_ref', 'chitinase_profile_ref', 'region1_profile_ref',
        'region2_profile_ref',
        'region3_profile_ref', 'region4_profile_ref', 'download')
    form_columns = ('name', 'profile')

    form_extra_fields = {'profile': models.BlobUploadField(
        label='File',
        allowed_extensions=['pdf', 'doc', 'docx', 'xls', 'xlsx', 'png', 'jpg', 'jpeg', 'gif', 'hmm', 'fa', 'fasta'],
        size_field='size',
        filename_field='filename',
        mimetype_field='mimetype'
    )}

    def _download_formatter(self, context, model, name):
        return Markup(
            "<a href='{url}' target='_blank'>Download</a>".format(url=self.get_url('download_blob', id=model.uid)))

    column_formatters = {
        'download': _download_formatter,
    }

    @action('item1_set_A1_reference', 'Set this profile as the A1 reference profile')
    def item1_set_a1_reference(self, ids):
        utilities.setProfileAsReference(ids, "a1")

    @action('item2_set_A2_reference', 'Set this profile as the A2 reference profile')
    def item1_set_a2_reference(self, ids):
        utilities.setProfileAsReference(ids, "a2")

    @action('item2_set_pore_reference', 'Set this profile as the pore reference profile')
    def item1_set_pore_reference(self, ids):
        utilities.setProfileAsReference(ids, "pore")

    @action('item2_set_chitinase_reference', 'Set this profile as the chitinase reference profile')
    def item1_set_chitinase_reference(self, ids):
        utilities.setProfileAsReference(ids, "chitinase")

    @action('item2_set_region1_reference', 'Set this profile as the region 1 reference profile')
    def item1_set_region1_reference(self, ids):
        utilities.setProfileAsReference(ids, "region1")

    @action('item2_set_region2_reference', 'Set this profile as the region 2 reference profile')
    def item1_set_region2_reference(self, ids):
        utilities.setProfileAsReference(ids, "region2")

    @action('item2_set_region3_reference', 'Set this profile as the region 3 reference profile')
    def item1_set_region3_reference(self, ids):
        utilities.setProfileAsReference(ids, "region3")

    @action('item2_set_region4_reference', 'Set this profile as the region 4 reference profile')
    def item1_set_region4_reference(self, ids):
        utilities.setProfileAsReference(ids, "region4")

class GenomeOverviewView(BaseView):
    @expose("/", methods=('GET', 'POST'))

    def genomeoverview(self):
        form = forms.GenomeOverviewSelctForm()
        form.genome.choices = [(genome.name, genome.species) for genome in models.GenomeRecords.query.all()]
        return self.render('genomeoverview.html', form=form)


class GenomeDetailView(BaseView):
    @expose("/", methods=('GET', 'POST'))

    def genomedetail(self):
        return self.render('genomedetail.html')



admin = Admin(app, 'Phylo Island', base_template='layout.html', url='/', template_mode='bootstrap3')
admin.add_view(UploadView(name='Upload', endpoint='upload_admin'))
admin.add_view(SetupView(name='Setup', endpoint='setup'))
admin.add_view(SequenceRecordsView(model=models.SequenceRecords, session=db.session, endpoint="sequence_records"))
admin.add_view(GenomeRecordsView(model=models.GenomeRecords, session=db.session, endpoint="genome_view"))
admin.add_view(ProfileView(model=models.Profile, session=db.session, name='Profiles'))
admin.add_view(GenomeOverviewView(name='Genome Overview', endpoint='genomeoverview'))
admin.add_view(GenomeDetailView(name='Genome Detail', endpoint='genomedetail'))
