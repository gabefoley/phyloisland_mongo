from flask_wtf import FlaskForm
from wtforms import FileField, StringField, SelectField, SelectMultipleField, BooleanField, PasswordField, SubmitField, validators
from wtforms.fields.html5 import IntegerField
from wtforms.validators import DataRequired


class LoginForm(FlaskForm):
    """Login form to access Phylo Island"""

    username = StringField('Username', validators=[DataRequired()])
    password = PasswordField('Password', validators=[DataRequired()])


class UploadForm(FlaskForm):
    """
    Form for uploading the initial files
    """
    file = FileField('Upload the file that contains the information we will map to the genome records.',
                     [validators.DataRequired()])
    # input_text = StringField(u'Input text', [validators.optional()])
    # type = SelectField('What type of file is this?', [validators.DataRequired()],
    #                    choices=[("protein", "FASTA (amino acids)"), ("nucleotide", "FASTA (nucleotides)"),
    #                             ("species", "Species list"), ("genome", "Genome ID list"), ("profile", "Profile")])
    type = SelectField('What type of file is this?', [validators.DataRequired()],
                       choices=[("protein", "FASTA (amino acids)"), ("species", "Species list"),  ("profile", "Profile")])
    add_sequence = BooleanField("Add sequences to sequence database?", default="checked")
    add_genome = BooleanField("Search for genomic records?", default="checked")
    single_genome = BooleanField("Retrieve just a single record for each genome?", default="checked")
    representative = BooleanField("If previous genome search fails, search RefSeq representative genomes",
                                  default="checked")
    assembly = BooleanField("If previous genome search fails, search RefSeq assembly genomes",
                            default="checked")
    genbank = BooleanField("If previous genome search fails, search GenBank assembly genomes",
                            default="checked")

    search_shotgun = BooleanField("Search for shotgun sequenced genomes if we can't find another "
                                  "genomic record?", default="checked")
    genome_type = SelectField('Which genome records should we return?', choices=[
        ('reference genome','RefSeq reference genome/s'),
        ('representative genome', 'RefSeq representative genome/s'),
        ('assembly', 'RefSeq assembly genome/s'),
        ('genbank', 'GenBank assembly genome/s')])

    upload_submit = SubmitField("Upload file")


class SetupForm(FlaskForm):
    """
    Form for customising user preferences
    """
    page_size = IntegerField(u'Number of records per page in Genome Records ', [validators.optional()])
    record_size = IntegerField(u'Number of records per page in Genome Detail', [validators.optional()])
    submit = SubmitField("Submit")

# class DeleteReferenceForm(FlaskForm):

class DeleteForm(FlaskForm):
    """
    Form for deleting references
    """
    del_references = SelectField('Reference to delete')
    submit = SubmitField("Delete reference")


class GenomeOverviewSelectForm(FlaskForm):
    """
    Form for selecting which genomes to look at in the Genome Overview
    """
    genome = SelectMultipleField('Genome', choices=[])
    submit = SubmitField("Submit")

class GenomeDiagramSelectForm(FlaskForm):
    """
    Form for selecting which genomes to look at in the Genome Diagram
    """
    genome = SelectField('Genome', choices=[])
    submit_diagram = SubmitField("Submit")
    tag_genome = SubmitField('Tag genome')
    clear_genome_tags = SubmitField('Clear genome tags')
    submit_hit = SubmitField("Add tag")
    hide_hit = SubmitField('Hide hit')

    delete_hit = SubmitField("Delete hit")

    associate_hits = SubmitField("Associate hits")

    hidden_hits = BooleanField("Hide hits marked 'hidden'")

    show_hits = SelectField('Which hits should we show?', choices=[('all', 'All hits'), ('initial', 'Just initial '
                                                                                                    'hits'),
                                                                   ('expanded', 'Just expanded hits')])

    show_promoters = BooleanField("Show promoters")
    show_stop_codons = BooleanField("Show stop codons")



class GenomeDiagramPageForm(FlaskForm):
    untagged = BooleanField("Limit selection to untagged genomes")
    limit_genomes = BooleanField("Limit selection to genomes tagged with:")
    genome_tagged = SelectField( choices=[])
    limit_selection = SubmitField('Limit selection')
    page = SelectField('Page to display', choices=[])
    select_page = SubmitField('Select page')


class ChartsForm(FlaskForm):
    select_tags = SelectMultipleField('Tags to include', choices=[])
    exclude_tags = SelectMultipleField('Tags to exclude', choices=[])

    update_chart = SubmitField('Update chart')


class GenomeDiagamShowRegions(FlaskForm):
    """
    Form for selecting which regions to look at in the Genome Diagram
    """

    showA1 = BooleanField("Show A1")
    showA2 = BooleanField("Show A2")
    showTcdA1 = BooleanField("Show TcdA1")
    showTcB = BooleanField("Show TcB")
    showTcC = BooleanField("Show TcC")
    showChitinase = BooleanField("Show Chitinase")

    show_regions = SubmitField("Update regions to show")


class GenomeHitForm(FlaskForm):
    """
    Form for selecting which genomes to look at in the Genome Diagram
    """
    tag = StringField('Add tag', [validators.optional()])
    submit_hit = SubmitField("Add tag")

    delete_hit = SubmitField("Delete hit")

class DownloadFastaForm(FlaskForm):
        region = SelectField('Which region should we download?', choices=[('A1', 'A1'), ('A2', 'A2'), ('TcdA1',
                                                                                                       'TcdA1'),
        ('Chitinase','Chitinase'), ('TcB', 'TcB'), ('TcC', 'TcC'), ('A1_expanded', 'A1_expanded'), ('A2_expanded',
                                                                                                    'A2_expanded'),
                                                                          ('TcdA1_expanded',
                                                                                                       'TcdA1_expanded'),
        ('Chitinase_expanded','Chitinase_expanded'), ('TcB_expanded', 'TcB_expanded'), ('TcC_expanded', 'TcC_expanded')])

        filename = StringField('Append extra text to filename?')

        include_genome = StringField('Only include genomes tagged with - ')

        exclude_genome = StringField('Exclude genomes tagged with - ')

        include_hits = StringField('Only include hits tagged with - ')

        exclude_hits = StringField('Exclude hits tagged with - ', default='hidden')



        translate = BooleanField("Translate nucleotides to proteins?",
                                default="checked")
        align = BooleanField("Also create an alignment?")


        submit = SubmitField("Submit")

class TempFixForm(FlaskForm):
    fix_assoc = SubmitField("Click this to fix the associated regions in the database")

class TagSimpleForm(FlaskForm):
    include_genome = StringField('Only include genomes tagged with - ')
    exclude_hits = StringField('Exclude hits tagged with - ', default='hidden')
    tag_simple = SubmitField("Tag genomes as simple")

class SearchForPromoters(FlaskForm):
    mismatch = IntegerField("Number of allowable mismatches")
    search_for_promoters = SubmitField("Search for promoters")

class UpdateTagsForm(FlaskForm):
     old_tag = SelectField('Tag to update', choices=[])
     new_tag = StringField('Change it to ', [validators.DataRequired()])
     update_tags = SubmitField("Update tags")


class UploadRegion(FlaskForm):
    name = StringField('Region name', [validators.DataRequired()])
    file = FileField('Upload the FASTA file.',
                     [validators.DataRequired()])
    upload_submit = SubmitField("Upload Region")

class RegionForm(FlaskForm):
    region = SelectField('Select region to use ', choices=[])
    profiles = SelectMultipleField('Search region using these profiles', choices=[])
    search_regions = SubmitField("Search regions")


class AlignmentForm(FlaskForm):
  name = StringField('Alignment name ', [validators.DataRequired()])
  region = SelectField('Make alignment based on ', choices=[])
  tool = SelectField('Select alignment tool - ', choices=[('MAFFT', 'MAFFT')])
  align = SubmitField('Make alignment')

class SelectAlignmentForm(FlaskForm):
    name = SelectField('Alignment name ', choices=[])
    submit = SubmitField("Select alignment")

class SelectRegionToProfilesForm(FlaskForm):
    name = SelectField('Region to profiles name ', choices=[])
    submit = SubmitField("Select region to profiles")

class MakeTreeForm(FlaskForm):
  name = StringField('Tree name ', [validators.DataRequired()])
  alignment = SelectField('Make tree based on ', choices=[])
  tool = SelectField('Select tree inference tool - ', choices=[('FastTree', 'FastTree')])
  make_tree = SubmitField('Make tree')

class TreeDownloadForm(FlaskForm):
    tree = SelectField('Select tree to download ', choices=[])
    download_tree = SubmitField('Download tree')

class RerootTreeForm(FlaskForm):
    tree = SelectField('Reroot this tree - ', choices=[])
    rerooted_tree_name = StringField('Tree name ', [validators.DataRequired()])
    seq = SelectField('Set this sequence as outgroup - ', choices=[])
    reroot_tree = SubmitField('Reroot tree')

class TrimToProfileForm(FlaskForm):
    trim_to_region = SelectField('Regions to trim to', choices=[])
    trim_to_name = StringField('Name of new Regions file ', [validators.DataRequired()])
    trim_to_profile = SelectField('Choose profile to trim Regions to :', choices=[])
    trim_to = SubmitField("Trim to profile")

class TrimAroundProfileForm(FlaskForm):
    trim_around_region = SelectField('Regions to trim around', choices=[])
    trim_around_name = StringField('Name of new Regions file ', [validators.DataRequired()])
    trim_around_profile = SelectMultipleField('Place content where you would like to trim to :', choices=[])
    section1 = BooleanField("Include the content from the profile labelled - ")
    section2 = BooleanField("Include the content from the profile labelled - ")
    trim_around_submit = SubmitField("Trim around profile")

class TreeSelectForm(FlaskForm):
    tree_select_name = SelectField('Tree name ', choices=[])
    profiles = SelectField('Add profile search results from', choices=[])
    region_order = SelectField('Add region order search results from', choices=[])
    full_names = BooleanField("Display full names on tree", default='unchecked')
    collapse_on_genome_tags = BooleanField("Collapse tree based on genome tags", default='unchecked')

    submit = SubmitField("Select tree")


class BatchDeleteForm(FlaskForm):
    delete_all_tags = SubmitField('Delete all tags')
    delete_all_hits = SubmitField('Delete all hits')
    # delete_all_region_tags = SubmitField('Delete all tags at the region level')


class DownloadAssociatedRegions(FlaskForm):

    associated_regions = SubmitField("Download Associated Regions")


class DownloadTags(FlaskForm):

    tags = SubmitField("Download Tags")



class DownloadRegionOrder(FlaskForm):
    include_genome = StringField('Only include genomes tagged with - ')

    exclude_genome = StringField('Exclude genomes tagged with - ')

    include_hits = StringField('Only include hits tagged with - ')

    exclude_hits = StringField('Exclude hits tagged with - ', default='hidden')

    save_to_db = BooleanField('Save to database as well', default='checked')

    submit = SubmitField("Submit")

