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
    page_size = IntegerField(u'Preferred page size', [validators.optional()])
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
    genome = SelectMultipleField('Genome', choices=[])
    submit_diagram = SubmitField("Submit")
    tag_genome = SubmitField('Tag genome')
    submit_hit = SubmitField("Add tag")
    hide_hit = SubmitField('Hide hit')

    delete_hit = SubmitField("Delete hit")

    show_hits = SelectField('Which hits should we show?', choices=[('all', 'All hits'), ('initial', 'Just initial '
                                                                                                    'hits'),
                                                                   ('expanded', 'Just expanded hits')])


class GenomeHitForm(FlaskForm):
    """
    Form for selecting which genomes to look at in the Genome Diagram
    """
    tag = StringField('Add tag', [validators.optional()])
    submit_hit = SubmitField("Add tag")

    delete_hit = SubmitField("Delete hit")

class DownloadFastaForm(FlaskForm):
        region = SelectField('Which region should we download for?', choices=[
        ('Chitinase','Chitinase'),
        ('TcC', 'TcC')])
        filename = StringField('Append extra text to filename?')

        translate = BooleanField("Translate nucleotides to proteins?",
                                default="checked")
        align = BooleanField("Also create an alignment?",
                               default="checked")
        submit = SubmitField("Submit")
