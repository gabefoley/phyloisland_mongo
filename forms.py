from flask_wtf import FlaskForm
from wtforms import FileField, StringField, SelectField, SelectMultipleField, BooleanField, SubmitField, validators
from wtforms.fields.html5 import IntegerField


class UploadForm(FlaskForm):
    """
    Form for uploading the initial files
    """
    file = FileField('Upload the file that contains the information we will map to the genome records.',
                     [validators.DataRequired()])
    input_text = StringField(u'Input text', [validators.optional()])
    type = SelectField('What type of file is this?', [validators.DataRequired()],
                       choices=[("protein", "FASTA (amino acids)"), ("nucleotide", "FASTA (nucleotides)"),
                                ("species", "Species list"), ("genome", "Genome ID list"), ("profile", "Profile")])
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
        ('reference genome','Retrieve RefSeq reference genome/s'),
        ('representative genome', 'Retrieve RefSeq representative genome/s'),
        ('assembly', 'Retrieve RefSeq assembly genome/s'),
        ('genbank', 'Retrieve GenBank assembly genome/s')])

    upload_submit = SubmitField("Upload file")

class SetupForm(FlaskForm):
    name = StringField(u'User name', [validators.optional()])
    page_size = IntegerField(u'Preferred page size', [validators.optional()])
    region_1_name = StringField(u'Region 1 name', [validators.optional()])
    submit = SubmitField("Submit")



class GenomeOverviewSelctForm(FlaskForm):
    genome = SelectMultipleField('Genome', choices=[])
    submit = SubmitField("Submit")
