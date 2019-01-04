from phyloisland import db
from wtforms import ValidationError, fields
from wtforms.widgets import FileInput
from werkzeug.datastructures import FileStorage
from wtforms.validators import required
from gettext import gettext


class User(db.Model):
    """
    Class to hold information about a User, such as their password and preferences for naming fields
    """
    __tablename__ = 'user'
    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text)
    page_size = db.Column(db.Integer)
    region_1_name = db.Column(db.Text)

    def __init__(self, name="", page_size=50, region_1_name='region_1'):
        self.name = name
        self.page_size = page_size
        self.region_1_name = region_1_name



class SequenceRecords(db.Model):
    """
    Class to hold information about a Sequence that is separate from a Genome. Generally will be uploaded by the user to
    create Profiles with
    """
    __tablename__ = 'seqrecord'
    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text)
    species = db.Column(db.String(255))
    description = db.Column(db.Text)
    sequence = db.Column(db.Text)
    a1_ref = db.Column(db.Boolean)
    a2_ref = db.Column(db.Boolean)
    pore_ref = db.Column(db.Boolean)
    chitinase_ref = db.Column(db.Boolean)
    region1_ref = db.Column(db.Boolean)
    region2_ref = db.Column(db.Boolean)
    region3_ref = db.Column(db.Boolean)
    region4_ref = db.Column(db.Boolean)

    def __init__(self, name="", species="", description="", sequence="", a1_ref=0, a2_ref=0, pore_ref=0,
                 chitinase_ref=0, region1_ref=0,
                 region2_ref=0, region3_ref=0, region4_ref=0, ):
        self.name = name
        self.species = species
        self.description = description
        self.sequence = sequence
        self.a1_ref = a1_ref
        self.a2_ref = a2_ref
        self.pore_ref = pore_ref
        self.chitinase_ref = chitinase_ref
        self.region1_ref = region1_ref
        self.region2_ref = region2_ref
        self.region3_ref = region3_ref
        self.region4_ref = region4_ref

class GenomeRecords(db.Model):
    __tablename__ = 'genomerecord'
    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text)
    species = db.Column(db.String(255))
    strain = db.Column(db.String(255))
    description = db.Column(db.Text)
    a1 = db.Column(db.Text)
    a1_length = db.Column(db.Integer)
    a1_loc = db.Column(db.Text)
    a2 = db.Column(db.Text)
    a2_length = db.Column(db.Integer)
    a2_loc = db.Column(db.Text)
    overlap = db.Column(db.String(255))
    distance = db.Column(db.VARCHAR(255))
    sequence = db.Column(db.Text)
    a1_ref = db.Column(db.Boolean)
    a2_ref = db.Column(db.Boolean)
    pore = db.Column(db.Text)
    pore_length = db.Column(db.Integer)
    pore_loc = db.Column(db.Text)
    pore_within_a2 = db.Column(db.String(255))
    pore_ref = db.Column(db.Boolean)
    chitinase = db.Column(db.Text)
    chitinase_length = db.Column(db.Integer)
    chitinase_loc = db.Column(db.Text)
    chitinase_distance_from_a2 = db.Column(db.VARCHAR(255))
    chitinase_ref = db.Column(db.Boolean)
    region1 = db.Column(db.Text)
    region1_length = db.Column(db.Integer)
    region1_loc = db.Column(db.Text)
    region2 = db.Column(db.Text)
    region2_length = db.Column(db.Integer)
    region2_loc = db.Column(db.Text)
    region3 = db.Column(db.Text)
    region3_length = db.Column(db.Integer)
    region3_loc = db.Column(db.Text)
    region4 = db.Column(db.Text)
    region4_length = db.Column(db.Integer)
    region4_loc = db.Column(db.Text)
    region1_ref = db.Column(db.Boolean)
    region2_ref = db.Column(db.Boolean)
    region3_ref = db.Column(db.Boolean)
    region4_ref = db.Column(db.Boolean)

    def __init__(self, name="", species="", strain="", description="", a1="", a1_length="", a1_loc="",
                 a2="", a2_length="", a2_loc="", overlap="", distance="", sequence="", a1_ref=0, a2_ref=0,
                 pore="", pore_length=None, pore_loc="", pore_within_a2="", pore_ref=0,
                 chitinase="", chitinase_length=None, chitinase_loc="", chitinase_distance_from_a2="", chitinase_ref=0,
                 region1="", region1_length=None, region1_loc="", region2="", region2_length=None, region2_loc="",
                 region3="", region3_length=None, region3_loc="", region4="", region4_length=None, region4_loc="",
                 region1_ref="", region2_ref="", region3_ref="", region4_ref=""):
        self.name = name
        self.species = species
        self.strain = strain
        self.description = description
        self.a1 = a1
        self.a1_length = a1_length
        self.a1_loc = a1_loc
        self.a2 = a2
        self.a2_length = a2_length
        self.a2_loc = a2_loc
        self.overlap = overlap
        self.distance = distance
        self.sequence = sequence
        self.a1_ref = a1_ref
        self.a2_ref = a2_ref
        self.pore = pore
        self.pore_length = pore_length
        self.pore_loc = pore_loc
        self.pore_within_a2 = pore_within_a2
        self.pore_ref = pore_ref
        self.chitinase = chitinase
        self.chitinase_length = chitinase_length
        self.chitinase_loc = chitinase_loc
        self.chitinase_distance_from_a2 = chitinase_distance_from_a2
        self.chitinase_ref = chitinase_ref
        self.region1 = region1
        self.region1_length = region1_length
        self.region1_loc = region1_loc
        self.region2 = region2
        self.region2_length = region2_length
        self.region2_loc = region2_loc
        self.region3 = region3
        self.region3_length = region3_length
        self.region3_loc = region3_loc
        self.region4 = region4
        self.region4_length = region4_length
        self.region4_loc = region4_loc
        self.region1_ref = region1_ref
        self.region2_ref = region2_ref
        self.region3_ref = region3_ref
        self.region4_ref = region4_ref

class BlobMixin(object):
    mimetype = db.Column(db.Unicode(length=255), nullable=False)
    filename = db.Column(db.Unicode(length=255), nullable=False)
    profile = db.Column(db.BLOB, nullable=False)
    size = db.Column(db.Integer, nullable=False)

    def __init__(self, mimetype, filename, profile, size):
        self.mimetype = mimetype
        self.filename = filename
        self.profile = profile
        self.size = size


class Profile(db.Model, BlobMixin):
    __tablename__ = 'profile_blob'

    uid = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Unicode(length=255), nullable=False, unique=True)
    a1_profile_ref = db.Column(db.Boolean)
    a2_profile_ref = db.Column(db.Boolean)
    pore_profile_ref = db.Column(db.Boolean)
    chitinase_profile_ref = db.Column(db.Boolean)
    region1_profile_ref = db.Column(db.Boolean)
    region2_profile_ref = db.Column(db.Boolean)
    region3_profile_ref = db.Column(db.Boolean)
    region4_profile_ref = db.Column(db.Boolean)

    def __init__(self, name="", blob_mix="", a1_profile_ref=0, a2_profile_ref=0, pore_profile_ref=0,
                 chitinase_profile_ref=0, region1_profile_ref=0,
                 region2_profile_ref=0, region3_profile_ref=0, region4_profile_ref=0, ):
        self.name = name
        self.a1_profile_ref = a1_profile_ref
        self.a2_profile_ref = a2_profile_ref
        self.pore_profile_ref = pore_profile_ref
        self.chitinase_profile_ref = chitinase_profile_ref
        self.region1_profile_ref = region1_profile_ref
        self.region2_profile_ref = region2_profile_ref
        self.region3_profile_ref = region3_profile_ref
        self.region4_profile_ref = region4_profile_ref

        self.blob_mix = blob_mix

    def set_blobMix(self, blob_mix):
        self.blob_mix = blob_mix
        self.mimetype = blob_mix.mimetype
        self.filename = blob_mix.filename
        self.profile = blob_mix.profile
        self.size = blob_mix.size

    def __unicode__(self):
        return u"name : {name}; filename : {filename})".format(name=self.name, filename=self.filename)

class BlobUploadField(fields.StringField):
    widget = FileInput()

    def __init__(self, label=None, allowed_extensions=None, size_field=None, filename_field=None, mimetype_field=None,
                 **kwargs):

        self.allowed_extensions = allowed_extensions
        self.size_field = size_field
        self.filename_field = filename_field
        self.mimetype_field = mimetype_field
        validators = [required()]

        super(BlobUploadField, self).__init__(label, validators, **kwargs)

    def is_file_allowed(self, filename):
        """
            Check if file extension is allowed.

            :param filename:
                File name to check
        """
        if not self.allowed_extensions:
            return True

        return ('.' in filename and
                filename.rsplit('.', 1)[1].lower() in
                map(lambda x: x.lower(), self.allowed_extensions))

    def _is_uploaded_file(self, data):
        return (data and isinstance(data, FileStorage) and data.filename)

    def pre_validate(self, form):
        super(BlobUploadField, self).pre_validate(form)
        if self._is_uploaded_file(self.data) and not self.is_file_allowed(self.data.filename):
            raise ValidationError(gettext('Invalid file extension'))

    def process_formdata(self, valuelist):
        if valuelist:
            data = valuelist[0]
            self.data = data

    def populate_obj(self, obj, name):

        if self._is_uploaded_file(self.data):

            _profile = self.data.read()
            setattr(obj, name, _profile)

            if self.size_field:
                setattr(obj, self.size_field, len(_profile))

            if self.filename_field:
                setattr(obj, self.filename_field, self.data.filename)

            if self.mimetype_field:
                setattr(obj, self.mimetype_field, self.data.content_type)