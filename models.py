from phyloisland import app, db
from flask_user import UserMixin, UserManager
from wtforms import ValidationError, fields
from wtforms.validators import required
from wtforms.widgets import FileInput
from werkzeug.datastructures import FileStorage
from gettext import gettext


class User(db.DynamicDocument, UserMixin):
    """
    User model
    """
    active = db.BooleanField(default=True)

    # User authentication information
    username = db.StringField(default='', unique=True)
    email = db.StringField(max_length=30)
    password = db.StringField()
    email_confirmed_at = db.DateTimeField()

    # User information
    first_name = db.StringField(default='')
    last_name = db.StringField(default='')

    # Customisable information
    page_size = db.IntField()
    record_size = db.IntField()
    references = db.ListField(db.StringField(), default=[])

    # Relationships
    roles = db.ListField(db.StringField(), default=[])

# Setup Flask-User and specify the User data-model
user_manager = UserManager(app, db, User)


class SequenceRecords(db.DynamicDocument):
    """
    Class for storing Sequence records
    """
    name = db.StringField()
    species = db.StringField()
    description = db.StringField()
    sequence = db.StringField()
    references = db.ListField(db.StringField(), default=list)



class Hits(db.EmbeddedDocument):
    """
    Class for storing a set of HMMER hits for a given region for a Genome Record or Region Record
    """
    id = db.ObjectIdField()
    name = db.StringField()
    region = db.StringField()
    score = db.StringField()
    start = db.StringField()
    end = db.StringField()
    expand = db.BooleanField()
    strand = db.StringField()
    sequence = db.StringField()
    tags = db.ListField(db.StringField(), default=list)

class AssociatedHits(db.DynamicDocument):
    genome_id = db.StringField()
    region1 = db.StringField()
    region2 = db.StringField()

class GenomeTags(db.DynamicDocument):
    tag_id = db.StringField()
    tags = db.ListField(db.StringField(), default=list)

class RegionRecords(db.DynamicDocument):
    name = db.StringField(unique=True)
    region_tags = db.ListField()
    regions = db.BinaryField()
    hits = db.EmbeddedDocumentListField(Hits)

class RegionToProfileRecords(db.DynamicDocument):
    rtp_id = db.StringField()
    region = db.StringField()
    profiles = db.ListField()
    region_dict = db.DictField()
    domain_dict = db.DictField()

class AlignmentRecords(db.DynamicDocument):
    name = db.StringField()
    alignment = db.BinaryField()

class TreeRecords(db.DynamicDocument):
    name = db.StringField()
    alignment = db.IntField()
    tree = db.BinaryField()







class GenomeRecords(db.DynamicDocument):
    """
    Class for storing Genome records
    """

    name = db.StringField()
    species = db.StringField()
    strain = db.StringField()
    description = db.StringField()
    sequence = db.StringField()
    present = db.DictField()
    hits = db.EmbeddedDocumentListField(Hits)
    references = db.ListField(db.StringField(), default=list)
    genome_overview = db.ImageField()
    genome_expanded_overview = db.FileField()
    tags = db.ListField(db.StringField(), default=list)

    meta = {
	'indexes': [
	{
	     'fields': ['+name']
	},
	{
	     'fields': ['#name']
	}]
    }






# class BlobMixin(db.DynamicDocument):
#     mimetype = db.StringField()
#     filename = db.StringField()
#     profile = db.BinaryField()
#     size = db.IntField()



# class BlobMixin(object):
#     mimetype = db.Column(db.Unicode(length=255), nullable=False)
#     filename = db.Column(db.Unicode(length=255), nullable=False)
#     profile = db.Column(db.BLOB, nullable=False)
#     size = db.Column(db.Integer, nullable=False)
#
#     def __init__(self, mimetype, filename, profile, size):
#         self.mimetype = mimetype
#         self.filename = filename
#         self.profile = profile
#         self.size = size

class Profile(db.DynamicDocument):

    name = db.StringField()
    profile = db.FileField()
    references = db.DictField(db.StringField(), default=list)


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