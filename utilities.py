import phyloisland
import models
import checkForFeature
import os
from flask import flash
import random
import time
import subprocess
import sys
from collections import defaultdict
import json


from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

region_name_mapper = {"A1": "track1", "A2": "track2", "Chitinase": "track1", "TcdA1": "track4",
                  "TcB": "track5", "TcC": "track5", "A1_expanded" : "track1_expanded",  "Chitinase_expanded":
                          "track3_expanded",
                      "TcC_expanded": "track5_expanded"}
def read_fasta(filename):
    """
    Read in a FASTA file
    :param filename:
    :return: Dictionary object containing a SeqRecord
    """
    return SeqIO.to_dict(SeqIO.parse(filename, "fasta"))


def readLinesFromFile(filepath):
    """
    Takes a file and reads each individual line into a set
    :param filepath: Path of the file
    :return: Set containing lines from the file
    """

    content = set()

    with open(filepath, 'r') as query_file:
        for line in query_file:
            if len(line) > 1:
                content.add(line.strip())
    return content


def remove_file(*args):
    """
    Remove files in the list from the directory

    :param args: Files to remove
    :return:
    """
    for arg in args:
        os.remove(arg)


def add_genome(genome_results):
    """
    Add a genome into the database
    :param genome_results:
    """
    for record in genome_results:

        current = genome_results[record]
        if type(current) == SeqRecord:
            name = current.id
            species = " ".join(current.annotations.get('organism').split()[0:2])
            strain = current.annotations['source']
            sequence = str(current.seq)
            description = current.description

            # Check to see if the genome record already exists
            if models.GenomeRecords.objects(name=name):
                print("The genome record - %s from species - %s already exists in the database" % (name, species))
                continue

            else:
                print("Adding the genome record - %s from species - %s to the genome database" % (name, species))

                genome = models.GenomeRecords(name=name, species=species, strain=strain, description=description,
                                              sequence=sequence, tags = ['first', 'second'])
                genome.save()



def addSequence(seq_records):
    """
    Add a sequence into the database
    :param seq_records:
    """
    for record in seq_records.values():
        seq_name = record.id
        seq_description = record.description.split(">")[0]
        seq_species = seq_description.split("[")[1].split("]")[0]
        seq_sequence = str(record.seq)

        # Check if the sequence record already exists
        if models.SequenceRecords.objects(name=seq_name):
            print('Sequence with ID - %s from species - %s already exists in the sequence database' % (
                seq_name, seq_species) + "\n")
            flash('Sequence with ID - %s from species - %s already exists in the sequence database' % (
                seq_name, seq_species) + "\n")


        else:
            print('Adding sequence with ID - %s from species - %s to the sequence database' % (
                seq_name, seq_species) + "\n")

            sequence = models.SequenceRecords(name=seq_name, species=seq_species, description=seq_description,
                                              sequence = seq_sequence)
            sequence.save()

def save_profile(profile):
    """
    Save a profile into the database
    :param profile: Profile to save
    """
    print ('in save profile')
    name = randstring(5)

    # print (str(profile))
    # print (profile)
    print (profile.read())
    profile_entry = models.Profile(name, profile, {})
    profile_entry.save()


    # print ('in save profile')
    # name = randstring(5)
    # blobMix = models.BlobMixin("application/octet-stream", name, profile.read(), '566666')
    # profileEntry = models.Profile(name)
    # profileEntry.set_blobMix(blobMix)
    # phyloisland.db.session.add(profileEntry)
    # phyloisland.db.session.commit()

def set_profile_as_reference(profile_ids, region):
    """

    :param ids:
    :param region:
    :return:
    """
    if len(profile_ids) > 1:
        flash('Only select a single record', category='error')
    else:

        # Check for a previous Profile set as this reference
        prev = models.Profile.objects(__raw__={"references.%s" % (region):{'$exists': True}})

        print ('prev is - ', prev)

        if prev:
            print ('found prev')
            prev.update(**{"unset__references__%s" % (region): "*"})

        profile_id = profile_ids[0]

        curr = models.Profile.objects().get(id=profile_id)

        # curr = models.Profile.objects(id=profile_id)


        # curr.update(references= {region: "*"})

        curr.update(**{"set__references__%s" % (region) : "*"})

        curr.save()

        print (curr.profile)

        print ('*****')

        print (curr.profile.read())


        # Write the new profile to the tmp folder ready to be used
        # with open("tmp/" + region + "_profile.hmm", 'w') as profile_path:
        #     profile_path.write(curr.profile.decode('utf-8'))

        flash("The profile named %s has been set as the reference profile for %s" % (curr.name, region), category='success')

def createAlignment(input, output):
    muscle_cline = MuscleCommandline(input=input)
    # result = subprocess.run(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) )
    child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, shell=(sys.platform != "win32"))
    child.wait()

    alignment = AlignIO.read(child.stdout, "fasta")
    AlignIO.write(alignment, output, "fasta")



def createProfile(align_list):

    SeqIO.write(align_list, "tmp/align.fasta", "fasta")

    hmm_path = "tmp/profile3.hmm"

    createAlignment('tmp/align.fasta', 'tmp/align.aln')

    outfile = open(hmm_path, "w")
    result = subprocess.call(["hmmbuild", hmm_path, "tmp/align.aln"], stdout=subprocess.PIPE)

    while not os.path.exists(hmm_path):
        time.sleep(1)

    if os.path.isfile(hmm_path):
        file = open(hmm_path, 'rb')

        save_profile(file)
        remove_file(hmm_path, "tmp/align.fasta", "tmp/align.aln")


def createFasta(fasta_list, name, align):

    print ('hee ')
    print (fasta_list)
    SeqIO.write(fasta_list,  name + ".fasta", "fasta")

    if align:
        createAlignment(name + ".fasta", name + ".aln")


def check_with_profile(ids, region):
    # Check if a reference profile for this region exists
    profile_reference = models.Profile.objects(__raw__={"references.%s" % (region):{'$exists': True}})
    if (profile_reference):
        for profile in profile_reference:
            print("Using the %s profile named %s to check for %s regions" % (region, profile.name, region))

            eval(
                'checkForFeature.get_feature_location_with_profile(ids, "hmm_outputs' + '", "' + profile.name + '", "' + region + '", "' + region + '_loc' + '","' + region + '")')
    else:
        flash("Please set a profile as the %s reference profile first" % (region), "error")

def get_genome_items(genome, hits='all'):
    """
    Format the items in a genome correctly for a Genome Detail view
    :param self:
    :return:
    """

    items = defaultdict(list)
    region_list = []

    for count, hit in enumerate(genome.hits):

        if ((hits == 'all') or ((hits == 'initial') and ('expanded' not in hit.region)) or ((hits == 'expanded') and
         'expanded' in hit.region)):

            print (hit.id)
            hit_details = dict()
            hit_details['id'] = count
            hit_details['hit_id'] = str(hit.id)
            hit_details['start'] = hit.start
            hit_details['end'] = hit.end
            hit_details['name'] = hit.region
            hit_details['strand'] = 1 if count % 2 == 0 else -1

            region_list.append(hit_details)


            items[hit.region].append(hit_details)

    print (items)

    tracks = build_tracks(items)

    return tracks

def build_tracks(items):

    tracks = []

    for region_name in items:

        print ('region name')
        print (region_name)

        regions = []

        #NOTE: I'm forcing the strands all to be 1 here to visualise on the same line in the linear genome

        for region in items[region_name]:
            print ('region')
            print (region)
            region_dict = {'id' : region['hit_id'], 'start' : int(region['start']), 'end' : int(region['end']),
                           'name' : region[
                'name'],
             'strand' : 1}
            regions.append(region_dict)

        track = { 'trackName': region_name_mapper[region_name],
		'trackType': "stranded",
		'visible': 'true',
		'inner_radius': 130,
		'outer_radius': 185,
		'trackFeatures': "complex",
		'featureThreshold': 7000000,
		'mouseclick': 'linearClick',
		'mouseover_callback': 'islandPopup',
		'mouseout_callback': 'islandPopupClear',
		'linear_mouseclick': 'linearPopup',
		'showLabels': 'true',
		'showTooltip': 'true',
		'linear_mouseclick': 'linearClick',
        'items' : regions }

        tracks.append(track)
    print ('and tracks are ')
    print (tracks)

    json_tracks = json.dumps(tracks)

    print (json_tracks)

    return json_tracks

def randstring(length=10):
    valid_letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return ''.join((random.choice(valid_letters) for i in range(length)))