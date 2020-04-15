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
import regex as re
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide


from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

region_name_mapper = {"A1": "track1", "A2": "track2", "Chitinase": "track3", "TcdA1": "track4",
                  "TcB": "track5", "TcC": "track6", "A1_expanded" : "track1_expanded", 'A2_expanded' : "track2",
                      "Chitinase_expanded": "track3_expanded", "TcdA1_expanded": "track4_expanded","TcB_expanded" :
                          "track5_expanded",
                      "TcC_expanded": "track6_expanded"}
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

def save_profile(profile, name=None):
    """
    Save a profile into the database
    :param profile: Profile to save
    """
    print ('in save profile')
    if not name:
        name = randstring(5)

    print (name)

    # print (str(profile))
    # print (profile)
    # print (profile.read())
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


        if prev:
            prev.update(**{"unset__references__%s" % (region): "*"})

        profile_id = profile_ids[0]


        curr = models.Profile.objects().get(id=profile_id)



        # curr.update(references= {region: "*"})

        curr.update(**{"set__references__%s" % (region) : "*"})

        curr.save()

        # print (curr.profile.read())
        #
        # print ('*****')
        #
        # # print (curr.profile.file.read())
        #
        # print (type(curr.profile.read()))
        #
        # print (curr.profile.decode("utf-8"))
        #
        # print ('decoder')
        #
        # print (str(curr.profile.read(), 'utf-8'))



        # Write the new profile to the tmp folder ready to be used
        with open("tmp/" + region + "_profile.hmm", 'wb') as profile_path:

            profile_path.write(curr.profile.read())

        # flash("The profile named %s has been set as the reference profile for %s" % (curr.name, region), category='success')

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

    SeqIO.write(fasta_list,  name + ".fasta", "fasta")

    if align:
        print ("And now an aln")
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

def get_genome_items(genome, hits='all', hidden_type=True, checked_regions=None):
    """
    Format the items in a genome correctly for a Genome Detail view
    :param self:
    :return:
    """

    items = defaultdict(list)
    glyphs = {}
    region_list = []
    genomesize = len(genome.sequence)

    # print ('CHECKED REGIONS IS ')
    # print (checked_regions)

    # Add in the expanded versions of the checked regions
    # if checked_regions != None:
    #     for x in range(len(checked_regions)):
    #         print (x)
    #         checked_regions.append (checked_regions[x] + "_expanded")

    # print (checked_regions)

    for count, hit in enumerate(genome.hits):

        # print ('count is ')
        #
        # print (count)

        if ((hits == 'all') or ((hits == 'initial') and ('expanded' not in hit.region)) or ((hits == 'expanded') and
         'expanded' in hit.region)):

            #
            # print ('hidden type is ' + str(hidden_type))
            # print (hit.tags)



            if (checked_regions == None) or hit.region in checked_regions or hit.region.split("_expanded")[0] in checked_regions:

                if ((hidden_type == False) or (hidden_type == True and 'hidden' not in hit.tags)):

                    # print (hit.id)
                    # print (hit.name)

                    if hit.name != None:
                        hit_details = dict()
                        hit_details['id'] = count
                        hit_details['hit_id'] = str(hit.id)
                        hit_details['start'] = hit.start
                        hit_details['end'] = hit.end
                        hit_details['name'] = hit.region
                        # hit_details['strand'] = 1 if count % 2 == 0 else -1
                        if hit.region == 'TcdA1' or hit.region == 'TcdA1_expanded' or hit.region == 'TcC' or \
                                hit.region == 'TcC_expanded':
                            hit_details['strand'] = -1
                        else:
                            hit_details['strand'] = 1



                        #
                        # print ('genome length')




                        idx1 = 0
                        idx2 = idx1 + 3

                        stop_codons = ["TAG", "TAA", "TGA"]

                        if hit.strand == 'backward':
                            sequence = Seq(hit.sequence, generic_nucleotide)

                            hit_sequence = sequence.reverse_complement()
                        else:
                            hit_sequence = hit.sequence
                        #
                        # print ('flipped seq')
                        #
                        # print (hit_sequence)
                        # print ('get the sequence')

                        # print (hit_sequence)
                        #
                        # print (hit.strand)
                        #
                        # print (hit.start)
                        #
                        # print (hit.end)


                        while idx2 <= len(hit.sequence):
                            if hit_sequence[idx1:idx2] in stop_codons:
                                # print('found', idx1)
                                # print (hit_sequence)
                                # print (hit.region)
                                # print (hit_sequence[idx1:idx2 + 20])
                                #
                                # print (hit.start)
                                # print (hit.end)
                                # print (idx1)

                                if hit.strand == 'backward':
                                    pos = int(hit.end) - idx1
                                else:
                                    pos = int(hit.start) + idx1

                                # print (pos)

                                if pos in glyphs:
                                    glyphs[pos].append(hit.region)
                                else:
                                    glyphs[pos] = [hit.region]
                            # print (hit_sequence[idx1:idx2])
                            idx1 += 3
                            idx2 += 3


                        region_list.append(hit_details)



                        items[hit.region].append(hit_details)

    # print (items)

    # print ('glyphs')
    # print (glyphs)
    tracks = build_tracks(items, glyphs)

    return tracks, genomesize

def build_tracks(items, glyphs):

    tracks = []

    for region_name in items:

        # print ('region name')
        # print (region_name)

        regions = []

        #NOTE: I'm forcing the strands all to be 1 here to visualise on the same line in the linear genome

        for region in items[region_name]:
            # print ('region')
            # print (region)
            region_dict = {'id' : region['hit_id'], 'start' : int(region['start']), 'end' : int(region['end']),
                           'name' : region[
                'name'],
             'strand' : region['strand']}
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
    # print ('and tracks are ')
    # print (tracks)

    glyph_regions = []

    count = 0

    for loc, names in glyphs.items():
        count +=1
        name = " ".join(names)
        glyph_dict = {'id' : count, 'bp' : loc, 'type': name, 'name': name}

        glyph_regions.append(glyph_dict)



    if glyphs:
        glyph_track = {'trackName': "track1",
        'trackType': 'glyph',
        'glyphType': 'circle',
        'radius': 155,
        'pixel_spacing': 5,
        'linear_pixel_spacing': 5,
        'glyph_buffer': 5,
        'linear_glyph_buffer': 5,
        'glyphSize': 20,
        'linear_glyphSize': 20,
        'linear_height': 100,
        'showTooltip': 'true',
        'items': glyph_regions
        }

    # print (glyph_track)

        tracks.append(glyph_track)





    json_tracks = json.dumps(tracks)

    # print (json_tracks)

    return json_tracks

def randstring(length=10):
    valid_letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return ''.join((random.choice(valid_letters) for i in range(length)))

def sort_func(elem):
    return int(elem.split("_position=_")[1].split("_")[0])


def rename_duplicates(genome_name, old):
    seen = {}
    index = {}
    for pos in range(len(old)):
        check = old[pos]

        if check in seen:
            seen[check] += 1
            if check in index:
                index.pop(check)
        else:
            seen[check] = 1
            index[check] = pos
        old[pos] = check + "_" + str(seen[check])

    print(index)

    # Either add the genome name, or remove the temporarily annotated number from the end

    for pos in range(len(old)):
        if pos in index.values():
            old[pos] = "_".join(old[pos].split("_")[0:-1])
        else:
            old[pos] = "_".join(old[pos].split("_")[0:-1]) + "_" + genome_name + "_" + old[pos].split("_")[-1]

    return old