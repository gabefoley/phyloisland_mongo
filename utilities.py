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
import shutil
import glob
import tree_code
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO, SearchIO
from Bio.Align.Applications import MuscleCommandline
import pickle

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

def remove_folder(*args):

    for arg in args:
        shutil.rmtree(arg)

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
                                              sequence=sequence, tags = [])
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
    # Version that gets called from Download FASTA
    muscle_cline = MuscleCommandline(input=input)
    # result = subprocess.run(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) )
    child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, shell=(sys.platform != "win32"))
    child.wait()

    alignment = AlignIO.read(child.stdout, "fasta")
    AlignIO.write(alignment, output, "fasta")


def search_regions_with_profiles(region_to_search, profile_ids):

    domScore = 1
    region_dict = {}

    regions = models.RegionRecords.objects.get(name=region_to_search)


    profile_folder = f"./tmp/profiles_{regions.name}/"

    os.mkdir(profile_folder)

    fasta_path = f"{profile_folder}{regions.name}.fasta"

    with open (fasta_path, "w+") as fasta_file:
        fasta_file.write(regions.regions.decode())

    while not os.path.exists(fasta_path):
        time.sleep(1)

    if os.path.isfile(fasta_path):

        seqs = read_fasta(fasta_path)

        for name in seqs.keys():
            region_dict[name.replace(".", "***")] = {}


    for profile_id in profile_ids:
        profile = models.Profile.objects().get(id=profile_id)


        with open(profile_folder + profile.name + "_profile.hmm", 'wb') as profile_path:
            profile_path.write(profile.profile.read())

        while not os.path.exists(profile_folder + profile.name + "_profile.hmm"):
            time.sleep(1)

        subprocess.getoutput(
            f'hmmsearch -o {profile_folder}{profile.name}.output --domT {domScore} {profile_folder}{profile.name}_profile.hmm {fasta_path}')

        while not os.path.exists(profile_folder + profile.name + ".output"):
            time.sleep(1)

    region_dict, domain_dict = process_hmmer_results(region_dict, profile_folder)

    # Delete the profile folder we used
    remove_folder(profile_folder)


    return region_dict, domain_dict

def process_hmmer_results(region_dict, profile_folder):

    domain_dict = {}

    print ('in process hmmer')


    print (region_dict)



    for infile in glob.glob(profile_folder + '*.output'):

        # seq_name = infile.split(output_folder)[1].split("_profile=")[0]

        qresult = SearchIO.read(infile, 'hmmer3-text')

        print ('charlie')

        print (region_dict)

        if len(qresult.hits) > 0:

            for hsp in qresult.hsps:
                # hsp = qresult[0][i]


                print (hsp.query.id)
                start = hsp.hit_start
                end = hsp.hit_end
                print ('pinto')
                print (qresult.id)
                print (hsp.hit_id)
                print (start)
                print (end)
                region_dict[hsp.hit_id.replace(".", "***")][hsp.query.id] =  (start, end)

                if hsp.query.id not in domain_dict:
                    domain_dict[hsp.query.id] = []
                domain_dict[hsp.query.id].append(hsp.hit_id)

    print (region_dict)

    print (domain_dict)

    return region_dict, domain_dict








def make_alignment_from_regions(aln_name, region_data, tool="MAFFT"):

    # Write out region data to a FASTA file

    fasta_path = "./tmp/" + aln_name + ".fasta"
    aln_path = "./tmp/" + aln_name + ".aln"

    with open (fasta_path, "w+") as fasta_file:
        fasta_file.write(region_data)

    while not os.path.exists(fasta_path):
        time.sleep(1)

    if os.path.isfile(fasta_path):
        if tool=="MAFFT":
            print ("Aligning with MAFFT")
            stdoutdata = subprocess.getoutput(f'mafft  --reorder {fasta_path} > {aln_path}')

    if os.path.isfile(aln_path):
        return aln_path

def make_tree(alignment_name, alignment, tool):

    aln_path = "./tmp/" + alignment_name + ".aln"
    tree_path = "./tmp/" + alignment_name + ".nwk"

    with open (aln_path, "w+") as fasta_file:
        fasta_file.write(alignment)

    while not os.path.exists(aln_path):
        time.sleep(1)

    if os.path.isfile(aln_path):
        if tool=="FastTree":
            print ("Making tree with FastTree")
            stdoutdata = subprocess.getoutput(f'fasttree {aln_path} > {tree_path}')

    if os.path.isfile(tree_path):
        return tree_path


def get_tree_image(tree, tree_name, tag_dict, region_dict, region_order_dict, colour_dict, full_names,
                   collapse_on_genome_tags):

    tree_path = f"./tmp/{tree_name}.nwk"
    img_path = f"static/img/trees/{tree_name}{'_full' if full_names else ''}{'_collapse' if collapse_on_genome_tags else ''}.png"
    tag_dict_path = f"./tmp/{tree_name}_tagdict.p"
    region_dict_path = f"./tmp/{tree_name}_regiondict.p"
    region_order_dict_path = f"./tmp/{tree_name}_regionorderdict.p"

    colour_dict_path = f"./tmp/{tree_name}_colourdict.p"

    print (img_path)

    # Write out tree to file
    with open (tree_path, "w+") as tree_file:
        tree_file.write(tree)

    while not os.path.exists(tree_path):
        time.sleep(1)

    pickle_dict(tag_dict, tag_dict_path)
    pickle_dict(region_dict, region_dict_path)
    pickle_dict(region_order_dict, region_order_dict_path)
    pickle_dict(colour_dict, colour_dict_path)

    if os.path.isfile(tree_path):

        # remove_file(img_path)

        print (tree_path)

        loaded_tree = tree_code.load_tree(tree_path)

        stdoutdata = subprocess.getoutput(f"python tree_code.py -t {tree_path} -o {img_path} -td {tag_dict_path} -rd {region_dict_path} -rod {region_order_dict_path} -cd {colour_dict_path} -fn {full_names} -cgt {collapse_on_genome_tags}")

        print (stdoutdata)

        # tree_code.colour_tips(loaded_tree, tag_dict, colour_dict, region_dict, outpath=img_path,
        #                                         custom_layout=False)

        if os.path.isfile(img_path):
            return img_path

def trim_to_profile(regions, profile, trimmed_name):

    fasta_path =  "./tmp/" + trimmed_name + ".fasta"
    profile_path = "./tmp/" + trimmed_name + ".hmm"
    results_outpath = "./tmp/" + trimmed_name + "_results.txt"
    trimmed_path =  "./tmp/" + trimmed_name  + "_trimmed"

    # Write out regions


    with open(fasta_path, 'w+') as regions_out:
        regions_out.write(regions.decode())
    # Write out profile

    with open(profile_path, 'wb') as profile_out:
        profile_out.write(profile.read())

    seqs = read_fasta(fasta_path)

    # Perform the hmmsearch on the regions file
    os.system('hmmsearch -o' + results_outpath + ' --domT 1 ' + profile_path + " " + fasta_path)

    pos_dict = get_pos_dict_from_hmm(results_outpath, seqs)

    trimmed = []

    failed_seqs = []

    for name, seq in seqs.items():
        if name in pos_dict:
            print ('trimmo seq')
            trimmed_seq = seq.seq.tomutable()
            print (trimmed_seq)
            print (pos_dict[name][0])
            print (pos_dict[name][1])
            trimmed_seq = trimmed_seq[int(pos_dict[name][0]):int(pos_dict[name][1])]
            trimmed.append(SeqRecord(trimmed_seq, seq.name))

            print (trimmed_seq)
        else:
            failed_seqs.append(name)

    print ('trimmed')
    print (trimmed)

    print ('failed')
    print (failed_seqs)

    # Write the newly trimmed file to disk
    createFasta(trimmed, trimmed_path, align=False)

    # Load and save the newly trimmed file
    with open(trimmed_path + ".fasta", 'rb') as trimmed_seqs:
        # Load files into database
        region = models.RegionRecords(name=trimmed_name, regions=trimmed_seqs.read())
        region.save()

    # Return sequences that failed to have a hmmer match
    return failed_seqs

def get_pos_dict_from_hmm(path, seqs):
    qresult = SearchIO.read(path, 'hmmer3-text')

    pos_dict = {}

    print(len(qresult.hsps))
    print(len(seqs))

    if len(qresult.hsps) > len(seqs):
        print("ERROR: More HSPs than expected")

    for hsp in qresult.hsps:

        pos_dict[hsp.hit.id] = (hsp.hit_start, hsp.hit_end)

    return pos_dict

def trim_sequence(seqs, pos_dict1=None, pos_dict2=None, pos1='start', pos2='start'):

    failed_seqs = []

    if pos_dict1==None:
        if pos_dict2 == None:
            raise NameError("Must provide at least one position dictionary")
        else: # From the start of a sequence to a profile match

            if pos2 == 'start': # Don't include the profile
                pass
            elif pos2 == 'end': # Include the profile
                pass
    elif pos_dict2 == None: # From a profile match to the end of a sequence
        pass

    else: # Between two profile matches
        pass


def search_for_promoters(mismatch):
    genomes = models.GenomeRecords.objects()

    regions = []

    prom_regex = "(TTGACA.{15,25}TATAAT){s<=" + str(mismatch) + "}"

    for g in genomes:
        # if g.name == "NZ_CP010029.1":
        for hit in g.hits:
            if "expanded" in hit.region:
                print(hit.region)
                print(hit.start)
                print(hit.end)
                if hit.strand == "forward":
                    seq_content = g.sequence[int(hit.start) - 50: int(hit.start)]

                    print(seq_content)
                    regions.append(seq_content)
                    match = re.findall(prom_regex, seq_content)
                    print(match)
                    if match:
                        hit.promoter = True
                    else:
                        hit.promoter = False


                elif hit.strand == "backward":
                    seq_content = g.sequence[int(hit.end): int(hit.end) + 50]

                    rev_content = str(Seq(str(seq_content)).reverse_complement())
                    print(rev_content)
                    regions.append(rev_content)
                    match = re.findall(prom_regex, rev_content)
                    print(match)
                    if match:
                        hit.promoter = True
                    else:
                        hit.promoter = False

        g.save()


def clear_all_promoters():
    """
    Clear all the promoters for all the hits
    :return:
    """
    genomes = models.GenomeRecords.objects()

    for g in genomes:
        for hit in g.hits:
            hit.promoter = False
        g.save()



def pickle_dict(dict, outpath):
    with open(outpath, 'wb') as handle:
        pickle.dump(dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


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

def get_genome_items(genome, hits='all', hidden_type=True, show_promoters=False, show_stop_codons=False,
checked_regions=None):
    """
    Format the items in a genome correctly for a Genome Detail view
    :param self:
    :return:
    """

    items = defaultdict(list)
    stop_codon_glyphs = {}
    promoter_glyphs = {}
    hit_tags = {}
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

                    if hit.tags:
                        hit_tags[str(hit.id)] = (str(hit.region) + " " + str(hit.start) + ":" + str(hit.end) + " " + \
                                         hit.strand, hit.tags)

                    if hit.name != None:
                        hit_details = dict()
                        hit_details['id'] = count
                        hit_details['hit_id'] = str(hit.id)
                        hit_details['start'] = hit.start
                        hit_details['end'] = hit.end
                        hit_details['name'] = hit.region
                        # hit_details['strand'] = 1 if count % 2 == 0 else -1

                        # We set a fake strand here just to force TcC and TcdA1 to a different track
                        if hit.region == 'TcdA1' or hit.region == 'TcdA1_expanded' or hit.region == 'TcC' or \
                                hit.region == 'TcC_expanded':
                            hit_details['strand'] = -1
                        else:
                            hit_details['strand'] = 1


                        # But store the real strand here so we can annotate the hits correctly
                        hit_details['actual_strand'] = hit.strand

                        # Add the promoters:

                        if show_promoters:

                            if hit.promoter:

                                promoter_pos = hit.start if hit.strand == 'forward' else hit.end

                                if promoter_pos in promoter_glyphs:

                                    promoter_glyphs[promoter_pos].append( hit.region + ' promoter')

                                else:
                                    promoter_glyphs[promoter_pos] = [hit.region + ' promoter']


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

                                if show_stop_codons:

                                    if pos in stop_codons:
                                        stop_codon_glyphs[pos].append(hit.region)
                                    else:
                                        stop_codon_glyphs[pos] = [hit.region]

                            idx1 += 3
                            idx2 += 3


                        region_list.append(hit_details)



                        items[hit.region].append(hit_details)

    # print (items)

    print ('stop_codons')
    print (stop_codon_glyphs)
    tracks = build_tracks(items, stop_codon_glyphs, promoter_glyphs)

    return tracks, hit_tags, genomesize
# def get_hit_tags((hits, hits='all', hidden_type=True, checked_regions=None):):


def build_tracks(items, stop_codons, promoters):

    tracks = []

    for region_name in items:

        regions = []

        #NOTE: I'm forcing the strands all to be 1 here to visualise on the same line in the linear genome

        for region in items[region_name]:

            region_dict = {'id' : region['hit_id'], 'start' : int(region['start']), 'end' : int(region['end']),
                           'name' : region[
                'name'],
             'strand' : region['strand'], 'actual_strand' : region['actual_strand']}
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





    if stop_codons:

        stop_codon_regions = []

        count = 0

        for loc, names in stop_codons.items():
            count += 1
            name = " ".join(names)
            stop_codon_dict = {'id': count, 'bp': loc, 'type': name, 'name': name}

            stop_codon_regions.append(stop_codon_dict)



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
        'items': stop_codon_regions
        }


        tracks.append(glyph_track)

    if promoters:

        promoter_regions = []

        count = 0

        for loc, names in promoters.items():
            count += 1
            name = names
            promoter_dict = {'id': count, 'bp': loc, 'type': 'promoter', 'name': name}

            promoter_regions.append(promoter_dict)

        promoter_track = {'trackName': "track3",
                           'trackType': 'glyph',
                           'glyphType': 'diamond',
                           'radius': 155,
                           'pixel_spacing': 5,
                           'linear_pixel_spacing': 5,
                           'glyph_buffer': 5,
                           'linear_glyph_buffer': 5,
                           'glyphSize': 40,
                           'linear_glyphSize': 40,
                           'linear_height': 100,
                           'showTooltip': 'true',
                           'items': promoter_regions
                           }

        tracks.append(promoter_track)

    json_tracks = json.dumps(tracks)

    return json_tracks

def get_associated_dict(genome):

    associated_dict = {}

    print ('assoc')

    print (genome.id)

    # assoc_hits = models.AssociatedHits.objects().get(genome_id=str(genome.id))


    aggregate = models.AssociatedHits._get_collection().aggregate([
        {"$match": {
            "genome_id": str(genome.id)
        }},
    ])

    assoc_hits = list(aggregate)


    # for as in aggregate:


    print (assoc_hits)

    for hit in assoc_hits:
        print (hit['_id'])
        print (hit['region1'])
        associated_dict[str(hit['_id'])] = (hit['region1'].split("region_")[1] + " and " + hit[
            'region2'].split("region_")[1])

    return associated_dict

def open_file(filename):

    print ('file name is ')
    print (filename)
    file_path = "static/uploads/" + filename

    print ('file path is')
    print(file_path)
    while not os.path.exists(file_path):
        time.sleep(1)
    if os.path.isfile(file_path):
        file = open(file_path, 'rb')
        return file


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

