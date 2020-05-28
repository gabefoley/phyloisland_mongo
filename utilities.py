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
import genome_overview
from genome_overview import models
import getGenomes



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

def create_pos_dict(regions, profile, trimmed_name, fasta_path):
    profile_path = "./tmp/" + trimmed_name + ".hmm"
    results_outpath = "./tmp/" + trimmed_name + "_results.txt"

    # Write out regions


    with open(fasta_path, 'w+') as regions_out:
        regions_out.write(regions.decode())
    # Write out profile

    with open(profile_path, 'wb') as profile_out:
        profile_out.write(profile.read())


    # Perform the hmmsearch on the regions file
    os.system('hmmsearch -o' + results_outpath + ' --domT 1 ' + profile_path + " " + fasta_path)

    seqs = read_fasta(fasta_path)

    pos_dict = get_pos_dict_from_hmm(results_outpath, seqs)

    return pos_dict

def get_pos_dict_from_hmm(path, seqs):
    """
    Read in a hmm output file and extract the positions of the hits for a given set of sequences
    :param path: path of the hmm output file
    :param seqs: SeqRecord of the sequences we want to search for
    :return: A dictionary mapping sequence name -> (start position, end position)
    """
    qresult = SearchIO.read(path, 'hmmer3-text')

    pos_dict = {}

    print(len(qresult.hsps))
    print(len(seqs))

    if len(qresult.hsps) > len(seqs):
        print("ERROR: More HSPs than expected")

    for hsp in qresult.hsps:

        pos_dict[hsp.hit.id] = (hsp.hit_start, hsp.hit_end)

    return pos_dict

def trim_to_profile(regions, profile, trimmed_name):
    fasta_path = "./tmp/" + trimmed_name + ".fasta"
    trimmed_path = "./tmp/" + trimmed_name + "_trimmed"

    pos_dict = create_pos_dict(regions, profile, trimmed_name, fasta_path)

    seqs = read_fasta(fasta_path)

    trimmed = []

    failed_seqs = []



    for name, seq in seqs.items():
        if name in pos_dict:
            trimmed_seq = seq.seq.tomutable()
            trimmed_seq = trimmed_seq[int(pos_dict[name][0]):int(pos_dict[name][1])]
            trimmed.append(SeqRecord(trimmed_seq, seq.name))

            print (trimmed_seq)
        else:
            failed_seqs.append(name)


    # Write the newly trimmed file to disk
    createFasta(trimmed, trimmed_path, align=False)

    # Load and save the newly trimmed file
    with open(trimmed_path + ".fasta", 'rb') as trimmed_seqs:
        # Load files into database
        region = models.RegionRecords(name=trimmed_name, regions=trimmed_seqs.read())
        region.save()

    # Return sequences that failed to have a hmmer match
    return failed_seqs

def trim_around_profile(regions, profile1, profile2, pos1, pos2, trimmed_name):
    fasta_path = "./tmp/" + trimmed_name + ".fasta"
    trimmed_path = "./tmp/" + trimmed_name + "_trimmed"

    if profile1:
        pos_dict1 = create_pos_dict(regions, profile1, trimmed_name, fasta_path)
    else:
        pos_dict1 = None

    if profile2:
        pos_dict2 = create_pos_dict(regions, profile2, trimmed_name, fasta_path)
    else:
        pos_dict2 = None

    seqs = read_fasta(fasta_path)
    trimmed = []
    failed_seqs = []

    trimmed, failed_seqs = trim_sequence(seqs, pos_dict1, pos_dict2, pos1, pos2)


    # Write the newly trimmed file to disk
    createFasta(trimmed, trimmed_path, align=False)

    # Load and save the newly trimmed file
    with open(trimmed_path + ".fasta", 'rb') as trimmed_seqs:
        # Load files into database
        region = models.RegionRecords(name=trimmed_name, regions=trimmed_seqs.read())
        region.save()

    # Return sequences that failed to have a hmmer match
    return failed_seqs




def trim_sequence(seqs, pos_dict1=None, pos_dict2=None, pos1='start', pos2='start'):

    trimmed = []
    failed_seqs = []

    for name, seq in seqs.items():

        if pos_dict1==None:
            if pos_dict2 == None:
                raise NameError("Must provide at least one position dictionary")

            else: # From the start of a sequence to a profile match



                if name in pos_dict2:

                    print("From the start of sequence to profile match")

                    print(pos2)


                    trimmed_seq = seq.seq.tomutable()

                    trimmed_seq = trimmed_seq[:int(pos_dict2[name][0 if
                    pos2 == 'start' else 1])]

                    if trimmed_seq:
                        trimmed.append(SeqRecord(trimmed_seq, seq.name))
                    else:
                        failed_seqs.append(name)

                else:
                    failed_seqs.append(name)


        elif pos_dict2 == None: # From a profile match to the end of a sequence
            if name in pos_dict1:
                trimmed_seq = seq.seq.tomutable()

                print ("From a profile match to the end of the sequence")

                print (pos1)

                trimmed_seq = trimmed_seq[int(pos_dict1[name][0 if pos1 == 'start' else 1]):]
                if trimmed_seq:
                    trimmed.append(SeqRecord(trimmed_seq, seq.name))
                else:
                    failed_seqs.append(name)
            else:
                failed_seqs.append(name)

        else: # Between two profile matches
            if name in pos_dict1 and name in pos_dict2:

                print ("Between two sequences")

                print (pos1)

                print (pos2)
                trimmed_seq = seq.seq.tomutable()

                trimmed_seq = trimmed_seq[int(pos_dict1[name][0 if pos1 == 'start' else 1]):int(pos_dict2[name][ 0 if
                pos2 == 'start' else 1])]
                if trimmed_seq:
                    trimmed.append(SeqRecord(trimmed_seq, seq.name))
                else:
                    failed_seqs.append(name)
            else:
                failed_seqs.append(name)

    return trimmed, failed_seqs




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

def test_auto_classify():
    original_classifications = {"NZ_WRXN01000057.1": "Type3", "NZ_QODI01000091.1": "Type1", "NZ_QAIK01000198.1":
        "Type2b",
                       "NZ_CP009451.1": "Type2b", "NC_015513.1": "Type3", "NC_013216.1": "Type3",
                       "NC_022663.1": "Type3", "NZ_CP004078.1": "Type3", "NZ_CP011809.2": "Type2b",
                       "NZ_VITD01000034.1": "Type3", "NZ_WSEM01000042.1": "Type3", "NZ_SOEK01000059.1": "Type2b",
                       "NZ_MQWC01000010.1": "Type3", "NZ_PVZA01000046.1": "Type2b", "NZ_BIFQ01000002.1": "Type3",
                       "NZ_QEOF01000027.1": "Type2b", "NZ_NJAK01000005.1": "Type2a", "NZ_FPBP01000034.1": "Type3",
                       "NZ_BBMZ01000055.1": "Type2b", "NZ_KI632511.1": "Type3", "NZ_CP027760.1": "Type2b",
                       "NZ_CP024793.1": "Type3", "NZ_FPJC01000063.1": "Type2b", "NZ_CP027734.1": "Type2b",
                       "NZ_SAXA01000055.1": "Type3", "NZ_NMRE01000216.1": "Type2b", "NZ_CP041692.1": "Type3",
                       "NC_008271.1": "Type3", "NZ_VCNA01000026.1": "Type3", "NZ_NBVR01000044.1": "Type1",
                       "NZ_QVIG01000004.1": "Type1", "NZ_BILZ01000167.1": "Type3", "NZ_FONE01000126.1": "Type3",
                       "NZ_WIVQ01000359.1": "Type2b", "NZ_UGTQ01000009.1": "Type1", "NZ_FOLC01000061.1": "Type1",
                       "NZ_SJOP01000106.1": "Type1", "JAAHIE010001451.1": "Type3", "NZ_KL647038.1": "Type3",
                       "NZ_WIBD01000081.1": "Type2b", "NZ_MWTQ01000158.1": "Type2b", "NZ_PHSU01000082.1": "Type2b",
                       "NZ_FUXU01000241.1": "Type1", "NZ_LT629762.1": "Type2b ", "NZ_SMSC01000216.1": "Type1",
                       "NZ_CP013429.1": "Type1", "NZ_AYLO01000184.1": "Type3", "NZ_KB905728.1": "Type3",
                       "NZNV01000041.1": "Type3", "NZ_POEF01000127.1": "Type3", "NZ_NVPT01000374.1": "Type3",
                       "NZ_KN173624.1": "Type3", "NZ_QTUF01000034.1": "Type2b", "NZ_KN266223.1": "Type3",
                       "NZ_PENW01000035.1": "Type1", "NZ_QAOO01000085.1": "Type3", "NZ_LT707062.1": "Type2b",
                       "NZ_SMKX01000348.1": "Type3", "NZ_LT855380.1": "Type1", "NZ_RQJP01000018.1": "Type3",
                       "NZ_JNYY01000082.1": "Type3", "NZ_CZQA01000015.1": "Type3", "NZ_AKJS01000210.1": "Type2b",
                       "NZ_LOYJ01000133.1": "Type2b", "NZ_BBCC01000558.1": "Type1", "NZ_OAOQ01000062.1": "Type3",
                       "MUGG01000223.1": "Type3", "NZ_CP014226.1": "Type3", "NZ_RZUM01000202.1": "Type3",
                       "NZ_CP027727.1": "Type2b", "NZ_NHML01000109.1": "Type3", "NZ_VCSG01000349.1": "Type3",
                       "NZ_AFWT01000141.1": "Type3", "NC_010830.1": "Type3", "RHGL01000187.1": "Type3",
                       "NZ_QUMQ01000001.1": "Type1", "NZ_VOBI01000057.1": "Type2b", "NZ_FNTY01000002.1": "Type2b",
                       "NZ_JYLF01000032.1": "Type2b", "NZ_CP029843.1": "Type3", "NZ_LT629732.1": "Type3",
                       "NZ_CP017687.1": "Type2b", "NZ_ONZJ01000003.1": "Type3", "NZ_AKJH01000183.1": "Type2b",
                       "NZ_SMKT01000421.1": "Type3", "NZ_RCOE01000307.1": "Type2b", "NZ_RCOD01000106.1": "Type2b",
                       "NZ_SSWO01000061.1": "Type1", "NZ_CP033700.1": "Type2b", "NZ_CP015381.1": "Type3",
                       "NZ_WIAO01000088.1": "Type3", "NZ_QARE02000057.1": "Type2b", "NZ_MASS01000135.1": "Type3",
                       "NC_020453.1": "Type1", "NC_021084.1": "Type3", "NZ_WTCR01000039.1": "Type3",
                       "NZ_CP029710.1": "Type3", "NZ_CP011521.2": "Type2b", "NZ_KI519434.1": "Type3",
                       "NZ_VHKL01000033.1": "Type3", "NZ_SNYU01000018.1": "Type2b", "LKCF01003901.1": "Type2b",
                       "NZ_JH941054.1": "Type3", "NZ_AALD02000179.1": "Type1", "NZ_KB903530.1": "Type3",
                       "NZ_FNNQ01000042.1": "Type3", "NZ_WIVR01000272.1": "Type2b", "NZ_LIUV01000061.1": "Type2b",
                       "NZ_AP020337.1": "Type2b", "NZ_KE386823.1": "Type3", "NZ_RAVY01000454.1": "Type3",
                       "NZ_RAVX01000101.1": "Type3", "NZ_CAEB01000078.1": "Type1", "NZ_MJMK01000061.1": "Type2b",
                       "NZ_BAXG01000116.1": "Type1", "NZ_RXLQ01000081.1": "Type3", "NZ_QFXK01000092.1": "Type3",
                       "MNDS01000223.1": "Type3", "NZ_LECZ01000027.1": "Type1", "NZ_CP024866.1": "Type2b",
                       "NZ_RCNX01000105.1": "Type2b", "NZ_FNBR01000014.1": "Type2b", "DLUT01000591.1": "Type1",
                       "NZ_QAIL01000624.1": "Type3", "NZ_AKJU01000280.1": "Type3", "NZ_KQ058904.1": "Type3",
                       "LADO01000066.1": "Type3", "NZ_RHHV01000044.1": "Type3", "NZ_JAABNE010000136.1": "Type1",
                       "NZ_MKMC01000100.1": "Type3", "NZ_NJAJ01000180.1": "Type2b", "NZ_FZPH01000044.1": "Type3",
                       "NZ_PPRY01000171.1": "Type2b", "NZ_CP022411.1": "Type2b", "NZ_AMZY02000039.1": "Type3",
                       "NZ_KK211074.1": "Type2b", "NZ_VEJO01000043.1": "Type2b", "NZ_CP018049.1": "Type2b",
                       "NZ_LGRC01000069.1": "Type3", "NZ_FMCR01000011.1": "Type3", "DPJM01000861.1": "Type3",
                       "NC_008027.1": "Type2b", "NC_017448.1": "Type3", "NZ_AZNV01000101.1": "Type3",
                       "NZ_CP028042.1": "Type2b", "NZ_AP017313.1": "Type3", "NZ_QUMR01000030.1": "Type2b",
                       "NZ_FODL01000038.1": "Type2b", "QJPH01000597.1": "Type3", "NZ_FORG01000094.1": "Type2b",
                       "NZ_CP007231.1": "Type2b", "NZ_NEVM01000005.1": "Type1", "NZ_CP009289.1": "Type3",
                       "NZ_VDFY01000321.1": "Type3", "NZ_BJMN01000147.1": "Type3", "NZ_NIBS01000173.1": "Type2a",
                       "JMHR01000441.1": "Type3", "JAABOU010001898.1": "Type3", "NZ_VIUK01000070.1": "Type1",
                       "JAAHFV010000687.1": "Type3", "NZ_PVZG01000092.1": "Type3", "NZ_KI911557.1": "Type3",
                       "NZ_JOIX01000228.1": "Type3", "JEMY01000085.1": "Type3", "NZ_KI421497.1": "Type3",
                       "NZ_CP045011.1": "Type1", "NZ_SSNI01000125.1": "Type3", "NZ_SZWE01000003.1": "Type3",
                       "NZ_CDSC02000515.1": "Type3", "NZ_FMYF01000036.1": "Type3", "NZ_RCWL01000032.1": "Type3",
                       "NZ_PHHE01000001.1": "Type2b", "NZ_LKBY01000178.1": "Type3", "NZ_AP018150.1": "Type2b",
                       "NZ_FMXV01000075.1": "Type2b", "NZ_KB897775.1": "Type3", "NZ_PENZ01000052.1": "Type1",
                       "NZ_NIRH01000051.1": "Type1", "NZ_PENX01000027.1": "Type3", "NZ_CP047651.1": "Type2b",
                       "NZ_SNXZ01000022.1": "Type3", "NC_017565.1": "Type1", "NZ_AP018449.1": "Type3",
                       "NZ_MSSW01000136.1": "Type3", "NZ_LR134373.1": "Type2b", "NZ_CP038274.1": "Type3",
                       "NZ_ASSC01000896.1": "Type3", "NZ_FOSU01000047.1": "Type3", "NZ_LAIJ01000019.1": "Type3",
                       "NZ_FUYT01000034.1": "Type2b", "NZ_MBLO01000280.1": "Type3", "NZ_KE384514.1": "Type3",
                       "NZ_CP009747.1": "Type2b", "NZ_QGSY01000345.1": "Type3", "NZ_QGSZ01000408.1": "Type3",
                       "NZ_JH725405.1": "Type3", "NZ_OUNR01000022.1": "Type3", "NZ_CP005927.1": "Type1",
                       "NZ_CP043925.1": "Type1", "NZ_ASRX01000182.1": "Type1", "NZ_BAHC01000261.1": "Type3",
                       "NZ_PVTJ01000022.1": "Type3", "NZ_LR590468.1": "Type3", "NZ_LLWH01000241.1": "Type2b",
                       "NC_017447.1": "Type1", "NZ_MWPQ01000095.1": "Type3", "NZ_RHQN01000027.1": "Type2b",
                       "NZ_QAOQ01000022.1": "Type3", "NZ_VUAZ01000259.1": "Type2b", "NZ_CP033931.1": "Type3",
                       "NZ_LS999839.1": "Type3", "NZ_KI421431.1": "Type3", "NZ_CP027759.1": "Type2b",
                       "NZ_QTUH01000043.1": "Type2b", "NZ_BCBA01000109.1": "Type2b", "NC_015559.1": "Type3",
                       "NZ_AKJQ01000091.1": "Type2b", "NZ_QTPO01000204.1": "Type3", "NZ_PDUD01000143.1": "Type3",
                       "JAAAKW010000055.1": "Type2b", "NZ_WSTC01000100.1": "Type3", "NZ_SOCG01000010.1": "Type2b",
                       "NZ_QTPW01000119.1": "Type3", "NZ_JAABMA010000050.1": "Type1", "NZ_LT629778.1": "Type2b",
                       "NZ_PJBP01000186.1": "Type3", "NZ_QFRW01000331.1": "Type3", "NZ_MKQR01000032.1": "Type3",
                       "NZ_CP046054.1": "Type3", "NZ_UPHT01000230.1": "Type3", "NZ_MCHY01000013.1": "Type3",
                       "NZ_MUNY01000094.1": "Type3", "NZ_RAWG01000802.1": "Type3", "NZ_JYLO01000042.1": "Type2b",
                       "NZ_WSQA01000032.1": "Type3", "NZ_CP028923.1": "Type3", "NZ_MUBJ01000149.1": "Type2a",
                       "NZ_JXRA01000201.1": "Type2b", "NZ_JYLD01000037.1": "Type3", "NZ_BDBY01000492.1": "Type3",
                       "NZ_LIPP01000561.1": "Type3", "NZ_AHAM01000375.1": "Type3", "NZ_BILY01000094.1": "Type3",
                       "NZ_SODV01000004.1": "Type3", "NZ_FYEA01000033.1": "Type2b", "NZ_WMBA01000144.1": "Type3",
                       "NZ_FNCO01000054.1": "Type2b", "NZ_FMUL01000047.1": "Type2b", "NZ_FCON02000657.1": "Type3",
                       "NZ_CP023969.1": "Type2b", "NZ_JJML01000118.1": "Type3", "NZ_JAABLV010000025.1": "Type1",
                       "NZ_JAABLU010000039.1": "Type1", "NZ_FORB01000034.1": "Type3", "NZ_JYLH01000043.1": "Type2b",
                       "NZ_PGGO01000087.1": "Type3", "NZ_LMGQ01000029.1": "Type2b", "NZ_JAABNK010000025.1": "Type1",
                       "NZ_KZ679081.1": "Type3", "NKIG01000124.1": "Type3", "NZ_LMGK01000026.1": "Type2b",
                       "WASQ01000153.1": "Type3", "NZ_BAOS01000047.1": "Type3", "NZ_BCQP01000133.1": "Type3",
                       "NZ_CP010898.2": "Type2b", "NC_021184.1": "Type3", "NZ_FOZR01000085.1": "Type3",
                       "NC_009253.1": "Type3", "NZ_QKLY01000024.1": "Type2b", "NZ_LVYD01000134.1": "Type3",
                       "NZ_VFIO01000040.1": "Type2b", "QQTZ01000066.1": "Type3", "NZ_FNJL01000093.1": "Type3",
                       "NZ_MVHE01000525.1": "Type3", "NZ_FMWY01000062.1": "Type3", "NZ_CP010408.1": "Type3",
                       "NZ_LT605205.1": "Type3", "LKBL01002861.1": "Type1", "NZ_KB944506.1": "Type3",
                       "NZ_CP029618.1": "Type3", "NZ_FPBO01000103.1": "Type3", "NZ_QJUG01000493.1": "Type3",
                       "NZ_QAJM01000070.1": "Type2b", "LGGF01000107.1": "Type3", "NZ_WUNA01000095.1": "Type3",
                       "NZ_MKCS01000005.1": "Type3", "NZ_CM002331.1": "Type2b", "NZ_CP023695.1": "Type3",
                       "NZ_SMFY01000011.1": "Type3", "NZ_QSNX01000060.1": "Type2b", "NZ_LMXH01000018.1": "Type3",
                       "NZ_CP014135.1": "Type2b", "NZ_JXDG01000126.1": "Type2b", "NZ_PIQI01000031.1": "Type1",
                       "NZ_QKTW01000033.1": "Type3", "NZ_LOWA01000060.1": "Type2b", "NZ_CP042382.1": "Type3",
                       "NZ_CP013461.1": "Type3", "NZ_LWBP01000264.1": "Type3", "NZ_LYRP01000050.1": "Type1",
                       "NZ_SEIT01000119.1": "Type2b", "NZ_JYLN01000037.1": "Type2b", "NZ_QTTH01000050.1": "Type2b",
                       "NZ_NITZ01000138.1": "Type1", "NZ_VOLC01000068.1": "Type3", "NZ_LFCV01000266.1": "Type1",
                       "NZ_MVIF01000387.1": "Type3", "NZ_VRLS01000101.1": "Type3", "NZ_MULM01000185.1": "Type3",
                       "JAABLX010000056.1": "Type1", "NC_016905.1": "Type1", "NZ_SOZA01000107.1": "Type2b",
                       "NZ_MJML01000057.1": "Type2b", "NZ_WTYM01000063.1": "Type3", "NZ_QOIO01000184.1": "Type3",
                       "NZ_FQUQ01000022.1": "Type3", "NZ_FAOZ01000083.1": "Type3", "NZ_JNZS01000088.1": "Type3",
                       "NZ_KQ257877.1": "Type3", "NZ_KB892704.1": "Type3", "AP018273.1": "Type3",
                       "NZ_LT985385.1": "Type3", "NZ_PYAC01000043.1": "Type3", "NZ_FOBB01000023.1": "Type3",
                       "NZ_JPMW01000007.1": "Type3", "NZ_CPYD01000045.1": "Type2a", "NZ_CP021135.1": "Type2b",
                       "MNDA01000671.1": "Type3", "NZ_QXQA01000045.1": "Type3", "NZ_OCSV01000008.1": "Type3",
                       "NZ_FXWP01000029.1": "Type1", "NZ_AWXZ01000044.1": "Type3", "NZ_UYJA01000022.1": "Type2b",
                       "NZ_LNTU01000041.1": "Type3", "NZ_QNVV01000061.1": "Type3", "NZ_CP022121.1": "Type3",
                       "NZ_SAIQ01000015.1": "Type3", "NZ_VCRA01000094.1": "Type3", "NZ_CP029197.1": "Type3",
                       "NZ_PODL01000171.1": "Type2b", "NZ_FOAF01000023.1": "Type3", "NZ_QKWJ01000248.1": "Type3",
                       "NZ_CP029608.1": "Type2b", "NZ_JFHN01000075.1": "Type1", "NZ_FXBM01000005.1": "Type3",
                       "NZ_CP048209.1": "Type3", "NZ_VJZE01000939.1": "Type3", "NC_013131.1": "Type3",
                       "NZ_JH651384.1": "Type3", "NZ_PYAW01000034.1": "Type3", "NZ_WMJZ01000174.1": "Type1",
                       "NZ_SNZP01000034.1": "Type3", "NZ_CP010896.1": "Type2b", "NZ_SMJU01000044.1": "Type3",
                       "NZ_FAOS01000004.1": "Type3", "NZ_RHLK01000044.1": "Type3", "NZ_VSFF01000027.1": "Type3",
                       "NZ_RQPI01000039.1": "Type3", "NZ_FSRS01000002.1": "Type3", "NZ_QLTF01000036.1": "Type2b",
                       "NZ_FNKR01000003.1": "Type3", "NZ_SMSL01000028.1": "Type3", "NZ_VZZK01000111.1": "Type3",
                       "NZ_CP028272.1": "Type1", "NZ_VDCQ01000167.1": "Type3", "NZ_OGTP01000072.1": "Type3",
                       "NZ_MPIN01000042.1": "Type3", "NZ_CDPK01000072.1": "Type1", "NZ_CP026364.1": "Type1",
                       "NZ_LXEN01000293.1": "Type3", "NZ_CABPSP010000048.1": "Type2b", "NZ_CP019686.1": "Type1",
                       "NZ_SJSL01000015.1": "Type3", "CABPSQ010000036.1": "Type2b", "JAAHFO010001461.1": "Type3",
                       "NZ_SOCQ01000042.1": "Type2b", "MNJJ01000259.1": "Type3", "NZ_CP012159.1": "Type3",
                       "NZ_QLTJ01000050.1": "Type2b", "NZ_JNWO01000211.1": "Type3", "NZ_CP013341.1": "Type3",
                       "NC_017807.1": "Type2b", "NZ_PYBV01000203.1": "Type3", "NZ_KE332397.1": "Type3",
                       "NZ_RCSU01000043.1": "Type3", "NZ_QLLL01000025.1": "Type3", "NZ_QTUB01000001.1": "Type2a",
                       "NZ_JAAGLX010000268.1": "Type3", "NZ_FODH01000041.1": "Type3", "NZ_FNVU01000044.1": "Type3",
                       "NZ_FXAS01000142.1": "Type1", "NZ_CP038255.1": "Type3", "NZ_QVNU01000027.1": "Type3",
                       "NZ_VOIW01000021.1": "Type2b", "NZ_LT629795.1": "Type1", "NZ_CP026110.1": "Type3",
                       "NZ_AEDD01000040.1": "Type3", "NZ_FNAD01000036.1": "Type3", "NZ_PVZV01000026.1": "Type3",
                       "NZ_PVNL01000192.1": "Type3", "NZ_LXYR01000213.1": "Type3", "NZ_LN623556.1": "Type2b",
                       "NZ_FNGF01000016.1": "Type3", "NZ_CP022961.1": "Type3", "NZ_CXPG01000027.1": "Type3",
                       "NZ_CP017236.1": "Type1", "NZ_CP012332.1": "Type3", "NZ_NIRS01000013.1": "Type2b",
                       "NZ_FMDM01000037.1": "Type3", "NZ_PTJB01000052.1": "Type3", "NZ_JAAFZB010000096.1": "Type3",
                       "NZ_CP016211.1": "Type3", "NZ_PQKR01000051.1": "Type2b", "NZ_NCXP01000142.1": "Type3",
                       "NZ_ANMG01000154.1": "Type3", "NC_020209.1": "Type2b", "NZ_JZSQ01000140.1": "Type3",
                       "NZ_LT629705.1": "Type2b", "NZ_PYAL01000013.1": "Type2b", "MKSF01000039.1": "Type3",
                       "LAQJ01000315.1": "Type3", "NZ_SMKU01000805.1": "Type3", "NZ_LT828648.1": "Type3",
                       "NZ_CP007215.2": "Type3", "NZ_WNKZ01000359.1": "Type3", "NZ_LR590482.1": "Type2b",
                       "NZ_LT907981.1": "Type3", "NZ_QAIP01000588.1": "Type1", "NZ_LNCD01000152.1": "Type3",
                       "NZ_KE384562.1": "Type3", "NZ_ATXB01000005.1": "Type3", "NZ_SMKK01000563.1": "Type3",
                       "NC_019762.1": "Type3", "NZ_JOGP01000180.1": "Type3", "KZ266893.1": "Type3",
                       "NZ_FNON01000025.1": "Type3", "NZ_SHKK01000001.1": "Type3", "NZ_FNUD01000002.1": "Type1",
                       "NZ_FQYP01000028.1": "Type3", "NZ_QGTQ01000078.1": "Type1", "NZ_JFJW01000247.1": "Type1",
                       "NZ_FOVS01000095.1": "Type2b", "NZ_CP012540.1": "Type3", "NZ_JUHO01000001.1": "Type2b",
                       "DNUG01000139.1": "Type3", "NZ_CP038630.1": "Type3", "NC_013954.1": "Type1",
                       "NZ_FXYF01000056.1": "Type3"}



    reclassify = True


    if reclassify:
        print ("Deleting all tags")
        delete_all_tags()
        all_genomes = models.GenomeRecords.objects()
        getGenomes.tag_as_simple(all_genomes, "hidden")
        queries = models.GenomeRecords.objects(tags="Simple").timeout(False)
        print ("Classifying the genomes")
        genome_overview.classify_genomes(queries)

    else:

        queries = models.GenomeRecords.objects(tags="Simple").timeout(False)


    count = 0
    diff_count = 0
    for query in queries:
        if query.name in original_classifications:
            count += 1

            new_classification = models.GenomeTags.objects().get(tag_id=query.name)


            if original_classifications[query.name] != new_classification.tags[0].split("Auto_")[1]:
                diff_count += 1
                print ("DIFFERENT")
                print(query.name)
                print("Original was " + original_classifications[query.name])
                print("Automatic classification was " + new_classification.tags[0].split("Auto_")[1])

    print ("Wrong: " + str(diff_count))
    print ("Correct " + str(count - diff_count))
    print ("Total " + str(count))




def delete_all_tags():
    queries = models.GenomeRecords.objects().all()

    queries.update(tags=[])

    for query in queries:

        for hit in query.hits:
            hit.tags = []

        query.save()

    models.GenomeTags.objects().all().delete()
