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
import ete3
import genome_overview
# from genome_overview import models

import getGenomes

region_name_mapper = {"A1": "track1", "A2": "track2", "Chitinase": "track3", "TcdA1": "track4",
                      "TcB": "track5", "TcC": "track6", "region1" : "track7", "A1_expanded": "track1_expanded", \
                                                                                       'A2_expanded':
                          "track2",
                      "Chitinase_expanded": "track3_expanded", "TcdA1_expanded": "track4_expanded", "TcB_expanded":
                          "track5_expanded",
                      "TcC_expanded": "track6_expanded", "region1_expanded" : "track7_expanded", "EXISTING:": "track8"}


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
            plasmid = current.annotations['plasmid']
            organism = current.annotations['organism']
            assembly_name = current.annotations['assembly_name']
            biosample = current.annotations['biosample']
            bioproject = current.annotations['bioproject']
            date = current.annotations['date']
            wgs_project = current.annotations['wgs_project']
            genome_coverage = current.annotations['genome_coverage']

            taxid = current.annotations['taxid']
            assembly_type = current.annotations['assembly_type']
            release_type = current.annotations['release_type']
            assembly_level = current.annotations['assembly_level']
            genome_representation = current.annotations['genome_representation']
            expected_final_version = current.annotations['expected_final_version']
            excluded = current.annotations['excluded']
            genbank_accession_id = current.annotations['genbank_accession_id']
            refseq_accession_id = current.annotations['refseq_accession_id']
            r_g_identical = current.annotations['r_g_identical']
            sequence = str(current.seq)
            description = current.description

            # Check to see if the genome record already exists
            if models.GenomeRecords.objects(name=name):
                print("The genome record - %s from species - %s already exists in the database" % (name, species))
                continue

            else:
                print("Adding the genome record - %s from species - %s to the genome database" % (name, species))

                genome = models.GenomeRecords(name=name, species=species, organism=organism,

                                              assembly_name=assembly_name,
                                              biosample = biosample,
                                              bioproject = bioproject,
                                              date = date,
                                              wgs_project = wgs_project,
                                              genome_coverage = genome_coverage,
                                              taxid=taxid,
                                              assembly_type=assembly_type,
                                              release_type=release_type,
                                              assembly_level=assembly_level,
                                              genome_representation=genome_representation,
                                              expected_final_version=expected_final_version,
                                              excluded=excluded,
                                              genbank_accession_id=genbank_accession_id,
                                              refseq_accession_id=refseq_accession_id,
                                              r_g_identical=r_g_identical,
                                              plasmid=plasmid, description=description,
                                              sequence=sequence, tags=[])
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
                                              sequence=seq_sequence)
            sequence.save()


def save_profile(profile, name=None):
    """
    Save a profile into the database
    :param profile: Profile to save
    """
    print('in save profile')
    if not name:
        name = randstring(5)

    print(name)

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
        prev = models.Profile.objects(__raw__={"references.%s" % (region): {'$exists': True}})

        if prev:
            prev.update(**{"unset__references__%s" % (region): "*"})

        profile_id = profile_ids[0]

        curr = models.Profile.objects().get(id=profile_id)

        # curr.update(references= {region: "*"})

        curr.update(**{"set__references__%s" % (region): "*"})

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

    # print (input)
    # Version that gets called from Download FASTA
    # muscle_cline = MuscleCommandline(input=input)
    # # result = subprocess.run(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) )
    # child = subprocess.Popen(str(muscle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    #                          universal_newlines=True, shell=(sys.platform != "win32"))

    subprocess.getoutput(f'muscle -in {input} -out {output}')
    # child.wait()


    # alignment = AlignIO.read(child.stdout, "fasta")
    #
    #
    #
    # print (alignment)
    # AlignIO.write(alignment, output, "fasta")
    #
    # while not os.path.exists(output):
    #     time.sleep(1)


def search_regions_with_profiles(region_to_search, profile_ids):
    domScore = 1
    region_dict = {}

    regions = models.RegionRecords.objects.get(name=region_to_search)

    profile_folder = f"./tmp/profiles_{regions.name}/"

    os.mkdir(profile_folder)

    fasta_path = f"{profile_folder}{regions.name}.fasta"

    with open(fasta_path, "w+") as fasta_file:
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

    print('in process hmmer')

    print(region_dict)

    for infile in glob.glob(profile_folder + '*.output'):

        # seq_name = infile.split(output_folder)[1].split("_profile=")[0]

        qresult = SearchIO.read(infile, 'hmmer3-text')

        print('charlie')

        print(region_dict)

        if len(qresult.hits) > 0:

            for hsp in qresult.hsps:
                # hsp = qresult[0][i]


                print(hsp.query.id)
                start = hsp.hit_start
                end = hsp.hit_end
                print('pinto')
                print(qresult.id)
                print(hsp.hit_id)
                print(start)
                print(end)
                if hsp.query.id in region_dict[hsp.hit_id.replace(".", "***")]:
                    region_dict[hsp.hit_id.replace(".", "***")][hsp.query.id + "_multiple_" + randstring(5)] = (start,
                                                                                                               end)

                else:
                    region_dict[hsp.hit_id.replace(".", "***")][hsp.query.id] = (start, end)

                if hsp.query.id not in domain_dict:
                    domain_dict[hsp.query.id] = []
                domain_dict[hsp.query.id].append(hsp.hit_id)

    print(region_dict)

    print(domain_dict)

    return region_dict, domain_dict


def make_alignment_from_regions(aln_name, region_data, tool="MAFFT"):
    # Write out region data to a FASTA file

    fasta_path = "./tmp/" + aln_name + ".fasta"
    aln_path = "./tmp/" + aln_name + ".aln"

    with open(fasta_path, "w+") as fasta_file:
        fasta_file.write(region_data)

    while not os.path.exists(fasta_path):
        time.sleep(1)

    if os.path.isfile(fasta_path):
        if tool == "MAFFT":
            print("Aligning with MAFFT")
            stdoutdata = subprocess.getoutput(f'mafft  --reorder {fasta_path} > {aln_path}')

    if os.path.isfile(aln_path):
        return aln_path

def load_list(*args):
    return_list = []
    for x in args:
        with open(x) as f:
            curr = [line.strip() for line in f]
        return_list += curr
    return return_list


def get_sequence_content_dict(region):
    fasta_path = "./tmp/tmp_regions.fasta"

    alns = models.AlignmentRecords.objects().get(name=region)

    # Write out the regions

    with open(fasta_path, "w+") as fasta_file:
        fasta_file.write(alns.alignment.read().decode())

    while not os.path.exists(fasta_path):
        time.sleep(1)

    # Load the regions in as a dictionary

    if os.path.isfile(fasta_path):
        seqs = read_fasta(fasta_path)

    print('pnky')
    print(seqs)

    seq_content_dict = {k: v.seq for k, v in seqs.items()}

    print(seq_content_dict)

    return seq_content_dict


def make_tree(alignment_name, alignment, tool):
    aln_path = "./tmp/" + alignment_name + ".aln"
    tree_path = "./tmp/" + alignment_name + ".nwk"

    with open(aln_path, "w+") as fasta_file:
        fasta_file.write(alignment)

    while not os.path.exists(aln_path):
        time.sleep(1)

    if os.path.isfile(aln_path):
        if tool == "FastTree":
            print("Making tree with FastTree")
            stdoutdata = subprocess.getoutput(f'fasttree -nosupport {aln_path} > {tree_path}')

    if os.path.isfile(tree_path):
        return tree_path


def get_tree_image(tree, tree_name, tag_dict, region_dict, region_order_dict, sequence_content_dict, colour_dict,
                   full_names,
                   collapse_on_genome_tags, display_circular, display_circular_180):
    tree_path = f"./tmp/{tree_name}.nwk"
    img_path = f"static/img/trees/{tree_name}{'_full' if full_names else ''}{'_collapse' if collapse_on_genome_tags else ''}{'_rd' if region_dict else ''}{'_ro' if region_order_dict else ''}{'_sc' if sequence_content_dict else ''}{'_circ' if display_circular else ''}{'_circ180' if display_circular_180 else ''}.png"
    tag_dict_path = f"./tmp/{tree_name}_tagdict.p"
    region_dict_path = f"./tmp/{tree_name}_regiondict.p"
    region_order_dict_path = f"./tmp/{tree_name}_regionorderdict.p"
    sequence_content_dict_path = f"./tmp/{tree_name}_seqcontentdict.p"
    colour_dict_path = f"./tmp/{tree_name}_colourdict.p"

    print(img_path)

    # Write out tree to file
    with open(tree_path, "w+") as tree_file:
        tree_file.write(tree)

    while not os.path.exists(tree_path):
        time.sleep(1)

    print ('here is the tag dict')
    print (tag_dict)

    # print (region_order_dict.keys())
    #
    # print (region_order_dict['QGTQ01'])
    # return

    # Override QGTQ01 here



    pickle_dict(tag_dict, tag_dict_path)
    pickle_dict(region_dict, region_dict_path)
    pickle_dict(region_order_dict, region_order_dict_path)
    pickle_dict(sequence_content_dict, sequence_content_dict_path)

    print ('here is the colour dict')
    print (colour_dict)

    pickle_dict(colour_dict, colour_dict_path)

    if os.path.isfile(tree_path):

        # remove_file(img_path)

        print(tree_path)

        loaded_tree = tree_code.load_tree(tree_path)

        stdoutdata = subprocess.getoutput(
            f'python tree_code.py -t {tree_path} -o {img_path} -td {tag_dict_path} -rd {region_dict_path} -rod {region_order_dict_path} -scd {sequence_content_dict_path} -cd {colour_dict_path} -fn {full_names} -cgt {collapse_on_genome_tags} {" -dc" if display_circular else ""} {" -dco" if display_circular_180 else ""}')

        print(stdoutdata)

        # tree_code.colour_tips(loaded_tree, tag_dict, colour_dict, region_dict, outpath=img_path,
        #                                         custom_layout=False)

        if os.path.isfile(img_path):
            return img_path


def highlight_taxonomy(tree):
    ts = ete3.TreeStyle()
    # disable default PhyloTree Layout
    ts.layout_fn = lambda x: True

    for n in tree.traverse():

        if not n.is_leaf():
            N = ete3.TextFace(n.sci_name + " (" + n.rank + ")", fsize=14, fgcolor="black")
            n.add_face(N, 12, position="branch-top")

            nstyle = ete3.NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = 'blue'
            nstyle["size"] = 10
            nstyle["hz_line_type"] = 1
            n.set_style(nstyle)

        else:
            pass

    return tree, ts


def get_species_tree_image(tree, tree_name):
    img_path = f"static/img/trees/{tree_name}.png"
    ncbi_img_path = f"static/img/trees/{tree_name}_ncbi.png"

    stree = ete3.PhyloTree(tree, sp_naming_function=lambda name: name.split('_taxid_')[1].split("_")[0])


    # Create annotated species tree image
    tax2names, tax2lineages, tax2rank = stree.annotate_ncbi_taxa()
    tcb, ts = highlight_taxonomy(stree)
    stree.render(img_path, tree_style=ts, dpi=300)



    # Create NCBI species tree image
    ncbi = ete3.NCBITaxa()
    taxids = get_taxids(stree)
    ncbi_tree = ncbi.get_topology(taxids)
    ncbi_tree, ts = highlight_taxonomy(ncbi_tree)
    ncbi_tree.render(ncbi_img_path, tree_style=ts, dpi=300)

    if os.path.isfile(img_path):
        return img_path


# Get out a list of the taxonomic IDs stored in the _taxids_ annotation on a tree
def get_taxids(tree):
    taxids = []
    for n in tree.traverse():
        taxids.append(n.taxid)
    return taxids

def get_ml_go_tree_image(tree, name, ancestral_orders, ref_ml_go_dict):
    tree_path = f"./tmp/{name}_ml.nwk"
    img_path = f"static/img/trees/{name}_ml.png"
    ancestral_orders_path = f"./tmp/{name}_ancestralorders.p"
    ref_ml_go_dict_path = f"./tmp/{name}_ref_ml_go_dict.p"

    print('hops')
    print(tree)

    # Write out tree to file
    with open(tree_path, "w+") as tree_file:
        tree_file.write(tree)

    while not os.path.exists(tree_path):
        time.sleep(1)

    pickle_dict(ancestral_orders, ancestral_orders_path)
    pickle_dict(ref_ml_go_dict, ref_ml_go_dict_path)

    if os.path.isfile(tree_path):

        print('call it')
        print(tree_path)
        print(img_path)
        print(ancestral_orders_path)
        stdoutdata = subprocess.getoutput(
            f'python tree_code.py -t {tree_path} -o {img_path} -ao {ancestral_orders_path} -mlgo {ref_ml_go_dict_path}')

        print(stdoutdata)

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

            print(trimmed_seq)
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

        if pos_dict1 == None:
            if pos_dict2 == None:
                raise NameError("Must provide at least one position dictionary")

            else:  # From the start of a sequence to a profile match



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


        elif pos_dict2 == None:  # From a profile match to the end of a sequence
            if name in pos_dict1:
                trimmed_seq = seq.seq.tomutable()

                print("From a profile match to the end of the sequence")

                print(pos1)

                trimmed_seq = trimmed_seq[int(pos_dict1[name][0 if pos1 == 'start' else 1]):]
                if trimmed_seq:
                    trimmed.append(SeqRecord(trimmed_seq, seq.name))
                else:
                    failed_seqs.append(name)
            else:
                failed_seqs.append(name)

        else:  # Between two profile matches
            if name in pos_dict1 and name in pos_dict2:

                print("Between two sequences")

                print(pos1)

                print(pos2)
                trimmed_seq = seq.seq.tomutable()

                trimmed_seq = trimmed_seq[int(pos_dict1[name][0 if pos1 == 'start' else 1]):int(pos_dict2[name][0 if
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
    SeqIO.write(fasta_list, name + ".fasta", "fasta")

    while not os.path.exists(name + ".fasta"):
        time.sleep(1)

    if align:
        # print("And now an aln")
        createAlignment(name + ".fasta", name + ".aln")


def check_with_profile(ids, region):
    # Check if a reference profile for this region exists
    profile_reference = models.Profile.objects(__raw__={"references.%s" % (region): {'$exists': True}})
    if (profile_reference):
        for profile in profile_reference:
            print("Using the %s profile named %s to check for %s regions" % (region, profile.name, region))

            eval(
                'checkForFeature.get_feature_location_with_profile(ids, "hmm_outputs' + '", "' + profile.name + '", "' + region + '", "' + region + '_loc' + '","' + region + '")')
    else:
        flash("Please set a profile as the %s reference profile first" % (region), "error")


def get_genome_items(genome, hits='all', hidden_type=True, show_promoters=False, show_stop_codons=False,
                     show_existing_features=False,
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

    # show_existing_features == bool(show_existing_features)


    # print ('CHECKED REGIONS IS ')
    # print (checked_regions)

    # Add in the expanded versions of the checked regions
    # if checked_regions != None:
    #     for x in range(len(checked_regions)):
    #         print (x)
    #         checked_regions.append (checked_regions[x] + "_expanded")

    # print (checked_regions)

    for count, hit in enumerate(genome.hits):

        if ((hits == 'all') or ((hits == 'initial') and ('expanded' not in hit.region)) or ((hits == 'expanded') and
                                                                                                    'expanded' in
                                                                                                    hit.region) or
                (show_existing_features and 'EXISTING' in hit.region)):

            #
            # print ('hidden type is ' + str(hidden_type))
            # print (hit.tags)



            if (checked_regions == None) or hit.region in checked_regions or hit.region.split("_expanded")[
                0] in checked_regions or (show_existing_features and 'EXISTING' in hit.region):

                if ((hidden_type == False) or (hidden_type == True and 'hidden' not in hit.tags)):

                    # Here is a place to update FASTA ID headers


                    if hit.tags:
                        hit_tags[str(hit.id)] = (str(hit.region) + " [" + str(hit.score) + "] " + str(hit.start) +
                                                 ":" + str(
                            hit.end) + " " + \
                                                 hit.strand, hit.tags)


                    if hit.name != None:



                        hit_details = dict()
                        hit_details['id'] = count
                        hit_details['hit_id'] = str(hit.id)
                        hit_details['start'] = hit.start
                        hit_details['end'] = hit.end
                        hit_details['name'] = hit.region
                        hit_details['score'] = hit.score
                        # hit_details['strand'] = 1 if count % 2 == 0 else -1



                        # We set a fake strand here just to force TcC and TcdA1 to a different track
                        if hit.region == 'TcdA1' or hit.region == 'TcdA1_expanded' or hit.region == 'TcC' or \
                                        hit.region == 'TcC_expanded' or hit.region.startswith("EXISTING:"):
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

                                    promoter_glyphs[promoter_pos].append(hit.region + ' promoter')

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

    print('stop_codons')
    print(stop_codon_glyphs)
    tracks = build_tracks(items, stop_codon_glyphs, promoter_glyphs, show_existing_features)

    return tracks, hit_tags, genomesize


# def get_hit_tags((hits, hits='all', hidden_type=True, checked_regions=None):):


def build_tracks(items, stop_codons, promoters, show_existing_features):
    tracks = []

    for region_name in items:

        regions = []

        # NOTE: I'm forcing the strands all to be 1 here to visualise on the same line in the linear genome


        # Here is a place to update FASTA ID headers

        for region in items[region_name]:
            region_dict = {'id': region['hit_id'], 'start': int(region['start']), 'end': int(region['end']),
                           'name': region[
                                       'name'] + " [" + region['score'] + "]",
                           'strand': region['strand'], 'actual_strand': region['actual_strand']}
            regions.append(region_dict)

        # if show_existing_features:
        #     print ('show features')
        #     print (region_name)

        # if show_existing_features and region_name.startswith('EXISTING'):
        #     print('&&&')
        #     print(region_name)
            region_name = region_name.split(" ")[0]



        track = {'trackName': region_name_mapper[region_name],
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
                 'items': regions}

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

    print('assoc')

    print(genome.id)

    # assoc_hits = models.AssociatedHits.objects().get(genome_id=str(genome.id))


    aggregate = models.AssociatedHits._get_collection().aggregate([
        {"$match": {
            "genome_id": str(genome.id)
        }},
    ])

    assoc_hits = list(aggregate)

    # for as in aggregate:


    print(assoc_hits)

    for hit in assoc_hits:
        print(hit['_id'])
        print(hit['region1'])
        associated_dict[str(hit['_id'])] = (hit['region1'].split("region_")[1] + " and " + hit[
            'region2'].split("region_")[1])

    return associated_dict


def open_file(filename):
    print('file name is ')
    print(filename)
    file_path = "static/uploads/" + filename

    print('file path is')
    print(file_path)
    while not os.path.exists(file_path):
        time.sleep(1)
    if os.path.isfile(file_path):
        file = open(file_path, 'rb')
        return file


def randstring(length=10):
    valid_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
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


def test_auto_classify(queries, skip_tags):
    original_classifications = {"AP018269.1": ["Incomplete"], "AP018270.1": ["Incomplete"],
                                "AP018271.1": ["Simple", "Type3"], "AP018272.1": ["Incomplete"],
                                "AP018273.1": ["Incomplete"], "CP034538.1": ["Incomplete"],
                                "DAATWL01": ["Single", "Type2b"], "OOHJ01": ["Incomplete"],
                                "NZ_CP022961.1": ["Multiple", "Type1", "Type3"], "JPPA01": ["Simple", "Type3"],
                                "NUCS01": ["Incomplete"], "LIVZ01": ["Simple", "Type3"],
                                "FTPJ01": ["Multiple", "Type1", "Type2b"], "LAVL01": ["Incomplete"],
                                "NC_013892.1": ["Multiple", "Type1", "Type2a"], "PCQL01": ["Simple", "Type2b"],
                                "JAAMRA01": ["Multiple", "Type2b"], "MDEO01": ["Simple", "Type3"],
                                "WIVP01": ["Simple", "Type2b"], "VHLG01": ["Incomplete"],
                                "MIFU01": ["Simple", "Type2b"], "NZ_AP018449.1": ["Simple", "Type3"],
                                "CABPSQ01": ["Simple", "Type2b"], "CP034537.1": ["Incomplete"],
                                "PHFJ01": ["Single", "Type1"], "JAABMF01": ["Simple", "Type1"],
                                "RBQC01": ["Single", "Type2b"], "AGJN02": ["Simple", "Type3"],
                                "ONZJ01": ["Simple", "Type3"], "WWHE01": ["Simple", "Type3"],
                                "FMXW01": ["Multiple", "Type2b"], "NGVR01": ["Simple", "Type1"],
                                "NZ_CP020038.1": ["Incomplete"], "NZ_CP021745.1": ["Incomplete"],
                                "NZ_CP021746.1": ["Incomplete"], "NZ_CP021747.1": ["Incomplete"],
                                "QGAC01": ["Incomplete"], "NZ_CP024081.1": ["Simple", "Type2b"],
                                "NZ_CP015613.1": ["Simple", "Type1"], "VIWL01": ["Single", "Type2b"],
                                "NJAK01": ["Simple", "Type2a"], "UUIW01": ["Single", "Type2b"],
                                "NZ_CP012672.1": ["Multiple", "Type1", "Type3"], "NZ_CP047267.1": ["Single", "Type2b"],
                                "PCQC01": ["Single", "Type2b"], "NHML01": ["Incomplete"], "FNJL01": ["Simple", "Type3"],
                                "NZ_CP012533.1": ["Single", "Type3"], "NZ_CP012534.1": ["Incomplete"],
                                "NZ_CP012535.1": ["Incomplete"], "NZ_CP012536.1": ["Incomplete"],
                                "NZ_CP012537.1": ["Incomplete"], "NZ_CP012538.1": ["Incomplete"],
                                "NZ_CP012539.1": ["Incomplete"], "NZ_CP012540.1": ["Incomplete"],
                                "RBSM01": ["Single", "Type2b"], "NZ_CP010029.1": ["Single", "Type2a"],
                                "PHHE01": ["Simple", "Type2b"], "NEJT01": ["Single", "Type2b"],
                                "NZ_CP031450.1": ["Simple", "Type2b"], "NZ_CP017708.1": ["Incomplete"],
                                "NZ_CP017709.1": ["Incomplete"], "NZ_CP017710.1": ["Incomplete"],
                                "FYEE01": ["Multiple", "Type2b"], "ATXB01": ["Simple", "Type3"],
                                "APLI01": ["Incomplete"], "AWQP01": ["Single", "Type2b"], "LMTZ01": ["Incomplete"],
                                "NCXP01": ["Simple", "Type3"], "NZ_CP045799.1": ["Single", "Type2b"],
                                "NZ_CP045800.1": ["Incomplete"], "NZ_CP045801.1": ["Incomplete"],
                                "AJXJ01": ["Single", "Type2b"], "NZ_CP047073.1": ["Single", "Type2b"],
                                "NVXX01": ["Single", "Type2b"], "MUIN01": ["Multiple", "Type2b"],
                                "VDNE01": ["Incomplete"], "LXYR01": ["Incomplete"],
                                "NZ_CP022411.1": ["Single", "Type2b"], "WIVR01": ["Multiple", "Type2b"],
                                "WIWH01": ["Single", "Type2b"], "FONE01": ["Incomplete"], "MKZS01": ["Incomplete"],
                                "QJUG01": ["Simple", "Type3"], "LWBP01": ["Simple", "Type3"],
                                "QUOK01": ["Multiple", "Type3"], "FQUQ01": ["Simple", "Type3"],
                                "VIWO01": ["Multiple", "Type3"], "NITZ01": ["Simple", "Type1"],
                                "CBLV01": ["Multiple", "Type2b"], "MASH01": ["Incomplete"], "LQOW01": ["Incomplete"],
                                "NEHI01": ["Simple", "Type2b"], "NZ_CP038254.1": ["Simple", "Type3"],
                                "JTBY01": ["Incomplete"], "FNTY01": ["Single", "Type2b"],
                                "NZ_CP028826.1": ["Single", "Type2b"], "NIRH01": ["Simple", "Type1"],
                                "LVYD01": ["Simple", "Type3"], "NZ_CP025800.1": ["Simple", "Type1"],
                                "NZ_CP025801.1": ["Incomplete"], "NZ_CP025802.1": ["Incomplete"],
                                "QLKY01": ["Incomplete"], "RCFQ01": ["Simple", "Type3"], "SMCH01": ["Simple", "Type2b"],
                                "NZ_CP015381.1": ["Incomplete"], "NMRE01": ["Single", "Type2b"],
                                "QSNX01": ["Single", "Type2b"], "NZ_CM001558.1": ["Simple", "Type2b"],
                                "FNNQ01": ["Incomplete"], "FQUS01": ["Multiple", "Type3"],
                                "NZ_LT629762.1": ["Single", "Type2b"], "NZ_CP031065.1": ["Incomplete"],
                                "NZ_CP031066.1": ["Simple", "Type3"], "NZ_CP031067.1": ["Incomplete"],
                                "VSJH01": ["Simple", "Type2b"], "FNYO01": ["Incomplete"],
                                "NZ_CP036313.1": ["Simple", "Type3"], "NZ_CP036314.1": ["Incomplete"],
                                "NZ_CP036315.1": ["Incomplete"], "NZ_CP041668.1": ["Simple", "Type3"],
                                "NZ_CP041669.1": ["Incomplete"], "FXYF01": ["Simple", "Type3"],
                                "NIRS01": ["Single", "Type2b"], "FNCO01": ["Single", "Type2b"],
                                "RHQN01": ["Multiple", "Type1", "Type2b"], "NZ_CP021983.2": ["Incomplete"],
                                "RCWL01": ["Simple", "Type3"], "QUMQ01": ["Incomplete"], "QAJM01": ["Single", "Type2b"],
                                "NBRZ01": ["Single", "Type2b"], "BJLR01": ["Incomplete"],
                                "NZ_CP039291.1": ["Incomplete"], "NC_013947.1": ["Multiple", "Type3"],
                                "JMCC02": ["Simple", "Type3"], "NBVR01": ["Single", "Type1"],
                                "QAIL01": ["Simple", "Type3"], "QWFB01": ["Single", "Type2b"],
                                "NZ_CP031648.1": ["Simple", "Type2b"], "MCHY01": ["Simple", "Type3"],
                                "NZ_CP041186.1": ["Multiple", "Type3"], "NC_013216.1": ["Simple", "Type3"],
                                "JYHW01": ["Single", "Type2b"], "WIVY01": ["Multiple", "Type1", "Type2b"],
                                "FYDX01": ["Single", "Type2b"], "MUGY01": ["Incomplete"], "FNYJ01": ["Incomplete"],
                                "JJML01": ["Simple", "Type3"], "FNTZ01": ["Multiple", "Type2b"],
                                "NZ_CP029064.1": ["Simple", "Type3"], "LRUN01": ["Single", "Type2b"],
                                "VIUF01": ["Simple", "Type2b"], "VZZK01": ["Simple", "Type3"],
                                "AJLJ01": ["Simple", "Type3"], "CAADIW01": ["Single", "Type1"],
                                "AXVJ01": ["Incomplete"], "VIUC01": ["Multiple", "Type2b"],
                                "AMBZ01": ["Multiple", "Type2b"], "QGGJ01": ["Incomplete"],
                                "VUOC01": ["Simple", "Type3"], "QAAE01": ["Simple", "Type3"], "NCWQ01": ["Incomplete"],
                                "PVTU01": ["Incomplete"], "BBMZ01": ["Single", "Type2b"],
                                "NZ_CP054043.1": ["Simple", "Type1"], "SJSL01": ["Simple", "Type3"],
                                "FQYP01": ["Simple", "Incomplete"], "NZ_CP011129.1": ["Simple", "Type3"],
                                "NC_012961.1": ["Incomplete"], "NC_012962.1": ["Multiple", "Type1", "Type2b"],
                                "BBXD01": ["Incomplete"], "NZ_CP029196.1": ["Simple", "Type3"],
                                "AKJT01": ["Multiple", "Type2b"], "NVPT01": ["Simple", "Type3"],
                                "BBXG01": ["Incomplete"], "ALVN01": ["Simple", "Type3"], "NJFA02": ["Single", "Type1"],
                                "NC_019738.1": ["Incomplete"], "NC_019739.1": ["Incomplete"],
                                "NC_019740.1": ["Incomplete"], "NC_019741.1": ["Incomplete"],
                                "NC_019742.1": ["Incomplete"], "NC_019743.1": ["Incomplete"],
                                "NC_019760.1": ["Incomplete"], "NC_019761.1": ["Incomplete"],
                                "NC_019762.1": ["Incomplete"], "QPCD01": ["Incomplete"], "QTPO01": ["Simple", "Type3"],
                                "FOEO01": ["Single", "Type2b"], "QWLL01": ["Incomplete"],
                                "QOVA01": ["Single", "Type2b"], "NZ_CP014262.1": ["Multiple", "Type2b"],
                                "FNDJ01": ["Incomplete"], "NZ_AP017422.1": ["Incomplete"],
                                "SNXZ01": ["Simple", "Type3"], "FXWP01": ["Single", "Type1"],
                                "UTBZ01": ["Multiple", "Type2b"], "BCBA01": ["Single", "Type2b"],
                                "VSRQ01": ["Incomplete"], "LFWB01": ["Incomplete"], "QTUB01": ["Simple", "Type2a"],
                                "NZ_CP053584.1": ["Single", "Type2b"], "NZ_CP010897.2": ["Simple", "Type2b"],
                                "NZ_CP010898.2": ["Incomplete"], "WIVZ01": ["Multiple", "Type1", "Type2b"],
                                "NZ_CP013341.1": ["Simple", "Type3"], "JACAQG01": ["Simple", "Type2b"],
                                "FNKR01": ["Simple", "Type3"], "NZ_CP027723.1": ["Single", "Type2b"],
                                "MDEN01": ["Incomplete"], "CVRZ01": ["Simple", "Type1"],
                                "NZ_CP038033.1": ["Incomplete"], "NZ_CP044217.1": ["Incomplete"],
                                "NZ_CP044218.1": ["Simple", "Type3"], "PENV01": ["Simple", "Type1"],
                                "NRQY01": ["Simple", "Type1"], "SISB01": ["Multiple", "Type2b"],
                                "NZ_LT629732.1": ["Simple", "Type3"], "AOCZ01": ["Simple", "Type1"],
                                "NZ_CP039371.1": ["Simple", "Type1"], "NZ_CP039372.1": ["Incomplete"],
                                "JAFA01": ["Incomplete"], "FNOY01": ["Incomplete"], "CABPSP01": ["Simple", "Type2b"],
                                "LGSI01": ["Single", "Type2b"], "VZRB01": ["Simple", "Type3"],
                                "MKWS01": ["Multiple", "Type2b"], "VIUI01": ["Multiple", "Type2b"],
                                "RXOM01": ["Incomplete"], "BCQP01": ["Incomplete"], "SMTE01": ["Simple", "Type1"],
                                "QMEY01": ["Incomplete"], "MBDT01": ["Simple", "Type2b"], "LKPJ01": ["Simple", "Type3"],
                                "OGTP01": ["Simple", "Type3"], "QKTW01": ["Simple", "Type3"],
                                "NC_005773.3": ["Multiple", "Type1", "Type2b"], "NC_007274.1": ["Incomplete"],
                                "NC_007275.1": ["Incomplete"], "NZ_CP048835.1": ["Simple", "Type3"],
                                "NC_010162.1": ["Multiple", "Type3"], "NEVM01": ["Single", "Type1"],
                                "FOUX01": ["Simple", "Type3"], "NZ_CP023526.1": ["Incomplete"],
                                "NZ_CP054422.1": ["Single", "Type2b"], "VOIX01": ["Single", "Type2b"],
                                "VIWA01": ["Simple", "Type3"], "VEBC01": ["Incomplete"],
                                "WIWK01": ["Multiple", "Type1", "Type2b"], "QREK01": ["Simple", "Type3"],
                                "NZ_CM002330.1": ["Single", "Type2b"], "NZ_CM002331.1": ["Incomplete"],
                                "BAHC01": ["Incomplete"], "NZ_CP042968.1": ["Incomplete"],
                                "NZ_CP018049.1": ["Single", "Type2b"], "VZPM01": ["Simple", "Type2b"],
                                "QLIN01": ["Single", "Type2b"], "AUYR01": ["Incomplete"], "NTYK01": ["Incomplete"],
                                "VSFF01": ["Simple", "Type3"], "LRTK01": ["Incomplete"], "ARBP01": ["Simple", "Type3"],
                                "ABCS01": ["Multiple", "Type3"], "BJNF01": ["Simple", "Type3"],
                                "VOQD01": ["Simple", "Type3"], "VIUL01": ["Multiple", "Type2b"],
                                "WHJD01": ["Multiple", "Type1", "Type2b"], "MLFS01": ["Simple", "Type1"],
                                "NZ_CP024900.1": ["Multiple", "Type1", "Type2a", "Type2b"],
                                "NZ_CP009555.1": ["Incomplete"], "NZ_CP009556.1": ["Incomplete"],
                                "NZ_CP013426.1": ["Simple", "Type3"], "NZ_CP013427.1": ["Incomplete"],
                                "NZ_CP013428.1": ["Incomplete"], "NZ_CP013429.1": ["Incomplete"],
                                "POUA01": ["Incomplete"], "AJUL01": ["Simple", "Type3"], "PCOS01": ["Simple", "Type2b"],
                                "QKZA01": ["Simple", "Type1"], "FNQW01": ["Incomplete"], "JADL01": ["Incomplete"],
                                "CABHXE01": ["Simple", "Type1"], "VIKS01": ["Incomplete"],
                                "MOBX01": ["Single", "Type2b"], "QKLR01": ["Simple", "Type3"],
                                "JCLE01": ["Simple", "Type1"], "FSRS01": ["Simple", "Type3"],
                                "NZ_LR134159.1": ["Multiple", "Type2b"], "VCKW01": ["Multiple", "Type3"],
                                "WTCR01": ["Simple", "Type3"], "LLWH01": ["Single", "Type2b"],
                                "NZ_CP027738.1": ["Simple", "Type3"], "QKVL01": ["Single", "Type2b"],
                                "NZ_CP033932.1": ["Simple", "Type3"], "NZ_CM001441.1": ["Incomplete"],
                                "QGTQ01": ["Simple", "Type3"], "RCZD01": ["Simple", "Type1"],
                                "PYLU01": ["Simple", "Type1"], "NZ_CP011288.1": ["Single", "Type2b"],
                                "FPLG01": ["Simple", "Type1"], "NZ_CP012371.1": ["Incomplete"],
                                "NZ_CP022478.1": ["Incomplete"], "NMQR01": ["Multiple", "Type1", "Type2a", "Type2b"],
                                "CTIO01": ["Simple", "Type1"], "VCNG01": ["Multiple", "Type2b"],
                                "NZ_CP007410.1": ["Multiple", "Type2b"], "NKHL01": ["Incomplete"],
                                "MVGR01": ["Simple", "Type3"], "NZ_CP056779.1": ["Incomplete"],
                                "NZ_CP056780.1": ["Multiple", "Type1"], "NZ_CP056781.1": ["Incomplete"],
                                "NZ_CP056782.1": ["Single", "Type2b"], "NQKQ01": ["Single", "Type2b"],
                                "JOGE01": ["Simple", "Type3"], "NZ_CP009533.1": ["Multiple", "Type2b"],
                                "NQKJ01": ["Multiple", "Type1", "Type2b"], "NETK01": ["Simple", "Type3"],
                                "NZ_CP031062.1": ["Incomplete"], "NZ_CP031063.1": ["Simple", "Type3"],
                                "NZ_CP031064.1": ["Incomplete"], "NZ_CP004078.1": ["Simple", "Type3"],
                                "PJZH01": ["Incomplete"], "FNPW01": ["Multiple", "Type1"],
                                "SEUB01": ["Multiple", "Type2b"], "UPHP01": ["Simple", "Type3"],
                                "JNGI01": ["Simple", "Type1"], "UUFD01": ["Incomplete"], "AAWS01": ["Incomplete"],
                                "NZ_CP021659.1": ["Multiple", "Type1"], "NZ_CP021660.1": ["Incomplete"],
                                "NZ_CP021661.1": ["Incomplete"], "NZ_CP021662.1": ["Incomplete"],
                                "MOBP01": ["Single", "Type2b"], "OIFR01": ["Simple", "Type3"],
                                "JSAL01": ["Multiple", "Type2b"], "NZ_CP011104.1": ["Multiple", "Type1", "Type2b"],
                                "MOBI01": ["Simple", "Type2b"], "PUJU01": ["Multiple", "Type1", "Type2a", "Type2b"],
                                "BIFQ01": ["Simple", "Type3"], "NZ_CP025035.2": ["Incomplete"],
                                "LIUV01": ["Single", "Type2b"], "NC_010830.1": ["Simple", "Type3"],
                                "CABPSR01": ["Single", "Type2b"], "CVTM01": ["Single", "Type2b"],
                                "RQJP01": ["Simple", "Type3"], "NZ_CP009288.1": ["Simple", "Type3"],
                                "NZ_CM001025.1": ["Simple", "Type2b"], "MOBT01": ["Single", "Type2b"],
                                "NZ_LR134318.1": ["Multiple", "Type2b"], "ABBN01": ["Incomplete"],
                                "NZ_CP039287.1": ["Incomplete"], "NZ_CP039288.1": ["Simple", "Type3"],
                                "NZ_CP039289.1": ["Incomplete"], "LAIJ01": ["Simple", "Type3"],
                                "LFCV01": ["Simple", "Type3"], "WWJN01": ["Simple", "Type3"], "VZPQ01": ["Incomplete"],
                                "VOBN01": ["Single", "Type2b"], "QGTG01": ["Simple", "Type3"],
                                "AYLO01": ["Simple", "Type3"], "NZ_LT707064.1": ["Single", "Type2b"],
                                "NZ_CP020076.1": ["Incomplete"], "NZ_CP020077.1": ["Incomplete"],
                                "NZ_CP020078.1": ["Incomplete"], "NZ_CP020079.1": ["Incomplete"],
                                "NZ_CP020080.1": ["Simple", "Type1"], "NZ_CP020081.1": ["Incomplete"],
                                "AMRI01": ["Incomplete"], "NZ_LT629705.1": ["Incomplete"],
                                "NRST01": ["Single", "Type2b"], "NZ_CP050291.1": ["Simple", "Type1"],
                                "NZ_CP025263.1": ["Simple", "Type2b"], "FXWM01": ["Simple", "Type3"],
                                "NZ_CP034725.1": ["Multiple", "Type2b"], "MKQR01": ["Incomplete"],
                                "FOCT01": ["Incomplete"], "NUVY01": ["Incomplete"], "MRCJ01": ["Simple", "Type3"],
                                "JUQG01": ["Simple", "Type1"], "LECZ01": ["Single", "Type1"],
                                "MTHI01": ["Single", "Type2b"], "NZ_CP022121.1": ["Simple", "Type3"],
                                "NZ_CM001561.1": ["Single", "Type2b"], "NZ_CP017141.1": ["Multiple", "Type3"],
                                "AZAN01": ["Simple", "Type3"], "AGFX01": ["Incomplete"], "VDCQ01": ["Simple", "Type3"],
                                "QHJL01": ["Simple", "Type3"], "QWEX01": ["Multiple", "Type3"],
                                "LMCT01": ["Incomplete"], "NTTM01": ["Incomplete"],
                                "VZZS01": ["Multiple", "Type1", "Type2b"], "SMFW01": ["Multiple", "Type2b"],
                                "UEXE01": ["Single", "Type1"], "NZ_CP013046.2": ["Incomplete"],
                                "NZ_CP013047.2": ["Simple", "Type1"], "FRBZ01": ["Single", "Type2b"],
                                "AKJH01": ["Single", "Type2b"], "BBLT01": ["Simple", "Type3"], "NBWC01": ["Incomplete"],
                                "NZ_CP007039.1": ["Multiple", "Type2b"], "FMCR01": ["Simple", "Type3"],
                                "VIUK01": ["Single", "Type2b"], "MVHE01": ["Incomplete"],
                                "RCOE01": ["Single", "Type2b"], "QGSY01": ["Simple", "Type3"], "AXWS01": ["Incomplete"],
                                "AYMJ01": ["Single", "Type2b"], "VOBI01": ["Single", "Type2b"],
                                "AKJK01": ["Single", "Type2b"], "FNUD01": ["Single", "Type2b"],
                                "MOBJ01": ["Multiple", "Type2b"], "CAAKGZ01": ["Single", "Type2b"],
                                "FOUB01": ["Simple", "Type3"], "MUXN01": ["Simple", "Type3"],
                                "LKBR01": ["Simple", "Type2b"], "UTVB01": ["Multiple", "Type1", "Type2b"],
                                "PENZ01": ["Simple", "Type1"], "NZ_CP009451.1": ["Single", "Type2b"],
                                "NZ_CP034148.1": ["Incomplete"], "NZ_CP034149.1": ["Incomplete"],
                                "NZ_CP034150.1": ["Incomplete"], "NZ_CP034151.1": ["Simple", "Type1"],
                                "NZ_AP017313.1": ["Simple", "Type3"], "FAOS01": ["Incomplete"],
                                "NZ_CP027727.1": ["Single", "Type2b"], "NZ_CP035319.1": ["Incomplete"],
                                "QAIP01": ["Simple", "Type2b"], "FNBR01": ["Multiple", "Type3"],
                                "AXDH01": ["Single", "Type2b"], "FMVH01": ["Simple", "Type3"],
                                "NZ_CP036488.1": ["Incomplete"], "NZ_CP036489.1": ["Incomplete"],
                                "NZ_CP036490.1": ["Single", "Type2b"], "CAQM01": ["Incomplete"],
                                "LOWA01": ["Simple", "Type2b"], "NZ_CP049044.1": ["Single", "Type2b"],
                                "NZ_CP010896.1": ["Single", "Type2b"], "NC_017168.1": ["Single", "Type1"],
                                "NC_017169.1": ["Incomplete"], "NC_017170.1": ["Incomplete"],
                                "NZ_CP031641.1": ["Single", "Type2b"], "VKDC01": ["Multiple", "Type3"],
                                "JOAG01": ["Incomplete"], "MWQG01": ["Incomplete"], "VDFY01": ["Simple", "Type3"],
                                "ALVK01": ["Multiple", "Type3"], "QFRW01": ["Simple", "Type3"],
                                "BILZ01": ["Simple", "Type3"], "BAXG01": ["Multiple", "Type1"],
                                "MWPQ01": ["Simple", "Type3"], "WIWM01": ["Single", "Type2b"],
                                "FOCU01": ["Single", "Type2b"], "MQZX01": ["Simple", "Type1"],
                                "RKHS01": ["Single", "Type1"], "QHHZ01": ["Incomplete"],
                                "MYFJ01": ["Multiple", "Type1", "Type3", "Type2a", "Type2b"],
                                "NC_016901.1": ["Single", "Type1"], "NC_016905.1": ["Incomplete"],
                                "PEIB01": ["Incomplete"], "MOBQ01": ["Single", "Type2b"],
                                "NXNJ01": ["Single", "Type2b"], "NZ_CP044407.1": ["Incomplete"],
                                "PYBV01": ["Incomplete"], "JABTYG01": ["Multiple", "Type2b"],
                                "NZ_CP042468.1": ["Multiple", "Type1", "Type3"], "NZ_CP014135.1": ["Simple", "Type2b"],
                                "NC_016818.1": ["Multiple", "Type2b"], "NC_016819.1": ["Incomplete"],
                                "NC_016835.1": ["Incomplete"], "NC_017092.1": ["Incomplete"],
                                "MTAX01": ["Simple", "Type3"], "NC_015559.1": ["Simple", "Type3"],
                                "LQRT01": ["Incomplete"], "NZ_LS999839.1": ["Simple", "Type3"],
                                "SOCV01": ["Single", "Type2b"], "ASRX01": ["Single", "Type3"],
                                "NZ_CP044064.1": ["Simple", "Type2b"], "AKJM01": ["Single", "Type2b"],
                                "SMKX01": ["Simple", "Type3"], "CAAJVF01": ["Simple", "Type3"],
                                "VIUJ01": ["Simple", "Type2b"], "LGTC01": ["Incomplete"],
                                "NZ_CP033893.1": ["Incomplete"], "NZ_CP033894.1": ["Single", "Type1"],
                                "NZ_CP033895.1": ["Incomplete"], "JXRA01": ["Simple", "Type3"],
                                "RQPI01": ["Simple", "Type3"], "NZ_CP023695.1": ["Simple", "Type3"],
                                "NZ_LR134335.1": ["Single", "Type1"], "SMJU01": ["Simple", "Type3"],
                                "LMCV01": ["Single", "Type2b"], "PKNM01": ["Simple", "Type1"],
                                "PIQI01": ["Simple", "Type1"], "FZPH01": ["Simple", "Type3"],
                                "WIWB01": ["Single", "Type2b"], "NC_009253.1": ["Simple", "Type3"],
                                "SOZA01": ["Single", "Type2b"], "NZ_LT855380.1": ["Incomplete"],
                                "NZ_CP014947.1": ["Multiple", "Type2b"], "ALVJ01": ["Incomplete"],
                                "NZ_CP013459.1": ["Incomplete"], "NZ_CP013460.1": ["Incomplete"],
                                "NZ_CP013461.1": ["Simple", "Type3"], "NZ_CP048408.1": ["Multiple", "Type2b"],
                                "NZ_CP003181.2": ["Simple", "Type2b"], "VFIO01": ["Single", "Type2b"],
                                "MASS01": ["Incomplete"], "NC_020453.1": ["Simple", "Type3"],
                                "PYUC01": ["Simple", "Type1"], "VEGT01": ["Simple", "Type3"],
                                "MKZO01": ["Single", "Type2b"], "WIWE01": ["Single", "Type2b"],
                                "FMWY01": ["Simple", "Type3"], "MWQL01": ["Incomplete"], "FMVD01": ["Incomplete"],
                                "NZ_CP023969.1": ["Single", "Type2b"], "NZ_CP029608.1": ["Single", "Type2b"],
                                "SMKU01": ["Incomplete"], "FUKJ01": ["Single", "Type3"],
                                "JONO01": ["Multiple", "Type1", "Type2a"], "RAVW01": ["Incomplete"],
                                "PDUD01": ["Incomplete"], "MKMC01": ["Simple", "Type3"],
                                "NC_017448.1": ["Single", "Type3"], "PVZV01": ["Incomplete"],
                                "NZ_CP031069.1": ["Incomplete"], "NZ_CP031070.1": ["Simple", "Type2b"],
                                "NZ_CP023269.1": ["Multiple", "Type2b"], "VLLP01": ["Incomplete"],
                                "NZ_CM001559.1": ["Simple", "Type3"], "NZ_CP029983.1": ["Simple", "Type2b"],
                                "VHKL01": ["Simple", "Type3"], "NZ_CP027218.1": ["Single", "Type2b"],
                                "JPPZ01": ["Simple", "Type3"], "AKJD01": ["Incomplete"],
                                "VCNJ01": ["Multiple", "Type2b"], "NZ_CP013423.1": ["Incomplete"],
                                "NZ_CP013424.1": ["Incomplete"], "NZ_CP013425.1": ["Incomplete"],
                                "SMKK01": ["Simple", "Type3"], "SODH01": ["Incomplete"],
                                "AZSS01": ["Multiple", "Type3"], "JFHN01": ["Simple", "Type1"],
                                "MUNM01": ["Single", "Type2b"], "NC_021492.1": ["Incomplete"],
                                "NC_021500.1": ["Single", "Type1"], "RCSU01": ["Simple", "Type3"],
                                "SMOD01": ["Simple", "Type3"], "NZ_CP042382.1": ["Simple", "Type3"],
                                "NC_008268.1": ["Incomplete"], "NC_008269.1": ["Incomplete"],
                                "NC_008270.1": ["Incomplete"], "NC_008271.1": ["Incomplete"],
                                "JYLE01": ["Multiple", "Type1", "Type2b"], "PYMM01": ["Single", "Type1"],
                                "NZ_CP007699.2": ["Simple", "Type3"], "QAOU01": ["Incomplete"],
                                "WBOI01": ["Multiple", "Type1", "Type2b"], "CAACVJ01": ["Incomplete"],
                                "BJMN01": ["Simple", "Type3"], "SMDG01": ["Simple", "Type1"],
                                "CABIWI01": ["Multiple", "Type1", "Type2b"],
                                "WHZZ01": ["Multiple", "Type1", "Type2a", "Type2b"], "QAOQ01": ["Simple", "Type3"],
                                "RCBZ01": ["Incomplete"], "NZ_CP022303.1": ["Incomplete"],
                                "NZ_CP022304.1": ["Incomplete"], "NZ_CP022305.1": ["Incomplete"],
                                "NZ_CP022306.1": ["Incomplete"], "AQRJ01": ["Simple", "Type3"],
                                "FNGP01": ["Incomplete"], "RJKM01": ["Incomplete"], "PKND01": ["Simple", "Type1"],
                                "FOSU01": ["Incomplete"], "AWZT01": ["Simple", "Type3"],
                                "NZ_CP009458.1": ["Single", "Type1"], "WEGH01": ["Incomplete"],
                                "VZZZ01": ["Simple", "Type3"], "NZ_CP043060.1": ["Multiple", "Type2b"],
                                "VZZR01": ["Incomplete"], "JAAQYP01": ["Simple", "Type3"],
                                "MTBD01": ["Simple", "Type1"], "NZ_CP019686.1": ["Single", "Type1"],
                                "VUAZ01": ["Simple", "Type2b"], "AZXK01": ["Incomplete"],
                                "BBIR01": ["Multiple", "Type2b"], "MWLO01": ["Incomplete"],
                                "QYZD01": ["Simple", "Type2b"], "MOBZ01": ["Single", "Type2b"],
                                "FOBB01": ["Simple", "Type3"], "FMDM01": ["Simple", "Type3"],
                                "NZ_CP019888.1": ["Simple", "Type1"], "AUAX01": ["Simple", "Type3"],
                                "AVEF02": ["Simple", "Type1"], "FNAD01": ["Simple", "Type3"],
                                "BBCC01": ["Simple", "Type1"], "QAOV01": ["Single", "Type2b"], "BAZX01": ["Incomplete"],
                                "NKQZ01": ["Simple", "Type3"], "NZ_CP022960.1": ["Simple", "Type2b"],
                                "SHKK01": ["Simple", "Type3"], "NZ_CP012673.1": ["Multiple", "Type1", "Type3"],
                                "WBKQ01": ["Simple", "Type3"], "LQAL01": ["SImple", "Type2b"],
                                "FCNY02": ["Simple", "Type3"], "VFEU01": ["Simple", "Type2b"], "SHKT01": ["Incomplete"],
                                "BAVR01": ["Incomplete"], "OAOQ01": ["Simple", "Type3"], "MSSW01": ["Simple", "Type3"],
                                "JOGR01": ["Simple", "Type3"], "NZ_CP034337.1": ["Incomplete"],
                                "CAACUY01": ["Incomplete"], "QVNU01": ["Simple", "Type3"],
                                "ASSC01": ["Simple", "Type3"], "NZ_CP028035.1": ["Single", "Type2b"],
                                "NZ_CP028036.1": ["Incomplete"], "NZ_CP028037.1": ["Incomplete"],
                                "NZ_CP028038.1": ["Incomplete"], "NZ_CP028039.1": ["Incomplete"],
                                "NZ_CP028040.1": ["Incomplete"], "NZ_CP028041.1": ["Incomplete"],
                                "NZ_CP028042.1": ["Incomplete"], "NZ_CP057330.1": ["Simple", "Type2b"],
                                "NZ_CP057331.1": ["Incomplete"], "NZ_CP057332.1": ["Incomplete"],
                                "NZ_CP057333.1": ["Incomplete"], "NZ_CP046052.1": ["Incomplete"],
                                "NZ_CP046053.1": ["Incomplete"], "NZ_CP046054.1": ["Incomplete"],
                                "BFCB01": ["Simple", "Type3"], "NZ_CM001489.1": ["Incomplete"],
                                "RBOV01": ["Single", "Type2b"], "NZ_CP030750.1": ["Multiple", "Type1", "Type2b"],
                                "SODV01": ["Simple", "Type3"], "QEKL01": ["Multiple", "Type2b"],
                                "QJRT01": ["Single", "Type2b"], "NZ_LT629778.1": ["Single", "Type2b"],
                                "NZ_CP024634.1": ["Incomplete"], "LVTS01": ["Simple", "Type1"],
                                "NZ_AP018150.1": ["Single", "Type2b"], "CPYD01": ["Single", "Type2a"],
                                "RIAR02": ["Multiple", "Type3"], "NZ_CP010407.1": ["Simple", "Type3"],
                                "NZ_CP010408.1": ["Incomplete"], "AUKP01": ["Simple", "Type3"],
                                "MTSA01": ["Single", "Type2b"], "MOBL01": ["Single", "Type2b"],
                                "NZ_CP028158.1": ["Incomplete"], "MOBY01": ["Multiple", "Type2b"],
                                "NIBV01": ["Multiple", "Type1", "Type2a"], "SJOP01": ["Incomplete"],
                                "VEJO01": ["Single", "Type2b"], "FWXB01": ["Incomplete"], "SMKT01": ["Single", "Type3"],
                                "FPBO01": ["Incomplete"], "NZ_CP012831.1": ["Multiple", "Type2b"],
                                "JAAQYX01": ["Multiple", "Type2b"], "QGSZ01": ["Simple", "Type3"],
                                "QARA01": ["Single", "Type2b"], "RXLQ01": ["Simple", "Type3"],
                                "ANOR01": ["Single", "Type2b"], "PJBP01": ["Simple", "Type3"],
                                "NQKI01": ["Single", "Type2b"], "CBSW01": ["Multiple", "Type1", "Type2a", "Type2b"],
                                "VCKX01": ["Incomplete"], "FNGF01": ["Simple", "Type3"],
                                "NZ_CP010310.2": ["Incomplete"], "WJIE01": ["Multiple", "Type1"],
                                "VAUR01": ["Multiple", "Type2b"], "NC_017956.1": ["Incomplete"],
                                "NC_017957.2": ["Incomplete"], "NC_017958.1": ["Incomplete"],
                                "NC_017959.1": ["Incomplete"], "NC_017966.1": ["Incomplete"],
                                "NZ_CP013368.1": ["Incomplete"], "NZ_CP013369.1": ["Single", "Type2b"],
                                "NZ_CP013370.1": ["Incomplete"], "NZ_CP013371.1": ["Incomplete"],
                                "VFVY01": ["Incomplete"], "NZ_CP023567.1": ["Single", "Type2b"],
                                "NZ_CP023568.1": ["Incomplete"], "RJKL01": ["Incomplete"],
                                "NZ_CP045767.1": ["Simple", "Type2b"], "JOGP01": ["Simple", "Type3"],
                                "CABHPT01": ["Simple", "Type1"], "UWFA01": ["Incomplete"],
                                "QVNQ01": ["Simple", "Type3"], "PJBC01": ["Simple", "Type2b"],
                                "SHKZ01": ["Simple", "Type3"], "JQFM01": ["Incomplete"], "SLVP01": ["Simple", "Type3"],
                                "NZ_CP005969.1": ["Multiple", "Type1", "Type2b"], "NZ_CP021106.3": ["Incomplete"],
                                "JYLH01": ["Single", "Type2b"], "LJCS01": ["Multiple", "Type1"],
                                "NZ_CP017606.1": ["Incomplete"], "NZ_CP017607.1": ["Single", "Type1"],
                                "NZ_CP017608.1": ["Incomplete"], "NZ_CP017609.1": ["Incomplete"],
                                "SOAB01": ["Multiple", "Type3"], "QKLY01": ["Simple", "Type2b"],
                                "JAAQWH01": ["Single", "Type2b"], "AKJS01": ["Single", "Type2b"],
                                "QPDS01": ["Multiple", "Type2b"], "NZ_CP023965.1": ["Single", "Type1"],
                                "NZ_CP029231.1": ["Incomplete"], "NZ_CP029232.1": ["Incomplete"],
                                "NZ_CP029233.1": ["Incomplete"], "NZ_CP029234.1": ["Incomplete"],
                                "NZ_CP029235.1": ["Simple", "Type3"], "LNCD01": ["Simple", "Type3"],
                                "NVDM01": ["Single", "Type1"], "NZ_CP011020.1": ["Multiple", "Type2b"],
                                "PZZQ01": ["Multiple", "Type1"], "NZ_CP011807.3": ["Simple", "Type2b"],
                                "NZ_CP011808.2": ["Incomplete"], "NZ_CP011809.2": ["Incomplete"],
                                "ALJC01": ["Single", "Type2b"], "QJRQ01": ["Incomplete"], "QEOK01": ["Simple", "Type3"],
                                "NZ_CP029618.1": ["Simple", "Type3"], "NZ_CP010016.1": ["Incomplete"],
                                "NZ_CP010017.1": ["Incomplete"], "NZ_CP010018.1": ["Incomplete"],
                                "SOCQ01": ["Single", "Type2b"], "RJKE01": ["Multiple", "Type3"],
                                "QAOO01": ["Incomplete"], "JMCL01": ["Multiple", "Type2b"], "QBKC01": ["Incomplete"],
                                "NZ_CP034335.1": ["Incomplete"], "VDLX02": ["Incomplete"],
                                "SSMR01": ["Multiple", "Type3"], "NSCM01": ["Multiple", "Type1", "Type2a", "Type2b"],
                                "VMSG01": ["Multiple", "Type2b"], "ABCQ01": ["Simple", "Type1"],
                                "OUND01": ["Incomplete"], "QAJI01": ["Simple", "Type3"],
                                "NZ_CP045761.1": ["Single", "Type2b"], "SUNB01": ["Incomplete"],
                                "LJSD01": ["Incomplete"], "NZ_CP041692.1": ["Simple", "Type3"],
                                "PZZR01": ["Incomplete"], "JPMW01": ["Incomplete"], "QPIJ01": ["Simple", "Type3"],
                                "LYUY01": ["Incomplete"], "SMJX01": ["Incomplete"], "VATL01": ["Single", "Type1"],
                                "FMYF01": ["Simple", "Type3"],
                                "PUWV01": ["Multiple", "Type1", "Type2a", "Type2b", "Type3"], "FQYM01": ["Incomplete"],
                                "PJRP01": ["Incomplete"], "QRBE01": ["Incomplete"], "FOVJ01": ["Multiple", "Type3"],
                                "SOCT01": ["Incomplete"], "CABMLW01": ["Multiple", "Type1"],
                                "BDBY01": ["Simple", "Type3"], "PYGV01": ["Incomplete"], "VRLS01": ["Single", "Type1"],
                                "ASTJ01": ["Simple", "Type3"], "LVEJ01": ["Single", "Type2b"],
                                "OUNR01": ["Simple", "Type3"], "FPBP01": ["Single", "Type3"],
                                "FSRU01": ["Simple", "Type2b"], "SMKN01": ["Incomplete"], "ASJB01": ["Simple", "Type3"],
                                "VIYH01": ["Single", "Type2b"], "SNZP01": ["Simple", "Type3"],
                                "NZ_CP014847.1": ["Incomplete"], "NZ_CP014848.1": ["Incomplete"],
                                "NZ_CP014849.1": ["Incomplete"], "NZ_CP014850.1": ["Incomplete"],
                                "NZ_CP014851.1": ["Incomplete"], "NZ_CP014852.1": ["Incomplete"],
                                "NZ_CP014853.1": ["Simple", "Type2b"], "NZ_CP034780.1": ["Incomplete"],
                                "NZ_CP034781.1": ["Incomplete"], "NZ_CP034782.1": ["Incomplete"],
                                "NZ_CP034783.1": ["Simple", "Type3"], "PYBJ01": ["Simple", "Type3"],
                                "PTJB01": ["Incomplete"], "NZ_CP024159.1": ["Incomplete"],
                                "JNYY01": ["Simple", "Type3"], "NZ_CP027756.1": ["Simple", "Type2b"],
                                "SSNZ01": ["Simple", "Incomplete"], "NZ_CP046874.1": ["Multiple", "Type2b"],
                                "WIBD01": ["Single", "Type2b"], "NZ_CP029710.1": ["Simple", "Type3"],
                                "RBRE01": ["Multiple", "Type1", "Type2b"], "NZ_CP024866.1": ["Single", "Type2b"],
                                "JAAAGD01": ["Single", "Type2b"], "JAAEFD01": ["Simple", "Type1"],
                                "RBUY01": ["Incomplete"], "QXQA01": ["Single", "Type3"], "QJRP01": ["Incomplete"],
                                "AXBA01": ["Simple", "Type3"], "OMPE01": ["Incomplete"],
                                "NZ_LT629790.1": ["Multiple", "Type2b"], "LLWI01": ["Multiple", "Type1", "Type2b"],
                                "NZ_LT629746.1": ["Multiple", "Type2b"], "BAOS01": ["Simple", "Type3"],
                                "VLPL01": ["Simple", "Type3"], "LYRP01": ["Simple", "Type1"],
                                "JXDG01": ["Multiple", "Type2b"], "LIPP01": ["Simple", "Type3"],
                                "JAAQWI01": ["Multiple", "Type2b"], "NZ_LT629795.1": ["Single", "Type2b"],
                                "LXEN01": ["Simple", "Type1"], "NZ_CM001514.1": ["Multiple", "Type2b"],
                                "NC_015731.1": ["Simple", "Type3"], "NZ_CP041754.1": ["Multiple", "Type1", "Type2b"],
                                "NZ_CP041755.1": ["Incomplete"], "NZ_CP041756.1": ["Incomplete"],
                                "NZ_CP045572.1": ["Incomplete"], "FOVO01": ["Incomplete"], "JMKF01": ["Incomplete"],
                                "JACCAY01": ["Simple", "Type3"], "MAID01": ["Incomplete"],
                                "NZ_CP009727.1": ["Incomplete"], "NZ_CP009728.1": ["Incomplete"],
                                "AKXV01": ["Simple", "Type3"], "CZQA01": ["Simple", "Type3"],
                                "NZ_CP054609.1": ["Incomplete"], "NZ_CP054610.1": ["Incomplete"],
                                "NZ_CP054611.1": ["Incomplete"], "NZ_CP054612.1": ["Incomplete"],
                                "NZ_CP054613.1": ["Simple", "Type3"], "NZ_CP014226.1": ["Simple", "Type3"],
                                "CABIVM01": ["Single", "Type2b"], "CDPK01": ["Single", "Type1"],
                                "WIAO01": ["Simple", "Type3"], "WSFG01": ["Multiple", "Type1", "Type2a", "Type2b"],
                                "WIVV01": ["Multiple", "Type1", "Type2b"], "NQKN01": ["Single", "Type2b"],
                                "FUWU01": ["Incomplete"], "UTLV01": ["Single", "Type2b"], "QOIO01": ["Simple", "Type3"],
                                "QGLF01": ["Incomplete"], "VSRO01": ["Multiple", "Type2b"],
                                "WIVQ01": ["Single", "Type2b"], "OBDY01": ["Incomplete"], "SJZK02": ["Single", "Type1"],
                                "NZ_CP029482.1": ["Single", "Type2b"], "LHVN01": ["Single", "Type2b"],
                                "UEXF01": ["Single", "Type1"], "VIUR01": ["Multiple", "Type2b"],
                                "QBJA02": ["Multiple", "Type1", "Type2b"], "LZEX01": ["Multiple", "Type1", "Type2b"],
                                "JYLN01": ["Simple", "Type2b"], "NZ_CP022504.1": ["Single", "Type1"],
                                "LFQK01": ["Single", "Type2b"], "NZ_CP038255.1": ["Simple", "Type3"],
                                "NZ_CP024646.1": ["Simple", "Type2b"], "CVJX01": ["Single", "Type1"],
                                "NZ_CP031146.1": ["Simple", "Type1"], "LACH01": ["Simple", "Type2b"],
                                "NZ_CP011253.3": ["Single", "Type2b"], "NZ_CP011518.2": ["Incomplete"],
                                "NZ_CP011519.2": ["Incomplete"], "VJZE01": ["Incomplete"], "QKYK01": ["Incomplete"],
                                "NZ_CP023525.1": ["Single", "Type1"], "AZUB01": ["Simple", "Type2b"],
                                "JRYA01": ["Multiple", "Type1", "Type2b"], "AUEZ01": ["Simple", "Type3"],
                                "VSFG01": ["Simple", "Type3"], "QOIL01": ["Incomplete"], "SJZD01": ["Incomplete"],
                                "VJWF01": ["Simple", "Type3"], "FNSJ01": ["Multiple", "Type3"],
                                "AWNH01": ["Incomplete"], "NZ_CP016211.1": ["Incomplete"], "SMKO01": ["Incomplete"],
                                "VXLC01": ["Simple", "Type3"], "JFCA01": ["Simple", "Type2b"],
                                "PPRY01": ["Simple", "Type2b"], "NZ_LT629708.1": ["Simple", "Type2b"],
                                "NZ_CP026364.1": ["Single", "Type1"], "NZ_CP015225.1": ["Simple", "Type2b"],
                                "MUNY01": ["Simple", "Type3"], "MBLO01": ["Multiple", "Type3"],
                                "NZ_CP012159.1": ["Single", "Type3"], "ALBT01": ["Simple", "Type3"],
                                "RPND01": ["Simple", "Type3"], "LMXH01": ["Simple", "Type3"],
                                "CPZI01": ["Simple", "Type2b"], "SAXA01": ["Simple", "Incomplete"],
                                "QPAO01": ["Single", "Type1"], "SSMQ01": ["Multiple", "Type3"],
                                "PQKR01": ["Single", "Type2b"], "NC_015510.1": ["Simple", "Type3"],
                                "CWJL01": ["Single", "Type1"], "BJMH01": ["Simple", "Type3"],
                                "JABWHT01": ["Simple", "Type1"], "NZ_CP024085.1": ["Simple", "Type3"],
                                "SMFY01": ["Incomplete"], "NZ_CP027732.1": ["Single", "Type2b"],
                                "NQKV01": ["Multiple", "Type2b"], "PYGG01": ["Incomplete"],
                                "AWXZ01": ["Simple", "Type3"],
                                "NZ_LN681227.1": ["Multiple", "Type1", "Type2a", "Type2b"],
                                "PQMB01": ["Simple", "Type1"], "FODH01": ["Simple", "Type3"],
                                "FOIG01": ["Simple", "Type3"], "ARBJ01": ["Simple", "Type3"],
                                "NZ_CP029843.1": ["Simple", "Type3"], "FOLC01": ["Simple", "Type3"],
                                "NZ_LT629691.1": ["Multiple", "Type1", "Type2b"], "FORG01": ["Multiple", "Type1"],
                                "WIWO01": ["Single", "Type2b"], "MUBJ01": ["Multiple", "Type2a", "Type2b"],
                                "QPIP01": ["Simple", "Type1"], "QNVV01": ["Simple", "Type3"],
                                "LHVM01": ["Single", "Type2b"], "PVZG01": ["Simple", "Type3"],
                                "NZ_CP017687.1": ["Single", "Type2b"], "FNQB01": ["Multiple", "Type3"],
                                "NQYG01": ["Single", "Type2b"], "JQNB01": ["Incomplete"], "FNVO01": ["Incomplete"],
                                "VIVV01": ["Simple", "Type3"], "CABLBP01": ["Simple", "Type3"],
                                "NZ_CP043422.1": ["Simple", "Type1"], "QEOF01": ["Simple", "Type2b"],
                                "PVNL01": ["Single", "Type3"], "CQAW01": ["Simple", "Type1"], "AUAF01": ["Incomplete"],
                                "SZQA01": ["Incomplete"], "NC_015379.1": ["Multiple", "Type2b"],
                                "NZ_CP026106.1": ["Simple", "Type3"], "VIUG01": ["Multiple", "Type2b"],
                                "SMSC01": ["Simple", "Type1"], "MQWC01": ["Simple", "Type3"],
                                "NZ_CP042804.1": ["Incomplete"], "PENT01": ["Simple", "Type1"],
                                "NC_014718.1": ["Multiple", "Type2b"], "NC_014722.1": ["Incomplete"],
                                "PENW01": ["Simple", "Type1"], "QRAV01": ["Single", "Type2b"],
                                "QEQQ01": ["Single", "Type2b"], "MKCS01": ["Simple", "Type3"],
                                "JPIX01": ["Single", "Type1"], "NJAH01": ["Incomplete"], "JPLA01": ["Incomplete"],
                                "VFPP01": ["Incomplete"], "BBXC01": ["Incomplete"], "JPOH01": ["Simple", "Type3"],
                                "VWSH01": ["Simple", "Type3"], "FNON01": ["Simple", "Type3"],
                                "FAOZ01": ["Simple", "Type3"], "NC_013131.1": ["Simple", "Type3"],
                                "FOOS01": ["Simple", "Type3"], "NZ_CP053682.1": ["Simple", "Type1"],
                                "NZ_CP013949.1": ["Simple", "Type3"], "VZPJ01": ["Multiple", "Type2b"],
                                "MIIW01": ["Single", "Type1"], "JXDI01": ["Single", "Type2b"],
                                "NZ_LT629704.1": ["Single", "Type2b"], "NZ_CP016634.1": ["Multiple", "Type1", "Type2b"],
                                "SLYK01": ["Incomplete"], "NZ_AP020337.1": ["Multiple", "Type1", "Type2b"],
                                "PYAL01": ["Single", "Type2b"], "FYDV01": ["Simple", "Type2b"],
                                "VCBA01": ["Incomplete"], "QFZQ01": ["Single", "Type3"], "SMSL01": ["Simple", "Type3"],
                                "FODZ01": ["Simple", "Type3"], "PYGD01": ["Multiple", "Type3"],
                                "LHVL01": ["Multiple", "Type2b"], "QAIK01": ["Single", "Type2b"],
                                "FODL01": ["Single", "Type2b"], "NZ_LT828648.1": ["Simple", "Type3"],
                                "FUXU01": ["Simple", "Type1"], "VZII01": ["Single", "Type2b"],
                                "QLLL01": ["Simple", "Type3"], "NZ_CP016176.1": ["Incomplete"],
                                "FQWR01": ["Multiple", "Type3"], "NZ_CP029197.1": ["Simple", "Type3"],
                                "VIUE01": ["Multiple", "Type2b"], "VIVC01": ["Single", "Type2b"],
                                "NEJJ01": ["Single", "Type2b"], "VOQB01": ["Multiple", "Type1"],
                                "NC_007954.1": ["Incomplete"], "AEDD01": ["Simple", "Type3"],
                                "NC_020209.1": ["Simple", "Type2b"], "LKBT01": ["Single", "Type2b"],
                                "NZ_LT629801.1": ["Simple", "Type2b"], "JHVK01": ["Single", "Type2b"],
                                "PVNG01": ["Incomplete"], "RBIS01": ["Incomplete"], "FNTT01": ["Simple", "Type2b"],
                                "VCRA01": ["Incomplete"], "ABXF01": ["Simple", "Type3"], "PJCP01": ["Single", "Type2b"],
                                "NZ_CP045158.1": ["Single", "Type1"], "BBXE01": ["Incomplete"],
                                "RKHU01": ["Simple", "Incomplete"], "LIPN01": ["Simple", "Type3"],
                                "NZ_LT629788.1": ["Single", "Type2b"], "NIBS01": ["Simple", "Type2a"],
                                "JAABNH01": ["Simple", "Type1"], "NZ_CP017599.1": ["Multiple", "Type3"],
                                "QKWJ01": ["Simple", "Type3"], "FNVU01": ["Simple", "Type3"],
                                "RAVX01": ["Simple", "Type3"], "NZ_CP024923.1": ["Simple", "Type3"],
                                "NKFP01": ["Incomplete"], "AEDB02": ["Incomplete"], "NZ_CP038613.1": ["Incomplete"],
                                "NZ_CP018319.1": ["Multiple", "Type2b"], "SLZA01": ["Simple", "Type3"],
                                "WIWL01": ["Single", "Type2b"],
                                "SMFY01_information_Ancylobacter_aquaticus_region_TcB_expanded_[640.9]_4636162_4644109_backward": [
                                    ""],
                                "SMFY01_information_Ancylobacter_aquaticus_region_A2_expanded_[100.2]_4628824_4636126_backward": [
                                    ""],
                                "SMFY01_information_Ancylobacter_aquaticus_region_TcC_expanded_[251.6]_4636162_4644109_backward": [
                                    ""], "QTSW01": ["Single", "Type2b"], "VIVY01": ["Incomplete"],
                                "MOBO01": ["Multiple", "Type2b"], "NZ_CP013997.1": ["Incomplete", "Simple"],
                                "RBIO01": ["Incomplete"], "SNYE01": ["Simple", "Type3"],
                                "NZ_CP038438.1": ["Single", "Type2b"], "SMAS01": ["Single", "Type1"],
                                "VCNA01": ["Simple", "Type3"]}

    # original_classifications = {"AP018269.1": ["Incomplete"], "AP018270.1": ["Incomplete"], "AP018271.1": ["Simple", "Type3"], "AP018272.1": ["Incomplete"], "AP018273.1": ["Incomplete"], "CP034538.1": ["Incomplete"], "DAATWL01": ["Single", "Type2b"], "OOHJ01": ["Incomplete"], "NZ_CP022961.1": ["Multiple", "Type1", "Type3"], "JPPA01": ["Simple", "Type3"], "NUCS01": ["Incomplete"], "LIVZ01": ["Simple", "Type3"], "FTPJ01": ["Multiple", "Type1", "Type2b"], "LAVL01": ["Incomplete"], "NC_013892.1": ["Multiple", "Type1", "Type2a"], "PCQL01": ["Simple", "Type2b"], "JAAMRA01": ["Multiple", "Type2b"], "MDEO01": ["Simple", "Type3"], "WIVP01": ["Simple", "Type2b"], "VHLG01": ["Incomplete"], "MIFU01": ["Simple", "Type2b"], "NZ_AP018449.1": ["Simple", "Type3"], "CABPSQ01": ["Simple", "Type2b"], "CP034537.1": ["Incomplete"], "PHFJ01": ["Single", "Type1"], "JAABMF01": ["Simple", "Type1"], "RBQC01": ["Single", "Type2b"], "AGJN02": ["Simple", "Type3"], "ONZJ01": ["Simple", "Type3"], "WWHE01": ["Simple", "Type3"], "FMXW01": ["Multiple", "Type2b"], "NGVR01": ["Simple", "Type1"], "NZ_CP020038.1": ["Incomplete"], "NZ_CP021745.1": ["Incomplete"], "NZ_CP021746.1": ["Incomplete"], "NZ_CP021747.1": ["Incomplete"], "QGAC01": ["Incomplete"], "NZ_CP024081.1": ["Simple", "Type2b"], "NZ_CP015613.1": ["Simple", "Type1"], "VIWL01": ["Single", "Type2b"], "NJAK01": ["Simple", "Type2a"], "UUIW01": ["Single", "Type2b"], "NZ_CP012672.1": ["Multiple", "Type1", "Type3"], "NZ_CP047267.1": ["Single", "Type2b"], "PCQC01": ["Single", "Type2b"], "NHML01": ["Incomplete"], "FNJL01": ["Simple", "Type3"], "NZ_CP012533.1": ["Single", "Type3"], "NZ_CP012534.1": ["Incomplete"], "NZ_CP012535.1": ["Incomplete"], "NZ_CP012536.1": ["Incomplete"], "NZ_CP012537.1": ["Incomplete"], "NZ_CP012538.1": ["Incomplete"], "NZ_CP012539.1": ["Incomplete"], "NZ_CP012540.1": ["Incomplete"], "RBSM01": ["Single", "Type2b"], "NZ_CP010029.1": ["Single", "Type2a"], "PHHE01": ["Simple", "Type2b"], "NEJT01": ["Single", "Type2b"], "NZ_CP031450.1": ["Simple", "Type2b"], "NZ_CP017708.1": ["Incomplete"], "NZ_CP017709.1": ["Incomplete"], "NZ_CP017710.1": ["Incomplete"], "FYEE01": ["Multiple", "Type2b"], "ATXB01": ["Simple", "Type3"], "APLI01": ["Incomplete"], "AWQP01": ["Single", "Type2b"], "LMTZ01": ["Incomplete"], "NCXP01": ["Simple", "Type3"], "NZ_CP045799.1": ["Single", "Type2b"], "NZ_CP045800.1": ["Incomplete"], "NZ_CP045801.1": ["Incomplete"], "AJXJ01": ["Single", "Type2b"], "NZ_CP047073.1": ["Single", "Type2b"], "NVXX01": ["Single", "Type2b"], "MUIN01": ["Multiple", "Type2b"], "VDNE01": ["Incomplete"], "LXYR01": ["Incomplete"], "NZ_CP022411.1": ["Single", "Type2b"], "WIVR01": ["Multiple", "Type2b"], "WIWH01": ["Single", "Type2b"], "FONE01": ["Incomplete"], "MKZS01": ["Incomplete"], "QJUG01": ["Simple", "Type3"], "LWBP01": ["Simple", "Type3"], "QUOK01": ["Multiple", "Type3"], "FQUQ01": ["Simple", "Type3"], "VIWO01": ["Multiple", "Type3"], "NITZ01": ["Simple", "Type1"], "CBLV01": ["Multiple", "Type2b"], "MASH01": ["Incomplete"], "LQOW01": ["Incomplete"], "NEHI01": ["Simple", "Type2b"], "NZ_CP038254.1": ["Simple", "Type3"], "JTBY01": ["Incomplete"], "FNTY01": ["Single", "Type2b"], "NZ_CP028826.1": ["Single", "Type2b"], "NIRH01": ["Simple", "Type1"], "LVYD01": ["Simple", "Type3"], "NZ_CP025800.1": ["Simple", "Type1"], "NZ_CP025801.1": ["Incomplete"], "NZ_CP025802.1": ["Incomplete"], "QLKY01": ["Incomplete"], "RCFQ01": ["Simple", "Type3"], "SMCH01": ["Simple", "Type2b"], "NZ_CP015381.1": ["Incomplete"], "NMRE01": ["Single", "Type2b"], "QSNX01": ["Single", "Type2b"], "NZ_CM001558.1": ["Simple", "Type2b"], "FNNQ01": ["Incomplete"], "FQUS01": ["Multiple", "Type3"], "NZ_LT629762.1": ["Single", "Type2b"], "NZ_CP031065.1": ["Incomplete"], "NZ_CP031066.1": ["Simple", "Type3"], "NZ_CP031067.1": ["Incomplete"], "VSJH01": ["Simple", "Type2b"], "FNYO01": ["Incomplete"], "NZ_CP036313.1": ["Simple", "Type3"], "NZ_CP036314.1": ["Incomplete"], "NZ_CP036315.1": ["Incomplete"], "NZ_CP041668.1": ["Simple", "Type3"], "NZ_CP041669.1": ["Incomplete"], "FXYF01": ["Simple", "Type3"], "NIRS01": ["Single", "Type2b"], "FNCO01": ["Single", "Type2b"], "RHQN01": ["Multiple", "Type1", "Type2b"], "NZ_CP021983.2": ["Incomplete"], "RCWL01": ["Simple", "Type3"], "QUMQ01": ["Incomplete"], "QAJM01": ["Single", "Type2b"], "NBRZ01": ["Single", "Type2b"], "BJLR01": ["Incomplete"], "NZ_CP039291.1": ["Incomplete"], "NC_013947.1": ["Multiple", "Type3"], "JMCC02": ["Simple", "Type3"], "NBVR01": ["Single", "Type1"], "QAIL01": ["Simple", "Type3"], "QWFB01": ["Single", "Type2b"], "NZ_CP031648.1": ["Simple", "Type2b"], "MCHY01": ["Simple", "Type3"], "NZ_CP041186.1": ["Multiple", "Type3"], "NC_013216.1": ["Simple", "Type3"], "JYHW01": ["Single", "Type2b"], "WIVY01": ["Multiple", "Type1", "Type2b"], "FYDX01": ["Single", "Type2b"], "MUGY01": ["Incomplete"], "FNYJ01": ["Incomplete"], "JJML01": ["Simple", "Type3"], "FNTZ01": ["Multiple", "Type2b"], "NZ_CP029064.1": ["Simple", "Type3"], "LRUN01": ["Single", "Type2b"], "VIUF01": ["Simple", "Type2b"], "VZZK01": ["Simple", "Type3"], "AJLJ01": ["Simple", "Type3"], "CAADIW01": ["Single", "Type1"], "AXVJ01": ["Incomplete"], "VIUC01": ["Multiple", "Type2b"], "AMBZ01": ["Multiple", "Type2b"], "QGGJ01": ["Incomplete"], "VUOC01": ["Simple", "Type3"], "QAAE01": ["Simple", "Type3"], "NCWQ01": ["Incomplete"], "PVTU01": ["Incomplete"], "BBMZ01": ["Single", "Type2b"], "NZ_CP054043.1": ["Simple", "Type1"], "SJSL01": ["Simple", "Type3"], "FQYP01": ["Simple", "Type3"], "NZ_CP011129.1": ["Simple", "Type3"], "NC_012961.1": ["Incomplete"], "NC_012962.1": ["Multiple", "Type1", "Type2b"], "BBXD01": ["Incomplete"], "NZ_CP029196.1": ["Simple", "Type3"], "AKJT01": ["Multiple", "Type2b"], "NVPT01": ["Simple", "Type3"], "BBXG01": ["Incomplete"], "ALVN01": ["Simple", "Type3"], "NJFA02": ["Single", "Type1"], "NC_019738.1": ["Incomplete"], "NC_019739.1": ["Incomplete"], "NC_019740.1": ["Incomplete"], "NC_019741.1": ["Incomplete"], "NC_019742.1": ["Incomplete"], "NC_019743.1": ["Incomplete"], "NC_019760.1": ["Incomplete"], "NC_019761.1": ["Incomplete"], "NC_019762.1": ["Incomplete"], "QPCD01": ["Incomplete"], "QTPO01": ["Simple", "Type3"], "FOEO01": ["Single", "Type2b"], "QWLL01": ["Incomplete"], "QOVA01": ["Single", "Type2b"], "NZ_CP014262.1": ["Multiple", "Type2b"], "FNDJ01": ["Incomplete"], "NZ_AP017422.1": ["Incomplete"], "SNXZ01": ["Simple", "Type3"], "FXWP01": ["Single", "Type1"], "UTBZ01": ["Multiple", "Type2b"], "BCBA01": ["Single", "Type2b"], "VSRQ01": ["Incomplete"], "LFWB01": ["Incomplete"], "QTUB01": ["Simple", "Type2a"], "NZ_CP053584.1": ["Single", "Type2b"], "NZ_CP010897.2": ["Simple", "Type2b"], "NZ_CP010898.2": ["Incomplete"], "WIVZ01": ["Multiple", "Type1", "Type2b"], "NZ_CP013341.1": ["Simple", "Type3"], "JACAQG01": ["Simple", "Type2b"], "FNKR01": ["Simple", "Type3"], "NZ_CP027723.1": ["Single", "Type2b"], "MDEN01": ["Incomplete"], "CVRZ01": ["Simple", "Type1"], "NZ_CP038033.1": ["Incomplete"], "NZ_CP044217.1": ["Incomplete"], "NZ_CP044218.1": ["Simple", "Type3"], "PENV01": ["Simple", "Type1"], "NRQY01": ["Simple", "Type1"], "SISB01": ["Multiple", "Type2b"], "NZ_LT629732.1": ["Simple", "Type3"], "AOCZ01": ["Simple", "Type1"], "NZ_CP039371.1": ["Simple", "Type1"], "NZ_CP039372.1": ["Incomplete"], "JAFA01": ["Incomplete"], "FNOY01": ["Incomplete"], "CABPSP01": ["Simple", "Type2b"], "LGSI01": ["Single", "Type2b"], "VZRB01": ["Simple", "Type3"], "MKWS01": ["Multiple", "Type2b"], "VIUI01": ["Multiple", "Type2b"], "RXOM01": ["Incomplete"], "BCQP01": ["Incomplete"], "SMTE01": ["Simple", "Type1"], "QMEY01": ["Incomplete"], "MBDT01": ["Simple", "Type2b"], "LKPJ01": ["Simple", "Type3"], "OGTP01": ["Simple", "Type3"], "QKTW01": ["Simple", "Type3"], "NC_005773.3": ["Multiple", "Type1", "Type2b"], "NC_007274.1": ["Incomplete"], "NC_007275.1": ["Incomplete"], "NZ_CP048835.1": ["Simple", "Type3"], "NC_010162.1": ["Multiple", "Type3"], "NEVM01": ["Single", "Type1"], "FOUX01": ["Simple", "Type3"], "NZ_CP023526.1": ["Incomplete"], "NZ_CP054422.1": ["Single", "Type2b"], "VOIX01": ["Single", "Type2b"], "VIWA01": ["Simple", "Type3"], "VEBC01": ["Incomplete"], "WIWK01": ["Multiple", "Type1", "Type2b"], "QREK01": ["Simple", "Type3"], "NZ_CM002330.1": ["Single", "Type2b"], "NZ_CM002331.1": ["Incomplete"], "BAHC01": ["Incomplete"], "NZ_CP042968.1": ["Incomplete"], "NZ_CP018049.1": ["Single", "Type2b"], "VZPM01": ["Simple", "Type2b"], "QLIN01": ["Single", "Type2b"], "AUYR01": ["Incomplete"], "NTYK01": ["Incomplete"], "VSFF01": ["Simple", "Type3"], "LRTK01": ["Incomplete"], "ARBP01": ["Simple", "Type3"], "ABCS01": ["Multiple", "Type3"], "BJNF01": ["Simple", "Type3"], "VOQD01": ["Simple", "Type3"], "VIUL01": ["Multiple", "Type2b"], "WHJD01": ["Multiple", "Type1", "Type2b"], "MLFS01": ["Simple", "Type1"], "NZ_CP024900.1": ["Multiple", "Type1", "Type2a", "Type2b"], "NZ_CP009555.1": ["Incomplete"], "NZ_CP009556.1": ["Incomplete"], "NZ_CP013426.1": ["Simple", "Type3"], "NZ_CP013427.1": ["Incomplete"], "NZ_CP013428.1": ["Incomplete"], "NZ_CP013429.1": ["Incomplete"], "POUA01": ["Incomplete"], "AJUL01": ["Simple", "Type3"], "PCOS01": ["Simple", "Type2b"], "QKZA01": ["Simple", "Type1"], "FNQW01": ["Incomplete"], "JADL01": ["Incomplete"], "CABHXE01": ["Simple", "Type1"], "VIKS01": ["Incomplete"], "MOBX01": ["Single", "Type2b"], "QKLR01": ["Simple", "Type3"], "JCLE01": ["Simple", "Type1"], "FSRS01": ["Simple", "Type3"], "NZ_LR134159.1": ["Multiple", "Type2b"], "VCKW01": ["Multiple", "Type3"], "WTCR01": ["Simple", "Type3"], "LLWH01": ["Single", "Type2b"], "NZ_CP027738.1": ["Simple", "Type3"], "QKVL01": ["Single", "Type2b"], "NZ_CP033932.1": ["Simple", "Type3"], "NZ_CM001441.1": ["Incomplete"], "QGTQ01": ["Simple", "Type3"], "RCZD01": ["Simple", "Type1"], "PYLU01": ["Simple", "Type1"], "NZ_CP011288.1": ["Single", "Type2b"], "FPLG01": ["Simple", "Type1"], "NZ_CP012371.1": ["Incomplete"], "NZ_CP022478.1": ["Incomplete"], "NMQR01": ["Multiple", "Type1", "Type2a", "Type2b"], "CTIO01": ["Simple", "Type1"], "VCNG01": ["Multiple", "Type2b"], "NZ_CP007410.1": ["Multiple", "Type2b"], "NKHL01": ["Incomplete"], "MVGR01": ["Simple", "Type3"], "NZ_CP056779.1": ["Incomplete"], "NZ_CP056780.1": ["Multiple", "Type1"], "NZ_CP056781.1": ["Incomplete"], "NZ_CP056782.1": ["Single", "Type2b"], "NQKQ01": ["Single", "Type2b"], "JOGE01": ["Simple", "Type3"], "NZ_CP009533.1": ["Multiple", "Type2b"], "NQKJ01": ["Multiple", "Type1", "Type2b"], "NETK01": ["Simple", "Type3"], "NZ_CP031062.1": ["Incomplete"], "NZ_CP031063.1": ["Simple", "Type3"], "NZ_CP031064.1": ["Incomplete"], "NZ_CP004078.1": ["Simple", "Type3"], "PJZH01": ["Incomplete"], "FNPW01": ["Multiple", "Type1"], "SEUB01": ["Multiple", "Type2b"], "UPHP01": ["Simple", "Type3"], "JNGI01": ["Simple", "Type1"], "UUFD01": ["Incomplete"], "AAWS01": ["Incomplete"], "NZ_CP021659.1": ["Multiple", "Type1"], "NZ_CP021660.1": ["Incomplete"], "NZ_CP021661.1": ["Incomplete"], "NZ_CP021662.1": ["Incomplete"], "MOBP01": ["Single", "Type2b"], "OIFR01": ["Simple", "Type3"], "JSAL01": ["Multiple", "Type2b"], "NZ_CP011104.1": ["Multiple", "Type1", "Type2b"], "MOBI01": ["Simple", "Type2b"], "PUJU01": ["Multiple", "Type1", "Type2a", "Type2b"], "BIFQ01": ["Simple", "Type3"], "NZ_CP025035.2": ["Incomplete"], "LIUV01": ["Single", "Type2b"], "NC_010830.1": ["Simple", "Type3"], "CABPSR01": ["Single", "Type2b"], "CVTM01": ["Single", "Type2b"], "RQJP01": ["Simple", "Type3"], "NZ_CP009288.1": ["Simple", "Type3"], "NZ_CM001025.1": ["Simple", "Type2b"], "MOBT01": ["Single", "Type2b"], "NZ_LR134318.1": ["Multiple", "Type2b"], "ABBN01": ["Incomplete"], "NZ_CP039287.1": ["Incomplete"], "NZ_CP039288.1": ["Simple", "Type3"], "NZ_CP039289.1": ["Incomplete"], "LAIJ01": ["Simple", "Type3"], "LFCV01": ["Simple", "Type3"], "WWJN01": ["Simple", "Type3"], "VZPQ01": ["Incomplete"], "VOBN01": ["Single", "Type2b"], "QGTG01": ["Simple", "Type3"], "AYLO01": ["Simple", "Type3"], "NZ_LT707064.1": ["Single", "Type2b"], "NZ_CP020076.1": ["Incomplete"], "NZ_CP020077.1": ["Incomplete"], "NZ_CP020078.1": ["Incomplete"], "NZ_CP020079.1": ["Incomplete"], "NZ_CP020080.1": ["Simple", "Type1"], "NZ_CP020081.1": ["Incomplete"], "AMRI01": ["Incomplete"], "NZ_LT629705.1": ["Incomplete"], "NRST01": ["Single", "Type2b"], "NZ_CP050291.1": ["Simple", "Type1"], "NZ_CP025263.1": ["Simple", "Type2b"], "FXWM01": ["Simple", "Type3"], "NZ_CP034725.1": ["Multiple", "Type2b"], "MKQR01": ["Incomplete"], "FOCT01": ["Incomplete"], "NUVY01": ["Incomplete"], "MRCJ01": ["Simple", "Type3"], "JUQG01": ["Simple", "Type1"], "LECZ01": ["Single", "Type1"], "MTHI01": ["Single", "Type2b"], "NZ_CP022121.1": ["Simple", "Type3"], "NZ_CM001561.1": ["Single", "Type2b"], "NZ_CP017141.1": ["Multiple", "Type3"], "AZAN01": ["Simple", "Type3"], "AGFX01": ["Incomplete"], "VDCQ01": ["Simple", "Type3"], "QHJL01": ["Simple", "Type3"], "QWEX01": ["Multiple", "Type3"], "LMCT01": ["Incomplete"], "NTTM01": ["Incomplete"], "VZZS01": ["Multiple", "Type1", "Type2b"], "SMFW01": ["Multiple", "Type2b"], "UEXE01": ["Single", "Type1"], "NZ_CP013046.2": ["Incomplete"], "NZ_CP013047.2": ["Simple", "Type1"], "FRBZ01": ["Single", "Type2b"], "AKJH01": ["Single", "Type2b"], "BBLT01": ["Simple", "Type3"], "NBWC01": ["Incomplete"], "NZ_CP007039.1": ["Multiple", "Type2b"], "FMCR01": ["Simple", "Type3"], "VIUK01": ["Single", "Type2b"], "MVHE01": ["Incomplete"], "RCOE01": ["Single", "Type2b"], "QGSY01": ["Simple", "Type3"], "AXWS01": ["Incomplete"], "AYMJ01": ["Single", "Type2b"], "VOBI01": ["Single", "Type2b"], "AKJK01": ["Single", "Type2b"], "FNUD01": ["Single", "Type2b"], "MOBJ01": ["Multiple", "Type2b"], "CAAKGZ01": ["Single", "Type2b"], "FOUB01": ["Simple", "Type3"], "MUXN01": ["Simple", "Type3"], "LKBR01": ["Simple", "Type2b"], "UTVB01": ["Multiple", "Type1", "Type2b"], "PENZ01": ["Simple", "Type1"], "NZ_CP009451.1": ["Single", "Type2b"], "NZ_CP034148.1": ["Incomplete"], "NZ_CP034149.1": ["Incomplete"], "NZ_CP034150.1": ["Incomplete"], "NZ_CP034151.1": ["Simple", "Type1"], "NZ_AP017313.1": ["Simple", "Type3"], "FAOS01": ["Incomplete"], "NZ_CP027727.1": ["Single", "Type2b"], "NZ_CP035319.1": ["Incomplete"], "QAIP01": ["Simple", "Type2b"], "FNBR01": ["Multiple", "Type3"], "AXDH01": ["Single", "Type2b"], "FMVH01": ["Simple", "Type1"], "NZ_CP036488.1": ["Incomplete"], "NZ_CP036489.1": ["Incomplete"], "NZ_CP036490.1": ["Single", "Type2b"], "CAQM01": ["Incomplete"], "LOWA01": ["Simple", "Type2b"], "NZ_CP049044.1": ["Single", "Type2b"], "NZ_CP010896.1": ["Single", "Type2b"], "NC_017168.1": ["Single", "Type1"], "NC_017169.1": ["Incomplete"], "NC_017170.1": ["Incomplete"], "NZ_CP031641.1": ["Single", "Type2b"], "VKDC01": ["Multiple", "Type3"], "JOAG01": ["Incomplete"], "MWQG01": ["Incomplete"], "VDFY01": ["Simple", "Type3"], "ALVK01": ["Multiple", "Type3"], "QFRW01": ["Simple", "Type3"], "BILZ01": ["Simple", "Type3"], "BAXG01": ["Multiple", "Type1"], "MWPQ01": ["Simple", "Type3"], "WIWM01": ["Single", "Type2b"], "FOCU01": ["Single", "Type2b"], "MQZX01": ["Simple", "Type1"], "RKHS01": ["Single", "Type1"], "QHHZ01": ["Incomplete"], "MYFJ01": ["Multiple", "Type1", "Type3", "Type2a", "Type2b"], "NC_016901.1": ["Single", "Type1"], "NC_016905.1": ["Incomplete"], "PEIB01": ["Incomplete"], "MOBQ01": ["Single", "Type2b"], "NXNJ01": ["Single", "Type2b"], "NZ_CP044407.1": ["Incomplete"], "PYBV01": ["Incomplete"], "JABTYG01": ["Multiple", "Type2b"], "NZ_CP042468.1": ["Multiple", "Type1", "Type3"], "NZ_CP014135.1": ["Simple", "Type2b"], "NC_016818.1": ["Multiple", "Type2b"], "NC_016819.1": ["Incomplete"], "NC_016835.1": ["Incomplete"], "NC_017092.1": ["Incomplete"], "MTAX01": ["Simple", "Type3"], "NC_015559.1": ["Simple", "Type3"], "LQRT01": ["Incomplete"], "NZ_LS999839.1": ["Simple", "Type3"], "SOCV01": ["Single", "Type2b"], "ASRX01": ["Single", "Type3"], "NZ_CP044064.1": ["Simple", "Type2b"], "AKJM01": ["Single", "Type2b"], "SMKX01": ["Simple", "Type3"], "CAAJVF01": ["Simple", "Type3"], "VIUJ01": ["Simple", "Type2b"], "LGTC01": ["Incomplete"], "NZ_CP033893.1": ["Incomplete"], "NZ_CP033894.1": ["Single", "Type1"], "NZ_CP033895.1": ["Incomplete"], "JXRA01": ["Simple", "Type3"], "RQPI01": ["Simple", "Type3"], "NZ_CP023695.1": ["Simple", "Type3"], "NZ_LR134335.1": ["Single", "Type1"], "SMJU01": ["Simple", "Type3"], "LMCV01": ["Single", "Type2b"], "PKNM01": ["Simple", "Type1"], "PIQI01": ["Simple", "Type1"], "FZPH01": ["Simple", "Type3"], "WIWB01": ["Single", "Type2b"], "NC_009253.1": ["Simple", "Type3"], "SOZA01": ["Single", "Type2b"], "NZ_LT855380.1": ["Incomplete"], "NZ_CP014947.1": ["Multiple", "Type2b"], "ALVJ01": ["Incomplete"], "NZ_CP013459.1": ["Incomplete"], "NZ_CP013460.1": ["Incomplete"], "NZ_CP013461.1": ["Simple", "Type3"], "NZ_CP048408.1": ["Multiple", "Type2b"], "NZ_CP003181.2": ["Simple", "Type2b"], "VFIO01": ["Single", "Type2b"], "MASS01": ["Incomplete"], "NC_020453.1": ["Simple", "Type3"], "PYUC01": ["Simple", "Type1"], "VEGT01": ["Simple", "Type3"], "MKZO01": ["Single", "Type2b"], "WIWE01": ["Single", "Type2b"], "FMWY01": ["Simple", "Type3"], "MWQL01": ["Incomplete"], "FMVD01": ["Incomplete"], "NZ_CP023969.1": ["Single", "Type2b"], "NZ_CP029608.1": ["Single", "Type2b"], "SMKU01": ["Incomplete"], "FUKJ01": ["Single", "Type3"], "JONO01": ["Multiple", "Type1", "Type2a"], "RAVW01": ["Incomplete"], "PDUD01": ["Incomplete"], "MKMC01": ["Simple", "Type3"], "NC_017448.1": ["Single", "Type3"], "PVZV01": ["Incomplete"], "NZ_CP031069.1": ["Incomplete"], "NZ_CP031070.1": ["Simple", "Type2b"], "NZ_CP023269.1": ["Multiple", "Type2b"], "VLLP01": ["Incomplete"], "NZ_CM001559.1": ["Simple", "Type3"], "NZ_CP029983.1": ["Simple", "Type2b"], "VHKL01": ["Simple", "Type3"], "NZ_CP027218.1": ["Single", "Type2b"], "JPPZ01": ["Simple", "Type3"], "AKJD01": ["Incomplete"], "VCNJ01": ["Multiple", "Type2b"], "NZ_CP013423.1": ["Incomplete"], "NZ_CP013424.1": ["Incomplete"], "NZ_CP013425.1": ["Incomplete"], "SMKK01": ["Simple", "Type3"], "SODH01": ["Incomplete"], "AZSS01": ["Multiple", "Type3"], "JFHN01": ["Simple", "Type1"], "MUNM01": ["Single", "Type2b"], "NC_021492.1": ["Incomplete"], "NC_021500.1": ["Single", "Type1"], "RCSU01": ["Simple", "Type3"], "SMOD01": ["Simple", "Type3"], "NZ_CP042382.1": ["Simple", "Type3"], "NC_008268.1": ["Incomplete"], "NC_008269.1": ["Incomplete"], "NC_008270.1": ["Incomplete"], "NC_008271.1": ["Incomplete"], "JYLE01": ["Multiple", "Type1", "Type2b"], "PYMM01": ["Single", "Type1"], "NZ_CP007699.2": ["Simple", "Type3"], "QAOU01": ["Incomplete"], "WBOI01": ["Multiple", "Type1", "Type2b"], "CAACVJ01": ["Incomplete"], "BJMN01": ["Simple", "Type3"], "SMDG01": ["Simple", "Type1"], "CABIWI01": ["Multiple", "Type1", "Type2b"], "WHZZ01": ["Multiple", "Type1", "Type2a", "Type2b"], "QAOQ01": ["Simple", "Type3"], "RCBZ01": ["Incomplete"], "NZ_CP022303.1": ["Incomplete"], "NZ_CP022304.1": ["Incomplete"], "NZ_CP022305.1": ["Incomplete"], "NZ_CP022306.1": ["Incomplete"], "AQRJ01": ["Simple", "Type3"], "FNGP01": ["Incomplete"], "RJKM01": ["Incomplete"], "PKND01": ["Simple", "Type1"], "FOSU01": ["Incomplete"], "AWZT01": ["Simple", "Type3"], "NZ_CP009458.1": ["Single", "Type1"], "WEGH01": ["Incomplete"], "VZZZ01": ["Simple", "Type3"], "NZ_CP043060.1": ["Multiple", "Type2b"], "VZZR01": ["Incomplete"], "JAAQYP01": ["Simple", "Type3"], "MTBD01": ["Simple", "Type1"], "NZ_CP019686.1": ["Single", "Type1"], "VUAZ01": ["Simple", "Type2b"], "AZXK01": ["Incomplete"], "BBIR01": ["Multiple", "Type2b"], "MWLO01": ["Incomplete"], "QYZD01": ["Simple", "Type2b"], "MOBZ01": ["Single", "Type2b"], "FOBB01": ["Simple", "Type3"], "FMDM01": ["Simple", "Type3"], "NZ_CP019888.1": ["Simple", "Type1"], "AUAX01": ["Simple", "Type3"], "AVEF02": ["Simple", "Type1"], "FNAD01": ["Simple", "Type3"], "BBCC01": ["Simple", "Type1"], "QAOV01": ["Single", "Type2b"], "BAZX01": ["Incomplete"], "NKQZ01": ["Simple", "Type3"], "NZ_CP022960.1": ["Simple", "Type2b"], "SHKK01": ["Simple", "Type3"], "NZ_CP012673.1": ["Multiple", "Type1", "Type3"], "WBKQ01": ["Simple", "Type3"], "LQAL01": ["SImple", "Type2b"], "FCNY02": ["Simple", "Type3"], "VFEU01": ["Simple", "Type2b"], "SHKT01": ["Incomplete"], "BAVR01": ["Incomplete"], "OAOQ01": ["Simple", "Type3"], "MSSW01": ["Simple", "Type3"], "JOGR01": ["Simple", "Type3"], "NZ_CP034337.1": ["Incomplete"], "CAACUY01": ["Incomplete"], "QVNU01": ["Simple", "Type3"], "ASSC01": ["Simple", "Type3"], "NZ_CP028035.1": ["Single", "Type2b"], "NZ_CP028036.1": ["Incomplete"], "NZ_CP028037.1": ["Incomplete"], "NZ_CP028038.1": ["Incomplete"], "NZ_CP028039.1": ["Incomplete"], "NZ_CP028040.1": ["Incomplete"], "NZ_CP028041.1": ["Incomplete"], "NZ_CP028042.1": ["Incomplete"], "NZ_CP057330.1": ["Simple", "Type2b"], "NZ_CP057331.1": ["Incomplete"], "NZ_CP057332.1": ["Incomplete"], "NZ_CP057333.1": ["Incomplete"], "NZ_CP046052.1": ["Incomplete"], "NZ_CP046053.1": ["Incomplete"], "NZ_CP046054.1": ["Incomplete"], "BFCB01": ["Simple", "Type3"], "NZ_CM001489.1": ["Incomplete"], "RBOV01": ["Single", "Type2b"], "NZ_CP030750.1": ["Multiple", "Type1", "Type2b"], "SODV01": ["Simple", "Type3"], "QEKL01": ["Multiple", "Type2b"], "QJRT01": ["Single", "Type2b"], "NZ_LT629778.1": ["Single", "Type2b"], "NZ_CP024634.1": ["Incomplete"], "LVTS01": ["Simple", "Type1"], "NZ_AP018150.1": ["Single", "Type2b"], "CPYD01": ["Single", "Type2a"], "RIAR02": ["Multiple", "Type3"], "NZ_CP010407.1": ["Simple", "Type3"], "NZ_CP010408.1": ["Incomplete"], "AUKP01": ["Simple", "Type3"], "MTSA01": ["Single", "Type2b"], "MOBL01": ["Single", "Type2b"], "NZ_CP028158.1": ["Incomplete"], "MOBY01": ["Multiple", "Type2b"], "NIBV01": ["Multiple", "Type1", "Type2a"], "SJOP01": ["Incomplete"], "VEJO01": ["Single", "Type2b"], "FWXB01": ["Incomplete"], "SMKT01": ["Single", "Type3"], "FPBO01": ["Incomplete"], "NZ_CP012831.1": ["Multiple", "Type2b"], "JAAQYX01": ["Multiple", "Type2b"], "QGSZ01": ["Simple", "Type3"], "QARA01": ["Single", "Type2b"], "RXLQ01": ["Simple", "Type3"], "ANOR01": ["Single", "Type2b"], "PJBP01": ["Simple", "Type3"], "NQKI01": ["Single", "Type2b"], "CBSW01": ["Multiple", "Type1", "Type2a", "Type2b"], "VCKX01": ["Incomplete"], "FNGF01": ["Simple", "Type3"], "NZ_CP010310.2": ["Incomplete"], "WJIE01": ["Multiple", "Type1"], "VAUR01": ["Multiple", "Type2b"], "NC_017956.1": ["Incomplete"], "NC_017957.2": ["Incomplete"], "NC_017958.1": ["Incomplete"], "NC_017959.1": ["Incomplete"], "NC_017966.1": ["Incomplete"], "NZ_CP013368.1": ["Incomplete"], "NZ_CP013369.1": ["Single", "Type2b"], "NZ_CP013370.1": ["Incomplete"], "NZ_CP013371.1": ["Incomplete"], "VFVY01": ["Incomplete"], "NZ_CP023567.1": ["Single", "Type2b"], "NZ_CP023568.1": ["Incomplete"], "RJKL01": ["Incomplete"], "NZ_CP045767.1": ["Simple", "Type2b"], "JOGP01": ["Simple", "Type3"], "CABHPT01": ["Simple", "Type1"], "UWFA01": ["Incomplete"], "QVNQ01": ["Simple", "Type3"], "PJBC01": ["Simple", "Type2b"], "SHKZ01": ["Simple", "Type3"], "JQFM01": ["Incomplete"], "SLVP01": ["Simple", "Type3"], "NZ_CP005969.1": ["Multiple", "Type1", "Type2b"], "NZ_CP021106.3": ["Incomplete"], "JYLH01": ["Single", "Type2b"], "LJCS01": ["Multiple", "Type1"], "NZ_CP017606.1": ["Incomplete"], "NZ_CP017607.1": ["Single", "Type1"], "NZ_CP017608.1": ["Incomplete"], "NZ_CP017609.1": ["Incomplete"], "SOAB01": ["Multiple", "Type3"], "QKLY01": ["Simple", "Type2b"], "JAAQWH01": ["Single", "Type2b"], "AKJS01": ["Single", "Type2b"], "QPDS01": ["Multiple", "Type2b"], "NZ_CP023965.1": ["Single", "Type1"], "NZ_CP029231.1": ["Incomplete"], "NZ_CP029232.1": ["Incomplete"], "NZ_CP029233.1": ["Incomplete"], "NZ_CP029234.1": ["Incomplete"], "NZ_CP029235.1": ["Simple", "Type3"], "LNCD01": ["Simple", "Type3"], "NVDM01": ["Single", "Type1"], "NZ_CP011020.1": ["Multiple", "Type2b"], "PZZQ01": ["Multiple", "Type1"], "NZ_CP011807.3": ["Simple", "Type2b"], "NZ_CP011808.2": ["Incomplete"], "NZ_CP011809.2": ["Incomplete"], "ALJC01": ["Single", "Type2b"], "QJRQ01": ["Incomplete"], "QEOK01": ["Simple", "Type3"], "NZ_CP029618.1": ["Simple", "Type3"], "NZ_CP010016.1": ["Incomplete"], "NZ_CP010017.1": ["Incomplete"], "NZ_CP010018.1": ["Incomplete"], "SOCQ01": ["Single", "Type2b"], "RJKE01": ["Multiple", "Type3"], "QAOO01": ["Incomplete"], "JMCL01": ["Multiple", "Type2b"], "QBKC01": ["Incomplete"], "NZ_CP034335.1": ["Incomplete"], "VDLX02": ["Incomplete"], "SSMR01": ["Multiple", "Type3"], "NSCM01": ["Multiple", "Type1", "Type2a", "Type2b"], "VMSG01": ["Multiple", "Type2b"], "ABCQ01": ["Simple", "Type1"], "OUND01": ["Incomplete"], "QAJI01": ["Simple", "Type3"], "NZ_CP045761.1": ["Single", "Type2b"], "SUNB01": ["Incomplete"], "LJSD01": ["Incomplete"], "NZ_CP041692.1": ["Simple", "Type3"], "PZZR01": ["Incomplete"], "JPMW01": ["Incomplete"], "QPIJ01": ["Simple", "Type3"], "LYUY01": ["Incomplete"], "SMJX01": ["Incomplete"], "VATL01": ["Single", "Type1"], "FMYF01": ["Simple", "Type3"], "PUWV01": ["Multiple", "Type1", "Type2a", "Type2b", "Type3"], "FQYM01": ["Incomplete"], "PJRP01": ["Incomplete"], "QRBE01": ["Incomplete"], "FOVJ01": ["Multiple", "Type3"], "SOCT01": ["Incomplete"], "CABMLW01": ["Multiple", "Type1"], "BDBY01": ["Simple", "Type3"], "PYGV01": ["Incomplete"], "VRLS01": ["Single", "Type1"], "ASTJ01": ["Simple", "Type3"], "LVEJ01": ["Single", "Type2b"], "OUNR01": ["Simple", "Type3"], "FPBP01": ["Single", "Type3"], "FSRU01": ["Simple", "Type2b"], "SMKN01": ["Incomplete"], "ASJB01": ["Simple", "Type3"], "VIYH01": ["Single", "Type2b"], "SNZP01": ["Simple", "Type3"], "NZ_CP014847.1": ["Incomplete"], "NZ_CP014848.1": ["Incomplete"], "NZ_CP014849.1": ["Incomplete"], "NZ_CP014850.1": ["Incomplete"], "NZ_CP014851.1": ["Incomplete"], "NZ_CP014852.1": ["Incomplete"], "NZ_CP014853.1": ["Simple", "Type2b"], "NZ_CP034780.1": ["Incomplete"], "NZ_CP034781.1": ["Incomplete"], "NZ_CP034782.1": ["Incomplete"], "NZ_CP034783.1": ["Simple", "Type3"], "PYBJ01": ["Simple", "Type3"], "PTJB01": ["Incomplete"], "NZ_CP024159.1": ["Incomplete"], "JNYY01": ["Simple", "Type3"], "NZ_CP027756.1": ["Simple", "Type2b"], "SSNZ01": ["Simple", "Type3"], "NZ_CP046874.1": ["Multiple", "Type2b"], "WIBD01": ["Single", "Type2b"], "NZ_CP029710.1": ["Simple", "Type3"], "RBRE01": ["Multiple", "Type1", "Type2b"], "NZ_CP024866.1": ["Single", "Type2b"], "JAAAGD01": ["Single", "Type2b"], "JAAEFD01": ["Simple", "Type1"], "RBUY01": ["Incomplete"], "QXQA01": ["Single", "Type3"], "QJRP01": ["Incomplete"], "AXBA01": ["Simple", "Type3"], "OMPE01": ["Incomplete"], "NZ_LT629790.1": ["Multiple", "Type2b"], "LLWI01": ["Multiple", "Type1", "Type2b"], "NZ_LT629746.1": ["Multiple", "Type2b"], "BAOS01": ["Simple", "Type3"], "VLPL01": ["Simple", "Type3"], "LYRP01": ["Simple", "Type1"], "JXDG01": ["Multiple", "Type2b"], "LIPP01": ["Simple", "Type3"], "JAAQWI01": ["Multiple", "Type2b"], "NZ_LT629795.1": ["Single", "Type2b"], "LXEN01": ["Simple", "Type1"], "NZ_CM001514.1": ["Multiple", "Type2b"], "NC_015731.1": ["Simple", "Type3"], "NZ_CP041754.1": ["Multiple", "Type1", "Type2b"], "NZ_CP041755.1": ["Incomplete"], "NZ_CP041756.1": ["Incomplete"], "NZ_CP045572.1": ["Incomplete"], "FOVO01": ["Incomplete"], "JMKF01": ["Incomplete"], "JACCAY01": ["Simple", "Type3"], "MAID01": ["Incomplete"], "NZ_CP009727.1": ["Incomplete"], "NZ_CP009728.1": ["Incomplete"], "AKXV01": ["Simple", "Type3"], "CZQA01": ["Simple", "Type3"], "NZ_CP054609.1": ["Incomplete"], "NZ_CP054610.1": ["Incomplete"], "NZ_CP054611.1": ["Incomplete"], "NZ_CP054612.1": ["Incomplete"], "NZ_CP054613.1": ["Simple", "Type3"], "NZ_CP014226.1": ["Simple", "Type3"], "CABIVM01": ["Single", "Type2b"], "CDPK01": ["Single", "Type1"], "WIAO01": ["Simple", "Type3"], "WSFG01": ["Multiple", "Type1", "Type2a", "Type2b"], "WIVV01": ["Multiple", "Type1", "Type2b"], "NQKN01": ["Single", "Type2b"], "FUWU01": ["Incomplete"], "UTLV01": ["Single", "Type2b"], "QOIO01": ["Simple", "Type3"], "QGLF01": ["Incomplete"], "VSRO01": ["Multiple", "Type2b"], "WIVQ01": ["Single", "Type2b"], "OBDY01": ["Incomplete"], "SJZK02": ["Single", "Type1"], "NZ_CP029482.1": ["Single", "Type2b"], "LHVN01": ["Single", "Type2b"], "UEXF01": ["Single", "Type1"], "VIUR01": ["Multiple", "Type2b"], "QBJA02": ["Multiple", "Type1", "Type2b"], "LZEX01": ["Multiple", "Type1", "Type2b"], "JYLN01": ["Simple", "Type2b"], "NZ_CP022504.1": ["Single", "Type1"], "LFQK01": ["Single", "Type2b"], "NZ_CP038255.1": ["Simple", "Type3"], "NZ_CP024646.1": ["Simple", "Type2b"], "CVJX01": ["Single", "Type1"], "NZ_CP031146.1": ["Simple", "Type1"], "LACH01": ["Simple", "Type2b"], "NZ_CP011253.3": ["Single", "Type2b"], "NZ_CP011518.2": ["Incomplete"], "NZ_CP011519.2": ["Incomplete"], "VJZE01": ["Incomplete"], "QKYK01": ["Incomplete"], "NZ_CP023525.1": ["Single", "Type1"], "AZUB01": ["Simple", "Type2b"], "JRYA01": ["Multiple", "Type1", "Type2b"], "AUEZ01": ["Simple", "Type3"], "VSFG01": ["Simple", "Type3"], "QOIL01": ["Incomplete"], "SJZD01": ["Incomplete"], "VJWF01": ["Simple", "Type3"], "FNSJ01": ["Multiple", "Type3"], "AWNH01": ["Incomplete"], "NZ_CP016211.1": ["Incomplete"], "SMKO01": ["Incomplete"], "VXLC01": ["Simple", "Type3"], "JFCA01": ["Simple", "Type2b"], "PPRY01": ["Simple", "Type2b"], "NZ_LT629708.1": ["Simple", "Type2b"], "NZ_CP026364.1": ["Single", "Type1"], "NZ_CP015225.1": ["Simple", "Type2b"], "MUNY01": ["Simple", "Type3"], "MBLO01": ["Multiple", "Type3"], "NZ_CP012159.1": ["Single", "Type3"], "ALBT01": ["Simple", "Type3"], "RPND01": ["Simple", "Type3"], "LMXH01": ["Simple", "Type3"], "CPZI01": ["Simple", "Type2b"], "SAXA01": ["Simple", "Type3"], "QPAO01": ["Single", "Type1"], "SSMQ01": ["Multiple", "Type3"], "PQKR01": ["Single", "Type2b"], "NC_015510.1": ["Simple", "Type3"], "CWJL01": ["Single", "Type1"], "BJMH01": ["Simple", "Type3"], "JABWHT01": ["Simple", "Type1"], "NZ_CP024085.1": ["Simple", "Type3"], "SMFY01": ["Incomplete"], "NZ_CP027732.1": ["Single", "Type2b"], "NQKV01": ["Multiple", "Type2b"], "PYGG01": ["Incomplete"], "AWXZ01": ["Simple", "Type3"], "NZ_LN681227.1": ["Multiple", "Type1", "Type2a", "Type2b"], "PQMB01": ["Simple", "Type1"], "FODH01": ["Simple", "Type3"], "FOIG01": ["Simple", "Type3"], "ARBJ01": ["Simple", "Type3"], "NZ_CP029843.1": ["Simple", "Type3"], "FOLC01": ["Simple", "Type3"], "NZ_LT629691.1": ["Multiple", "Type1", "Type2b"], "FORG01": ["Multiple", "Type1"], "WIWO01": ["Single", "Type2b"], "MUBJ01": ["Multiple", "Type2a", "Type2b"], "QPIP01": ["Simple", "Type1"], "QNVV01": ["Simple", "Type3"], "LHVM01": ["Single", "Type2b"], "PVZG01": ["Simple", "Type3"], "NZ_CP017687.1": ["Single", "Type2b"], "FNQB01": ["Multiple", "Type3"], "NQYG01": ["Single", "Type2b"], "JQNB01": ["Incomplete"], "FNVO01": ["Incomplete"], "VIVV01": ["Simple", "Type3"], "CABLBP01": ["Simple", "Type3"], "NZ_CP043422.1": ["Simple", "Type1"], "QEOF01": ["Simple", "Type2b"], "PVNL01": ["Single", "Type3"], "CQAW01": ["Simple", "Type1"], "AUAF01": ["Incomplete"], "SZQA01": ["Incomplete"], "NC_015379.1": ["Multiple", "Type2b"], "NZ_CP026106.1": ["Simple", "Type3"], "VIUG01": ["Multiple", "Type2b"], "SMSC01": ["Simple", "Type1"], "MQWC01": ["Simple", "Type3"], "NZ_CP042804.1": ["Incomplete"], "PENT01": ["Simple", "Type1"], "NC_014718.1": ["Multiple", "Type2b"], "NC_014722.1": ["Incomplete"], "PENW01": ["Simple", "Type1"], "QRAV01": ["Single", "Type2b"], "QEQQ01": ["Single", "Type2b"], "MKCS01": ["Simple", "Type3"], "JPIX01": ["Single", "Type1"], "NJAH01": ["Incomplete"], "JPLA01": ["Incomplete"], "VFPP01": ["Incomplete"], "BBXC01": ["Incomplete"], "JPOH01": ["Simple", "Type3"], "VWSH01": ["Simple", "Type3"], "FNON01": ["Simple", "Type3"], "FAOZ01": ["Simple", "Type3"], "NC_013131.1": ["Simple", "Type3"], "FOOS01": ["Simple", "Type3"], "NZ_CP053682.1": ["Simple", "Type1"], "NZ_CP013949.1": ["Simple", "Type3"], "VZPJ01": ["Multiple", "Type2b"], "MIIW01": ["Single", "Type1"], "JXDI01": ["Single", "Type2b"], "NZ_LT629704.1": ["Single", "Type2b"], "NZ_CP016634.1": ["Multiple", "Type1", "Type2b"], "SLYK01": ["Incomplete"], "NZ_AP020337.1": ["Multiple", "Type1", "Type2b"], "PYAL01": ["Single", "Type2b"], "FYDV01": ["Simple", "Type2b"], "VCBA01": ["Incomplete"], "QFZQ01": ["Single", "Type3"], "SMSL01": ["Simple", "Type3"], "FODZ01": ["Simple", "Type3"], "PYGD01": ["Multiple", "Type3"], "LHVL01": ["Multiple", "Type2b"], "QAIK01": ["Single", "Type2b"], "FODL01": ["Single", "Type2b"], "NZ_LT828648.1": ["Simple", "Type3"], "FUXU01": ["Simple", "Type1"], "VZII01": ["Single", "Type2b"], "QLLL01": ["Simple", "Type3"], "NZ_CP016176.1": ["Incomplete"], "FQWR01": ["Multiple", "Type3"], "NZ_CP029197.1": ["Simple", "Type3"], "VIUE01": ["Multiple", "Type2b"], "VIVC01": ["Single", "Type2b"], "NEJJ01": ["Single", "Type2b"], "VOQB01": ["Multiple", "Type1"], "NC_007954.1": ["Incomplete"], "AEDD01": ["Simple", "Type3"], "NC_020209.1": ["Simple", "Type2b"], "LKBT01": ["Single", "Type2b"], "NZ_LT629801.1": ["Simple", "Type2b"], "JHVK01": ["Single", "Type2b"], "PVNG01": ["Incomplete"], "RBIS01": ["Incomplete"], "FNTT01": ["Simple", "Type2b"], "VCRA01": ["Incomplete"], "ABXF01": ["Simple", "Type3"], "PJCP01": ["Single", "Type2b"], "NZ_CP045158.1": ["Single", "Type1"], "BBXE01": ["Incomplete"], "RKHU01": ["Simple", "Type3"], "LIPN01": ["Simple", "Type3"], "NZ_LT629788.1": ["Single", "Type2b"], "NIBS01": ["Simple", "Type2a"], "JAABNH01": ["Simple", "Type1"], "NZ_CP017599.1": ["Multiple", "Type3"], "QKWJ01": ["Simple", "Type3"], "FNVU01": ["Simple", "Type3"], "RAVX01": ["Simple", "Type3"], "NZ_CP024923.1": ["Simple", "Type3"], "NKFP01": ["Incomplete"], "AEDB02": ["Incomplete"], "NZ_CP038613.1": ["Incomplete"], "NZ_CP018319.1": ["Multiple", "Type2b"], "SLZA01": ["Simple", "Type3"], "WIWL01": ["Single", "Type2b"], "SMFY01_information_Ancylobacter_aquaticus_region_TcB_expanded_[640.9]_4636162_4644109_backward": [""], "SMFY01_information_Ancylobacter_aquaticus_region_A2_expanded_[100.2]_4628824_4636126_backward": [""], "SMFY01_information_Ancylobacter_aquaticus_region_TcC_expanded_[251.6]_4636162_4644109_backward": [""], "QTSW01": ["Single", "Type2b"], "VIVY01": ["Incomplete"], "MOBO01": ["Multiple", "Type2b"], "NZ_CP013997.1": ["Incomplete", "Simple"], "RBIO01": ["Incomplete"], "SNYE01": ["Simple", "Type3"], "NZ_CP038438.1": ["Single", "Type2b"], "SMAS01": ["Single", "Type1"], "VCNA01": ["Simple", "Type3"]}
    # original_classifications = {"NZ_WRXN01000057.1": ["Single", "Simple", "Type3"],
    #                             "NZ_QODI01000091.1": ["Single", "Type1"], "NZ_QAIK01000198.1": ["Single", "Type2b"],
    #                             "NZ_CP009451.1": ["Single", "Type2b"], "NC_015513.1": ["Single", "Type3"],
    #                             "NC_013216.1": ["Single", "Type3"], "NC_022663.1": ["Single", "Type3"],
    #                             "NZ_CP004078.1": ["Single", "Type3", "Simple"], "NZ_CP011809.2": ["Single", "Type2b"],
    #                             "NZ_VITD01000034.1": ["Single", "Type3"], "NZ_WSEM01000042.1": ["Single", "Type3"],
    #                             "NZ_SOEK01000059.1": ["Single", "Type2b"],
    #                             "NZ_MQWC01000010.1": ["Single", "Type3", "Simple"],
    #                             "NZ_PVZA01000046.1": ["Single", "Type2b"], "NZ_BIFQ01000002.1": ["Single", "Type3"],
    #                             "NZ_SSMR01000050.1": ["Multiple", "Type3"], "NZ_QEOF01000027.1": ["Single", "Type2b"],
    #                             "NZ_CP041186.1": ["Multiple", "Type3", ],
    #                             "NZ_NJAK01000005.1": ["Single", "Type2a"], "NZ_FPBP01000034.1": ["Single", "Type3"],
    #                             "NZ_BBMZ01000055.1": ["Single", "Type2b"], "NZ_KI632511.1": ["Single", "Type3"],
    #                             "NZ_CP027760.1": ["Single", "Type2b"], "NZ_CP024793.1": ["Single", "Type3"],
    #                             "NC_005126.1": ["Multiple", "Type2b", "Type1"],
    #                             "NZ_AYSJ01000017.1": ["Multiple", "Type2b", "Type2a", "Type1"],
    #                             "NZ_FPJC01000063.1": ["Single", "Type2b"], "NZ_CP027734.1": ["Single", "Type2b"],
    #                             "NZ_SAXA01000055.1": ["Single", "Type3"],
    #                             "NZ_CP024901.1": ["Multiple", "Type1", "Type2b", "Type2a"],
    #                             "NZ_NMRE01000216.1": ["Single", "Type2b"], "NZ_CP041692.1": ["Single", "Type3"],
    #                             "NC_008271.1": ["Single", "Type3"], "NZ_VCNA01000026.1": ["Single", "Type3"],
    #                             "NZ_NBVR01000044.1": ["Single", "Type1"], "NZ_QVIG01000004.1": ["Single", "Type1"],
    #                             "NZ_BILZ01000167.1": ["Single", "Type3"], "NZ_FONE01000126.1": ["Single", "Type3"],
    #                             "NZ_WIVQ01000359.1": ["Single", "Type2b"], "NZ_UGTQ01000009.1": ["Single", "Type1"],
    #                             "NZ_LN681228.1": ["Multiple", "Type1", "Type2b", "Type2a"],
    #                             "NZ_FOLC01000061.1": ["Single", "Type1"], "NZ_SJOP01000106.1": ["Single", "Type1"],
    #                             "NZ_PHHS01000064.1": ["Multiple", "Type2b"], "JAAHIE010001451.1": ["Single", "Type3"],
    #                             "NZ_KL647038.1": ["Single", "Type3"],
    #                             "NZ_FQUS01000064.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_WIBD01000081.1": ["Single", "Type2b"], "NZ_MWTQ01000158.1": ["Single", "Type2b"],
    #                             "NC_015379.1": ["Multiple", "Type2b"], "NZ_PHSU01000082.1": ["Single", "Type2b"],
    #                             "NZ_FUXU01000241.1": ["Single", "Type1"], "NZ_AKJE01000126.1": ["Multiple", "Type2b"],
    #                             "NZ_LT629762.1": ["Type2b", "Single"], "NZ_SMSC01000216.1": ["Single", "Type1"],
    #                             "NZ_CP013429.1": ["Single", "Type1"], "NZ_AYLO01000184.1": ["Single", "Type3"],
    #                             "NZ_KB905728.1": ["Single", "Type3"], "NZNV01000041.1": ["Single", "Type3"],
    #                             "NZ_POEF01000127.1": ["Single", "Type3"], "NZ_NVPT01000374.1": ["Single", "Type3"],
    #                             "NZ_KN173624.1": ["Single", "Type3"], "NZ_QTUF01000034.1": ["Single", "Type2b"],
    #                             "NZ_KN266223.1": ["Single", "Type3"], "NZ_PENW01000035.1": ["Single", "Type1"],
    #                             "NZ_CBLI010000886.1": [""], "NZ_QAOO01000085.1": ["Single", "Type3"],
    #                             "NZ_LT707062.1": ["Single", "Type2b"], "NZ_SMKX01000348.1": ["Single", "Type3"],
    #                             "NZ_LT855380.1": ["Single", "Type1"], "NZ_RQJP01000018.1": ["Single", "Type3"],
    #                             "NZ_JNYY01000082.1": ["Single", "Type3"], "NZ_CZQA01000015.1": ["Single", "Type3"],
    #                             "NZ_AKJS01000210.1": ["Single", "Type2b"], "NZ_LOYJ01000133.1": ["Single", "Type2b"],
    #                             "NZ_BBCC01000558.1": ["Single", "Type1"],
    #                             "QHVH01000106.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_OAOQ01000062.1": ["Single", "Type3"], "MUGG01000223.1": ["Single", "Type3"],
    #                             "NZ_CP014226.1": ["Single", "Type3"], "NZ_RZUM01000202.1": ["Single", "Type3"],
    #                             "NZ_CP027727.1": ["Single", "Type2b"], "NZ_NHML01000109.1": ["Single", "Type3"],
    #                             "NZ_VCSG01000349.1": ["Single", "Type3"], "NZ_AFWT01000141.1": ["Single", "Type3"],
    #                             "NC_010830.1": ["Single", "Type3"], "RHGL01000187.1": ["Single", "Type3"],
    #                             "NC_013892.1": ["Multiple", "Type2a", "Type1"],
    #                             "NZ_QUMQ01000001.1": ["Single", "Type1"], "NZ_VOBI01000057.1": ["Single", "Type2b"],
    #                             "NZ_FNTY01000002.1": ["Single", "Type2b"], "NZ_SSMQ01000201.1": ["Multiple", "Type3"],
    #                             "NZ_JYLF01000032.1": ["Single", "Type2b"], "NZ_CP029843.1": ["Single", "Type3"],
    #                             "NZ_LT629732.1": ["Single", "Type3"], "NZ_CP017687.1": ["Single", "Type2b"],
    #                             "NZ_ONZJ01000003.1": ["Single", "Type3"], "NZ_AKJH01000183.1": ["Single", "Type2b"],
    #                             "NZ_SMKT01000421.1": ["Single", "Type3"],
    #                             "NZ_QUOK01000045.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_RCOE01000307.1": ["Single", "Type2b"], "NZ_RCOD01000106.1": ["Single", "Type2b"],
    #                             "NZ_SSWO01000061.1": ["Single", "Type1"], "NZ_CP033700.1": ["Single", "Type2b"],
    #                             "NZ_CP015381.1": ["Single", "Type3"], "NZ_WIAO01000088.1": ["Single", "Type3"],
    #                             "NZ_QARE02000057.1": ["Single", "Type2b"], "NZ_MASS01000135.1": ["Single", "Type3"],
    #                             "NC_020453.1": ["Single", "Type1"], "NC_021084.1": ["Single", "Type3"],
    #                             "NC_020418.1": ["Multiple", "Type2b", "Type1"],
    #                             "NZ_WTCR01000039.1": ["Single", "Type3"],
    #                             "NZ_PGEO01000001.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_CP029710.1": ["Single", "Type3"], "NZ_CP011521.2": ["Single", "Type2b"],
    #                             "NZ_KI519434.1": ["Single", "Type3"], "NZ_VHKL01000033.1": ["Single", "Type3"],
    #                             "NZ_SNYU01000018.1": ["Single", "Type2b"], "LKCF01003901.1": ["Single", "Type2b"],
    #                             "NZ_JH941054.1": ["Single", "Type3"], "NZ_AALD02000179.1": ["Single", "Type1"],
    #                             "NZ_KB903530.1": ["Single", "Type3"], "NZ_FNNQ01000042.1": ["Single", "Type3"],
    #                             "NZ_WIVR01000272.1": ["Single", "Type2b"], "NZ_LIUV01000061.1": ["Single", "Type2b"],
    #                             "NZ_CQAZ01000234.1": ["Multiple", "Type3", "Type1"],
    #                             "NZ_VOBM01000105.1": ["Multiple", "Type2b"], "NZ_AP020337.1": ["Single", "Type2b"],
    #                             "NZ_KE386823.1": ["Single", "Type3"], "NZ_RAVY01000454.1": ["Single", "Type3"],
    #                             "NZ_RAVX01000101.1": ["Single", "Type3"], "NZ_CAEB01000078.1": ["Single", "Type1"],
    #                             "NZ_MJMK01000061.1": ["Single", "Type2b"], "NZ_BAXG01000116.1": ["Single", "Type1"],
    #                             "NZ_RXLQ01000081.1": ["Single", "Type3"], "NZ_QFXK01000092.1": ["Single", "Type3"],
    #                             "NZ_FNPW01000042.1": ["Multiple", "Type1", "Type3"],
    #                             "MNDS01000223.1": ["Single", "Type3"], "NZ_LECZ01000027.1": ["Single", "Type1"],
    #                             "NZ_CP024866.1": ["Single", "Type2b"], "NZ_RCNX01000105.1": ["Single", "Type2b"],
    #                             "NZ_FNBR01000014.1": ["Single", "Type2b"], "DLUT01000591.1": ["Single", "Type1"],
    #                             "NZ_QAIL01000624.1": ["Single", "Type3"], "NZ_AKJU01000280.1": ["Single", "Type3"],
    #                             "NZ_KQ058904.1": ["Single", "Type3"], "LADO01000066.1": ["Single", "Type3"],
    #                             "NZ_RHHV01000044.1": ["Single", "Type3"], "NZ_JAABNE010000136.1": ["Single", "Type1"],
    #                             "NZ_AKJT01000140.1": ["Multiple", "Type2b"], "NZ_MKMC01000100.1": ["Single", "Type3"],
    #                             "NZ_CP017600.1": ["Multiple", "Type3"], "NZ_NJAJ01000180.1": ["Single", "Type2b"],
    #                             "NZ_FZPH01000044.1": ["Single", "Type3"], "NZ_PPRY01000171.1": ["Single", "Type2b"],
    #                             "NZ_CP022411.1": ["Single", "Type2b"], "NZ_AMZY02000039.1": ["Single", "Type3"],
    #                             "NZ_KK211074.1": ["Single", "Type2b"], "NZ_RPSA01000026.1": ["Multiple", "Type3"],
    #                             "NZ_VEJO01000043.1": ["Single", "Type2b"], "NZ_CP018049.1": ["Single", "Type2b"],
    #                             "NZ_LGRC01000069.1": ["Single", "Type3"], "NZ_FMCR01000011.1": ["Single", "Type3"],
    #                             "DPJM01000861.1": ["Single", "Type3"], "NC_008027.1": ["Single", "Type2b"],
    #                             "NC_017448.1": ["Single", "Type3"], "NZ_AZNV01000101.1": ["Single", "Type3"],
    #                             "NZ_CP028042.1": ["Single", "Type2b"],
    #                             "NZ_NSCD01000185.1": ["Multiple", "Type1", "Type2b", "Type2a"],
    #                             "NZ_AP017313.1": ["Single", "Type3"], "NZ_QUMR01000030.1": ["Single", "Type2b"],
    #                             "NZ_FODL01000038.1": ["Single", "Type2b"], "QJPH01000597.1": ["Single", "Type3"],
    #                             "NZ_FORG01000094.1": ["Single", "Type2b"], "NZ_CP007231.1": ["Single", "Type2b"],
    #                             "NZ_NEVM01000005.1": ["Single", "Type1"], "NZ_KB235915.1": ["Multiple", "Type3"],
    #                             "NZ_CP009289.1": ["Single", "Type3"], "NZ_VDFY01000321.1": ["Single", "Type3"],
    #                             "NZ_BJMN01000147.1": ["Single", "Type3"],
    #                             "NZ_LOIC01000105.1": ["Multiple", "Type2a", "Type2b", "Type1"],
    #                             "NZ_NIBS01000173.1": ["Single", "Type2a"], "JMHR01000441.1": ["Single", "Type3"],
    #                             "JAABOU010001898.1": ["Single", "Type3"], "MSXM01000117.1": ["Multiple", "Type3"],
    #                             "NZ_VIUK01000070.1": ["Single", "Type1"], "JAAHFV010000687.1": ["Single", "Type3"],
    #                             "NZ_PVZG01000092.1": ["Single", "Type3"], "NZ_KI911557.1": ["Single", "Type3"],
    #                             "NZ_JOIX01000228.1": ["Single", "Type3"], "JEMY01000085.1": ["Single", "Type3"],
    #                             "NZ_CP011104.1": ["Multiple", "Type1", "Type2b"], "NZ_KI421497.1": ["Single", "Type3"],
    #                             "NZ_CP045011.1": ["Single", "Type1"], "NZ_SSNI01000125.1": ["Single", "Type3"],
    #                             "NZ_SZWE01000003.1": ["Single", "Type3"], "NZ_CDSC02000515.1": ["Single", "Type3"],
    #                             "NZ_FMYF01000036.1": ["Single", "Type3"], "NZ_RCWL01000032.1": ["Single", "Type3"],
    #                             "NZ_PHHE01000001.1": ["Single", "Type2b"], "NZ_LKBY01000178.1": ["Single", "Type3"],
    #                             "NZ_AP018150.1": ["Single", "Type2b"], "NZ_SSBS01000012.1": ["Multiple", "Type2b"],
    #                             "NZ_FMXV01000075.1": ["Single", "Type2b"], "NZ_KB897775.1": ["Single", "Type3"],
    #                             "NZ_PENZ01000052.1": ["Single", "Type1"], "NZ_NIRH01000051.1": ["Single", "Type1"],
    #                             "NZ_PENX01000027.1": ["Single", "Type3"], "NZ_CP047651.1": ["Single", "Type2b"],
    #                             "NZ_SNXZ01000022.1": ["Single", "Type3"], "NC_017565.1": ["Single", "Type1"],
    #                             "NZ_SOBU01000009.1": ["Multiple", "Type2a", "Type2b", "Type1"],
    #                             "NZ_FPIZ01000088.1": ["Multiple", "Type3"], "NZ_AP018449.1": ["Single", "Type3"],
    #                             "NZ_MSSW01000136.1": ["Single", "Type3"], "NZ_LR134373.1": ["Single", "Type2b"],
    #                             "NZ_CP038274.1": ["Single", "Type3"], "NZ_ASSC01000896.1": ["Single", "Type3"],
    #                             "NZ_FOSU01000047.1": ["Single", "Type3"], "NZ_LAIJ01000019.1": ["Single", "Type3"],
    #                             "NZ_FUYT01000034.1": ["Single", "Type2b"], "NZ_MBLO01000280.1": ["Single", "Type3"],
    #                             "NZ_KE384514.1": ["Single", "Type3"], "NZ_CP009747.1": ["Single", "Type2b"],
    #                             "NZ_QGSY01000345.1": ["Single", "Type3"], "NZ_QGSZ01000408.1": ["Single", "Type3"],
    #                             "NZ_JH725405.1": ["Single", "Type3"], "NZ_OUNR01000022.1": ["Single", "Type3"],
    #                             "NZ_CP005927.1": ["Single", "Type1"], "NZ_CP043925.1": ["Single", "Type1"],
    #                             "NZ_ASRX01000182.1": ["Single", "Type1"], "NZ_BAHC01000261.1": ["Single", "Type3"],
    #                             "NZ_PVTJ01000022.1": ["Single", "Type3"], "NZ_LR590468.1": ["Single", "Type3"],
    #                             "NZ_CABIVL010000085.1": ["Multiple", "Type2b", "Type1"],
    #                             "NZ_LLWH01000241.1": ["Single", "Type2b"], "NC_017447.1": ["Single", "Type1"],
    #                             "NZ_CP014262.1": ["Multiple", "Type2b"], "NZ_MWPQ01000095.1": ["Single", "Type3"],
    #                             "NZ_RHQN01000027.1": ["Single", "Type2b"], "NZ_CP009533.1": ["Multiple", "Type2b"],
    #                             "NZ_VRLV01000055.1": ["Multiple", "Type3"], "NZ_QAOQ01000022.1": ["Single", "Type3"],
    #                             "NZ_VUAZ01000259.1": ["Single", "Type2b"], "NZ_CP033931.1": ["Single", "Type3"],
    #                             "NZ_CP014947.1": ["Multiple", "Type2b"], "NZ_LS999839.1": ["Single", "Type3"],
    #                             "NZ_RAZO01000515.1": ["Multiple", "Type2b"], "NZ_KI421431.1": ["Single", "Type3"],
    #                             "NZ_CP017141.1": ["Multiple", "Type3"], "NZ_CP027759.1": ["Single", "Type2b"],
    #                             "NZ_QTUH01000043.1": ["Single", "Type2b"], "NZ_PYGD01000024.1": ["Multiple", "Type3"],
    #                             "NZ_BCBA01000109.1": ["Single", "Type2b"], "NC_015559.1": ["Single", "Type3"],
    #                             "NZ_AKJQ01000091.1": ["Single", "Type2b"],
    #                             "NZ_AKJR01000399.1": ["Multiple", "Type1", "Type2b"],
    #                             "NZ_QTPO01000204.1": ["Single", "Type3"], "NZ_PDUD01000143.1": ["Single", "Type3"],
    #                             "JAAAKW010000055.1": ["Single", "Type2b"], "NZ_AMBZ01000025.1": ["Multiple", "Type3"],
    #                             "NZ_WSTC01000100.1": ["Single", "Type3"], "NZ_SOCG01000010.1": ["Single", "Type2b"],
    #                             "NZ_QTPW01000119.1": ["Single", "Type3"],
    #                             "NZ_NMQR01000210.1": ["Multiple", "Type2a", "Type1", "Type2b"],
    #                             "NZ_JAABMA010000050.1": ["Single", "Type1"], "NZ_LT629778.1": ["Single", "Type2b"],
    #                             "NZ_PJBP01000186.1": ["Single", "Type3"], "NZ_QFRW01000331.1": ["Single", "Type3"],
    #                             "NZ_KV791721.1": ["Multiple", "Type1"], "NZ_MKQR01000032.1": ["Single", "Type3"],
    #                             "NZ_CP046054.1": ["Single", "Type3"], "NZ_KL543992.1": ["Multiple", "Type1", "Type2a"],
    #                             "NZ_JXSK01000016.1": ["Multiple", "Type2a", "Type1", "Type2b"],
    #                             "NZ_UPHT01000230.1": ["Single", "Type3", "Single", "Type3"],
    #                             "NZ_MCHY01000013.1": ["Single", "Type3"], "NZ_MUNY01000094.1": ["Single", "Type3"],
    #                             "NZ_NSCM01000192.1": ["Multiple", "Type2a", "Type1", "Type2b"],
    #                             "NZ_RAWG01000802.1": ["Single", "Type3"], "NZ_JYLO01000042.1": ["Single", "Type2b"],
    #                             "NZ_WSQA01000032.1": ["Single", "Type3"], "NZ_CP028923.1": ["Single", "Type3"],
    #                             "NZ_MUBJ01000149.1": ["Single", "Type2a"], "NZ_JXRA01000201.1": ["Single", "Type2b"],
    #                             "NZ_JYLD01000037.1": ["Single", "Type3"], "NZ_BDBY01000492.1": ["Single", "Type3"],
    #                             "NZ_LIPP01000561.1": ["Single", "Type3"], "NZ_AHAM01000375.1": ["Single", "Type3"],
    #                             "NZ_BILY01000094.1": ["Single", "Type3"],
    #                             "NZ_VKDC01000148.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_FPBA01000069.1": ["Multiple", "Type3"], "NZ_BCBN01000125.1": ["Multiple", "Type2b"],
    #                             "NZ_SODV01000004.1": ["Single", "Type3"], "NZ_PPRZ01000380.1": ["Multiple", "Type2b"],
    #                             "NZ_FYEA01000033.1": ["Single", "Type2b"], "NZ_WMBA01000144.1": ["Single", "Type3"],
    #                             "NZ_FNCO01000054.1": ["Single", "Type2b"], "NZ_FMUL01000047.1": ["Single", "Type2b"],
    #                             "NZ_FCON02000657.1": ["Single", "Type3"], "NZ_CP023969.1": ["Single", "Type2b"],
    #                             "NZ_JJML01000118.1": ["Single", "Type3"], "NZ_JAABLV010000025.1": ["Single", "Type1"],
    #                             "NZ_FQWR01000065.1": ["Multiple", "Type3"], "NZ_JAABLU010000039.1": ["Single", "Type1"],
    #                             "NZ_BCBO01000195.1": ["Multiple", "Type1", "Type2b"],
    #                             "NZ_FORB01000034.1": ["Single", "Type3"], "NZ_JYLH01000043.1": ["Single", "Type2b"],
    #                             "NZ_PGGO01000087.1": ["Single", "Type3"], "NZ_LMGQ01000029.1": ["Single", "Type2b"],
    #                             "NZ_JAABNK010000025.1": ["Single", "Type1"], "NZ_KZ679081.1": ["Single", "Type3"],
    #                             "NKIG01000124.1": ["Single", "Type3"], "NZ_LMGK01000026.1": ["Single", "Type2b"],
    #                             "WASQ01000153.1": ["Single", "Type3"], "NZ_BAOS01000047.1": ["Single", "Type3"],
    #                             "NZ_BCQP01000133.1": ["Single", "Type3"], "NZ_CP010898.2": ["Single", "Type2b"],
    #                             "NC_021184.1": ["Single", "Type3"], "NZ_FOZR01000085.1": ["Single", "Type3"],
    #                             "NC_009253.1": ["Single", "Type3"], "NZ_QKLY01000024.1": ["Single", "Type2b"],
    #                             "NZ_LVYD01000134.1": ["Single", "Type3"], "NZ_VFIO01000040.1": ["Single", "Type2b"],
    #                             "QQTZ01000066.1": ["Single", "Type3"], "NC_013947.1": ["Multiple", "Type3"],
    #                             "NZ_VZZS01000049.1": ["Multiple", "Type2b", "Type1"],
    #                             "NZ_FNJL01000093.1": ["Single", "Type3"], "NZ_MVHE01000525.1": ["Single", "Type3"],
    #                             "NZ_FMWY01000062.1": ["Single", "Type3"], "NZ_CP010408.1": ["Single", "Type3"],
    #                             "NZ_LT605205.1": ["Single", "Type3"], "LKBL01002861.1": ["Single", "Type1"],
    #                             "NZ_KB944506.1": ["Single", "Type3"], "NZ_CP029618.1": ["Single", "Type3"],
    #                             "NZ_FPBO01000103.1": ["Single", "Type3"], "NZ_QJUG01000493.1": ["Single", "Type3"],
    #                             "NZ_QAJM01000070.1": ["Single", "Type2b"], "LGGF01000107.1": ["Single", "Type3"],
    #                             "NZ_WUNA01000095.1": ["Single", "Type3"], "NZ_MKCS01000005.1": ["Single", "Type3"],
    #                             "NZ_VCNG01000086.1": ["Multiple", "Type2b"],
    #                             "NZ_ABCS01000237.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_CM002331.1": ["Single", "Type2b"], "NZ_CP023695.1": ["Single", "Type3"],
    #                             "NZ_SMFY01000011.1": ["Single", "Type3"], "NZ_QSNX01000060.1": ["Single", "Type2b"],
    #                             "NZ_LMXH01000018.1": ["Single", "Type3"], "NZ_CP014135.1": ["Single", "Type2b"],
    #                             "NZ_JXDG01000126.1": ["Single", "Type2b"], "NZ_PIQI01000031.1": ["Single", "Type1"],
    #                             "NZ_JYLB01000026.1": ["Multiple", "Type2b"], "NZ_QKTW01000033.1": ["Single", "Type3"],
    #                             "NZ_LOWA01000060.1": ["Single", "Type2b"],
    #                             "JAAHFU010000658.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_CP042382.1": ["Single", "Type3"], "NZ_CP013461.1": ["Single", "Type3"],
    #                             "NZ_LWBP01000264.1": ["Single", "Type3"], "NZ_LYRP01000050.1": ["Single", "Type1"],
    #                             "NZ_SEIT01000119.1": ["Single", "Type2b"], "NZ_JYLN01000037.1": ["Single", "Type2b"],
    #                             "NZ_QTTH01000050.1": ["Single", "Type2b"], "NZ_NITZ01000138.1": ["Single", "Type1"],
    #                             "NZ_RJKE01000001.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_VOLC01000068.1": ["Single", "Type3"], "NZ_LFCV01000266.1": ["Single", "Type1"],
    #                             "NZ_MVIF01000387.1": ["Single", "Type3"], "NZ_VRLS01000101.1": ["Single", "Type3"],
    #                             "NZ_MULM01000185.1": ["Single", "Type3"], "JAABLX010000056.1": ["Single", "Type1"],
    #                             "NC_014723.1": ["Multiple", "Type2b", "Type3"], "NC_016905.1": ["Single", "Type1"],
    #                             "NZ_SOZA01000107.1": ["Single", "Type2b"], "NZ_MJML01000057.1": ["Single", "Type2b"],
    #                             "NZ_WTYM01000063.1": ["Single", "Type3"], "NZ_QOIO01000184.1": ["Single", "Type3"],
    #                             "NZ_FQUQ01000022.1": ["Single", "Type3"], "NZ_FAOZ01000083.1": ["Single", "Type3"],
    #                             "NZ_JNZS01000088.1": ["Single", "Type3"], "NZ_KQ257877.1": ["Single", "Type3"],
    #                             "NZ_KB892704.1": ["Single", "Type3"], "NZ_MUIN01000073.1": ["Multiple", "Type2b"],
    #                             "AP018273.1": ["Single", "Type3"], "NZ_LT985385.1": ["Single", "Type3"],
    #                             "NZ_PYAC01000043.1": ["Single", "Type3"],
    #                             "JAAHGQ010001185.1": ["Multiple", "Type1", "Type3"],
    #                             "NZ_FOBB01000023.1": ["Single", "Type3"], "NC_010162.1": ["Multiple", "Type3"],
    #                             "NZ_JPMW01000007.1": ["Single", "Type3"], "NZ_CPYD01000045.1": ["Single", "Type2a"],
    #                             "NZ_CP021135.1": ["Single", "Type2b"], "MNDA01000671.1": ["Single", "Type3"],
    #                             "NZ_QXQA01000045.1": ["Single", "Type3"], "NZ_OCSV01000008.1": ["Single", "Type3"],
    #                             "NZ_FXWP01000029.1": ["Single", "Type1"], "NZ_AWXZ01000044.1": ["Single", "Type3"],
    #                             "NZ_UYJA01000022.1": ["Single", "Type2b"], "NZ_LNTU01000041.1": ["Single", "Type3"],
    #                             "NZ_QNVV01000061.1": ["Single", "Type3"], "NZ_CP022121.1": ["Single", "Type3"],
    #                             "NZ_SAIQ01000015.1": ["Single", "Type3"], "NZ_VCRA01000094.1": ["Single", "Type3"],
    #                             "NZ_CP029197.1": ["Single", "Type3"], "NZ_PODL01000171.1": ["Single", "Type2b"],
    #                             "NZ_FOAF01000023.1": ["Single", "Type3"], "NZ_QKWJ01000248.1": ["Single", "Type3"],
    #                             "NZ_CP029608.1": ["Single", "Type2b"], "NZ_JFHN01000075.1": ["Single", "Type1"],
    #                             "NZ_FXBM01000005.1": ["Single", "Type3"], "NZ_CP048209.1": ["Single", "Type3"],
    #                             "NZ_VJZE01000939.1": ["Single", "Type3"], "NC_013131.1": ["Single", "Type3"],
    #                             "NZ_JH651384.1": ["Single", "Type3"], "NZ_PYAW01000034.1": ["Single", "Type3"],
    #                             "NZ_WMJZ01000174.1": ["Single", "Type1"], "NZ_SNZP01000034.1": ["Single", "Type3"],
    #                             "NZ_CP010896.1": ["Single", "Type2b"], "NZ_SMJU01000044.1": ["Single", "Type3"],
    #                             "NZ_FAOS01000004.1": ["Single", "Type3"], "NZ_RHLK01000044.1": ["Single", "Type3"],
    #                             "NZ_VSFF01000027.1": ["Single", "Type3"], "NZ_RQPI01000039.1": ["Single", "Type3"],
    #                             "NC_012962.1": ["Multiple", "Type1", "Type2a"],
    #                             "NZ_FSRS01000002.1": ["Single", "Type3"],
    #                             "NZ_CBXF010000164.1": ["Multiple", "Type2a", "Type1"],
    #                             "NZ_QLTF01000036.1": ["Single", "Type2b"],
    #                             "NZ_CACSJQ010000143.1": ["Multiple", "Type2b"],
    #                             "NZ_FNKR01000003.1": ["Single", "Type3"], "NZ_SMSL01000028.1": ["Single", "Type3"],
    #                             "NZ_VZZK01000111.1": ["Single", "Type3"],
    #                             "NZ_LRSO01000023.1": ["Multiple", "Type1", "Type2b"],
    #                             "NZ_CP028272.1": ["Single", "Type1"], "NZ_WJIE01000061.1": ["Multiple", "Type3"],
    #                             "NZ_VDCQ01000167.1": ["Single", "Type3"], "NZ_OGTP01000072.1": ["Single", "Type3"],
    #                             "NZ_MPIN01000042.1": ["Single", "Type3"], "NZ_CDPK01000072.1": ["Single", "Type1"],
    #                             "NZ_CP026364.1": ["Single", "Type1"], "NZ_LXEN01000293.1": ["Single", "Type3"],
    #                             "NZ_CABPSP010000048.1": ["Single", "Type2b"], "NZ_CP019686.1": ["Single", "Type1"],
    #                             "NZ_SJSL01000015.1": ["Single", "Type3"], "CABPSQ010000036.1": ["Single", "Type2b"],
    #                             "JAAHFO010001461.1": ["Single", "Type3"], "NZ_SOCQ01000042.1": ["Single", "Type2b"],
    #                             "MNJJ01000259.1": ["Single", "Type3"], "NZ_CP012159.1": ["Single", "Type3"],
    #                             "NZ_QLTJ01000050.1": ["Single", "Type2b"], "NZ_JNWO01000211.1": ["Single", "Type3"],
    #                             "NZ_CP013341.1": ["Single", "Type3"], "NC_017807.1": ["Single", "Type2b"],
    #                             "NZ_PYBV01000203.1": ["Single", "Type3"], "NZ_KE332397.1": ["Single", "Type3"],
    #                             "NZ_RCSU01000043.1": ["Single", "Type3"],
    #                             "NZ_FNQB01000011.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_QLLL01000025.1": ["Single", "Type3"], "NZ_QTUB01000001.1": ["Single", "Type2a"],
    #                             "NZ_NSCE01000184.1": ["Multiple", "Type1", "Type2b", "Type2a"],
    #                             "NZ_JAAGLX010000268.1": ["Single", "Type3"], "NZ_VOQB01000043.1": ["Multiple", "Type1"],
    #                             "NZ_FODH01000041.1": ["Single", "Type3"], "NZ_FNVU01000044.1": ["Single", "Type3"],
    #                             "NZ_FXAS01000142.1": ["Single", "Type1"], "NZ_CP038255.1": ["Single", "Type3"],
    #                             "NZ_QVNU01000027.1": ["Single", "Type3"], "NZ_VOIW01000021.1": ["Single", "Type2b"],
    #                             "NZ_LT629795.1": ["Single", "Type1"], "NZ_CP026110.1": ["Single", "Type3"],
    #                             "NZ_AEDD01000040.1": ["Single", "Type3"], "NZ_FNAD01000036.1": ["Single", "Type3"],
    #                             "NZ_PVZV01000026.1": ["Single", "Type3"], "NZ_PVNL01000192.1": ["Single", "Type3"],
    #                             "NZ_LXYR01000213.1": ["Single", "Type3"], "NZ_LN623556.1": ["Single", "Type2b"],
    #                             "NZ_FNGF01000016.1": ["Single", "Type3"], "NZ_CP022961.1": ["Single", "Type3"],
    #                             "NZ_CXPG01000027.1": ["Single", "Type3"], "NZ_CP017236.1": ["Single", "Type1"],
    #                             "NZ_CP012332.1": ["Single", "Type3"],
    #                             "NZ_FTPM01000002.1": ["Multiple", "Type2b", "Type1"],
    #                             "NZ_BLAD01000170.1": ["Multiple", "Type3", "Type3"],
    #                             "NZ_FTPJ01000003.1": ["Multiple", "Type1", "Type2b"],
    #                             "NZ_NIRS01000013.1": ["Single", "Type2b"], "NZ_FMDM01000037.1": ["Single", "Type3"],
    #                             "NZ_PTJB01000052.1": ["Single", "Type3"], "NZ_JAAFZB010000096.1": ["Single", "Type3"],
    #                             "NZ_CP016211.1": ["Single", "Type3"], "NZ_PQKR01000051.1": ["Single", "Type2b"],
    #                             "NZ_SOAB01000040.1": ["Multiple", "Type3"], "NZ_NCXP01000142.1": ["Single", "Type3"],
    #                             "NZ_ANMG01000154.1": ["Single", "Type3"], "NC_020209.1": ["Single", "Type2b"],
    #                             "NZ_JZSQ01000140.1": ["Single", "Type3"], "NZ_LT629705.1": ["Single", "Type2b"],
    #                             "NZ_PYAL01000013.1": ["Single", "Type2b"], "JAAHIG010000776.1": ["Multiple", "Type3"],
    #                             "MKSF01000039.1": ["Single", "Type3"], "LAQJ01000315.1": ["Single", "Type3"],
    #                             "NZ_SMKU01000805.1": ["Single", "Type3"], "NZ_LT828648.1": ["Single", "Type3"],
    #                             "NZ_CP007215.2": ["Single", "Type3"], "NZ_WNKZ01000359.1": ["Single", "Type3"],
    #                             "NZ_LR590482.1": ["Single", "Type2b"], "NZ_LT907981.1": ["Single", "Type3"],
    #                             "NZ_QAIP01000588.1": ["Single", "Type1"], "NZ_LNCD01000152.1": ["Single", "Type3"],
    #                             "NZ_KE384562.1": ["Single", "Type3"], "NZ_LJCS01000262.1": ["Multiple", "Type1"],
    #                             "NZ_ATXB01000005.1": ["Single", "Type3"], "NZ_SMKK01000563.1": ["Single", "Type3"],
    #                             "NC_019762.1": ["Single", "Type3"], "NZ_FNTZ01000002.1": ["Multiple", "Type2b"],
    #                             "NZ_QFZQ01000023.1": ["Multiple", "Type3"], "NZ_JOGP01000180.1": ["Single", "Type3"],
    #                             "KZ266893.1": ["Single", "Type3"], "NZ_FNON01000025.1": ["Single", "Type3"],
    #                             "NZ_SHKK01000001.1": ["Single", "Type3"], "NZ_FNUD01000002.1": ["Single", "Type1"],
    #                             "NZ_FQYP01000028.1": ["Single", "Type3"], "NZ_QGTQ01000078.1": ["Single", "Type1"],
    #                             "NZ_JFJW01000247.1": ["Single", "Type1"], "NZ_FOVS01000095.1": ["Single", "Type2b"],
    #                             "NZ_CP012540.1": ["Single", "Type3"], "NZ_JUHO01000001.1": ["Single", "Type2b"],
    #                             "NZ_CP007039.1": ["Multiple", "Type2b"], "DNUG01000139.1": ["Single", "Type3"],
    #                             "NZ_CP038630.1": ["Single", "Type3"], "NC_013954.1": ["Single", "Type1"],
    #                             "NZ_VCKW01000623.1": ["Multiple", "Type3"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_A1_expanded_930570_933549_forward": [
    #                                 "Type2b"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_A2_expanded_933600_938013_forward": [
    #                                 "Type2b"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_TcdA1_expanded_1144367_1151012_backward": [
    #                                 "Type1"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_TcdA1_expanded_1136560_1144054_backward": [
    #                                 "Type1"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_TcdA1_expanded_1119871_1127122_backward": [
    #                                 "Type1"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_TcdA1_expanded_1107514_1115125_backward": [
    #                                 "Type1"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_A2_expanded_2891462_2895557_backward": [
    #                                 "Type2a"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_A1_expanded_2895543_2899062_backward": [
    #                                 "Type2a"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_A2_expanded_4869364_4874056_backward": [
    #                                 "Type2b"],
    #                             "NC_005126.1_information_Photorhabdus_laumondii_region_A1_expanded_4874153_4877051_backward": [
    #                                 "Type2b"],
    #                             "NZ_AYSJ01000017.1_information_Photorhabdus_khanii_region_A2_expanded_1807955_1812674_backward": [
    #                                 "Type2b"],
    #                             "NZ_AYSJ01000017.1_information_Photorhabdus_khanii_region_A1_expanded_1812775_1815703_backward": [
    #                                 "Type2b"],
    #                             "NZ_AYSJ01000017.1_information_Photorhabdus_khanii_region_TcdA1_expanded_1746224_1753757_backward": [
    #                                 "Type1"],
    #                             "NZ_AYSJ01000017.1_information_Photorhabdus_khanii_region_A1_expanded_2879332_2882866_forward": [
    #                                 "Type2a"],
    #                             "NZ_AYSJ01000017.1_information_Photorhabdus_khanii_region_A2_expanded_2882852_2886944_forward": [
    #                                 "Type2a"],
    #                             "NZ_AYSJ01000017.1_information_Photorhabdus_khanii_region_TcdA1_expanded_4632407_4639025_backward": [
    #                                 "Type1"],
    #                             "NZ_AYSJ01000017.1_information_Photorhabdus_khanii_region_TcdA1_expanded_4624260_4631886_backward": [
    #                                 "Type1"],
    #                             "NZ_NJAK01000005.1_information_Xenorhabdus_ishibashii_region_A2_expanded_700848_705030_forward": [
    #                                 "Type2a"],
    #                             "NZ_NJAK01000005.1_information_Xenorhabdus_ishibashii_region_A1_expanded_697368_700839_forward": [
    #                                 "Type2a"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_A1_expanded_930569_933548_forward": [
    #                                 "Type2b"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_A2_expanded_933599_938012_forward": [
    #                                 "Type2b"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_TcdA1_expanded_1144366_1151011_backward": [
    #                                 "Type1"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_TcdA1_expanded_1136559_1144053_backward": [
    #                                 "Type1"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_TcdA1_expanded_1119870_1127121_backward": [
    #                                 "Type1"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_TcdA1_expanded_1107513_1115124_backward": [
    #                                 "Type1"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_A2_expanded_4868065_4872757_backward": [
    #                                 "Type2b"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_A1_expanded_4872854_4875752_backward": [
    #                                 "Type2b"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_A2_expanded_2891462_2895557_backward": [
    #                                 "Type2a"],
    #                             "NZ_CP024901.1_information_Photorhabdus_laumondii_region_A1_expanded_2895543_2899062_backward": [
    #                                 "Type2a"],
    #                             "NZ_LN681228.1_information_Xenorhabdus_nematophila_region_TcdA1_expanded_2110775_2118188_backward": [
    #                                 "Type1"],
    #                             "NZ_LN681228.1_information_Xenorhabdus_nematophila_region_TcdA1_expanded_2546301_2553873_backward": [
    #                                 "Type1"],
    #                             "NZ_LN681228.1_information_Xenorhabdus_nematophila_region_TcdA1_expanded_2530974_2538543_forward": [
    #                                 "Type1"],
    #                             "NZ_LN681228.1_information_Xenorhabdus_nematophila_region_A1_expanded_2518328_2521799_forward": [
    #                                 "Type2a"],
    #                             "NZ_LN681228.1_information_Xenorhabdus_nematophila_region_A2_expanded_2521808_2525981_forward": [
    #                                 "Type2a"],
    #                             "NZ_LN681228.1_information_Xenorhabdus_nematophila_region_A1_expanded_2272553_2275616_forward": [
    #                                 "Type2b"],
    #                             "NZ_LN681228.1_information_Xenorhabdus_nematophila_region_A2_expanded_2275690_2280313_forward": [
    #                                 "Type2b"],
    #                             "NC_013892.1_information_Xenorhabdus_bovienii_region_A2_expanded_576007_580111_backward": [
    #                                 "Type2a"],
    #                             "NC_013892.1_information_Xenorhabdus_bovienii_region_A1_expanded_580106_583658_backward": [
    #                                 "Type2a"],
    #                             "NC_013892.1_information_Xenorhabdus_bovienii_region_TcdA1_expanded_1524550_1532101_forward": [
    #                                 "Type1"],
    #                             "NZ_QUOK01000045.1_information_Rhodohalobacter_sp._region_TcdA1_expanded_2159005_2169691_forward": [
    #                                 "Type3"],
    #                             "NZ_QUOK01000045.1_information_Rhodohalobacter_sp._region_A2_expanded_2141911_2145157_forward": [
    #                                 "Type?"],
    #                             "NC_020418.1_information_Morganella_morganii_region_TcdA1_expanded_1872868_1880275_backward": [
    #                                 "Type1"],
    #                             "NC_020418.1_information_Morganella_morganii_region_A1_expanded_1374832_1378387_forward": [
    #                                 "Type2b"],
    #                             "NC_020418.1_information_Morganella_morganii_region_A2_expanded_1382402_1385276_forward": [
    #                                 "Type2b"],
    #                             "NC_020418.1_information_Morganella_morganii_region_TcdA1_expanded_2005504_2009848_backward": [
    #                                 "Type1"],
    #                             "NZ_PGEO01000001.1_information_Streptomyces_sp._region_A2_expanded_8169184_8172529_backward": [
    #                                 "Type?"],
    #                             "NZ_PGEO01000001.1_information_Streptomyces_sp._region_TcdA1_expanded_375542_385436_forward": [
    #                                 "Type3"],
    #                             "NZ_CQAZ01000234.1_information_Yersinia_pekkanenii_region_TcdA1_expanded_2059147_2062405_forward": [
    #                                 "Type1"],
    #                             "NZ_CQAZ01000234.1_information_Yersinia_pekkanenii_region_TcdA1_expanded_2064365_2070638_forward": [
    #                                 "Type3"],
    #                             "NZ_FNPW01000042.1_information_Pseudomonas_sp._region_TcdA1_expanded_6614186_6621635_backward": [
    #                                 "Type1"],
    #                             "NZ_FNPW01000042.1_information_Pseudomonas_sp._region_TcdA1_expanded_6628816_6636343_forward": [
    #                                 "Type3"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_TcdA1_expanded_843495_850140_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_TcdA1_expanded_835688_843182_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_TcdA1_expanded_818999_826250_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_TcdA1_expanded_806642_814253_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_A2_expanded_1481951_1486364_forward": [
    #                                 "Type2b"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_A1_expanded_1478921_1481900_forward": [
    #                                 "Type2b"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_A1_expanded_3354173_3357071_forward": [
    #                                 "Type2b"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_A2_expanded_3357168_3361860_forward": [
    #                                 "Type2b"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_A1_expanded_3627693_3631212_backward": [
    #                                 "Type2a"],
    #                             "NZ_NSCD01000185.1_information_Photorhabdus_sp._region_A2_expanded_3623612_3627707_backward": [
    #                                 "Type2a"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_TcdA1_expanded_36104_43622_forward": [
    #                                 "Type1"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_TcdA1_expanded_2349716_2356853_forward": [
    #                                 "Type1"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_TcdA1_expanded_3626745_3634320_backward": [
    #                                 "Type1"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_TcdA1_expanded_3617993_3624551_backward": [
    #                                 "Type1"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_TcdA1_expanded_3610392_3617688_backward": [
    #                                 "Type1"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_TcdA1_expanded_3588748_3595918_backward": [
    #                                 "Type1"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_TcdA1_expanded_3576411_3583962_backward": [
    #                                 "Type1"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_TcdA1_expanded_4750346_4757417_forward": [
    #                                 "Type1"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A1_expanded_1746473_1750022_backward": [
    #                                 "Type2a"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A2_expanded_1742386_1746487_backward": [
    #                                 "Type2a"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A2_expanded_2814033_2818173_forward": [
    #                                 "Type2a"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A1_expanded_2810554_2814010_forward": [
    #                                 "Type2a"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A1_expanded_2099000_2101886_forward": [
    #                                 "Type2b"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A2_expanded_2101936_2106331_forward": [
    #                                 "Type2b"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A1_expanded_4202590_4205512_forward": [
    #                                 "Type2b"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A2_expanded_4205609_4210301_forward": [
    #                                 "Type2b"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A1_expanded_4742654_4746626_forward": [
    #                                 "Type2b"],
    #                             "NZ_LOIC01000105.1_information_Photorhabdus_namnaonensis_region_A2_expanded_4746784_4749883_forward": [
    #                                 "Type2b"],
    #                             "NZ_CP011104.1_information_Photorhabdus_thracensis_region_TcdA1_expanded_2899156_2906767_backward": [
    #                                 "Type1"],
    #                             "NZ_CP011104.1_information_Photorhabdus_thracensis_region_A1_expanded_3227753_3231725_forward": [
    #                                 "Type2b"],
    #                             "NZ_CP011104.1_information_Photorhabdus_thracensis_region_A2_expanded_3231822_3236538_forward": [
    #                                 "Type2b"],
    #                             "NZ_CP011104.1_information_Photorhabdus_thracensis_region_TcdA1_expanded_3191309_3197951_forward": [
    #                                 "Type1"],
    #                             "NZ_CP011104.1_information_Photorhabdus_thracensis_region_TcdA1_expanded_3183237_3190896_forward": [
    #                                 "Type1"],
    #                             "NZ_CP011104.1_information_Photorhabdus_thracensis_region_TcdA1_expanded_3173661_3181236_forward": [
    #                                 "Type1"],
    #                             "NZ_SOBU01000009.1_information_Photorhabdus_temperata_region_A1_expanded_582667_586192_forward": [
    #                                 "Type2a"],
    #                             "NZ_SOBU01000009.1_information_Photorhabdus_temperata_region_A2_expanded_587016_590268_forward": [
    #                                 "Type2a"],
    #                             "NZ_SOBU01000009.1_information_Photorhabdus_temperata_region_A1_expanded_953038_955966_forward": [
    #                                 "Type2b"],
    #                             "NZ_SOBU01000009.1_information_Photorhabdus_temperata_region_A2_expanded_956061_960756_forward": [
    #                                 "Type2b"],
    #                             "NZ_SOBU01000009.1_information_Photorhabdus_temperata_region_TcdA1_expanded_1464525_1471155_backward": [
    #                                 "Type1"],
    #                             "NZ_SOBU01000009.1_information_Photorhabdus_temperata_region_TcdA1_expanded_1456693_1464112_backward": [
    #                                 "Type1"],
    #                             "NZ_CABIVL010000085.1_information_Pseudomonas_carnis_region_A1_expanded_2592488_2594525_forward": [
    #                                 "Type2b"],
    #                             "NZ_CABIVL010000085.1_information_Pseudomonas_carnis_region_A2_expanded_2594556_2599239_forward": [
    #                                 "Type2b"],
    #                             "NZ_CABIVL010000085.1_information_Pseudomonas_carnis_region_TcdA1_expanded_3215854_3218848_forward": [
    #                                 "Type1"],
    #                             "NZ_RHQN01000027.1_information_Pseudomonas_sp._region_A1_expanded_2865132_2868165_backward": [
    #                                 "Type2b"],
    #                             "NZ_RHQN01000027.1_information_Pseudomonas_sp._region_A2_expanded_2860228_2865118_backward": [
    #                                 "Type2b"],
    #                             "NZ_AKJR01000399.1_information_Pseudomonas_sp._region_TcdA1_expanded_267894_271566_forward": [
    #                                 "Type1"],
    #                             "NZ_AKJR01000399.1_information_Pseudomonas_sp._region_A1_expanded_857988_861717_forward": [
    #                                 "Type2b"],
    #                             "NZ_AKJR01000399.1_information_Pseudomonas_sp._region_A2_expanded_861716_865820_forward": [
    #                                 "Type2b"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_A1_expanded_90073_93523_forward": [
    #                                 "Type2a"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_A2_expanded_93547_97681_forward": [
    #                                 "Type2a"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_A1_expanded_2232559_2236093_forward": [
    #                                 "Type2a"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_A2_expanded_2236079_2240180_forward": [
    #                                 "Type2a"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_TcdA1_expanded_1581866_1589459_forward": [
    #                                 "Type1"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_TcdA1_expanded_1569951_1577154_forward": [
    #                                 "Type1"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_TcdA1_expanded_1553426_1560176_forward": [
    #                                 "Type1"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_TcdA1_expanded_1546321_1552942_forward": [
    #                                 "Type1"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_TcdA1_expanded_1538340_1545900_forward": [
    #                                 "Type1"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_TcdA1_expanded_4121789_4132319_backward": [
    #                                 "Type1"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_A2_expanded_2026472_2031185_backward": [
    #                                 "Type2b"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_A1_expanded_2031288_2034213_backward": [
    #                                 "Type2b"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_A1_expanded_2787806_2790692_forward": [
    #                                 "Type2b"],
    #                             "NZ_NMQR01000210.1_information_Photorhabdus_sp._region_A2_expanded_2790743_2795159_forward": [
    #                                 "Type2b"],
    #                             "NZ_KL543992.1_information_Photorhabdus_australis_region_A2_expanded_129145_133240_backward": [
    #                                 "Type2a"],
    #                             "NZ_KL543992.1_information_Photorhabdus_australis_region_A1_expanded_133229_136754_backward": [
    #                                 "Type2a"],
    #                             "NZ_KL543992.1_information_Photorhabdus_australis_region_TcdA1_expanded_672137_679775_forward": [
    #                                 "type1"],
    #                             "NZ_KL543992.1_information_Photorhabdus_australis_region_TcdA1_expanded_1950214_1957753_backward": [
    #                                 "type1"],
    #                             "NZ_KL543992.1_information_Photorhabdus_australis_region_TcdA1_expanded_1943255_1949780_backward": [
    #                                 "type1"],
    #                             "NZ_KL543992.1_information_Photorhabdus_australis_region_TcdA1_expanded_1935583_1942777_backward": [
    #                                 "type1"],
    #                             "NZ_KL543992.1_information_Photorhabdus_australis_region_TcdA1_expanded_3899266_3905800_backward": [
    #                                 "type1"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_A1_expanded_432044_435578_forward": [
    #                                 "Type2a"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_A2_expanded_435564_439665_forward": [
    #                                 "Type2a"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_TcdA1_expanded_557184_564402_backward": [
    #                                 "Type1"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_TcdA1_expanded_2390999_2398565_forward": [
    #                                 "Type1"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_TcdA1_expanded_2369647_2377138_forward": [
    #                                 "Type1"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_TcdA1_expanded_2362258_2369164_forward": [
    #                                 "Type1"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_A2_expanded_2524015_2528428_backward": [
    #                                 "Type2b"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_A1_expanded_2528479_2531365_backward": [
    #                                 "Type2b"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_A2_expanded_3937610_3942299_backward": [
    #                                 "Type2b"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_A1_expanded_3942402_3945327_backward": [
    #                                 "Type2b"],
    #                             "NZ_JXSK01000016.1_information_Photorhabdus_luminescens_region_A2_expanded_4049378_4052948_forward": [
    #                                 "Type1"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_A2_expanded_191734_195832_backward": [
    #                                 "Type2a"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_A1_expanded_195818_199337_backward": [
    #                                 "Type2a"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_TcdA1_expanded_2764973_2771600_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_TcdA1_expanded_2757244_2764660_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_TcdA1_expanded_2734525_2741752_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_TcdA1_expanded_2722185_2729760_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_A2_expanded_3540851_3545546_backward": [
    #                                 "Type2b"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_A1_expanded_3545643_3548565_backward": [
    #                                 "Type2b"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_A1_expanded_3817262_3820166_forward": [
    #                                 "Type2b"],
    #                             "NZ_NSCM01000192.1_information_Photorhabdus_bodei_region_A2_expanded_3820341_3824757_forward": [
    #                                 "Type2b"],
    #                             "NZ_VKDC01000148.1_information_Fulvivirga_sp._region_A2_expanded_7833204_7836423_backward": [
    #                                 "Type?", "Type3"],
    #                             "NZ_VKDC01000148.1_information_Fulvivirga_sp._region_TcdA1_expanded_2351126_2360498_forward": [
    #                                 "Type3"],
    #                             "NZ_VKDC01000148.1_information_Fulvivirga_sp._region_TcdA1_expanded_4853485_4862983_forward": [
    #                                 "Type3"],
    #                             "NZ_FPBA01000069.1_information_Geodermatophilus_amargosae_region_A2_expanded_4821940_4825159_backward": [
    #                                 "Type?"],
    #                             "NZ_FPBA01000069.1_information_Geodermatophilus_amargosae_region_A2_expanded_5630192_5640827_backward": [
    #                                 "Type3"],
    #                             "NZ_BCBO01000195.1_information_Pseudomonas_sp._region_A1_expanded_1940890_1943362_forward": [
    #                                 "Type2b"],
    #                             "NZ_BCBO01000195.1_information_Pseudomonas_sp._region_A2_expanded_1943554_1948237_forward": [
    #                                 "Type2b"],
    #                             "NZ_BCBO01000195.1_information_Pseudomonas_sp._region_TcdA1_expanded_6024212_6028535_forward": [
    #                                 "Type1"],
    #                             "NZ_VZZS01000049.1_information_Pseudomonas_sp._region_A1_expanded_1401257_1405886_forward": [
    #                                 "Type2b"],
    #                             "NZ_VZZS01000049.1_information_Pseudomonas_sp._region_A2_expanded_1405931_1410509_forward": [
    #                                 "Type2b"],
    #                             "NZ_VZZS01000049.1_information_Pseudomonas_sp._region_TcdA1_expanded_1382544_1391286_backward": [
    #                                 "Type1"],
    #                             "NC_014723.1_information_Paraburkholderia_rhizoxinica_region_TcdA1_expanded_263892_267498_forward": [
    #                                 "Type1"],
    #                             "NC_014723.1_information_Paraburkholderia_rhizoxinica_region_A2_expanded_279466_283009_forward": [
    #                                 "Type2b"],
    #                             "NC_014723.1_information_Paraburkholderia_rhizoxinica_region_A1_expanded_275920_279082_forward": [
    #                                 "Type2b"],
    #                             "NC_014723.1_information_Paraburkholderia_rhizoxinica_region_A2_expanded_2388313_2390908_backward": [
    #                                 "Type2b"],
    #                             "NC_014723.1_information_Paraburkholderia_rhizoxinica_region_A1_expanded_2393199_2397900_backward": [
    #                                 "Type2b"],
    #                             "JAAHGQ010001185.1_information_Symploca_sp._region_TcdA1_expanded_1614659_1629566_backward": [
    #                                 "Type1"],
    #                             "JAAHGQ010001185.1_information_Symploca_sp._region_TcdA1_expanded_2763629_2771996_backward": [
    #                                 "Type1"],
    #                             "JAAHGQ010001185.1_information_Symploca_sp._region_A2_expanded_622115_630707_forward": [
    #                                 "Type3", "Type3"],
    #                             "NC_012962.1_information_Photorhabdus_asymbiotica_region_TcdA1_expanded_1054728_1061256_backward": [
    #                                 "Type1"],
    #                             "NC_012962.1_information_Photorhabdus_asymbiotica_region_TcdA1_expanded_1045698_1052889_backward": [
    #                                 "Type1"],
    #                             "NC_012962.1_information_Photorhabdus_asymbiotica_region_A2_expanded_2339002_2343097_forward": [
    #                                 "Type2a"],
    #                             "NC_012962.1_information_Photorhabdus_asymbiotica_region_A1_expanded_2335464_2339013_forward": [
    #                                 "Type2a"],
    #                             "NZ_CBXF010000164.1_information_Xenorhabdus_szentirmaii_region_TcdA1_expanded_3670905_3678513_backward": [
    #                                 "Type1"],
    #                             "NZ_CBXF010000164.1_information_Xenorhabdus_szentirmaii_region_TcdA1_expanded_3655617_3663132_forward": [
    #                                 "Type1"],
    #                             "NZ_CBXF010000164.1_information_Xenorhabdus_szentirmaii_region_TcdA1_expanded_3599091_3606933_forward": [
    #                                 "Type1"],
    #                             "NZ_CBXF010000164.1_information_Xenorhabdus_szentirmaii_region_A1_expanded_3645585_3649056_forward": [
    #                                 "Type2a"],
    #                             "NZ_CBXF010000164.1_information_Xenorhabdus_szentirmaii_region_A2_expanded_3649065_3653265_forward": [
    #                                 "Type2a"], "NZ_FXYF01000056.1": ["Single", "Type3"],
    #                             "NZ_VCKW01000623.1_information_Actinomadura_sp._region_A2_expanded_747005_750221_forward": [
    #                                 "Type?"],
    #                             "NZ_VCKW01000623.1_information_Actinomadura_sp._region_A2_expanded_6267139_6277399_backward": [
    #                                 "Type3"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_TcdA1_expanded_556407_563052_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_TcdA1_expanded_548600_556094_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_TcdA1_expanded_531911_539162_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_TcdA1_expanded_519554_527165_backward": [
    #                                 "Type1"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_A2_expanded_1184450_1188863_forward": [
    #                                 "Type2b"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_A1_expanded_1181420_1184399_forward": [
    #                                 "Type2b"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_A2_expanded_3344032_3348724_forward": [
    #                                 "Type2b"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_A1_expanded_3341037_3343935_forward": [
    #                                 "Type2b"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_A1_expanded_3676266_3679785_backward": [
    #                                 "Type2a"],
    #                             "NZ_NSCE01000184.1_information_Photorhabdus_sp._region_A2_expanded_3672185_3676280_backward": [
    #                                 "Type2a"],
    #                             "NZ_FTPM01000002.1_information_Burkholderia_sp._region_A1_expanded_577789_581314_backward": [
    #                                 "Type2b"],
    #                             "NZ_FTPM01000002.1_information_Burkholderia_sp._region_A2_expanded_572951_575870_backward": [
    #                                 "Type2b"],
    #                             "NZ_FTPM01000002.1_information_Burkholderia_sp._region_A2_expanded_2487802_2491402_forward": [
    #                                 "Type2b", "Type3"],
    #                             "NZ_FTPM01000002.1_information_Burkholderia_sp._region_A2_expanded_2470787_2475539_forward": [
    #                                 "Type2b"],
    #                             "NZ_FTPM01000002.1_information_Burkholderia_sp._region_A1_expanded_2467223_2470688_forward": [
    #                                 "Type2b"],
    #                             "NZ_FTPM01000002.1_information_Burkholderia_sp._region_TcdA1_expanded_3110262_3117861_backward": [
    #                                 "Type1"],
    #                             "NZ_FTPJ01000003.1_information_Burkholderia_sp._region_TcdA1_expanded_21233_28058_forward": [
    #                                 "Type1"],
    #                             "NZ_FTPJ01000003.1_information_Burkholderia_sp._region_A1_expanded_1746754_1750222_backward": [
    #                                 "Type2b"],
    #                             "NZ_FTPJ01000003.1_information_Burkholderia_sp._region_A2_expanded_1741903_1746655_backward": [
    #                                 "Type2b"],
    #                             "NZ_CP041186.1_information_Bradymonas_sp._region_A2_expanded_4319135_4327433_backward": [
    #                                 "Type3"],
    #                             "NZ_ABCS01000237.1_information_Plesiocystis_pacifica_region_A2_expanded_2871032_2875556_forward": [
    #                                 "Type3"],
    #                             "JAAHFU010000658.1_information_Leptolyngbya_sp._region_A2_expanded_7706383_7715437_backward": [
    #                                 "Type3"],
    #                             "JAAHFU010000658.1_information_Leptolyngbya_sp._region_TcdA1_expanded_5813429_5827880_backward": [
    #                                 "Type3"],
    #                             "NZ_RJKE01000001.1_information_Actinocorallia_herbida_region_A2_expanded_4657113_4660347_forward": [
    #                                 "Type3"],
    #                             "NZ_LRSO01000023.1_information_Pseudomonas_thivervalensis_region_A2_expanded_2622756_2627292_backward": [
    #                                 "Type2b"],
    #                             "NZ_LRSO01000023.1_information_Pseudomonas_thivervalensis_region_A1_expanded_2627328_2629860_backward": [
    #                                 "Type2b"],
    #                             "NZ_BLAD01000170.1_information_Acrocarpospora_corrugata_region_A2_expanded_2785887_2789091_backward": [
    #                                 "Type3"],
    #                             "NZ_CP004078.1_information_Paenibacillus_sabinae_region_TcB_expanded_2796003_2803788_backward": [
    #                                 "hit_tag", "more"]}

    # original_classifications = {"NZ_WRXN01000057.1": "Type3", "NZ_QODI01000091.1": "Type1", "NZ_QAIK01000198.1":
    #     "Type2b",
    #                    "NZ_CP009451.1": "Type2b", "NC_015513.1": "Type3", "NC_013216.1": "Type3",
    #                    "NC_022663.1": "Type3", "NZ_CP004078.1": "Type3", "NZ_CP011809.2": "Type2b",
    #                    "NZ_VITD01000034.1": "Type3", "NZ_WSEM01000042.1": "Type3", "NZ_SOEK01000059.1": "Type2b",
    #                    "NZ_MQWC01000010.1": "Type3", "NZ_PVZA01000046.1": "Type2b", "NZ_BIFQ01000002.1": "Type3",
    #                    "NZ_QEOF01000027.1": "Type2b", "NZ_NJAK01000005.1": "Type2a", "NZ_FPBP01000034.1": "Type3",
    #                    "NZ_BBMZ01000055.1": "Type2b", "NZ_KI632511.1": "Type3", "NZ_CP027760.1": "Type2b",
    #                    "NZ_CP024793.1": "Type3", "NZ_FPJC01000063.1": "Type2b", "NZ_CP027734.1": "Type2b",
    #                    "NZ_SAXA01000055.1": "Type3", "NZ_NMRE01000216.1": "Type2b", "NZ_CP041692.1": "Type3",
    #                    "NC_008271.1": "Type3", "NZ_VCNA01000026.1": "Type3", "NZ_NBVR01000044.1": "Type1",
    #                    "NZ_QVIG01000004.1": "Type1", "NZ_BILZ01000167.1": "Type3", "NZ_FONE01000126.1": "Type3",
    #                    "NZ_WIVQ01000359.1": "Type2b", "NZ_UGTQ01000009.1": "Type1", "NZ_FOLC01000061.1": "Type1",
    #                    "NZ_SJOP01000106.1": "Type1", "JAAHIE010001451.1": "Type3", "NZ_KL647038.1": "Type3",
    #                    "NZ_WIBD01000081.1": "Type2b", "NZ_MWTQ01000158.1": "Type2b", "NZ_PHSU01000082.1": "Type2b",
    #                    "NZ_FUXU01000241.1": "Type1", "NZ_LT629762.1": "Type2b ", "NZ_SMSC01000216.1": "Type1",
    #                    "NZ_CP013429.1": "Type1", "NZ_AYLO01000184.1": "Type3", "NZ_KB905728.1": "Type3",
    #                    "NZNV01000041.1": "Type3", "NZ_POEF01000127.1": "Type3", "NZ_NVPT01000374.1": "Type3",
    #                    "NZ_KN173624.1": "Type3", "NZ_QTUF01000034.1": "Type2b", "NZ_KN266223.1": "Type3",
    #                    "NZ_PENW01000035.1": "Type1", "NZ_QAOO01000085.1": "Type3", "NZ_LT707062.1": "Type2b",
    #                    "NZ_SMKX01000348.1": "Type3", "NZ_LT855380.1": "Type1", "NZ_RQJP01000018.1": "Type3",
    #                    "NZ_JNYY01000082.1": "Type3", "NZ_CZQA01000015.1": "Type3", "NZ_AKJS01000210.1": "Type2b",
    #                    "NZ_LOYJ01000133.1": "Type2b", "NZ_BBCC01000558.1": "Type1", "NZ_OAOQ01000062.1": "Type3",
    #                    "MUGG01000223.1": "Type3", "NZ_CP014226.1": "Type3", "NZ_RZUM01000202.1": "Type3",
    #                    "NZ_CP027727.1": "Type2b", "NZ_NHML01000109.1": "Type3", "NZ_VCSG01000349.1": "Type3",
    #                    "NZ_AFWT01000141.1": "Type3", "NC_010830.1": "Type3", "RHGL01000187.1": "Type3",
    #                    "NZ_QUMQ01000001.1": "Type1", "NZ_VOBI01000057.1": "Type2b", "NZ_FNTY01000002.1": "Type2b",
    #                    "NZ_JYLF01000032.1": "Type2b", "NZ_CP029843.1": "Type3", "NZ_LT629732.1": "Type3",
    #                    "NZ_CP017687.1": "Type2b", "NZ_ONZJ01000003.1": "Type3", "NZ_AKJH01000183.1": "Type2b",
    #                    "NZ_SMKT01000421.1": "Type3", "NZ_RCOE01000307.1": "Type2b", "NZ_RCOD01000106.1": "Type2b",
    #                    "NZ_SSWO01000061.1": "Type1", "NZ_CP033700.1": "Type2b", "NZ_CP015381.1": "Type3",
    #                    "NZ_WIAO01000088.1": "Type3", "NZ_QARE02000057.1": "Type2b", "NZ_MASS01000135.1": "Type3",
    #                    "NC_020453.1": "Type1", "NC_021084.1": "Type3", "NZ_WTCR01000039.1": "Type3",
    #                    "NZ_CP029710.1": "Type3", "NZ_CP011521.2": "Type2b", "NZ_KI519434.1": "Type3",
    #                    "NZ_VHKL01000033.1": "Type3", "NZ_SNYU01000018.1": "Type2b", "LKCF01003901.1": "Type2b",
    #                    "NZ_JH941054.1": "Type3", "NZ_AALD02000179.1": "Type1", "NZ_KB903530.1": "Type3",
    #                    "NZ_FNNQ01000042.1": "Type3", "NZ_WIVR01000272.1": "Type2b", "NZ_LIUV01000061.1": "Type2b",
    #                    "NZ_AP020337.1": "Type2b", "NZ_KE386823.1": "Type3", "NZ_RAVY01000454.1": "Type3",
    #                    "NZ_RAVX01000101.1": "Type3", "NZ_CAEB01000078.1": "Type1", "NZ_MJMK01000061.1": "Type2b",
    #                    "NZ_BAXG01000116.1": "Type1", "NZ_RXLQ01000081.1": "Type3", "NZ_QFXK01000092.1": "Type3",
    #                    "MNDS01000223.1": "Type3", "NZ_LECZ01000027.1": "Type1", "NZ_CP024866.1": "Type2b",
    #                    "NZ_RCNX01000105.1": "Type2b", "NZ_FNBR01000014.1": "Type2b", "DLUT01000591.1": "Type1",
    #                    "NZ_QAIL01000624.1": "Type3", "NZ_AKJU01000280.1": "Type3", "NZ_KQ058904.1": "Type3",
    #                    "LADO01000066.1": "Type3", "NZ_RHHV01000044.1": "Type3", "NZ_JAABNE010000136.1": "Type1",
    #                    "NZ_MKMC01000100.1": "Type3", "NZ_NJAJ01000180.1": "Type2b", "NZ_FZPH01000044.1": "Type3",
    #                    "NZ_PPRY01000171.1": "Type2b", "NZ_CP022411.1": "Type2b", "NZ_AMZY02000039.1": "Type3",
    #                    "NZ_KK211074.1": "Type2b", "NZ_VEJO01000043.1": "Type2b", "NZ_CP018049.1": "Type2b",
    #                    "NZ_LGRC01000069.1": "Type3", "NZ_FMCR01000011.1": "Type3", "DPJM01000861.1": "Type3",
    #                    "NC_008027.1": "Type2b", "NC_017448.1": "Type3", "NZ_AZNV01000101.1": "Type3",
    #                    "NZ_CP028042.1": "Type2b", "NZ_AP017313.1": "Type3", "NZ_QUMR01000030.1": "Type2b",
    #                    "NZ_FODL01000038.1": "Type2b", "QJPH01000597.1": "Type3", "NZ_FORG01000094.1": "Type2b",
    #                    "NZ_CP007231.1": "Type2b", "NZ_NEVM01000005.1": "Type1", "NZ_CP009289.1": "Type3",
    #                    "NZ_VDFY01000321.1": "Type3", "NZ_BJMN01000147.1": "Type3", "NZ_NIBS01000173.1": "Type2a",
    #                    "JMHR01000441.1": "Type3", "JAABOU010001898.1": "Type3", "NZ_VIUK01000070.1": "Type1",
    #                    "JAAHFV010000687.1": "Type3", "NZ_PVZG01000092.1": "Type3", "NZ_KI911557.1": "Type3",
    #                    "NZ_JOIX01000228.1": "Type3", "JEMY01000085.1": "Type3", "NZ_KI421497.1": "Type3",
    #                    "NZ_CP045011.1": "Type1", "NZ_SSNI01000125.1": "Type3", "NZ_SZWE01000003.1": "Type3",
    #                    "NZ_CDSC02000515.1": "Type3", "NZ_FMYF01000036.1": "Type3", "NZ_RCWL01000032.1": "Type3",
    #                    "NZ_PHHE01000001.1": "Type2b", "NZ_LKBY01000178.1": "Type3", "NZ_AP018150.1": "Type2b",
    #                    "NZ_FMXV01000075.1": "Type2b", "NZ_KB897775.1": "Type3", "NZ_PENZ01000052.1": "Type1",
    #                    "NZ_NIRH01000051.1": "Type1", "NZ_PENX01000027.1": "Type3", "NZ_CP047651.1": "Type2b",
    #                    "NZ_SNXZ01000022.1": "Type3", "NC_017565.1": "Type1", "NZ_AP018449.1": "Type3",
    #                    "NZ_MSSW01000136.1": "Type3", "NZ_LR134373.1": "Type2b", "NZ_CP038274.1": "Type3",
    #                    "NZ_ASSC01000896.1": "Type3", "NZ_FOSU01000047.1": "Type3", "NZ_LAIJ01000019.1": "Type3",
    #                    "NZ_FUYT01000034.1": "Type2b", "NZ_MBLO01000280.1": "Type3", "NZ_KE384514.1": "Type3",
    #                    "NZ_CP009747.1": "Type2b", "NZ_QGSY01000345.1": "Type3", "NZ_QGSZ01000408.1": "Type3",
    #                    "NZ_JH725405.1": "Type3", "NZ_OUNR01000022.1": "Type3", "NZ_CP005927.1": "Type1",
    #                    "NZ_CP043925.1": "Type1", "NZ_ASRX01000182.1": "Type1", "NZ_BAHC01000261.1": "Type3",
    #                    "NZ_PVTJ01000022.1": "Type3", "NZ_LR590468.1": "Type3", "NZ_LLWH01000241.1": "Type2b",
    #                    "NC_017447.1": "Type1", "NZ_MWPQ01000095.1": "Type3", "NZ_RHQN01000027.1": "Type2b",
    #                    "NZ_QAOQ01000022.1": "Type3", "NZ_VUAZ01000259.1": "Type2b", "NZ_CP033931.1": "Type3",
    #                    "NZ_LS999839.1": "Type3", "NZ_KI421431.1": "Type3", "NZ_CP027759.1": "Type2b",
    #                    "NZ_QTUH01000043.1": "Type2b", "NZ_BCBA01000109.1": "Type2b", "NC_015559.1": "Type3",
    #                    "NZ_AKJQ01000091.1": "Type2b", "NZ_QTPO01000204.1": "Type3", "NZ_PDUD01000143.1": "Type3",
    #                    "JAAAKW010000055.1": "Type2b", "NZ_WSTC01000100.1": "Type3", "NZ_SOCG01000010.1": "Type2b",
    #                    "NZ_QTPW01000119.1": "Type3", "NZ_JAABMA010000050.1": "Type1", "NZ_LT629778.1": "Type2b",
    #                    "NZ_PJBP01000186.1": "Type3", "NZ_QFRW01000331.1": "Type3", "NZ_MKQR01000032.1": "Type3",
    #                    "NZ_CP046054.1": "Type3", "NZ_UPHT01000230.1": "Type3", "NZ_MCHY01000013.1": "Type3",
    #                    "NZ_MUNY01000094.1": "Type3", "NZ_RAWG01000802.1": "Type3", "NZ_JYLO01000042.1": "Type2b",
    #                    "NZ_WSQA01000032.1": "Type3", "NZ_CP028923.1": "Type3", "NZ_MUBJ01000149.1": "Type2a",
    #                    "NZ_JXRA01000201.1": "Type2b", "NZ_JYLD01000037.1": "Type3", "NZ_BDBY01000492.1": "Type3",
    #                    "NZ_LIPP01000561.1": "Type3", "NZ_AHAM01000375.1": "Type3", "NZ_BILY01000094.1": "Type3",
    #                    "NZ_SODV01000004.1": "Type3", "NZ_FYEA01000033.1": "Type2b", "NZ_WMBA01000144.1": "Type3",
    #                    "NZ_FNCO01000054.1": "Type2b", "NZ_FMUL01000047.1": "Type2b", "NZ_FCON02000657.1": "Type3",
    #                    "NZ_CP023969.1": "Type2b", "NZ_JJML01000118.1": "Type3", "NZ_JAABLV010000025.1": "Type1",
    #                    "NZ_JAABLU010000039.1": "Type1", "NZ_FORB01000034.1": "Type3", "NZ_JYLH01000043.1": "Type2b",
    #                    "NZ_PGGO01000087.1": "Type3", "NZ_LMGQ01000029.1": "Type2b", "NZ_JAABNK010000025.1": "Type1",
    #                    "NZ_KZ679081.1": "Type3", "NKIG01000124.1": "Type3", "NZ_LMGK01000026.1": "Type2b",
    #                    "WASQ01000153.1": "Type3", "NZ_BAOS01000047.1": "Type3", "NZ_BCQP01000133.1": "Type3",
    #                    "NZ_CP010898.2": "Type2b", "NC_021184.1": "Type3", "NZ_FOZR01000085.1": "Type3",
    #                    "NC_009253.1": "Type3", "NZ_QKLY01000024.1": "Type2b", "NZ_LVYD01000134.1": "Type3",
    #                    "NZ_VFIO01000040.1": "Type2b", "QQTZ01000066.1": "Type3", "NZ_FNJL01000093.1": "Type3",
    #                    "NZ_MVHE01000525.1": "Type3", "NZ_FMWY01000062.1": "Type3", "NZ_CP010408.1": "Type3",
    #                    "NZ_LT605205.1": "Type3", "LKBL01002861.1": "Type1", "NZ_KB944506.1": "Type3",
    #                    "NZ_CP029618.1": "Type3", "NZ_FPBO01000103.1": "Type3", "NZ_QJUG01000493.1": "Type3",
    #                    "NZ_QAJM01000070.1": "Type2b", "LGGF01000107.1": "Type3", "NZ_WUNA01000095.1": "Type3",
    #                    "NZ_MKCS01000005.1": "Type3", "NZ_CM002331.1": "Type2b", "NZ_CP023695.1": "Type3",
    #                    "NZ_SMFY01000011.1": "Type3", "NZ_QSNX01000060.1": "Type2b", "NZ_LMXH01000018.1": "Type3",
    #                    "NZ_CP014135.1": "Type2b", "NZ_JXDG01000126.1": "Type2b", "NZ_PIQI01000031.1": "Type1",
    #                    "NZ_QKTW01000033.1": "Type3", "NZ_LOWA01000060.1": "Type2b", "NZ_CP042382.1": "Type3",
    #                    "NZ_CP013461.1": "Type3", "NZ_LWBP01000264.1": "Type3", "NZ_LYRP01000050.1": "Type1",
    #                    "NZ_SEIT01000119.1": "Type2b", "NZ_JYLN01000037.1": "Type2b", "NZ_QTTH01000050.1": "Type2b",
    #                    "NZ_NITZ01000138.1": "Type1", "NZ_VOLC01000068.1": "Type3", "NZ_LFCV01000266.1": "Type1",
    #                    "NZ_MVIF01000387.1": "Type3", "NZ_VRLS01000101.1": "Type3", "NZ_MULM01000185.1": "Type3",
    #                    "JAABLX010000056.1": "Type1", "NC_016905.1": "Type1", "NZ_SOZA01000107.1": "Type2b",
    #                    "NZ_MJML01000057.1": "Type2b", "NZ_WTYM01000063.1": "Type3", "NZ_QOIO01000184.1": "Type3",
    #                    "NZ_FQUQ01000022.1": "Type3", "NZ_FAOZ01000083.1": "Type3", "NZ_JNZS01000088.1": "Type3",
    #                    "NZ_KQ257877.1": "Type3", "NZ_KB892704.1": "Type3", "AP018273.1": "Type3",
    #                    "NZ_LT985385.1": "Type3", "NZ_PYAC01000043.1": "Type3", "NZ_FOBB01000023.1": "Type3",
    #                    "NZ_JPMW01000007.1": "Type3", "NZ_CPYD01000045.1": "Type2a", "NZ_CP021135.1": "Type2b",
    #                    "MNDA01000671.1": "Type3", "NZ_QXQA01000045.1": "Type3", "NZ_OCSV01000008.1": "Type3",
    #                    "NZ_FXWP01000029.1": "Type1", "NZ_AWXZ01000044.1": "Type3", "NZ_UYJA01000022.1": "Type2b",
    #                    "NZ_LNTU01000041.1": "Type3", "NZ_QNVV01000061.1": "Type3", "NZ_CP022121.1": "Type3",
    #                    "NZ_SAIQ01000015.1": "Type3", "NZ_VCRA01000094.1": "Type3", "NZ_CP029197.1": "Type3",
    #                    "NZ_PODL01000171.1": "Type2b", "NZ_FOAF01000023.1": "Type3", "NZ_QKWJ01000248.1": "Type3",
    #                    "NZ_CP029608.1": "Type2b", "NZ_JFHN01000075.1": "Type1", "NZ_FXBM01000005.1": "Type3",
    #                    "NZ_CP048209.1": "Type3", "NZ_VJZE01000939.1": "Type3", "NC_013131.1": "Type3",
    #                    "NZ_JH651384.1": "Type3", "NZ_PYAW01000034.1": "Type3", "NZ_WMJZ01000174.1": "Type1",
    #                    "NZ_SNZP01000034.1": "Type3", "NZ_CP010896.1": "Type2b", "NZ_SMJU01000044.1": "Type3",
    #                    "NZ_FAOS01000004.1": "Type3", "NZ_RHLK01000044.1": "Type3", "NZ_VSFF01000027.1": "Type3",
    #                    "NZ_RQPI01000039.1": "Type3", "NZ_FSRS01000002.1": "Type3", "NZ_QLTF01000036.1": "Type2b",
    #                    "NZ_FNKR01000003.1": "Type3", "NZ_SMSL01000028.1": "Type3", "NZ_VZZK01000111.1": "Type3",
    #                    "NZ_CP028272.1": "Type1", "NZ_VDCQ01000167.1": "Type3", "NZ_OGTP01000072.1": "Type3",
    #                    "NZ_MPIN01000042.1": "Type3", "NZ_CDPK01000072.1": "Type1", "NZ_CP026364.1": "Type1",
    #                    "NZ_LXEN01000293.1": "Type3", "NZ_CABPSP010000048.1": "Type2b", "NZ_CP019686.1": "Type1",
    #                    "NZ_SJSL01000015.1": "Type3", "CABPSQ010000036.1": "Type2b", "JAAHFO010001461.1": "Type3",
    #                    "NZ_SOCQ01000042.1": "Type2b", "MNJJ01000259.1": "Type3", "NZ_CP012159.1": "Type3",
    #                    "NZ_QLTJ01000050.1": "Type2b", "NZ_JNWO01000211.1": "Type3", "NZ_CP013341.1": "Type3",
    #                    "NC_017807.1": "Type2b", "NZ_PYBV01000203.1": "Type3", "NZ_KE332397.1": "Type3",
    #                    "NZ_RCSU01000043.1": "Type3", "NZ_QLLL01000025.1": "Type3", "NZ_QTUB01000001.1": "Type2a",
    #                    "NZ_JAAGLX010000268.1": "Type3", "NZ_FODH01000041.1": "Type3", "NZ_FNVU01000044.1": "Type3",
    #                    "NZ_FXAS01000142.1": "Type1", "NZ_CP038255.1": "Type3", "NZ_QVNU01000027.1": "Type3",
    #                    "NZ_VOIW01000021.1": "Type2b", "NZ_LT629795.1": "Type1", "NZ_CP026110.1": "Type3",
    #                    "NZ_AEDD01000040.1": "Type3", "NZ_FNAD01000036.1": "Type3", "NZ_PVZV01000026.1": "Type3",
    #                    "NZ_PVNL01000192.1": "Type3", "NZ_LXYR01000213.1": "Type3", "NZ_LN623556.1": "Type2b",
    #                    "NZ_FNGF01000016.1": "Type3", "NZ_CP022961.1": "Type3", "NZ_CXPG01000027.1": "Type3",
    #                    "NZ_CP017236.1": "Type1", "NZ_CP012332.1": "Type3", "NZ_NIRS01000013.1": "Type2b",
    #                    "NZ_FMDM01000037.1": "Type3", "NZ_PTJB01000052.1": "Type3", "NZ_JAAFZB010000096.1": "Type3",
    #                    "NZ_CP016211.1": "Type3", "NZ_PQKR01000051.1": "Type2b", "NZ_NCXP01000142.1": "Type3",
    #                    "NZ_ANMG01000154.1": "Type3", "NC_020209.1": "Type2b", "NZ_JZSQ01000140.1": "Type3",
    #                    "NZ_LT629705.1": "Type2b", "NZ_PYAL01000013.1": "Type2b", "MKSF01000039.1": "Type3",
    #                    "LAQJ01000315.1": "Type3", "NZ_SMKU01000805.1": "Type3", "NZ_LT828648.1": "Type3",
    #                    "NZ_CP007215.2": "Type3", "NZ_WNKZ01000359.1": "Type3", "NZ_LR590482.1": "Type2b",
    #                    "NZ_LT907981.1": "Type3", "NZ_QAIP01000588.1": "Type1", "NZ_LNCD01000152.1": "Type3",
    #                    "NZ_KE384562.1": "Type3", "NZ_ATXB01000005.1": "Type3", "NZ_SMKK01000563.1": "Type3",
    #                    "NC_019762.1": "Type3", "NZ_JOGP01000180.1": "Type3", "KZ266893.1": "Type3",
    #                    "NZ_FNON01000025.1": "Type3", "NZ_SHKK01000001.1": "Type3", "NZ_FNUD01000002.1": "Type1",
    #                    "NZ_FQYP01000028.1": "Type3", "NZ_QGTQ01000078.1": "Type1", "NZ_JFJW01000247.1": "Type1",
    #                    "NZ_FOVS01000095.1": "Type2b", "NZ_CP012540.1": "Type3", "NZ_JUHO01000001.1": "Type2b",
    #                    "DNUG01000139.1": "Type3", "NZ_CP038630.1": "Type3", "NC_013954.1": "Type1",
    #                    "NZ_FXYF01000056.1": "Type3"}




    # skip_tags = ['Single', 'Multiple']
    #
    # skip_tags = ['Simple']

    test_results = ""

    count = 0
    diff_count = 0
    for query in queries:
        if query.name in original_classifications:
            count += 1

            new_classification = models.GenomeTags.objects.get(tag_id=query.name)

            set1 = set(original_classifications[query.name])

            check = list(set(new_classification.tags).difference(set(original_classifications[query.name]))) + list(set(
                original_classifications[query.name]).difference(set(new_classification.tags)))

            check_skip = [x for x in check if x not in skip_tags]

            if check_skip:
                diff_count += 1
                print("\n\nFound an automatically classified genome that was different\n")
                print(query.name)
                print("Original was ")
                print([x for x in original_classifications[query.name] if x not in skip_tags])
                print("Automatic classification was ")
                print(new_classification.tags)

                test_results += "\n\nFound an automatically classified genome that was different\n\n"
                test_results += query.name + "\n"
                test_results += "Original was \n"
                test_results +=  str([x for x in original_classifications[query.name] if x not in skip_tags])
                test_results += "\nAutomatic classification was \n"
                test_results += str([x for x in new_classification.tags if x not in skip_tags])




                # if original_classifications[query.name] != new_classification.tags[0].split("Auto_")[1]:

    print("\nWrong: " + str(diff_count))
    print("Correct " + str(count - diff_count))
    print("Total " + str(count))
    print()


    test_results = "Total " + str(count) + "\n" + test_results

    test_results = "Correct " + str(count - diff_count) + "\n" + test_results

    test_results = "\nWrong: " + str(diff_count) + "\n" + test_results

    return test_results


def delete_all_tags():

    queries = models.GenomeRecords.objects().all()

    queries.update(tags=[])

    for query in queries:

        # By default we keep the tag 'hidden' as this won't get GenomeTags out of sync

        for hit in query.hits:
            if 'hidden' in hit.tags:
                hit.tags = ['hidden']
            else:
                hit.tags = []

        query.save()

    models.GenomeTags.objects().all().delete()


def get_mlgo_dict(gene_orders):
    mlgo_dict = {}

    lines = gene_orders.split("\n")
    for i in range(0, len(lines)):
        line = lines[i]
        if line.startswith('>'):
            mlgo_dict[line.split(">")[1].strip()] = lines[i + 1].split("$")[0].strip().split(" ")

    print(mlgo_dict)

    return mlgo_dict


def colour_alignment_by_profiles(alignment, profiles):
    colour_dict = {'RBD_A': 'lightgreen',
                   'RBDA' : 'lightgreen',

                   'RBD_C': 'blue', 'RBD_B': 'orange', 'Neuraminidase': 'mediumPurple',
                   'TcA_RBD': 'grey',
                   'TcB_BD_seed': 'lawnGreen',
                   'PF18276_ncbi' : 'lawnGreen',
                   'VRP1_Full' : 'sandyBrown',
                   'Big_1_Full' : 'lightYellow',
                   'TcB_BLAST_500' : 'orange',
                   'TcC_BLAST_500': 'blue',
                   'Rhs_repeat' : 'green',

                   'Overlap' : 'pink'}

    split = [x for x in alignment.split(">") if len(x) > 0]

    print (split[0:5])


    size = len(split)

    output = {split[i].replace(" <unknown description", "").replace(".", "***"): split[i + 1].replace("\n", "") if size
                                                                                                                 > i
                                                                                                                   + 1
    else None
              for i in range(0, size, 2)}


    # print(profiles.region_dict.keys())
    #
    # print (output)

    # for s,d in profiles.region_dict.items():
    #     print (s)
    #     print (d)


    # for seqname in output.keys():
    #     if seqname in profiles.region_dict:
    #         print('found')
    #         for domain, pos in profiles.region_dict[seqname].items():
    #             print(domain)
    #             print(pos)
    #             print (output[seqname])
    #             print ('and now')
    #             output[seqname] = output[seqname][0:pos[0]] + '<span style = "background-color:' + colour_dict[
    #                 domain] + \
    #                               '">' + output[seqname][pos[0]: pos[1]] + '</span>' + output[seqname][pos[1]:]

    # for seqname, domains in output.items():
    #     if seqname in profiles.region_dict:
    #
    #         len_offset = 0
    #         furtherst_pos = -1
    #
    #         orig_seq = output[seqname]
    #
    #         for domain, pos in sorted(profiles.region_dict[seqname].items(), key=lambda k: k[1]):
    #             gap_offset = 0
    #
    #             if pos[0] < furtherst_pos:
    #                 print("WARNING: OVERLAP")
    #
    #             count = 0
    #             for aa in orig_seq:
    #                 if aa == "-":
    #                     gap_offset += 1
    #                 else:
    #                     count += 1
    #                     if count == pos[0] + 1:
    #                         first_gap_offset = gap_offset
    #
    #                     if count == pos[1]:
    #                         second_gap_offset = gap_offset
    #                         break
    #
    #             prev_len = len(output[seqname])
    #             output[seqname] = output[seqname][
    #                           0:pos[0] + len_offset + first_gap_offset] + '<span style = "background-color:' + \
    #                           colour_dict[
    #                               domain] + '">' + \
    #                               output[seqname][pos[0] + len_offset + first_gap_offset: pos[
    #                                                                                   1] + len_offset + second_gap_offset] + '</span>' + \
    #                               output[seqname][pos[1] +
    #                                    len_offset +
    #                                    second_gap_offset:]
    #             len_offset = len(output[seqname]) - prev_len
    #
    #             furtherst_pos = pos[1]

    for seqname, domains in profiles.region_dict.items():

        print (seqname)
        print (output.keys())

        if seqname not in output:
            print ('wowzers')
        # if seqname in profiles.region_dict:

        len_offset = 0
        furtherst_pos = -1


        # newlist = sorted(list_to_be_sorted, key=lambda k: k['name'])

        if seqname in output:
            orig_seq = output[seqname]

            sorted_list = sorted(profiles.region_dict[seqname].items(), key=lambda k: k[1])

            # print(sorted_list)

            list_w_overlaps = []

            overlap = False

            for idx, entry in enumerate(sorted_list):
                domain = entry[0]
                pos = entry[1]

                if idx + 1 < len(sorted_list):

                    next_entry = sorted_list[idx + 1]
                    next_domain = next_entry[0]
                    next_pos = next_entry[1]

                    if pos[1] > next_pos[0]:

                        overlap = True
                        print("WARNING: OVERLAP")
                        overlap_pos_1 = next_pos[0]
                        overlap_pos_2 = pos[1]

                        prev_entry = (domain, [pos[0], next_pos[0] ])
                        overlap_entry = ('Overlap', [overlap_pos_1, overlap_pos_2])
                        sorted_list[idx + 1] = (next_domain, [overlap_pos_2, next_pos[1]])

                        list_w_overlaps.append(prev_entry)
                        list_w_overlaps.append(overlap_entry)

                    else:
                        list_w_overlaps.append(entry)

                else:
                    list_w_overlaps.append(entry)


            if overlap:

                print('list with overlaps')
                print (seqname.replace("***", "."))
                print (list_w_overlaps)


            for entry in list_w_overlaps:

                if entry[1][1] < entry[1][0]:
                    print("WARNING: Serious error with overlap")
                    print(seqname.replace("***", "."))
                    print(list_w_overlaps)

            for domain, pos in list_w_overlaps:
                gap_offset = 0
                first_gap_offset = 0
                second_gap_offset = 0


                # if pos[0] < furtherst_pos:
                #     print("WARNING: OVERLAP")

                count = 0
                for aa in orig_seq:
                    if aa == "-":
                        gap_offset += 1
                    else:
                        count += 1
                        if count == pos[0] + 1:
                            first_gap_offset = gap_offset

                        if count == pos[1]:
                            second_gap_offset = gap_offset
                            break
                        else:
                            print ('count was ' + str(count))

                # print ('len offset')
                # print(len_offset)
                # print(gap_offset)

                # print(pos[0] + len_offset + first_gap_offset)
                #
                # print(pos[1] + len_offset + second_gap_offset)

                prev_len = len(output[seqname])

                output[seqname] = output[seqname][0:pos[0] + len_offset + first_gap_offset] + '<span style = "background-color:' + \
                              colour_dict[domain.split("_multiple_")[0]] + '">' \
                              + output[seqname][
                                pos[0] + len_offset + first_gap_offset: pos[1] + len_offset + second_gap_offset] + '</span>' \
                              + output[seqname][pos[1] + len_offset + second_gap_offset:]

                len_offset = len(output[seqname]) - len(orig_seq)

                # print ('prev len')
                # print (prev_len)
                #
                # print ('len aligns seq')
                # print (len(aligns[seq]))
                #
                # print ('len offset')
                #
                # print (len_offset)


                # print(aligns[seq])

                # furtherst_pos = pos[1]
                # print ('charlie')
                # print (output)


    return output
