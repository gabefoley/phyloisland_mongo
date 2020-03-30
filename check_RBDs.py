import os
import sequence
import subprocess
import glob
from Bio import SearchIO


# Change this base directory to your directory

# Original folder with profiles built from single P. luminescens or the Pfam seed alignments
# base_dir = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200309_RBD/20200324_RBD_profile_search"

# New folder with profiles built from all the Leidreiter profiles or the Pfam seed alignments
base_dir = '/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200309_RBD/20200326_RBD_Profile_HMMs'

input_folder = base_dir + "/input/"

output_folder = base_dir + "/output/"

profile_folder = base_dir + "/profiles/"

domScore = 1

region_dict = {}
domain_dict = {}

files = [x for x in os.listdir(input_folder) if x != ".DS_Store" and x != 'tmp']
profiles = [x for x in os.listdir(profile_folder) if x != ".DS_Store"]

for file in files:
    seqs = sequence.readFastaFile(input_folder + file)

    # Write out each region in the file to a separate file

    if not os.path.isdir(input_folder + "tmp"):
        os.mkdir(input_folder + "tmp")

    for seq in seqs:
        sequence.writeFastaFile(input_folder + "tmp/" + seq.name, [seq])

        if seq.name not in region_dict:
            region_dict[seq.name] = []


    # Check each region with each profile

        for profile in profiles:
            profile_name = profile.split(".hmm")[0]
            stdoutdata = subprocess.getoutput("hmmsearch -o %s --domT %s %s %s" % (output_folder + seq.name + "_" +
                                                                                   "profile=" +
            profile_name + ".output", domScore, profile_folder + profile, input_folder + "tmp/" + seq.name ))

for infile in glob.glob(output_folder + '/*.output'):

    seq_name = infile.split(output_folder)[1].split("_profile=")[0]

    # if seq_name not in region_dict:
    #     region_dict[seq_name] = []

    qresult = SearchIO.read(infile, 'hmmer3-text')
    if len(qresult.hits) > 0:
        domain = infile.split("profile=")[1].split(".")[0]
        region_dict[seq_name].append(domain)


# Remove the region files and profile search results

files = glob.glob(output_folder + '/*.output')
for f in files:
    os.remove(f)

files = glob.glob(input_folder + 'tmp/*')
for f in files:
    os.remove(f)

print (region_dict)




