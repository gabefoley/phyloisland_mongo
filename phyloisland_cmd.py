import argparse
import os
import mongoengine
import cmd_code
import utilities
import models
import genome_overview
from flask import Flask
from flask_mongoengine import MongoEngine
from bson.objectid import ObjectId
import getGenomes
import gzip
import refseq_code
import numpy
from collections import defaultdict
import alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import pandas as pd
import glob

from Bio.Alphabet import generic_nucleotide


parser = argparse.ArgumentParser()

parser.add_argument("-g", "--add_genomes", help="path to list of species")
# parser.add_argument("-i", "--input_file", help="path to list of species")
parser.add_argument("-p", "--add_profiles", help="path to profile folder")

parser.add_argument("-u", "--update_genomes", help="update genomes", action="store_true")
parser.add_argument("-f", "--fasta", help="save all regions to fasta files", action="store_true")
parser.add_argument("-r", "--region_order", help="write out order of regions in all genomes ", action="store_true")

parser.add_argument("-o", "--overview", help="get overview of database", action="store_true")
parser.add_argument("-dg", "--delete_genomes", help="delete genomes", action="store_true")
parser.add_argument("-dt", "--delete_genome_tags", help="delete current genome classifications", action="store_true")
parser.add_argument("-c", "--classify", help="classify genomes based on their regions ", action="store_true")
parser.add_argument("-m", "--mlgo", help="write out regions to mlgo format", action="store_true")
parser.add_argument("--refseq", help="run full refseq check", action="store_true")

parser.add_argument("--refseq_again", help="run full refseq check", action="store_true")


parser.add_argument("--qc", help="run full refseq check", action="store_true")
parser.add_argument("--count", help="run full refseq check", action="store_true")
parser.add_argument("--sc", help="run full refseq check", action="store_true")
parser.add_argument("--create_csv", help="create csv", action="store_true")

# Check the assembly just to check I didn't add back in things we don't need
parser.add_argument("--check_assembly", help="create csv", action="store_true")


parser.add_argument("--query_db", help="query db", action="store_true")
parser.add_argument("--load_genomes", help="load genomes")
parser.add_argument("--get_accession_ids", help="get accession ids", action="store_true")
parser.add_argument("--check_feature_table", help="check feature table")

# This is a specific one
parser.add_argument("--update_tags_for_list", help="temporary to change the tags in A2 only Type1s and Type3s",
                    action="store_true")

# This is the general one (but not fully realised)
parser.add_argument("--update_tags_via_list", action="store_true")

parser.add_argument("--delete_all_profiles", action="store_true")

parser.add_argument



args = parser.parse_args()

# if args.add_genomes and (args.input_file is None):
#     parser.error("-g requires that you provide -i the path to list of species ")
#
# if args.add_profiles and (args.profile_folder is None):
#     parser.error("-p requires that you provide -f the path to folder with profiles")
#


# Print out the submitted input
# print (f"Input file is {args.add_genomes}")
# print (f"Profile folder is {args.add_profiles}")
# print (f"Selected database is {args.database_name}")

# # Configure mongo database
# app = Flask(__name__)
# # app.config.from_pyfile('configs/mongoconfig.py')
# app.config.MONGO_DB = 'raptor'
#
# # Connect to mongo database
# MONGODB_DB = 'raptor'
#
# db = MongoEngine(app)
# mongoengine.connect(db=args.database_name)
#


if args.add_genomes:
    # Retrieve the genomes
    print ("Adding genomes")
    cmd_code.get_genomes(args)

if args.add_profiles:
    # Delete the existing profiles if there is a new one in the Profile folder

    print("Adding profiles")
    cmd_code.delete_profiles(args)

    # Update the profiles
    cmd_code.update_profiles(args)

if args.update_genomes:

    print ("Updating genomes")
    # Update the genomes



    cmd_code.update_genomes()

if args.overview:
    print ('Overview of database')
    cmd_code.get_overview()

if args.delete_genomes:
    queries = models.GenomeRecords.objects.all().timeout(False).delete()
    print("Deleting all genomes")

if args.delete_genome_tags:
    queries = models.GenomeRecords.objects.all().timeout(False)
    print ("Deleting all existing genome classification tags")
    genome_overview.delete_genome_tags(queries)


if args.classify:
    queries = models.GenomeRecords.objects.all().timeout(False)
    print ("Classifying the genomes")
    genome_overview.classify_genomes(queries)


if args.fasta:
    profile_names = models.Profile.objects().all()

    for profile in profile_names:
        getGenomes.download_fasta_regions(profile.name + "_expanded", filename='cmd', split_strands=False,
                                          align=False)

if args.region_order:
    genomes = models.GenomeRecords.objects.all().timeout(False)
    getGenomes.write_genome_order(genomes, split_strands=False, path ='./fasta_folder/genome_order_from_cmd.txt')

if args.mlgo:
    genomes = models.GenomeRecords.objects.all().timeout(False)
    getGenomes.write_mlgo_order(genomes, path ='./fasta_folder/mlgo.txt')

if args.refseq:

    # Add 100 in

    filepath = "./files/refseq/working/1/"

    genomes = [x for x in os.listdir(filepath) if x.endswith(".gz")]

    chunk = numpy.array_split(numpy.array(genomes), 100)

    print ('chunk')

    print (chunk)

    for genomes in chunk:

        for genome in genomes:

            print (genome)

            genome_path = filepath + genome


            file_from_zip = gzip.open(genome_path, mode="rb")

            outpath = ".".join(genome_path.split(".")[0:-1]) + ".fasta"

            with open(outpath, 'w') as query_file:
                for line in file_from_zip:
                    query_file.write(line.decode('utf-8'))

            file_from_zip.close()


            genome_list = []

            genome_dict = refseq_code.read_genome(outpath)

            for k, v in genome_dict.items():
                print (k)
                print (v)

            print ("Adding genome\n")

            print (genome_dict)

            utilities.add_genome(genome_dict)

            print ("Remove unzipped FASTA file from disk\n")
            utilities.remove_file(outpath)


            print ("Remove zip files from disk\n")

            if os.path.exists(genome_path):

                utilities.remove_file(genome_path)


        print ("Search for profile in new genomes\n")

        queries = models.GenomeRecords.objects(hits__size=0).timeout(False)

        profiles = models.Profile.objects(name='TcB_BD_NCBI')

        for profile in profiles:
            cmd_code.get_feature_location_with_profile_cmd(queries, "hmm_outputs", profile.name, "", "", profile.name)

        del profiles
        del queries

        print ("Remove genomes without a hit\n")
        missing_hit = models.GenomeRecords.objects(hits__size=0)

        with open(filepath + "removed_genomes.txt", "a") as removed_genomes:
            for genome in missing_hit:
                removed_genomes.write(genome.description + "\n")



        models.GenomeRecords.objects(hits__size=0).delete()







        # utilities.remove_file(outpath)
    # utilities.remove_file(filepath)

    # Get all genomes without a hit

    # Search for TcB_BD in genomes

    # Delete anything without a hit

if args.qc:
    queries = models.GenomeRecords.objects(hits__size=0).timeout(False)

    for query in queries:
        print (query.name)


# Count the species
if args.count:
    species_dict = defaultdict(list)

    genomes = models.GenomeRecords.objects()

    delete_list = []

    region_dict = defaultdict(list)

    for genome in genomes:
        species_dict[genome.species].append(genome)

    sorted_dict = sorted(species_dict.items(), key=lambda x: len(x[1]), reverse=True)

    for entry in sorted(species_dict.items(), key=lambda x: len(x[1]), reverse=True):
        print(entry[0], len(entry[1]))


if args.sc:

    # for species_name in ['Proteus vulgaris', 'Photorhabdus laumondii']:

    species_dict = defaultdict(list)

    # genomes = models.GenomeRecords.objects(species='Proteus vulgaris')
    # genomes = models.GenomeRecords.objects(species='Photorhabdus laumondii')
    # genomes = models.GenomeRecords.objects(species=species_name)
    # genomes = models.GenomeRecords.objects()
    genomes = models.GenomeRecords.objects(species='Proteus cibarius')




    delete_list = []

    region_dict = defaultdict(list)

    for genome in genomes:
        species_dict[genome.species].append(genome)

    sorted_dict = sorted(species_dict.items(), key=lambda x: len(x[1]), reverse=True)

    with open ("./pi_output.txt", "a+") as output:

        for entry in sorted(species_dict.items(), key=lambda x: len(x[1]), reverse=True):
            print (entry[0], len(entry[1]))
            output.write("\n\n" + entry[0])
            output.write("\n" + "Number of strains present: " + str(len(entry[1])))

            for genome in entry[1]:
                expanded_hits = [x for x in genome.hits if 'expanded' in
                                 x.region and
                                 'TcB_BD' not in x.region]
                # print (len(genome.hits))
                print (len(expanded_hits))
                print (" ".join([x.region for x in expanded_hits]))

                output.write("\n" + str(len(expanded_hits)) + "\n")
                output.write(" ".join([x.region for x in expanded_hits]))

                # print (genome.tags)
                for hit in expanded_hits:

                    sequence = Seq(hit['sequence'], generic_nucleotide)

                    if hit['strand'] == 'backward':
                        sequence = sequence.reverse_complement()


                    region_dict[hit.region].append(SeqRecord(sequence.translate(),
                                                   id = genome.name))


            for region, seqs in region_dict.items():
                # print (region)

                if region != 'A2_expanded':
                    print()
                    print (genome.name, region)
                    output.write("\n\n")
                    output.write(genome.name + " " + region)
                    rand_string = utilities.randstring(4)

                    print ('Writing to - tmp/refseq/' + region + rand_string + ".aln")

                    utilities.createFasta(seqs, 'tmp/refseq/' + region + rand_string, align=True)


                    while not os.path.exists('tmp/refseq/' + region + rand_string + ".aln"):
                        time.sleep(1)

                    align = alignment.read_alignment('tmp/refseq/' + region + rand_string + ".aln", 'fasta')
                    # print (region)
                    # for seq in seqs:
                    percent_identity = alignment.get_percent_identity_of_alignment(align)
                    # print (percent_identity)
                    mean = alignment.get_mean(percent_identity)
                    std = alignment.get_standard_deviation(percent_identity)
                    print (mean)
                    print (std)
                    output.write("\n" + " Mean: " + str(mean))
                    output.write("\n" + " Std: " + str(std))



        # Delete entries with a single hit
        # for entry in sorted(species_dict.items(), key=lambda x: len(x[1]), reverse=True):
        #     if len(entry[1]) == 1:
        #         delete_list.append(entry[0])
        #
        # for species in delete_list:
        #     print (species)
        #
        # models.GenomeRecords.objects(species__in=delete_list).delete()



# Used to pull out all the genomes with at least one strain with a hit (i.e. add back in deleted genomes if other
# strains had a hit).


hit_list = ['Candidatus Carsonella', 'Burkholderia pseudomallei', 'Yersinia pestis', 'Pseudomonas sp.', 'Yersinia pseudotuberculosis', 'Burkholderia mallei', 'Salmonella enterica', 'Burkholderia thailandensis', 'Morganella morganii', 'Pseudomonas syringae', 'Pseudomonas fluorescens', 'Yersinia enterocolitica', 'Pseudomonas putida', 'Proteus mirabilis', 'Burkholderia sp.', 'Proteus cibarius', 'Pseudomonas brassicacearum', 'Pseudomonas chlororaphis', 'Pseudomonas orientalis', 'Shewanella baltica', 'Yersinia ruckeri', 'Burkholderia oklahomensis', 'Proteus vulgaris', 'Pseudomonas synxantha', 'Photorhabdus laumondii', 'Sorangium cellulosum', 'Vibrio parahaemolyticus', 'Wolbachia endosymbiont', 'Legionella israelensis', 'Streptomyces venezuelae', 'Bacillus sp.', 'Cedecea neteri', 'Cupriavidus necator', 'Rahnella aquatilis', 'Xenorhabdus nematophila', 'Erwinia pyrifoliae', 'Nonomuraea sp.', 'Pseudomonas entomophila', 'Desulfococcus multivorans', 'Pseudomonas simiae', 'Pseudomonas aeruginosa', 'Plantactinospora sp.', 'Xenorhabdus bovienii', 'Pandoraea oxalativorans', 'Fibrobacter succinogenes', 'Pseudomonas mosselii', 'Pseudomonas parafulva', 'Chryseobacterium bernardetii', 'Halomicronema hongdechloris', 'Kosakonia sacchari', 'Paraburkholderia rhizoxinica', 'Vibrio campbellii', 'Psychrobacter sp.', 'Candidatus Amoebophilus', 'Terasakiella sp.', 'Arsenophonus nasoniae', 'Nitrosococcus wardiae', 'Cellulomonas shaoxiangyii', 'Rhodococcus sp.', 'Methylocystis heyeri', 'Nocardia sp.', 'Paenibacillus cellulosilyticus', 'Nitrosomonas sp.', 'Pantoea agglomerans', 'Sinorhizobium fredii', 'Proteus sp.', 'Serratia marcescens', 'Streptomyces lydicus', 'Bradyrhizobium sp.', 'Shewanella sp.', 'Pedobacter sp.', 'Proteus hauseri', 'Aquimarina sp.', 'Pseudomonas libanensis', 'Xylophilus sp.', 'Tistrella mobilis', 'Serratia liquefaciens', 'Phytohabitans flavus', 'Streptomyces spectabilis', 'Bathymodiolus thermophilus', 'Streptomyces sp.', 'Pseudomonas mediterranea', 'Chitinophaga sp.', 'Lysobacter antibioticus', 'Paenibacillus sp.', 'Bradymonadales bacterium', 'Agarilytica rhodophyticola', 'Cohnella sp.', 'Actinomadura sp.', 'Burkholderia stagnalis', 'Pseudomonas graminis', 'Pseudomonas amygdali', 'Pseudomonas savastanoi', 'Tateyamaria omphalii', 'Microcoleus sp.', 'Serratia plymuthica', 'Microbulbifer sp.', 'Serratia sp.', 'Bordetella hinzii', 'Bacillus thuringiensis', 'Enterobacter sp.', 'Pseudomonas rhodesiae', 'Microlunatus sp.', 'Pseudomonas frederiksbergensis', 'Yersinia mollaretii', 'Proteus columbae', 'Nitrospira japonica', 'Pistricoccus aurantiacus', 'Lysobacter maris', 'Blautia producta', 'Streptomyces vietnamensis', 'Desulfotomaculum reducens', 'Pedobacter steynii', 'Stackebrandtia nassauensis', 'Pseudomonas cichorii', 'Pseudomonas kribbensis', 'Paenibacillus sabinae', 'Paenibacillus durus', 'Mycoavidus cysteinexigens', 'Pandoraea faecigallinarum', 'Moorea producens', 'Pandoraea pulmonicola', 'Kordia antarctica', 'Streptomyces alboniger', 'Pseudomonas rhizosphaerae', 'Haliscomenobacter hydrossis', 'Filimonas lacunae', 'Pseudomonas corrugata', 'Pseudomonas plecoglossicida', 'Mesorhizobium sp.', 'Nitrosospira multiformis', 'Nitrosospira lacus', 'Proteiniphilum saccharofermentans', 'Methylomusa anaerophila', 'Bradyrhizobium oligotrophicum', 'Pseudomonas viridiflava', 'Pseudomonas poae', 'Paraburkholderia hospita', 'Dehalobacterium formicoaceticum', 'Desulfofarcimen acetoxidans', 'Yersinia similis', 'Mucilaginibacter gotjawali', 'Bradyrhizobium paxllaeri', 'Shewanella denitrificans', 'Pseudomonas koreensis', 'Marinomonas posidonica', 'Nostoc flagelliforme', 'Pandoraea vervacti', 'Nitrosomonas ureae', 'Mycobacterium kansasii', 'Rhodococcus jostii', 'Photorhabdus thracensis', 'Nitrosospira briensis', 'Vulgatibacter incomptus', 'Desulfallas gibsoniae', 'Thalassococcus sp.', 'Xenorhabdus hominickii', 'Pseudomonas fragi', 'Pseudomonas psychrophila', 'Fabibacter pacificus', 'Pseudomonas lundensis', 'Cellulophaga baltica', 'Photorhabdus asymbiotica', 'Mixta intestinalis', 'Sphingomonas sp.', 'Pseudomonas monteilii', 'Halomonas chromatireducens', 'Erwinia sp.', 'Microbacterium sp.', 'Candidatus Nitrospira', 'Burkholderia ubonensis', 'Granulibacter bethesdensis', 'Chlorobium phaeobacteroides', 'Candidatus Fukatsuia', 'Catenulispora acidiphila', 'Chondromyces crocatus', 'Minicystis rosea']

if args.refseq_again:

    # Add 100 in

    filepath = "./files/refseq/working/1/"


    genomes = [x for x in os.listdir(filepath) if x.endswith(".gz")]


    chunk = numpy.array_split(numpy.array(genomes), 1)

    # print ('chunk')
    #
    # print (chunk)

    for genomes in chunk:

        for genome in genomes:

            genome_path = filepath + genome


            file_from_zip = gzip.open(genome_path, mode="rb")

            outpath = ".".join(genome_path.split(".")[0:-1]) + ".fasta"

            with open(outpath, 'w') as query_file:
                for line in file_from_zip:
                    query_file.write(line.decode('utf-8'))

            file_from_zip.close()


            genome_list = []

            genome_dict = refseq_code.read_genome(outpath, hit_list)

            if len (genome_dict) > 0:
                print ('Found something!')

            # print (genome_dict)


        #     for k, v in genome_dict.items():
        #         print (k)
        #         print (v)
        #
        #     print ("Adding genome\n")
        #
        #     print (genome_dict)
        #
            utilities.add_genome(genome_dict)

            print ("Remove unzipped FASTA file from disk\n")
            utilities.remove_file(outpath)


            print ("Remove zip files from disk\n")

            if os.path.exists(genome_path):

                utilities.remove_file(genome_path)
        #
        #
        # print ("Search for profile in new genomes\n")
        #
        # queries = models.GenomeRecords.objects(hits__size=0).timeout(False)
        #
        # profiles = models.Profile.objects(name='TcB_BD_NCBI')
        #
        # for profile in profiles:
        #     cmd_code.get_feature_location_with_profile_cmd(queries, "hmm_outputs", profile.name, "", "", profile.name)
        #
        # del profiles
        # del queries
        #
        # print ("Remove genomes without a hit\n")
        # missing_hit = models.GenomeRecords.objects(hits__size=0)
        #
        # with open(filepath + "removed_genomes.txt", "a") as removed_genomes:
        #     for genome in missing_hit:
        #         removed_genomes.write(genome.description + "\n")
        #
        #
        #
        # models.GenomeRecords.objects(hits__size=0).delete()
        #

if args.create_csv:

    # for species_name in ['Proteus vulgaris', 'Photorhabdus laumondii']:

    species_dict = defaultdict(list)

    tag_dict = getGenomes.get_tags()

    for x, k in tag_dict.items():
        print (x, k)

    # genomes = models.GenomeRecords.objects(species='Proteus vulgaris')
    # genomes = models.GenomeRecords.objects(species='Photorhabdus laumondii')
    # genomes = models.GenomeRecords.objects(species=species_name)
    # genomes = models.GenomeRecords.objects()
    genomes = models.GenomeRecords.objects()

    delete_list = []

    region_dict = defaultdict(list)

    for genome in genomes:
        species_dict[genome.species].append(genome)

    sorted_dict = sorted(species_dict.items(), key=lambda x: len(x[1]), reverse=True)

    df = pd.DataFrame(columns=['species', 'refseq_id', 'name', 'description', 'assembly_name', 'assembly_level' ,
                               'plasmid', 'hit_num','Full'])



    for entry in sorted(species_dict.items(), key=lambda x: len(x[1]), reverse=True):

        entry_dict = {}

        # print(entry[0], len(entry[1]))

        for genome in entry[1]:
            entry_dict['plasmid'] = genome.plasmid

            entry_dict['species'] = genome.species
            entry_dict['name'] = genome.name
            entry_dict['description'] = genome.description
            entry_dict['assembly_level'] = genome.assembly_level
            entry_dict['assembly_name'] = genome.assembly_name
            entry_dict['refseq_id'] = genome.refseq_accession_id


            expanded_hits = [x for x in genome.hits if 'expanded' in
                             x.region and
                             'TcB_BD' not in x.region]
            TcC_hits = [x for x in genome.hits if 'TcC_expanded' in x.region]
            TcB_hits = [x for x in genome.hits if 'TcB_expanded' in x.region]
            TcdA1_hits = [x for x in genome.hits if 'TcdA1_expanded' in x.region]
            A1_hits = [x for x in genome.hits if 'A1_expanded' in x.region]
            A2_hits = [x for x in genome.hits if 'A2_expanded' in x.region]
            A2_hits = [x for x in genome.hits if 'A2_expanded' in x.region]
            Chitinase_hits = [x for x in genome.hits if 'Chitinase_expanded' in x.region]

            full = True if len(TcC_hits) > 0 and len(TcB_hits) > 0 and (len(TcdA1_hits) > 0 or len(A2_hits) > 0 or
                                                                        len(A1_hits) > 0) else \
                False

            # print(len(expanded_hits))
            # print(" ".join([x.region for x in expanded_hits]))

            entry_dict['hit_num'] = len(expanded_hits)
            entry_dict['all_hits'] = " ".join([x.region for x in expanded_hits])

            entry_dict['TcC_count'] = len(TcC_hits)
            entry_dict['TcB_count'] = len(TcB_hits)
            entry_dict['TcdA1_count'] = len(TcdA1_hits)
            entry_dict['A1_count'] = len(A1_hits)
            entry_dict['A2_count'] = len(A2_hits)
            entry_dict['Chitinase_count'] = len(Chitinase_hits)
            entry_dict['Full'] = full
            entry_dict['Classification'] = tag_dict[genome.name]




            print (entry_dict)

            df = df.append(entry_dict, ignore_index=True)
    df.to_csv("./database_dump.csv")

if args.check_assembly:
    genomes = models.GenomeRecords.objects(genbank_accession_id='GCA_900323885.1')
    for g in genomes:
        print (g)
        print (g.name)
        print (g.excluded)
        print (g.genbank_accession_id)
        print (g.refseq_accession_id)
        for hit in g.hits:
            print (hit)

    genomes = models.GenomeRecords.objects(genbank_accession_id='GCA_011110255.1')
    for g in genomes:
        print (g)
        print (g.name)
        print (g.excluded)
        print (g.genbank_accession_id)
        print (g.refseq_accession_id)
        for hit in g.hits:
            print (hit)

if args.query_db:
    genomes = models.GenomeRecords.objects(species='Halomicronema hongdechloris')

    # genomes = models.GenomeRecords.objects(species='Photorhabdus heterorhabditis')

    for genome in genomes:
        print (genome)
        for hit in genome.hits:
            print (hit)


if args.load_genomes:

    if not args.load_genomes[-1] == "/":
        filepath = args.load_genomes + "/"
    else:
        filepath = args.load_genomes




    # filepath = "./files/test_genomes/"
    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200817_Checking_full_9000_genomes/genbank" \
    #            "/bacteria/"

    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200831_Getting_just_matches/genbank/bacteria/"
    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200831_Getting_just_matches/refseq/bacteria/"

    # Just the YP
    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20201028_Testing_yersinia_pseudotuberculosis" \
    #            "/refseq/bacteria/"
    #
    # # Just the YP
    # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20201119_Testing_Feature_Tables" \
    #            "/refseq/bacteria/"

    # Octopus - eight genomes four each from genbank / refseq with different conditions for testing feature tables
    # filepath = '/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20201215_Checking_all_types_reading_features' \
    #            '/genbank/bacteria/'




    genome_name = [x for x in os.listdir(filepath) if x != '.DS_Store']


    chunk = numpy.array_split(numpy.array(genome_name), 5)

    print (chunk)


    for genomes in chunk:

        print (genomes)



        for genome in genomes:

            print ("NEW GENOME")

            print (filepath)
            print (genome)



            genome_path = glob.glob(filepath + genome + "/*_genomic.fna.gz")[0]

            report_path = glob.glob(filepath + genome + "/*_assembly_report.txt")[0]

            print (genome_path)

            print (report_path)


            # genome_path = filepath + genome + "/" + genome + '_genomic.fna.gz'

            file_from_zip = gzip.open(genome_path, mode="rb")

            outpath = ".".join(genome_path.split(".")[0:-1]) + ".fasta"

            with open(outpath, 'w') as query_file:
                for line in file_from_zip:
                    query_file.write(line.decode('utf-8'))

            file_from_zip.close()


            genome_list = []

            genome_dict = refseq_code.get_genome_dict(outpath, report_path)


            # genome_dict = refseq_code.read_genome(outpath, genome)

            if len (genome_dict) > 0:
                print ('Found something!')

            print (genome_dict)


        #     for k, v in genome_dict.items():
        #         print (k)
        #         print (v)
        #
        #     print ("Adding genome\n")
        #
        #     print (genome_dict)
        #
            utilities.add_genome(genome_dict)

            # print ("Remove unzipped FASTA file from disk\n")
            # utilities.remove_file(outpath)
            #
            #
            # print ("Remove zip files from disk\n")
            #
            # if os.path.exists(genome_path):
            #
            #     utilities.remove_file(genome_path)
        #
        #
        # print ("Search for profile in new genomes\n")
        #
        # queries = models.GenomeRecords.objects(hits__size=0).timeout(False)
        #
        # profiles = models.Profile.objects(name='TcB_BD_NCBI')
        #
        # for profile in profiles:
        #     cmd_code.get_feature_location_with_profile_cmd(queries, "hmm_outputs", profile.name, "", "", profile.name)
        #
        # del profiles
        # del queries
        #
        # print ("Remove genomes without a hit\n")
        # missing_hit = models.GenomeRecords.objects(hits__size=0)
        #
        # with open(filepath + "removed_genomes.txt", "a") as removed_genomes:
        #     for genome in missing_hit:
        #         removed_genomes.write(genome.description + "\n")
        #
        #
        #
        # models.GenomeRecords.objects(hits__size=0).delete()

if args.get_accession_ids:
    genomes = models.GenomeRecords.objects()

    expected_ids = set(utilities.load_list(
        '/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200831_Getting_just_matches/candidates_refseq_genbank_mapped_892.txt',
        '/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200831_Getting_just_matches'
        '/candidates_genbank_only_4.txt'))

    accession_ids = set()

    for g in genomes:
        if not g.refseq_accession_id:
            accession_ids.add(g.genbank_accession_id)
        else:
            accession_ids.add(g.refseq_accession_id)

    print (f"Looking for {len(expected_ids)} accessions ")

    missing = expected_ids.difference(accession_ids)

    if len(missing) > 0:
        print (f"We are missing {len(missing)}")

        print ("And they are")
        print (missing)

    else:
        print ("And they are all there")




if args.check_feature_table:
    genomes = models.GenomeRecords.objects()

    print (args.check_feature_table)

    if not args.check_feature_table[-1] == "/":
        path = args.check_feature_table + "/"
    else:
        path = args.check_feature_table

    genbank_path = path + "genbank/bacteria/"

    refseq_path = path + "refseq/bacteria/"

    print ('Base path is ' + path)
    print ('Genbank path is ' + genbank_path)
    print ('Refseq path is ' + refseq_path)
    print ('Starting')


    for g in genomes:

        print ('Checking this genome')
        print (f'{g.name}  {g.species}')

        # Get genbank or refseq accession ID
        if not g.refseq_accession_id:
            print ('It has a Genbank accession ID')
            accession_id = g.genbank_accession_id
            filepath = genbank_path
        else:
            print ('It has a RefSeq accession id')
            accession_id = g.refseq_accession_id
            filepath = refseq_path
            if not os.path.exists(filepath + accession_id):
                print('No entry for a RefSeq id in the feature table folder')
                accession_id = g.genbank_accession_id
                filepath = genbank_path
                if not os.path.exists(filepath + accession_id):
                    print ('No entry for either Refseq or Genbank in the feature table folder ')



            # filepath = "/Users/gabefoley/Dropbox/PhD/Projects/Phylo_Island/2020/20200831_Getting_just_matches/refseq" \
            #            "/bacteria/"

        # Load in the feature table
        print ("Feature table path is " + filepath + accession_id)

        feature_path = glob.glob(filepath + accession_id + "/*_feature_table.txt.gz")[0]
        df = pd.read_csv(feature_path, sep='\t', compression='gzip')
        df.rename(columns={'# feature': 'feature'}, inplace=True)

        for hit in g.hits:
            if 'expanded' in hit.region and not 'region1' in hit.region:
                # print (hit.region)
                # print (hit.start)
                # print (hit.end)

                genome_interval = pd.Interval(int(hit.start), int(hit.end), closed='both')

                filtered = df.loc[df['genomic_accession'] == g.name]

                # print (filtered.head(2))

                for x in filtered.itertuples():
                    if x.feature == 'CDS':
                        if not (x.start > int(hit.end)) and not (x.end < int(hit.start)):
                            check_interval = pd.Interval(x.start, x.end, closed='both')

                            if check_interval.overlaps(genome_interval):
                                print(x.name, x.product_accession, x.start, x.end)
                                object_id = ObjectId()

                            #
                                hit = models.Hits(object_id,
                                                  g.name,
                                                  "EXISTING:" + " " + x.name + " " + str(x.product_accession),
                                                  '-1',
                                                  str(x.start),
                                                  str(x.end),
                                                  True,
                                                  hit.strand, '')

                                g.hits.append(hit)
                                g.save()
                            #
                            # print('We are adding a refseq sequence in')
                            #
                            # print(hit.start)
                            # print(hit.end)
                            # curr.hits.append(hit)
                            # curr.save()




if args.update_tags_for_list:
    select_genomes = ['NZ_CP048835.1', 'JOGP01', 'NZ_CP031066.1', 'NZ_CP031063.1', 'FXWM01', 'NZ_CP038254.1', 'NZ_CP041668.1', 'FMVH01', 'NZ_CP036313.1', 'MKMC01', 'VWSH01', 'AGJN02', 'MTAX01', 'RCWL01', 'VZZK01', 'AJLJ01', 'NVPT01', 'ARBP01', 'NZ_CP033932.1', 'FPLG01', 'NC_010830.1', 'BBLT01', 'JXRA01', 'VEGT01', 'AUAX01', 'MSSW01', 'LNCD01', 'ASJB01', 'PYBJ01', 'VLPL01', 'CZQA01', 'NZ_CP014226.1', 'QOIO01', 'AWXZ01', 'FODH01', 'QNVV01', 'PVZG01', 'QLLL01', 'LIPN01']
    genomes = models.GenomeTags.objects(tag_id__in=select_genomes)
    for x in genomes:
        print (x)
        print (x.tags)
        for idx, tag in enumerate(x.tags):
            if tag == 'Type3':
                x.tags.pop(idx)
                x.tags.append('Type3_A2')
            if tag == 'Type1':
                x.tags.pop(idx)
                x.tags.append('Type1_A2')
        print (x.tags)
        x.save()

if args.update_tags_via_list:
    # OD_Skip
    # select_genomes = ['FUXU01', 'SMOD01', 'SODV01', 'VJWF01', 'NZ_CP013949.1', 'NZ_CP031070.1', 'NZ_CP036313.1', 'VZZK01', 'NZ_CP011129.1', 'OGTP01', 'NZ_CP048835.1', 'NC_020453.1', 'VZZZ01', 'AUAX01', 'FCNY02', 'ASJB01', 'PYBJ01', 'VLPL01', 'NZ_CP014226.1', 'FAOZ01', 'LIPN01', 'NZ_CP024923.1', 'NZ_CP031066.1', 'NZ_CP031063.1', 'VEGT01', 'AGJN02', 'NZ_CP036313.1', 'VZZK01', 'OGTP01', 'NZ_CP048835.1', 'MKMC01', 'BJMN01', 'VZZZ01', 'AUAX01', 'JOGP01', 'ASJB01', 'PYBJ01', 'VLPL01', 'CZQA01', 'NZ_CP014226.1', 'QOIO01', 'PVZG01', 'LIPN01', 'CABPSQ01']
    #
    # #OD_Skip ver 4
    # select_genomes = ['FXWM01', 'MTBD01','VUOC01', 'NC_015559.1' ]
    #
    # #OD_Skip ver 4_2 (A1 / TcdA1 alignment)
    # select_genomes = ['ALVN01', 'SMKK01', 'PJBP01' ]
    #
    # # OD_Skip ver 4_3 (From Cropped edges alignment)
    # select_genoems = ['QJUG01', 'VDCQ01', 'FNVU01']

    # All the 219 genomes used in the final ABC trees
    select_genomes = ['AP018271.1','LIVZ01','MDEO01','NZ_AP018449.1','ONZJ01','WWHE01','FNJL01','ATXB01','NCXP01','QJUG01','LWBP01','FQUQ01','NZ_CP038254.1','MTAX01','LVYD01','RCFQ01','NZ_CP041668.1','FXYF01','RCWL01','MCHY01','NC_013216.1','NZ_CP029064.1','AJLJ01','QAAE01','SJSL01','NZ_CP029196.1','NVPT01','QTPO01','SNXZ01','NZ_CP013341.1','FNKR01','NZ_CP044218.1','NZ_LT629732.1','VZRB01','LKPJ01','FOUX01','VIWA01','QREK01','VSFF01','ARBP01','BJNF01','VOQD01','NZ_CP013426.1','AJUL01','QKLR01','FSRS01','NZ_CP027738.1','NZ_CP033932.1','QGTQ01','MVGR01','JOGE01','NETK01','NZ_CP004078.1','UPHP01','OIFR01','BIFQ01','NC_010830.1','RQJP01','NZ_CP009288.1','NZ_CP039288.1','LAIJ01','WWJN01','QGTG01','AYLO01','MRCJ01','NZ_CP022121.1','AZAN01','VDCQ01','BBLT01','FMCR01','QGSY01','FOUB01','MUXN01','NZ_AP017313.1','FMVH01','BILZ01','MWPQ01','NZ_LS999839.1','SMKX01','JXRA01','RQPI01','NZ_CP023695.1','SMJU01','FZPH01','NC_009253.1','NZ_CP013461.1','NZ_CM001559.1','VHKL01','RCSU01','NZ_CP042382.1','NZ_CP007699.2','QAOQ01','AQRJ01','AWZT01','JAAQYP01','AUEZ01','FMDM01','FNAD01','NKQZ01','VSFG01','SHKK01','WBKQ01','OAOQ01','MSSW01','JOGR01','QVNU01','BFCB01','NZ_CP010407.1','AUKP01','QGSZ01','RXLQ01','FNGF01','QVNQ01','SHKZ01','SLVP01','NZ_CP029235.1','LNCD01','QEOK01','NZ_CP029618.1','SNYE01','QAJI01','NZ_CP041692.1','QPIJ01','FMYF01','BDBY01','ASTJ01','OUNR01','SNZP01','NZ_CP034783.1','JNYY01','NZ_CP029710.1','QXQA01','AXBA01','JACCAY01','AKXV01','NZ_CP054613.1','WIAO01','NZ_CP038255.1','VCNA01','VXLC01','MUNY01','ALBT01','RPND01','LMXH01','NC_015510.1','BJMH01','AWXZ01','FODH01','FOIG01','ARBJ01','NZ_CP029843.1','FOLC01','QNVV01','VIVV01','CABLBP01','NZ_CP026106.1','MQWC01','MKCS01','JPOH01','VWSH01','FNON01','NC_013131.1','FOOS01','SMSL01','QLLL01','NZ_CP029197.1','AEDD01','ABXF01','QKWJ01','FNVU01','RAVX01','SLZA01','PCQL01','PHHE01','NZ_CP031450.1','NITZ01','NEHI01','NZ_CP010897.2','CVRZ01','PENV01','AOCZ01','NZ_CP039371.1','CABPSP01','PYLU01','FPLG01','NZ_CP020080.1','NZ_CP050291.1','NZ_CP025263.1','NZ_CP013047.2','LKBR01','PENZ01','PKNM01','PIQI01','NZ_CP029983.1','JFHN01','PKND01','VUAZ01','NZ_CP022960.1','LVTS01','NZ_CP045767.1','NZ_CP011807.3','ABCQ01','LYRP01','JYLN01','NZ_CP024646.1','NZ_CP031146.1','LACH01','PPRY01','NZ_CP015225.1','PQMB01','CQAW01','PENT01','NZ_CP053682.1','NC_020209.1','FNTT01','JAABNH01','NJAK01','NZ_CP010029.1','CPYD01']

    # Missing NMD
    missing_nmd = ['QJUG01', 'MTAX01', 'FNVU01', 'VDCQ01', 'AJLJ01', 'AWXZ01', 'NZ_CP038254.1', 'NZ_CP041668.1',
                 'FMVH01', 'QTPO01', 'LKPJ01', 'NZ_CP013426.1', 'FOLC01', 'VIVV01', 'BIFQ01', 'SNZP01', 'MKCS01', 'VSFF01', 'VSFG01', 'QVNQ01', 'VOQD01', 'NZ_CP007699.2', 'JOGE01', 'MUXN01', 'JNYY01', 'WWJN01', 'AQRJ01', 'NZ_CP023695.1', 'BJNF01', 'MWPQ01', 'NZ_CP039288.1', 'RXLQ01', 'SLZA01', 'NC_013131.1', 'FSRS01']

    select_genomes = [ x for x in select_genomes if x not in missing_nmd]


    # Missing NMD to add missing NMD to the tags
    select_genomes = ['QJUG01', 'MTAX01', 'FNVU01', 'VDCQ01', 'AJLJ01', 'AWXZ01', 'NZ_CP038254.1', 'NZ_CP041668.1',
                 'FMVH01', 'QTPO01', 'LKPJ01', 'NZ_CP013426.1', 'FOLC01', 'VIVV01', 'BIFQ01', 'SNZP01', 'MKCS01', 'VSFF01', 'VSFG01', 'QVNQ01', 'VOQD01', 'NZ_CP007699.2', 'JOGE01', 'MUXN01', 'JNYY01', 'WWJN01', 'AQRJ01', 'NZ_CP023695.1', 'BJNF01', 'MWPQ01', 'NZ_CP039288.1', 'RXLQ01', 'SLZA01', 'NC_013131.1', 'FSRS01']


    genomes = models.GenomeRecords.objects(name__in=select_genomes)

    for x in genomes:
        x.tags.append('Missing_NMD')
        x.save()

if args.delete_all_profiles:
    # Delete the profiles from the database (so we can easily update them with new ones if need be)
    profiles_to_delete = models.Profile.objects()
    profiles_to_delete.delete()
