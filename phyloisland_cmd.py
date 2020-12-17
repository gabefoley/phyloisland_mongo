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

parser.add_argument("--query_db", help="query db", action="store_true")
parser.add_argument("--load_genomes", help="load genomes")
parser.add_argument("--get_accession_ids", help="get accession ids", action="store_true")
parser.add_argument("--check_feature_table", help="check feature table")

# parser.add_argument("-d", "--database_name", help="database name", required=True)


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

    df = pd.DataFrame()

    for entry in sorted(species_dict.items(), key=lambda x: len(x[1]), reverse=True):

        entry_dict = {}

        # print(entry[0], len(entry[1]))

        for genome in entry[1]:
            entry_dict['species'] = genome.species
            entry_dict['name'] = genome.name
            entry_dict['description'] = genome.description
            expanded_hits = [x for x in genome.hits if 'expanded' in
                             x.region and
                             'TcB_BD' not in x.region]

            # print(len(expanded_hits))
            # print(" ".join([x.region for x in expanded_hits]))

            entry_dict['hit_num'] = len(expanded_hits)
            entry_dict['all_hits'] = " ".join([x.region for x in expanded_hits])

            print (entry_dict)

            df = df.append(entry_dict, ignore_index=True)
    df.to_csv("./pi_output_gabepulse4.csv")



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

                print (filtered.head(2))

                for x in filtered.itertuples():
                    if x.feature == 'CDS':
                        if not (x.start > int(hit.end)) and not (x.end < int(hit.start)):
                            check_interval = pd.Interval(x.start, x.end, closed='both')

                            print ('checking')

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








