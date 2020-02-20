import argparse
import os
import mongoengine
import cmd_code
import utilities
import models
import genome_overview
from flask import Flask
from flask_mongoengine import MongoEngine
import getGenomes



parser = argparse.ArgumentParser()

parser.add_argument("-g", "--add_genomes", help="path to list of species")
# parser.add_argument("-i", "--input_file", help="path to list of species")
parser.add_argument("-p", "--add_profiles", help="path to profile folder")

parser.add_argument("-u", "--update_genomes", help="update genomes", action="store_true")
parser.add_argument("-f", "--fasta", help="save all regions to fasta files", action="store_true")
parser.add_argument("-r", "--region_order", help="write out order of regions in all genomes ", action="store_true")

parser.add_argument("-o", "--overview", help="get overview of database", action="store_true")
parser.add_argument("-d", "--delete_genome_tags", help="delete current genome classifications", action="store_true")
parser.add_argument("-c", "--classify", help="classify genomes based on their regions ", action="store_true")


# parser.add_argument("-d", "--database_name", help="database name", required=True)


args = parser.parse_args()

# if args.add_genomes and (args.input_file is None):
#     parser.error("-g requires that you provide -i the path to list of species ")
#
# if args.add_profiles and (args.profile_folder is None):
#     parser.error("-p requires that you provide -f the path to folder with profiles")
#


# Print out the submitted input
print (f"Input file is {args.add_genomes}")
print (f"Profile folder is {args.add_profiles}")
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
        getGenomes.download_fasta_regions(profile.name, "grobs")

if args.region_order:
    pass

