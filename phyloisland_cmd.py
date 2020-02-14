import argparse
import os
import mongoengine
import cmd_code
import utilities
import models
from flask import Flask
from flask_mongoengine import MongoEngine
import configs.mongoconfig




parser = argparse.ArgumentParser()

parser.add_argument("-g", "--add_genomes", help="add genomes to the database", action="store_true")
parser.add_argument("-i", "--input_file", help="path to list of species")
parser.add_argument("-p", "--add_profiles", help="add profiles to the database", action="store_true")
parser.add_argument("-f", "--profile_folder", help="path to profile folder")
parser.add_argument("-o", "--overview", help="get overview of database", action="store_true")
parser.add_argument("-d", "--database_name", help="database name", required=True)


args = parser.parse_args()

if args.add_genomes and (args.input_file is None):
    parser.error("-g requires that you provide -i the path to list of species ")

if args.add_profiles and (args.profile_folder is None):
    parser.error("-p requires that you provide -f the path to folder with profiles")



# Print out the submitted input
print (f"Input file is {args.input_file}")
print (f"Profile folder is {args.profile_folder}")
print (f"Selected database is {args.database_name}")

# Configure mongo database
# app = Flask(__name__)
# app.config.from_pyfile('configs/mongoconfig.py')
#
# # Connect to mongo database
# db = MongoEngine(app)
# mongoengine.connect(configs.mongoconfig.MONGODB_DB)
#
#
# # Delete the existing profiles
# cmd_code.delete_profiles()
#
# # Retrieve the genomes
# cmd_code.get_genomes(args)
#
# # Update the profiles
# cmd_code.update_profiles(args)
#
# # Update the genomes
# cmd_code.update_genomes()

if args.overview:
    print ('overview')




