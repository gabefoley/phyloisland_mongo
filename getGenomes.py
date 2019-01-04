import argparse, subprocess, sys, fnmatch
from ftplib import FTP
import gzip
from Bio import SeqIO
from datetime import datetime, MINYEAR
import pandas as pd
from collections import defaultdict
import operator
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import utilities


def read_genome(outpath, species_name):
    print ("Reading in genome")


    # Collate all of the nucleotide records together to make the genome
    concatenated_genome = ""

    my_dict = SeqIO.to_dict(SeqIO.parse(outpath, "fasta"))
    for r in sorted(my_dict.values(), key=operator.attrgetter('id')):

        if "plasmid" not in r.description:
            concatenated_genome += str(r.seq)
            description = r.description
            genome_id = r.id

        # Temporary measure to reduce the name so it can fit in the database. Edge case but occurs with
        # 'bacterium endosymbiont of Mortierella elongata FMR23-6', for example

        if len(species_name) > 40:
            species_name = species_name[0:40]

    return SeqRecord(Seq(concatenated_genome), id=genome_id, name= species_name, description=description,
                       annotations={"organism": species_name, "source": ""})


def retrieve_genome(records, species_name, category, database):

    print ("Got to retrieve genome")

    genome_dict = {}

    if category == "assembly" or category == "genbank":
        folder = "latest_assembly_versions"
    else:
        folder = category.split(" ")[0]

    for val in records.values():
        location = val[1]
        break

    file_type = "--exclude='*cds_from*' --exclude='*rna_from*' --include='*genomic.fna.gz' --exclude='*'"

    print ("rsync -Lrtv --chmod=+rwx -p %s rsync://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/%s/%s/*/%s %s" % (
                    file_type, database, species_name.replace(" ", "_"), folder, location, "./tmp"))

    try:
        process = subprocess.Popen(
                "rsync -Lrtv --chmod=+rwx -p %s rsync://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/%s/%s/*/%s %s" % (
                    file_type, database, species_name.replace(" ", "_"), folder, location, "./tmp"), stderr=subprocess.PIPE,
                stdout=subprocess.PIPE, shell=True)

        out, err = process.communicate()
        errcode = process.returncode

        print ("errcode ", errcode)
        print ("output ", out.decode("utf-8"))

        if errcode != 0:
            return

        out_decoded = out.decode("utf-8")

        if "incremental" in out_decoded:
            file_list = out_decoded.split("receiving incremental file list")[1].split("sent")[0]

        else:
            file_list = out_decoded.split("receiving file list ... done")[1].split("sent")[0]

        if not file_list:
            return


        for filename in file_list.split():

            if len(filename) > 2:
                record = "_".join(filename.split("_")[0:2])
                assembly_level = records[record][0]

                filepath = "./tmp/" + filename

                file_from_zip = gzip.open(filepath, mode="rb")

                outpath = ".".join(filepath.split(".")[0:-1]) + ".fasta"

                with open(outpath, 'w') as query_file:
                    for line in file_from_zip:
                            query_file.write(line.decode('utf-8'))

                file_from_zip.close()

                if assembly_level == "Contig":
                    print ("%s with id %s was a contig" % (species_name, record))
                elif assembly_level == "Chromosome" or assembly_level == "Complete Genome":
                    print ("%s with id %s was a chromosome or full sequence" % (species_name, record))

                elif assembly_level == "Scaffold":
                    print ("%s with id %s was a scaffold" % (species_name, record))

                else:
                    print ("%s with id %s had an unhandled assembly level is, which was %s - " % (species_name, record, assembly_level))

                genome = read_genome(outpath, species_name)

                genome_dict[record] = genome
                print ("Genome dict is ")
                print (genome_dict)

                utilities.removeFile(outpath)
                utilities.removeFile(filepath)



    except subprocess.CalledProcessError as exc:
        return

    return genome_dict


def get_record_list(summary, category, single):
    """
    Get the list of records that match the category and return the accession and assembly level
    :param summary: The summary data frame containing all records
    :param category: The specific type of refseq category we're searching for
    :param single: Whether or not to only return a single record
    :return: A dictionary mapping accession id to location, assembly level
    """
    ref_dict = defaultdict(list)

    if category == "assembly" or category == "genbank":
        category = "na"


    refs = summary.loc[summary['refseq_category'] == category]

    if refs.empty:
        return

    if len(refs) == 1 or len(refs) > 1 and not single:

        for ref in refs.itertuples():
            ref_dict[ref._1] = (ref.assembly_level, "")

    elif len(refs) > 1 and single:

        # Sort the records by release date
        summary['seq_rel_date'] = pd.to_datetime(summary.seq_rel_date)
        summary.sort_values(by=['seq_rel_date'], inplace=True, ascending=False)

        for ref in refs.itertuples():
            ref_dict[ref._1] = (ref.assembly_level, ref.ftp_path.split("/")[-1] + "*")
            break

    return ref_dict


def add_genome2(species_name, categories, single):

    for x in ["reference genome", "representative genome", "assembly"]:
        if x in categories:
            database = "refseq"
            break
        else:
            database = "genbank"
            # Need to reinstate all the categories because GenBank can have these even if it doesn't have a RefSeq record
            categories = ["reference genome", "representative genome", "assembly", "genbank"]
            break


    try:
        # Add a v to the end of -Lrt to get verbose print outs to the console
        process = subprocess.Popen(
                "rsync -Lrt --chmod=+rwx -p rsync://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/%s/assembly_summary.txt %s" % (
                    database, species_name.replace(" ", "_"), "./tmp"), stderr=subprocess.PIPE,
                stdout=subprocess.PIPE, shell=True)

        out, err = process.communicate()
        errcode = process.returncode

        # print ("errcode ", errcode)
        # print ("output ", out.decode("utf-8"))

        if errcode != 0:
            return

        summary = pd.read_csv("./tmp/assembly_summary.txt", sep='\t', header=1)

        # print (summary)

        for category in categories:
            records = get_record_list(summary, category, single)

            # print ('category and records found are - ')
            # print (category, records)

            if records:
                genome_dict = retrieve_genome(records, species_name, category, database)
                return genome_dict
        return

    except subprocess.CalledProcessError as exc:
        return
