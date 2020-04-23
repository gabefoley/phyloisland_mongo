import argparse, subprocess, sys, fnmatch
from ftplib import FTP
import gzip
from Bio import SeqIO
from Bio.Alphabet import generic_nucleotide
from datetime import datetime, MINYEAR
import pandas as pd
from collections import defaultdict
import operator
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import utilities
import models
import json


def read_genome(outpath, species_name):

    # Collate all of the nucleotide records together to make the genome
    concatenated_genome = ""

    my_dict = SeqIO.to_dict(SeqIO.parse(outpath, "fasta"))
    for r in sorted(my_dict.values(), key=operator.attrgetter('id')):

        concatenated_genome += str(r.seq)
        description = r.description
        genome_id = r.id

        # if r.id == "NZ_NIBS01000003.1":
        #     print (concatenated_genome)

        # if "plasmid" not in r.description:
        #     print ('Was not a plasmid')
        #     print (r.description)
        #     # concatenated_genome += str(r.seq)
        #     # description = r.description
        #     # genome_id = r.id
        #
        # else:
        #     # print ('Was a plasmid')
        #     # print (r.description)


        # Temporary measure to reduce the name so it can fit in the database. Edge case but occurs with
        # 'bacterium endosymbiont of Mortierella elongata FMR23-6', for example

        if len(species_name) > 40:
            species_name = species_name[0:40]

    return SeqRecord(Seq(concatenated_genome), id=genome_id, name= species_name, description=description,
                       annotations={"organism": species_name, "source": ""})


def retrieve_genome(records, species_name, category, database):

    print (f"Retrieving genome for {species_name}")

    genome_dict = {}

    if category == "assembly" or category == "genbank":
        folder = "latest_assembly_versions"
    else:
        folder = category.split(" ")[0]

    for val in records.values():
        location = val[1]
        break

    file_type = "--exclude='*cds_from*' --exclude='*rna_from*' --include='*genomic.fna.gz' --exclude='*'"

    print ("Genome retrieval called with following command - ")

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

                # utilities.remove_file(outpath)
                utilities.remove_file(filepath)



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


def add_genome(species_name, categories, single):

    print ('got here')

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

        if errcode != 0:
            print (errcode)
            print (err)
            return

        summary = pd.read_csv("./tmp/assembly_summary.txt", sep='\t', header=1)

        for category in categories:
            records = get_record_list(summary, category, single)

            if records:
                genome_dict = retrieve_genome(records, species_name, category, database)
                return genome_dict

        # Clean up the assembly summary
        utilities.remove_file('./tmp/assembly_summary.txt')

        return

    except subprocess.CalledProcessError as exc:
        return

def download_associated_regions():

    ar = models.AssociatedHits.objects()
    ar_dict = {}

    for x in ar:
        ar_dict[x['region1']] = x['region2']


    print (ar_dict)


    print (ar)

    out_path = "./fasta_folder/associated_regions.txt"

    with open (out_path, "w+") as outfile:
        outfile.write(json.dumps(ar_dict))

    return out_path

def download_tags():

    tags = models.GenomeTags.objects()

    tag_dict = {}

    for x in tags:
        tag_dict[x['tag_id']] = x['tags']

    out_path = "./fasta_folder/tags.txt"

    with open (out_path, "w+") as outfile:
        outfile.write(json.dumps(tag_dict))

    return out_path


def download_fasta_regions(region, filename="", include_genome=[], exclude_genome=[],
                           include_hits=[], \
                                                                                                     exclude_hits=[],
                           translate=True, align=True,  split_strands=False):
    fasta_dict = {}
    forward_dict = {}
    backward_dict = {}

    # Don't add an underscore if we're not also adding extra text to the filename
    filename = region + "_" + filename if filename else region
    count = 1



    aggregate = models.GenomeRecords._get_collection().aggregate([
        {"$match": {
            "hits.region": region
        }},
        {"$redact": {
            "$cond": {
                "if": {"$eq": [{"$ifNull": ["$region", region]}, region]},
                "then": "$$DESCEND",
                "else": "$$PRUNE"
            }
        }}
    ])


    for genome in aggregate:




        # seq_count = defaultdict(list)

        # Check that this genome should be included (if include_genome is non-empty) and that it
        # shouldn't be excluded
        if ((include_genome == [''] or bool(set(genome['tags']).intersection(set(include_genome))))) \
                and not bool(set(genome['tags']).intersection(set(exclude_genome))):

            print('got to here')

            print (genome['name'])


            # Check that this hit should be included (if include_hit is non-empty) and that it
            # shouldn't be excluded

            for hit in genome['hits']:
                # print ("*******" + hit)
                print (include_hits)
                if not include_hits:
                    print ('1')
                if bool(set(hit['tags']).intersection(set(include_hits))):
                    print ('2')
                if (include_hits == [""] or bool(set(hit['tags']).intersection(set(include_hits)))) and not bool(set(
                        hit['tags']).intersection(set(exclude_hits))):

                    print ('found hit')
                    #
                    # print()
                    print(hit['name'])
                    # print(hit['region'])
                    #
                    # print(hit['start'])
                    # print(hit['end'])
                    # print(hit['strand'])

                    sequence = Seq(hit['sequence'], generic_nucleotide)

                    if hit['strand'] == 'backward':
                        sequence = sequence.reverse_complement()

                    # Do we want to translate the sequences into protein?
                    if translate:
                        sequence = sequence.translate()

                    id_name = hit['name'] + "_information_" + genome['species'].replace(" ", "_") + '_region_' + hit[
                        'region'] + "_" + hit['start'] + "_" + hit['end'] + "_" + hit['strand']

                    # seq_count[hit['name']].append(id_name)

                    fasta_record = SeqRecord(sequence, id_name)

                    # We want to separate forward and backward strands
                    if split_strands:
                        print ('splitting strands')
                        if hit['strand'] == 'forward':
                            forward_dict[id_name] = (fasta_record)
                        else:
                            backward_dict[id_name] = (fasta_record)



                    else:

                        fasta_dict[id_name] = (fasta_record)



        # for hit_name, id_names in seq_count.items():
        #     if len(id_names) > 1:
        #         print('sorting')
        #         for id_name in sorted(id_names, key=utilities.sort_func):
        #             utilities.createFasta(fasta_dict[id_name], "./fasta_folder/" + filename + "_" + str(count),
        #                                   align)
        #
        #             print (id_name)
        #
        #             fasta_dict.pop(id_name)
        #             count += 1

    if fasta_dict:

        print (fasta_dict)


        print ("Writing out to " + filename)

        outpath = "./fasta_folder/" + filename + ".fasta"

        utilities.createFasta(fasta_dict.values(), "./fasta_folder/" + filename, align)

        return outpath

    else:

        forward_path = ""
        backward_path = ""


        if forward_dict:
            print("Writing out forward dict to ./fasta_folder/" + filename + "_forward")

            utilities.createFasta(forward_dict.values(), "./fasta_folder/" + filename + "_forward", align)

            forward_path = "./fasta_folder/" + filename + "_forward.fasta"

        if backward_dict:
            print("Writing out backward dict to ./fasta_folder/" + filename + "_backward")

            utilities.createFasta(backward_dict.values(), "./fasta_folder/" + filename + "_backward", align)

            backward_path = "./fasta_folder/" + filename + "_backward.fasta"

        return " and ".join(forward_path, backward_path)




def write_genome_order(genomes, split_strands=True, path="./fasta_outputs/genome_order.txt"):

    # Clear previous file if it exists
    open(path, 'w').close()

    for genome in genomes:
        print ('genome')
        print (genome.name)

        hits = sorted([(int(hit.start), hit.region + "_" + hit.strand if split_strands else hit.region) for hit in
                       genome.hits if 'expanded' in
                       hit.region])

        regions = [x[1] for x in hits]

        print(regions)

        renamed_regions = utilities.rename_duplicates(genome.name, regions)


        with open(path, "a") as genome_order:
            genome_order.write(">" + genome.name + "\n")
            genome_order.write(",".join(x for x in renamed_regions) + "\n")

def write_mlgo_order(genomes, split_strands=True, path="./fasta_outputs/genome_order.txt"):

    # Clear previous file if it exists
    open(path, 'w').close()

    region_id_count = 1
    seen_dict = {}

    for genome in genomes:
        print ('genome')
        print (genome.name)

        hits = sorted([(int(hit.start), hit.region +"_strand=" + hit.strand) for hit in
                       genome.hits if 'expanded' in
                       hit.region])

        regions = [x[1] for x in hits]

        print ('here come the regions')

        print(regions)


        # If we've already seen this region in another genome, give it that number, otherwise create new number
        for region in regions:
            print (region)
            if region.split("_strand=")[0] in seen_dict.keys():
                pass
            else:
                seen_dict[region.split("_strand=")[0]] = str(region_id_count)
                region_id_count += 1



        renamed_regions = ["-" + seen_dict[region.split("_strand=")[0]] if region.split("_strand=")[1] == 'backward'
                           else
                           seen_dict[region.split("_strand=")[0]] for region in regions]




        #

        with open(path, "a") as genome_order:
            genome_order.write(">" + genome.name + "\n")
            genome_order.write(" ".join(x for x in renamed_regions) + " $\n")