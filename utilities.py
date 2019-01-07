import phyloisland
import models
import os
from flask import flash

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def read_fasta(filename):
    """
    Read in a FASTA file
    :param filename:
    :return:
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

def removeFile(*args):
    """
    Remove files in the list from the directory

    :param args: Files to remove
    :return:
    """
    for arg in args:
        os.remove(arg)

def addGenome(genome_results):

    # print ("Here is where I would add a genome")

    # print (genome_results)
    for record in genome_results:

        current = genome_results[record]
        if type(current) == SeqRecord:
            name = current.id
            species = " ".join(current.annotations.get('organism').split()[0:2])
            strain = current.annotations['source']
            sequence = str(current.seq)
            description = current.description

            # Check to see if the genome record already exists
            if  models.GenomeRecords.objects(name=name):
                print("The genome record - %s from species - %s already exists in the database" % (name, species))
                continue

            else:
                print("Adding the genome record - %s from species - %s to the genome database" % (name, species))

                genome = models.GenomeRecords(name=name, species=species, strain=strain, description=description, sequence=sequence)
                genome.save()

def addSequence(seq_records):
    for record in seq_records.values():
        seq_name = record.id
        seq_description = record.description.split(">")[0]
        seq_species = seq_description.split("[")[1].split("]")[0]
        seq_sequence = str(record.seq)

        seq_entry = models.SequenceRecords(seq_name, seq_species, seq_description, seq_sequence, 0, 0)

        # Check if the sequence record already exists
        seq_check = models.SequenceRecords.query.filter_by(name=seq_name).first()
        if seq_check == None:
            print('Adding sequence with ID - %s from species - %s to the sequence database' % (
                seq_name, seq_species) + "\n")
            phyloisland.db.session.add(seq_entry)
            phyloisland.db.session.commit()
        else:
            print('Sequence with ID - %s from species - %s already exists in the sequence database' % (
                seq_name, seq_species) + "\n")

def setProfileAsReference(ids, region):
    if len(ids) > 1:
        flash('Only select a single record', category='error')
    else:
        query = models.Profile.query.filter(models.Profile.uid.in_(ids))
        for record in query.all():
            # Check for a previous reference profile
            old_profile_reference = eval("models.Profile.query.filter_by(" + region + "_profile_ref=1).first()")

            if (old_profile_reference):
                # Remove the previous reference profile
                setattr(old_profile_reference, region + "_profile_ref", 0)
                phyloisland.db.session.add(old_profile_reference)

            # Set the new reference profile
            setattr(record, region + "_profile_ref", 1)

            # Commit the changed record
            phyloisland.db.session.add(record)
            phyloisland.db.session.commit()

            # Write the new profile to the tmp folder ready to be used
            with open("tmp/" + region + "profile.hmm", 'w') as profile_path:
                profile_path.write(record.profile.decode('utf-8'))

            flash("The profile named %s has been set as the reference profile for %s" % (record.name, region), category='success')