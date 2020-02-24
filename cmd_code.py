import utilities
import models
import getGenomes
import resultread
import cmd_code
import timeit
import Bio
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import os
import time
import subprocess
import glob




# Get the genomes
def get_genomes(args):
    species_names = utilities.readLinesFromFile(args.add_genomes)
    get_genomes_cmd(species_names)


def update_profiles(args):

    # Ensure there is a trailing slash for the profile folder
    profile_folder = args.add_profiles + "/" if not args.add_profiles[-1] == "/" else args.add_profiles

    # For all the profiles in the profile folder, save them to the database and set them as the reference profile
    profiles = [x for x in os.listdir(profile_folder) if x != ".DS_Store"]


    for profile in profiles:
        profile_name = profile.split("_")[0]

        file = open(profile_folder + profile, 'rb')

        # Save the profile to database
        utilities.save_profile(file, profile_name)

        # Get the profile
        saved_profile = models.Profile.objects(name=profile_name).first()

        # Set the profile as the reference
        utilities.set_profile_as_reference([saved_profile.id], profile_name)


def update_genomes():
    # Check the genomes with the profiles
    queries = models.GenomeRecords.objects.all().timeout(False)


    profiles = models.Profile.objects.all()

    for profile in profiles:
        get_feature_location_with_profile_cmd(queries, "hmm_outputs", profile.name, "", "", profile.name)

    del profiles
    del queries


def delete_profiles(args):

    # Ensure there is a trailing slash for the profile folder
    profile_folder = args.add_profiles + "/" if not args.add_profiles[-1] == "/" else args.add_profiles

    # For all the profiles in the profile folder, save them to the database and set them as the reference profile
    profiles = [x for x in os.listdir(profile_folder) if x != ".DS_Store"]

    profile_names = [profile.split("_")[0] for profile in profiles]




    # Delete the profiles from the database (so we can easily update them with new ones if need be)
    profiles_to_delete = models.Profile.objects(name__in=profile_names)
    profiles_to_delete.delete()




def get_genomes_cmd(species_names):

    failed_genomes = []

    start_time = timeit.default_timer()
    for species_name in species_names:
        print("Species name is ", species_name)
        destinations = ['reference genome', "representative genome", "assembly", "genbank"]

        genome_results = getGenomes.add_genome(species_name, destinations, single=True)

        if genome_results and genome_results != "Fail":
            utilities.add_genome(genome_results)

        else:

            # Couldn't find it in RefSeq, let's try genbank
            destinations = ["genbank"]
            genome_results = getGenomes.add_genome(species_name, destinations, single=True)

            if genome_results and genome_results != "Fail":

                utilities.add_genome(genome_results)
            else:
                failed_genomes.append(species_name)


    seconds = timeit.default_timer() - start_time
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)

    periods = [('hours', hours), ('minutes', minutes), ('seconds', seconds)]
    time_string = ', '.join('{} {}'.format(value, name)
                            for name, value in periods
                            if value)

    print('\nFINISHED GETTING RECORDS: Time taken was {} \n'.format(time_string))
    if failed_genomes:
        print("The following genomes failed to map - " + " ".join(x for x in failed_genomes))

def get_feature_location_with_profile_cmd(queries, reference, profile_name, recordName, recordLocation, region):
    """
    Annotate a genome sequence with a feature location based on a profile
    :param ids: Genome sequences to annotate
    :param reference: Profile to annotate based on
    :param recordName: Which feature field to update
    :param recordLocation: Which feature location field to update
    :return:
    """


    for query in queries:

        seq_record = query.sequence

        # Create a path to write the translated genomic sequence to
        random_id = utilities.randstring(5)

        # Get the nucleotide sequence of the genome
        nuc_seq = Bio.Seq.Seq(str(query.sequence))

        seq_record = SeqRecord(Seq(query.sequence, generic_dna), '', '', '')

        # seq_record = SeqRecord(seq=Seq(nuc_seq, Alphabet()), id='', name='', description='', dbxrefs=[])



        outpath = reference + "/" + query.name + "/" + query.species + "/" + region + "/"
        outpath = outpath.replace(" ", "_")

        print ("outpath is " + outpath)

        # Check three forward reading frames
        if not os.path.exists(outpath):
            os.makedirs(outpath.replace(" ", "_"))

        for forward in [True, False]:
            for i in list(range(0, 3)):

                strand = "_forward_" + str(i) if forward else "_backward_" + str(i)
                sequence = nuc_seq[i:] if forward else nuc_seq.reverse_complement()[i:]

                cleaned_path = outpath + query.name.replace("/",
                                                               "_") + "_" + random_id + strand + "_translated_genome.fasta"
                hmmsearch_results = outpath + query.name.replace("/",
                                                                    "_") + "_" + random_id + strand + "_hmmsearch_results.fasta"
                domScore = 100

                cleaned_path = cleaned_path.replace(" ", "_")
                hmmsearch_results = hmmsearch_results.replace(" ", "_")

                # Translate the nucleotide genome sequence to a protein sequence
                with open(cleaned_path, 'w') as handle:

                    if query.name == "<unknown name>":
                        handle.write(">" + query.description + "\n" + str(
                            sequence.translate(stop_symbol="*")))
                    else:

                        handle.write(">" + query.description + "\n" + str(
                            sequence.translate(stop_symbol="*")))

                print("Writing the %s sequence with the species %s to %s" % (
                query.name, query.species, cleaned_path))

                while not os.path.exists(cleaned_path):
                    time.sleep(1)

                if os.path.isfile(cleaned_path):
                    stdoutdata = subprocess.getoutput("hmmsearch -o %s --domT %s %s %s" % (
                    hmmsearch_results, domScore, 'tmp/' + region + "_profile.hmm", cleaned_path))

                    print(stdoutdata)
                    # result = subprocess.call(["hmmsearch -o %s %s %s" % (hmmsearch_results, reference, cleaned_path)])

                    print("The results from the HMM search have been written to %s \n" % hmmsearch_results)
                     # read_hmmer_results(hmmsearch_results)
                    # result = subprocess.call(["hmmsearch", 'files/output.txt', reference, cleaned_path], stdout=subprocess.PIPE)
                    # for x in result:
                    #     print (x)

        hmmerout = []
        hmmerout_expanded = []

        reg = os.path.join(reference + "/" + query.name + "/" + query.species.replace(" ", "_") + "/"  + profile_name)

        hmmerout.append(resultread.HMMread(reg, query))
        hmmerout_expanded.append(resultread.HMMread(reg, query, expand=True))

        for infile in glob.glob(reg + '/*.fasta'):
            utilities.remove_file(infile)

def get_overview():
    queries = models.GenomeRecords.objects.all().timeout(False)
    print (f"There are {len(queries)} genomes stored in the database\n")

    print ("The genomes are - \n")


    for query in queries:
        print (query.description + "\n")

        print (query.hits)

        for hit in sorted(query.hits, key=lambda hit: int(hit.start)):
            if 'expanded' in hit.region:
                print (f"{hit.region} - {hit.start} : {hit.end}" )



        print ("And it has hits for " + " and ".join([hit.region for hit in query.hits if "expanded" not in
                                                      hit.region]))
        print ("And it has the following tags - " + " ".join([tag for tag in query.tags]))

        # print ("Gene order is " + " ".join(hit.region for hit ))
        print ()

