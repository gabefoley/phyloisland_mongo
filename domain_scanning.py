import sys
sys.path.append('/Users/gabefoley/Dropbox/Tutoring/public_binfpy/binfpy')
import subprocess
import sequence
import pickle
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import utilities


email = 'gabriel.foley@uqconnect.edu.au'
sequences = sequence.readFastaFile('fasta_outputs/TcB_test.fasta')


def parse_pfamscan(job_id):
    results = json.loads(open(job_id + '.out.txt').read())
    result_dict = {}

    for result in results:
        print (type(result))
        print (result)
        print (result['name'])
        result_dict[result['name']] = result

    return result_dict



def search_sequences_for_domains(filename):

    sequences = sequence.readFastaFile(filename)

    domain_names = set()

    sequence_dict = {}




    for seq in sequences:

        stdoutdata = subprocess.getoutput("python3 pfamscan.py --email %s --database pfam-a --sequence %s" % (email,
                                                                                                             seq))



        print(stdoutdata)

        job_id = stdoutdata.split("JobId:")[1].split("\n")[0].strip()
        print ('Job id is ' + job_id)

        result_dict = parse_pfamscan(job_id)

        for domain_name in result_dict.keys():
            domain_names.add(domain_name)

        sequence_dict[seq.info.split()[0]] = result_dict
        sequence_dict[seq.info.split()[0]]['sequence'] = seq.sequence

    print ("Found the following domain names")
    print (domain_names)
    print(sequence_dict)


    for domain in domain_names:
        fasta_list = []

        for hit_name, hit in sequence_dict.items():

            if domain in hit:

                print (hit[domain]['seq'])
                print(hit[domain]['seq']['to'])
                print(hit[domain]['seq']['from'])



                fasta_record = SeqRecord(Seq("".join(hit['sequence'][int(hit[domain]['seq']['from']):int(hit[domain][
                                                                                                          'seq'][
                                                                                                   'to'])])), id=
                                        hit_name +
                                         "_" +
                                         domain, description="")

                print (fasta_record)

                fasta_list.append(fasta_record)

        print('fasta list')
        print (fasta_list)

        if len(fasta_list) > 0:

            utilities.createFasta(fasta_list, './fasta_outputs/' + domain, align=True)






# job_id = 'pfamscan-R20200210-002446-0151-95364425-p2m'
#
# job_id = 'pfamscan-R20200210-002406-0465-37409833-p1m'
#
# parse_pfamscan(job_id)

# search_sequences_for_domains('fasta_outputs/TcB_expanded_.cutdown.fasta')

search_sequences_for_domains('fasta_outputs/TcC_expanded_.fasta')