from Bio import SeqIO
import operator
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def get_genome_dict(genome_path, report_path):
    genome_dict = {}
    report_dict = get_report_dict(report_path)

    print (report_dict['assembly_level'])

    if report_dict['assembly_level'] in ['Complete Genome', 'Chromosome']:
        seq_dict = get_complete_seq_dict(genome_path)

    elif report_dict['assembly_level'] in ['Contig', 'Scaffold']:
        seq_dict = {}
        seq = get_concatenated_seq(genome_path)
        seq_dict[report_dict['wgs_project']] = seq


    for genome_seq in seq_dict.keys():
        seq_record = SeqRecord(Seq(seq_dict[genome_seq]), id=genome_seq.split(" ")[0], name=report_dict['o_name'],
                               description=genome_seq,
                           annotations={"organism": report_dict['o_name'] or "",
                                        "assembly_name" :report_dict['assembly_name'] if "assembly_name" in report_dict else "",
                                        "biosample": report_dict['biosample'] if "biosample" in report_dict else "",
                                        "bioproject": report_dict['bioproject'] if "bioproject" in report_dict else "",
                                        "date": report_dict['date'] if "date" in report_dict else "",

                                        "taxid": report_dict['taxid'] if "taxid" in report_dict else "",
                                        "assembly_type": report_dict['assembly_type'] if "assembly_type" in report_dict else "",
                                        "release_type": report_dict['release_type'] if "release_type" in report_dict else "",
                                        "assembly_level": report_dict['assembly_level'] if "assembly_level" in report_dict else "",
                                        "genome_representation": report_dict['genome_representation'] if "genome_representation" in report_dict else "",
                                        "wgs_project": report_dict['wgs_project'] if
                                        "wgs_project" in report_dict else "",

                                        "expected_final_version": report_dict['expected_final_version'] if "expected_final_version" in report_dict else "",
                                        "genome_coverage": report_dict['genome_coverage'] if
                                        "genome_coverage" in report_dict else "",

                                        "excluded": report_dict['excluded'] if "excluded" in report_dict else "",
                                        "genbank_accession_id": report_dict['genbank_accession_id'] if "genbank_accession_id" in report_dict else "",
                                        "refseq_accession_id": report_dict['refseq_accession_id'] if "refseq_accession_id" in report_dict else "",
                                        "r_g_identical": report_dict['r_g_identical'] if "r_g_identical" in report_dict else "",

                                        "source": "",
                                        "plasmid": True if "plasmid" in genome_seq.lower() else False,

                                        })
        genome_dict[genome_seq] = seq_record

    return genome_dict




def get_complete_seq_dict(genome_path):
    seq_dict = {}

    my_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))

    for r in sorted(my_dict.values(), key=operator.attrgetter('id')):
        seq_dict[r.description] = str(r.seq)
    return seq_dict

def get_concatenated_seq(genome_path):
    concatenated_genome = ""

    my_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))

    for r in sorted(my_dict.values(), key=operator.attrgetter('id')):
        print (r.id)
        concatenated_genome += str(r.seq)

    return concatenated_genome

def get_report_dict(report_path):
    report_dict = {}
    for line in open(report_path):

        if line.startswith("# Assembly name:"):
            report_dict['assembly_name'] = line.split(":")[1].strip()

        if line.startswith('# Organism name:'):
            report_dict['o_name'] = line.split(":")[1].strip()

        if line.startswith('# Taxid:'):
            report_dict['taxid'] = line.split(":")[1].strip()

        if line.startswith('# BioSample:'):
            report_dict['biosample'] = line.split(":")[1].strip()

        if line.startswith('# BioProject:'):
            report_dict['bioproject'] = line.split(":")[1].strip()

        if line.startswith('# Date:'):
            report_dict['date'] = line.split(":")[1].strip()

        if line.startswith('# Expected final version:'):
            report_dict['expected_final_version'] = line.split(":")[1].strip()

        if line.startswith('# Assembly type:'):
            report_dict['assembly_type'] = line.split(":")[1].strip()

        if line.startswith('# Release type:'):
            report_dict['release_type'] = line.split(":")[1].strip()

        if line.startswith('# Assembly level:'):
            report_dict['assembly_level'] = line.split(":")[1].strip()

        if line.startswith('# Genome representation:'):
            report_dict['genome_representation'] = line.split(":")[1].strip()

        if line.startswith('# WGS project:'):
            report_dict['wgs_project'] = line.split(":")[1].strip()

        if line.startswith('# Expected final version:'):
            report_dict['expected_final_version'] = line.split(":")[1].strip()

        if line.startswith('# Genome coverage:'):
            report_dict['genome_coverage'] = line.split(":")[1].strip()

        if line.startswith('# Excluded from RefSeq:'):
            report_dict['excluded'] = line.split(":")[1].strip()

        if line.startswith('# GenBank assembly accession:'):
            report_dict['genbank_accession_id'] = line.split(":")[1].strip()

        if line.startswith('# RefSeq assembly accession:'):
            report_dict['refseq_accession_id'] = line.split(":")[1].strip()

        if line.startswith('# RefSeq assembly and GenBank assemblies identical:'):
            report_dict['r_g_identical'] = line.split(":")[1].strip()
    return report_dict


def read_genome(outpath, genome_id, hit_list=None):

    # Collate all of the nucleotide records together to make the genome
    concatenated_genome = ""
    concatenated_plasmid_genome = ""

    my_dict = SeqIO.to_dict(SeqIO.parse(outpath, "fasta"))

    genome_dict = {}

    for r in sorted(my_dict.values(), key=operator.attrgetter('id')):

        if 'plasmid' in r.description:

            description = r.description
            genome_id = genome_id
            species_name = " ".join(r.description.split(",")[0].split(" ")[1:])

            concatenated_plasmid_genome += str(r.seq)

            # Temporary measure to reduce the name so it can fit in the database. Edge case but occurs with
            # 'bacterium endosymbiont of Mortierella elongata FMR23-6', for example

            if len(species_name) > 40:
                species_name = species_name[0:40]

            seq_record = SeqRecord(Seq(concatenated_genome), id=genome_id, name=species_name, description=description,
                                   annotations={"organism": species_name, "source": "", "plasmid": True })

            print(species_name)
            print(seq_record)

            if hit_list and " ".join(species_name.split()[0:2]) in hit_list:

                genome_dict[description] = seq_record
            elif not hit_list:
                print('here')
                genome_dict[description] = seq_record


        else:

            genome = str(r.seq)
            description = r.description
            genome_id = genome_id
            species_name = " ".join(r.description.split(",")[0].split(" ")[1:])

            concatenated_genome += str(r.seq)


            # print ('here')
            #
            # print (r)
            # print (r.description)
            # print (r.description.split(",")[0])

            print (species_name)


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

        seq_record = SeqRecord(Seq(concatenated_genome), id=genome_id, name= species_name, description=description,
                       annotations={"organism": species_name, "source": "", "plasmid" : True if 'plasmid' in
                                                                                                description else False, })

    print (species_name)
    print (seq_record)


    if hit_list and " ".join(species_name.split()[0:2]) in hit_list:

        genome_dict[description] = seq_record
    elif not hit_list:
        print ('here')
        genome_dict[description] = seq_record

    print (len(concatenated_genome))

    print(genome_dict)
    return genome_dict

