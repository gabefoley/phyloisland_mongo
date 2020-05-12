import getGenomes
import utilities

seq_dict = {
    'NZ_CP014947.1_information_Pseudomonas_koreensis_region_A1_4428453_4431915_forward' : 1,
    'NZ_CP014947.1_information_Pseudomonas_koreensis_region_A1_5387500_5389606_backward' : 2
}

def check_output(generated_output, correct):

    generated = utilities.read_fasta(generated_output)

    check = sorted([seq_dict[x] for x in generated.keys()])

    assert(check == correct)


def generate_output(region, outpath, include_genome=[""], exclude_genome=[""], include_hits=[""], exclude_hits=[""]):

    outpath = "./files/fasta_output/" + outpath

    generated_output = getGenomes.download_fasta_regions(region, "", include_genome, exclude_genome, \
                                                include_hits,
                                                exclude_hits, True, False, False, outpath)
    return generated_output

def test_single_genome():

    generated_output = generate_output("A1", 'test_single', include_genome =['possum'])

    check_output(generated_output, [1, 2])