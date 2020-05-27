import Bio
from views import models
import regex as re


def search_for_promoters(mismatch):
    genomes = models.GenomeRecords.objects()

    regions = []

    prom_regex = "(TTGACA.{15,25}TATAAT){s<=" + mismatch + "}"

    for g in genomes:
        # if g.name == "NZ_CP010029.1":
        for hit in g.hits:
            if "expanded" in hit.region:
                print (hit.region)
                print (hit.start)
                print (hit.end)
                if hit.strand == "forward":
                    seq_content = g.sequence[int(hit.start) - 50: int(hit.start)]

                    print (seq_content)
                    regions.append(seq_content)
                    match = re.findall(prom_regex, seq_content)
                    print(match)
                    if match:
                        hit.promoter = True
                    else:
                        hit.promoter = False


                elif hit.strand == "backward":
                    seq_content = g.sequence[int(hit.end): int(hit.end) + 50]

                    rev_content = str(Bio.Seq.Seq(str(seq_content)).reverse_complement())
                    print (rev_content)
                    regions.append(rev_content)
                    match = re.findall(prom_regex, rev_content)
                    print(match)
                    if match:
                        hit.promoter = True
                    else:
                        hit.promoter = False

        g.save()


def clear_all_promoters():
    """
    Clear all the promoters for all the hits
    :return:
    """
    genomes = models.GenomeRecords.objects()

    for g in genomes:
        for hit in g.hits:
            hit.promoter = False

