from ete3 import PhyloTree, Tree, TreeStyle, TextFace, add_face_to_node, SeqMotifFace, NodeStyle, faces, ImgFace, \
    CircleFace, \
    AttrFace
import os
import sys
import argparse
import pickle


def load_phylo_tree(tree_path, aln_path=None):
    """
    Load a tree, associate an alignment with it if given
    """
    tree = PhyloTree(tree_path, alignment=aln_path, format=1, alg_format='fasta')
    return tree


def load_tree(tree_path, aln_path=None):
    """
    Load a tree, associate an alignment with it if given
    """
    tree = Tree(tree_path, format=1)
    return tree


###
###
###


def get_region_domains(regions):
    region_domains = []

    region_colour_dict = {"TcdA1": 'purple', 'A1': 'orange', 'A2': 'red', 'TcB': 'dodgerblue', 'TcC': 'pink',
                          'Chi': 'green', 'Chi/Chi': 'green', 'TcB/TcC': 'yellow', 'A2/TcdA1': 'red'}

    start = 20
    for region in regions:
        # print ("**")
        # print (region)
        region = region.replace("Chitinase", "Chi")
        region_split = region.split("_joined_")

        region_short = "/".join([region.split("_")[0] for region in region_split])
        if region_short == 'Chitinase':
            region_name = 'Chi'
        else:
            region_name = region_short

            if region_short in region_colour_dict:
                region_domain = [start, start + 120, ">" if 'forward' in region else "<", None, 40, "black",
                                 "rgradient:" + region_colour_dict[region_short], "arial|40|black|" + region_name]
            else:
                print("WARNING: Got a value for region colour dict that wasn't in there - " + region_short)
                region_domain = [start, start + 120, ">" if 'forward' in region else "<", None, 40, "black",
                                 "rgradient:black", "arial|40|black|" + region_name]

        region_domains.append(region_domain)
        start += 130

    return region_domains

def get_ancestor_domains(regions):
    ancestor_domains = []

    region_colour_dict = {"TcdA1": 'purple', 'A1': 'orange', 'A2': 'red', 'TcB': 'dodgerblue', 'TcC': 'pink',
                          'Chi': 'green', 'Chi/Chi': 'green', 'TcB/TcC': 'yellow', 'A2/TcdA1': \
        'red'}

    start = 20
    for region in regions:
        # print ("**")
        print (region)
        region = region.replace("Chitinase", "Chi")

        region_name = region.replace("-", "")
        # region_split = region.split("_joined_")

        # region_short = "/".join([region.split("_")[0] for region in region_split])
        # if region_short == 'Chitinase':
        # #     region_name = 'Chi'
        # else:
        #     region_name = region_short

        if region_name in region_colour_dict:
            region_domain = [start, start + 30, ">" if 'forward' in region else "<", None, 20, "black",
                             "rgradient:" + region_colour_dict[region_name], "arial|4|black|" + region_name]
        else:
            print("WARNING: Got a value for region colour dict that wasn't in there - " + region_name)
            region_domain = [start, start + 30, ">" if 'forward' in region else "<", None, 20, "black",
                             "rgradient:black", "arial|4|black|" + region_name]

        ancestor_domains.append(region_domain)
        start += 40

    return ancestor_domains


def get_leaf_tag(node, skip_tags=['Single', 'Multiple']):
    for descend in node.get_descendants():
        if descend.is_leaf():
            # descend.name

            long_name = descend.name.split("_joined")[0]
            short_name = descend.name.split("_information")[0]

            tag = []

            print (long_name)
            print (short_name)
            if long_name in tag_dict:
                tag = [x for x in tag_dict[long_name] if x not in skip_tags]
            elif short_name in tag_dict:
                tag = [x for x in tag_dict[short_name] if x not in skip_tags]

            if tag:
                return tag[0]
            else:
                return 'unknown'


def continue_wo_collapse(node, skip_tags=["Single", "Multiple"]):
    """
    If there is still a node in highlight nodes that is a descendant of current node, we want to keep going
    """
    check = False

    leaves = []

    for descend in node.get_descendants():
        if descend.is_leaf():
            leaves.append(descend)

    names = [x.name for x in leaves]

    tag_set = set()

    for name in names:
        long_name = name.split("_joined")[0]
        short_name = name.split("_information")[0]

        if long_name in tag_dict:
            tag = [x for x in tag_dict[long_name] if x not in skip_tags]
        elif short_name in tag_dict:
            tag = [x for x in tag_dict[short_name] if x not in skip_tags]
        else:
            tag = ['missing']

        if len(tag) > 1:
            print("Incorrectly tagged " + name)

        tag_set.add(tag[0])


    if len(tag_set) > 1:
        return True
    else:
        check = False

    return check

def get_domains(domains):
    pos_dict = {'Big_1_Full' : 0, 'VRP1_Full' : 1, 'RBD_A': 2, 'RBD_C': 3, 'RBD_B': 4, 'Neuraminidase': 5,
                'TcA_RBD': 6, 'TcB_BD_seed': 7}

    domain_list = [
        # seq.start, seq.end, shape, width, height, fgcolor, bgcolor
        [10, 70, "[]", None, 20, "black", "rgradient:red", "arial|40|black|Big_1"],
        [80, 140, "[]", None, 20, "black", "rgradient:pink", "arial|40|black|VRP1"],
        [150, 210, "[]", None, 20, "black", "rgradient:lightgreen", "arial|40|black|RBD_A"],
        [220, 280, "[]", None, 20, "black", "rgradient:blue", "arial|40|black|RBD_C"],
        [290, 350, "[]", None, 20, "black", "rgradient:orange", "arial|40|black|RBD_B"],
        [360, 430, "[]", None, 20, "black", "rgradient:purple", "arial|40|black|NMD"],
        [440, 500, "[]", None, 20, "black", "rgradient:gray", "arial|40|black|RBD_D"],
        [510, 560, "[]", None, 20, "black", "rgradient:darkgreen", "arial|40|black|TCB_BD"]]


    for k, v in pos_dict.items():
        if k not in domains:
            print(k + "was not there")
            domain_list[v][6] = 'white'
            domain_list[v][7] = "arial|3|black|"

    return domain_list

def setup_domain_list(region_dict):
    region_list = []

    print ('setup domain list')
    for entry in region_dict.keys():
        region_names = [x[0] for x in sorted(region_dict[entry].items(), key=lambda k: k[1])]
        region_list.append(region_names)

    print (region_list)




def get_example_tree(tree, tag_dict, colour_dict, region_dict, region_order_dict, sequence_content_dict, skip_list, \
                     outpath, \
                     full_names=False,
                     highlight_nodes=[], collapse_on_genome_tags=True):
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.branch_vertical_margin = 15
    ts.layout_fn = lambda x: None

    if region_dict:

        domain_list = setup_domain_list(region_dict)

    #


    # Get the colours for each extant genome
    for node in tree.iter_descendants("postorder"):
        node_style = ""
        if (collapse_on_genome_tags and continue_wo_collapse(node) == True) or node.is_leaf():
            if node.is_leaf():
                node.show_leaf_name = True

                long_name = node.name.split("_joined")[0]
                short_name = node.name.split("_information")[0]

                if long_name in tag_dict:

                    # print ('long name was in tag dict')
                    if tag_dict[long_name][0] in colour_dict:

                        colour = colour_dict[tag_dict[long_name][0]]
                    else:
                        colour = 'black'
                        # print(tag_dict[long_name][0] + " wasn't there")


                elif short_name in tag_dict:

                    tags = [x for x in tag_dict[short_name] if x not in skip_list]

                    if len(tags) > 1:
                        print('\nWARNING: The following genome had multiple tags associated with it')
                        print(long_name)
                        print("The tags were ")
                        print(tags)
                        colour = 'black'
                    else:
                        if tags[0] in colour_dict:
                            colour = colour_dict[tags[0]]

                        else:
                            "WARNING: Colour dict doesn't have an entry for " + tags[0]
                            colour = 'black'


                else:
                    print("\nWARNING: We couldn't find an entry for the following genome in the tag dictionary")
                    print(node.name)
                    colour = 'black'

                if full_names:
                    spaced_name = " ".join(long_name.split("_"))
                else:
                    spaced_name = " ".join(node.name.split("_")[3:5])

                nameFace = TextFace("  " + spaced_name, fsize=15, fgcolor='black')
                node.add_face(nameFace, column=0)

                region_names = []

                cleaned_name = node.name.replace(".", "***")

                if cleaned_name in region_dict:
                    # region_names = [x for x in region_dict[cleaned_name].keys()]
                    region_names = [x[0] for x in sorted(region_dict[cleaned_name].items(), key=lambda k: k[1])]



                    # region_dict = sorted(region_dict.items(), key=lambda k: k[1])

                if region_dict:

                    print ('horseradish')
                    print (region_names)

                    box_domains = get_domains(region_names)

                    # box_domains = get_domains([x for x in region_dict[node.name.replace("***", ".")].keys()])

                    seqFace = SeqMotifFace(seq=None, motifs=box_domains, gap_format="line")
                    node.add_face(seqFace, 0, "aligned")

                elif region_order_dict:

                    region_order_name = node.name.split("_information_")[0].replace(".", "***")

                    if region_order_name in region_order_dict:
                        regions = region_order_dict[region_order_name].split(",")

                        region_domains = get_region_domains(regions)

                        seqFace = SeqMotifFace(seq=None, motifs=region_domains, gap_format="line")
                        node.add_face(seqFace, 0, "aligned")
                elif sequence_content_dict:

                    if node.name in sequence_content_dict:
                        S = SeqMotifFace(sequence_content_dict[node.name], seq_format="seq", width=6)
                        node.add_face(S, 1, position="aligned")
                    else:
                        print("\nWARNING: " + node.name + " was not in the sequence content dict")

            else:
                colour = 'black'

            if colour == None:
                colour = 'black'

        else:
            #             print ('collapse')
            if node.name not in highlight_nodes:

                if collapse_on_genome_tags:
                    node_style = 'skip'
                    tag = get_leaf_tag(node)

                    format_text = " extant sequence " if len(node) == 1 else " extant sequences "
                    N = TextFace(" " + str(len(node)) + format_text, fsize=14, fgcolor="black")
                    node.add_face(N, 1, position="branch-right")

                    ceiling = 50
                    length = 100

                    if display_circular or display_circular_180:
                        ceiling = 100
                        length = 800


                    wid = min(len(node) * 10, ceiling)




                    #                 node.add_face(ImgFace("./triangle.png",width=wid), 0)
                    box_domains = [10, length, "<", None, wid, "black", "rgradient:" + colour_dict[tag],
                                   "arial|30|black|" + tag]

                    seqFace = SeqMotifFace(seq=None, motifs=[box_domains], gap_format="line")
                    seqFace.rotable = False
                    node.add_face(seqFace, 0, "branch-right")
                    node.img_style['draw_descendants'] = False

        if node_style != 'skip':
            nstyle = NodeStyle()
            nstyle["fgcolor"] = colour
            nstyle["size"] = 20
            node.set_style(nstyle)

    if display_circular:
        ts.mode = "c"

    if display_circular_180:
        ts.mode = "c"
        ts.arc_start = -180  # 0 degrees = 3 o'clock
        ts.arc_span = 180
    ts.root_opening_factor = 1

    return tree, ts

def get_ancestor_tree(tree, ancestral_order_dict, ref_ml_go_dict,
                     outpath,):

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.branch_vertical_margin = 15
    ts.layout_fn = lambda x: None
    #


    # Get the colours for each extant genome
    for node in tree.iter_descendants("postorder"):
        if not node.is_leaf():
            # node_style = ""
            print (node.name)
            print (ancestral_order_dict)
            ancestral_orders = [ref_ml_go_dict[x.replace("-", "")] for x in ancestral_order_dict[node.name]]
            print (ancestral_orders)

            region_domains = get_ancestor_domains(ancestral_orders)


            if region_domains:

                seqFace = SeqMotifFace(seq=None, motifs=region_domains, gap_format="line")
                node.add_face(seqFace, 0, "branch-right")






    ts.root_opening_factor = 1

    return tree, ts


def colour_tips(tree, tag_dict, colour_dict, region_dict=None, region_order_dict=None, sequence_content_dict=None,
                skip_list=[],
                display_circular=False, display_circular_180=False, outpath=None,
                full_names=False, highlight_nodes=[],
                custom_layout=False, collapse_on_genome_tags=False):
    tree, ts = get_example_tree(tree, tag_dict, colour_dict, region_dict=region_dict,
                                region_order_dict=region_order_dict, sequence_content_dict=sequence_content_dict,
                                skip_list=skip_list,
                                outpath=outpath,
                                full_names=full_names,
                                collapse_on_genome_tags=collapse_on_genome_tags)
    # if outpath:
    # ts.layout_fn = lambda x: None

    print('outpath is ')
    print(outpath)
    tree.render(outpath, dpi=30, tree_style=ts, w=1600)


def display_ancestors(tree, ancestral_order_dict, ref_ml_go_dict, outpath):
    tree, ts = get_ancestor_tree(tree, ancestral_order_dict=ancestral_order_dict, ref_ml_go_dict=ref_ml_go_dict,
                                 outpath=outpath,
                                 )

    tree.render(outpath, dpi=30, tree_style=ts, w=1600)


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", help="Path to tree", required=True)
    parser.add_argument("-s", "--seqs", help="Path to sequences")
    parser.add_argument("-td", "--tag_dict", help="Path to tag dict")
    parser.add_argument("-rd", "--region_dict", help="Path to region dict")
    parser.add_argument("-rod", "--region_order_dict", help="Path to region order dict")
    parser.add_argument("-scd", "--sequence_content_dict", help="Path to sequence content dict")
    parser.add_argument("-cd", "--colour_dict", help="Path to colour dict")
    parser.add_argument("-o", "--outpath", help="Outpath", default="treegaze.png")
    parser.add_argument("-fn", "--full_names", help="Add full name to leaf node", default=False)
    parser.add_argument("-cgt", "--collapse_on_genome_tags", help="Collapse tree based on genome tags", default=False)
    parser.add_argument("-dc", "--display_circular", help="Display tree as circular", action='store_true')
    parser.add_argument("-dco", "--display_circular_180", help="Display tree as circular (180 degrees)",
                        action='store_true')
    parser.add_argument("-ao", "--ancestral_order")
    parser.add_argument("-mlgo", "--ref_ml_go_dict")


    return parser.parse_args(args)


def pickle_open(filename):
    with open(filename, 'rb') as handle:
        pickle_dict = pickle.load(handle)
    return pickle_dict


if __name__ == "__main__":

    parser = parse_args(sys.argv[1:])

    loaded_tree = load_tree(parser.tree)
    outpath = parser.outpath

    if parser.ancestral_order != None:
        ancestral_order_dict = pickle_open(parser.ancestral_order)
        ref_ml_go_dict = pickle_open(parser.ref_ml_go_dict)

        display_ancestors(loaded_tree, ancestral_order_dict, ref_ml_go_dict, outpath)

    else:

        tag_dict = pickle_open(parser.tag_dict)
        region_dict = pickle_open(parser.region_dict)
        region_order_dict = pickle_open(parser.region_order_dict)
        sequence_content_dict = pickle_open(parser.sequence_content_dict)

        colour_dict = pickle_open(parser.colour_dict)

        full_names = True if parser.full_names == "True" else False

        collapse_on_genome_tags = True if parser.collapse_on_genome_tags == 'True' else False

        display_circular = True if parser.display_circular else False
        display_circular_180 = True if parser.display_circular_180 else False
        skip_list = ["Single", "Multiple"]  # Tags to skip when looking for tags to colour on

        print("\nMake tree called with the following dictionaries - ")

        print('tag dict')
        print(tag_dict)
        print('region dict')
        print(region_dict)


        print('region order dict')


        print(region_order_dict)
        print('sequence_content_dict')

        print(sequence_content_dict)
        print('colour dict')
        print(colour_dict)
        print('skip list')
        print(skip_list)

        colour_tips(loaded_tree, tag_dict, colour_dict, region_dict, region_order_dict, sequence_content_dict, skip_list, \
                    display_circular,
                    display_circular_180,
                    outpath, \
                    full_names=full_names,
                    collapse_on_genome_tags=collapse_on_genome_tags)

