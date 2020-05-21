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

def get_domains(domains):
    pos_dict = {'RBD_A': 0, 'RBD_C': 1, 'RBD_B': 2, 'Neuraminidase': 3, 'RBD_D': 4, 'TcB_BD_seed': 5}

    domain_list = [
        # seq.start, seq.end, shape, width, height, fgcolor, bgcolor
        [10, 70, "[]", None, 20, "black", "rgradient:lightgreen", "arial|3|black|RBD_A"],
        [80, 140, "[]", None, 20, "black", "rgradient:blue", "arial|3|black|RBD_C"],
        [150, 210, "[]", None, 20, "black", "rgradient:orange", "arial|3|black|RBD_B"],
        [220, 280, "[]", None, 20, "black", "rgradient:purple", "arial|3|black|NMD"],
        [290, 350, "[]", None, 20, "black", "rgradient:gray", "arial|3|black|RBD_D"],
        [360, 430, "[]", None, 20, "black", "rgradient:darkgreen", "arial|3|black|TCB_BD"]]

    print (domains)

    for k, v in pos_dict.items():
        if k not in domains:
            print (k + "was not there")
            domain_list[v][6] = 'white'
            domain_list[v][7] = "arial|3|black|"


    return domain_list


def get_region_domains(regions):
    region_domains = []

    region_colour_dict = {"TcdA1": 'purple', 'A1': 'orange', 'A2': 'red', 'TcB': 'dodgerblue', 'TcC': 'pink',
                          'Chi' : 'green', 'Chi/Chi': 'green', 'TcB/TcC' : 'yellow' }

    start = 20
    for region in regions:
        print ("**")
        print (region)
        region = region.replace("Chitinase", "Chi")
        region_split = region.split("_joined_")

        region_short = "/".join([region.split("_")[0] for region in region_split])
        if region_short == 'Chitinase':
            region_name = 'Chi'
        else:
            region_name = region_short
        region_domain = [start, start + 120, ">" if 'forward' in region else "<", None, 40, "black",
                         "rgradient:" + region_colour_dict[region_short], "arial|4|black|" + region_name]
        region_domains.append(region_domain)
        start += 130

    return region_domains


def get_leaf_tag(node, skip_tags=['Single', 'Multiple']):
    for descend in node.get_descendants():
        if descend.is_leaf():
            descend.name

            long_name = descend.name.split("_joined")[0]
            short_name = descend.name.split("_information")[0]

            if long_name in tag_dict:
                tag = [x for x in tag_dict[long_name] if x not in skip_tags]
            elif short_name in tag_dict:
                tag = [x for x in tag_dict[short_name] if x not in skip_tags]

            return tag[0]


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

    # print ('set was ')

    #     print (tag_set)

    if len(tag_set) > 1:
        return True
    else:
        check = False

    return check


def get_example_tree(tree, tag_dict, colour_dict, region_dict, region_order_dict, outpath, full_names=False,
                     highlight_nodes=[], collapse_on_genome_tags=True):
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.branch_vertical_margin = 15
    ts.layout_fn = lambda x: None
    #


    # Get the colours for each extant genome
    for node in tree.iter_descendants("postorder"):
        node_style = ""
        print (node.name)
        if (collapse_on_genome_tags and continue_wo_collapse(node) == True) or node.is_leaf():
            print ('here')
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

                    # print ('short name was in tag dict')

                    if len(tag_dict[short_name]) > 2:
                        if "Type?" in tag_dict[short_name]:
                            colour = 'pink'
                        else:
                            colour = 'red'

                    elif len(tag_dict[short_name]) == 2:
                        colour = colour_dict[tag_dict[short_name][1]]

                        # print ('colour was ')
                        #
                        # print (colour)

                    else:
                        colour = 'black'


                        #                         if tag_dict[short_name][0] in colour_dict:
                        #                             colour = colour_dict[tag_dict[short_name][0]]
                        #                         else:
                        #                             print (tag_dict[long_name][0] + " wasn't there")
                        #                             colour = 'black'
                else:
                    # print ('short and long both were not in tag dict')
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
                    region_names = [x for x in region_dict[cleaned_name].keys()]

                if region_dict:

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

                    wid = max(len(node) * 10, 10)

                    #                 node.add_face(ImgFace("./triangle.png",width=wid), 0)
                    box_domains = [10, 200, "<", None, wid, "black", "rgradient:" + colour_dict[tag],
                                   "arial|4|black|" + tag]

                    seqFace = SeqMotifFace(seq=None, motifs=[box_domains], gap_format="line")
                    node.add_face(seqFace, 0, "branch-right")
                    node.img_style['draw_descendants'] = False

        if node_style != 'skip':
            nstyle = NodeStyle()
            nstyle["fgcolor"] = colour
            nstyle["size"] = 20
            node.set_style(nstyle)

    # ts.mode = "c"
    # ts.arc_start = -180 # 0 degrees = 3 o'clock
    # ts.arc_span = 180
    ts.root_opening_factor = 1

    return tree, ts


def colour_tips(tree, tag_dict, colour_dict, region_dict=None, region_order_dict=None, outpath=None,
                full_names=False, highlight_nodes=[],
                custom_layout=False, collapse_on_genome_tags=False):
    tree, ts = get_example_tree(tree, tag_dict, colour_dict, region_dict=region_dict,
                                region_order_dict=region_order_dict,
                                outpath=outpath,
                                full_names=full_names,
                                collapse_on_genome_tags=collapse_on_genome_tags)
    # if outpath:
    # ts.layout_fn = lambda x: None

    print ('outpath is ')
    print (outpath)
    tree.render(outpath, dpi=300, tree_style=ts, w=600)

def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", help="Path to tree", required=True)
    parser.add_argument("-s", "--seqs", help="Path to sequences")
    parser.add_argument("-td", "--tag_dict", help="Path to tag dict")
    parser.add_argument("-rd", "--region_dict", help="Path to region dict")
    parser.add_argument("-rod", "--region_order_dict", help="Path to region order dict")

    parser.add_argument("-cd", "--colour_dict", help="Path to colour dict")

    parser.add_argument("-o", "--outpath", help="Outpath", default="treegaze.png")
    parser.add_argument("-fn", "--full_names", help="Add full name to leaf node", default=False)

    parser.add_argument("-cgt", "--collapse_on_genome_tags", help="Collapse tree based on genome tags", default=False)


    return parser.parse_args(args)

def pickle_open(filename):
    with open(filename, 'rb') as handle:
        pickle_dict = pickle.load(handle)
    return pickle_dict


if __name__ == "__main__":

    print ('called')

    parser = parse_args(sys.argv[1:])

    loaded_tree = load_tree(parser.tree)
    outpath = parser.outpath

    tag_dict = pickle_open(parser.tag_dict)
    region_dict = pickle_open(parser.region_dict)
    region_order_dict = pickle_open(parser.region_order_dict)

    colour_dict = pickle_open(parser.colour_dict)

    full_names = True if parser.full_names == "True" else False

    collapse_on_genome_tags = True if parser.collapse_on_genome_tags == 'True' else False

    print ('tag dict')
    print (tag_dict)
    print ('region dict')

    print (region_dict)
    print ('region order dict')

    print (region_order_dict)
    print ('colour dict')
    print (colour_dict)

    colour_tips(loaded_tree, tag_dict, colour_dict, region_dict, region_order_dict, outpath, full_names=full_names,
                collapse_on_genome_tags=collapse_on_genome_tags)





    # colour_tips(loaded_tree, tag_dict, colour_dict, region_dict, outpath=outpath,
    #                       custom_layout=True)


