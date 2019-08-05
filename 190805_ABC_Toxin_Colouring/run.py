from ete3 import Tree
import colour_tree
import attributes


# Define the tree to use
tree = Tree('./Files/abc_toxins.nwk')

attribute_dict = attributes.attribute_dict

# Define the colours to use for each type - list of colours are here -
colour_dict = {'type1_single' :'dodgerblue', 'type2b_single': 'gold', 'type2a_single':'forestgreen', 'multiple' : 'red', 'unknown' : 'black'}

output_tree, ts = colour_tree.colour_tips(tree, attribute_dict, colour_dict, outpath='./Output/coloured_tips22.png',
                                          custom_layout=False)