from ete3 import Tree, TreeStyle, TextFace, add_face_to_node, SeqMotifFace, NodeStyle, faces, ImgFace, CircleFace, AttrFace


# Add and style the node name text
def layout(node):
    print ('called')
    if node.is_leaf():
        N = AttrFace("name", fsize=60, colour='green')
        faces.add_face_to_node(N, node, 0, position="aligned", colour='blue')
        node.show_leaf_name = False

def colour_tips(tree, attribute_dict, colour_dict, outpath=None, custom_layout=False):

    def get_example_tree():

        # Label all internal nodes
        edge = 0
        for node in tree.traverse():
            if not node.is_leaf():
                node.name = "N%d" % edge
                edge += 1

        # Get the colours for each extant genome
        for node in tree.iter_descendants("postorder"):
            if node.is_leaf():
                spaced_name = node.name.replace("_", " ")
                if spaced_name in attribute_dict.keys():  # If we have a match in the attribute dict
                    colour = colour_dict[attribute_dict[spaced_name]]  # Get the associated colour for this attribute
                    
                if node.name == "Serratia_entomophila" or node.name == "Yersinia_entomophaga":
                    nameFace = TextFace("  " + spaced_name, fsize=15, fgcolor='red')
                    
                else:
                    nameFace = TextFace("  " + spaced_name, fsize=15, fgcolor='black')

                node.add_face(nameFace, column=0)

            else:
                colour = 'black'

            nstyle = NodeStyle()
            nstyle["fgcolor"] = colour
            nstyle["size"] = 20
            node.set_style(nstyle)
            


            ts = TreeStyle()
            ts.show_leaf_name = False

            if custom_layout:
                ts.layout_fn = layout
                ts.show_leaf_name = False

        # ts.mode = "c"
        ts.root_opening_factor = 1

        return tree, ts

    tree, ts = get_example_tree()
    if outpath:
        tree.render(outpath, dpi=300, tree_style=ts)

    return tree, ts




def colour_internal(tree, attribute_dict, colour_dict, first_attr, second_attr, outpath=None, custom_layout=False):

    def get_child_agreement(child1, child2):
        if child1 == child2:
            return child1
        else:
            if child1 == '2':
                return child2
            if child2 == '2':
                return child1
            else:
                return 'unknown'

    def get_example_tree():
        edge = 0
        for node in tree.traverse():
            if not node.is_leaf():
                node.name = "N%d" % edge
                edge += 1

        for node in tree.iter_descendants("postorder"):
            spaced_name = node.name.replace("_", " ")
            if node.is_leaf():

                if spaced_name in attribute_dict.keys() and (attribute_dict[spaced_name] == first_attr or
                                                                         attribute_dict[spaced_name] == second_attr):
                    colour = colour_dict[attribute_dict[spaced_name]]

                else:
                    colour = 'green'
                    
                nameFace = AttrFace("name", fsize=30, fgcolor='red')
                faces.add_face_to_node(nameFace, node, 0, position="branch-right")
                nameFace.background.color = '#eeeeee'
                nameFace.border.width = 1
            else:
                num = get_child_agreement(attribute_dict[node.children[0].name.replace("_", " ")],
                                          attribute_dict[node.children[1].name.replace("_", " ")])
                colour = colour_dict[num]

                attribute_dict[spaced_name] = num


            nstyle = NodeStyle()
            nstyle["fgcolor"] = colour
            nstyle["size"] = 15
            node.set_style(nstyle)

        ts = TreeStyle()
        if custom_layout:
            ts.layout_fn = layout
            ts.show_leaf_name = False

        # ts.mode = "c"
        ts.root_opening_factor = 1

        return tree, ts

    tree, ts = get_example_tree()
    if outpath:
        tree.render(outpath, dpi=300, tree_style=ts)

    return tree, ts
