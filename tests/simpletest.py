import graphbound.graphbound as gb


# # def test_get_positions_with_content_root(simple_tree_4, simple_aln_4):
# #     positions = ac.get_positions_with_content(simple_tree_4, simple_aln_4, node="root")
# #     assert positions == [0,2]
#
# def test_get_positions_with_content_root(simple_tree_6, simple_aln_6):
#     positions = ac.get_positions_with_content(simple_tree_6, simple_aln_6)
#     assert positions == [0,2,3]
#
# # def test_no_matching_label(simple_tree_6, simple_aln_6):
# #     with pytest.raises(NameError):
# #         positions = ac.get_positions_with_content(simple_tree_6, simple_aln_6, node_name="Won't find me")
#
# def test_get_positions_with_content_labelled_position1(simple_tree_6, simple_aln_6):
#     " Test re"
#     positions = ac.get_positions_with_content(simple_tree_6, simple_aln_6, node_name="Labelled_node_1")
#     assert positions == [0,1,2,3]
#
#
# def test_get_positions_with_content_labelled_position2(simple_tree_6, simple_aln_6):
#     " Test re"
#     positions = ac.get_positions_with_content(simple_tree_6, simple_aln_6, node_name="Labelled_node_2")
#     assert positions == [0,1,3]

def test_final_domain_is_at_end_of_sequence():
    pass

def test_sequence_not_in_order_file():
    pass

def test_sequence_from_order_file_not_in_list():
    pass

def test_folder_with_DS_store_works():
    pass

def test_order_list_with_spaces_around_names():
    pass

def test_ABC_simple(ABC_simple_3_region):






graphbound = gb.GraphBound(seq_region_dict)
graphbound.create_alignment()
graphbound.write_graph_to_image("/Users/gabefoley/Dropbox/wowzers.png")
