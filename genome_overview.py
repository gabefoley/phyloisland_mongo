from Bio.Graphics import GenomeDiagram
import utilities
from Bio.SeqFeature import SeqFeature, FeatureLocation



def writeHMMToImage(hmm_dict, reference, seq_record, name, query_id, species, expand=False):
    # Currently taking region to be run simultaneously with hmmer, will need to manipulate some stuff for later
    # Initiation of GenomeDiagram assets"

    # Convert the length of the genome into the length of the translated genome.
    genome_length = round(len(seq_record) / 3)

    # Parameters that can be set

    # The highest number of tracks we'll attempt to add (if a feature overlaps on all tracks it'll just be added to the bottom track)
    # Set max_tracks to None to just add as many as we need to get no overlap
    max_tracks = None

    # The amount we want to allow a feature to overlap with a neighbouring feature
    overlap_amount = 0

    name = reference + "_" + species + "_GenomeDiagram"
    name += "_expanded" if expand else ""



    gd_diagram = GenomeDiagram.Diagram(name)
    max_len = 0
    output_path = "static/user_images/%s_%s%s.png" % (query_id, species, "_expanded" if expand else "")

    output_path = output_path.replace(" ", "_")

    start = 0
    # For my work I was considering changing 'region1, 2, and 3' to a3, TcB, and TcC for convenience
    # Up to others though if I fully change that (is just a UI thing tbh)
    region_colours = {"A1": "orange", "A2": "red", "Chitinase": "green", "A3": "yellow",
                      "TcB": "blue", "TcC": "magenta", "pore": "grey", "region1": "lightblue", "region2": "pink",
                      "region3": "purple", "region4": "black"}
    locs = {}
    strand_dict = {}
    strandd = 1
    for result in hmm_dict:
        i = 0
        for reg in result:
            """ Create a dictionary for key = feature type -> value = location """
            if "forward" in reg:
                strandd = 1
                location = reg.split("/")[2] + utilities.randstring(5)
                locs[location] = result[reg].split(":")
                strand_dict[location] = strandd

            elif "backward" in reg:
                strandd = -1
                location = reg.split("/")[2] + utilities.randstring(5)
                locs[location] = result[reg].split(":")
                strand_dict[location] = strandd



            i += 1
            # Prepare for literally the worst code in existence
            """ We have to pull the features from location data in database
            instead of directly from database because of current limitations """
            # Brute force if statements to pull from database
            # TODO - Could likely turn this into a for loop in so

    """ Add the features from dictionary to the seq_record features """
    feature_locations = []
    # If want to add all features to image comment below line out
    seq_record.features = []
    i = 0
    for location in locs:
        """ Extract start and end values from each location, and add to independent lists """
        """ create and add features based on locations """
        feature = SeqFeature(
            location=FeatureLocation(int(locs[location][0]), int(locs[location][1]), strand=strand_dict[location]),
            type=location[0:-5])
        seq_record.features.append(feature)

    """ Set up the Genome Diagram """
    max_len = max(max_len, len(seq_record))

    gd_track1 = gd_diagram.new_track(0, name=name + " Track 1", greytrack=True, start=0,
                                     end=len(seq_record))
    gd_feature_set1 = gd_track1.new_set()

    for feature in seq_record.features:
        feature_locations.append((feature.location.start, feature.location.end))


    total_tracks = 1
    # Dictionary to keep track of which locations are at which track
    forward_tracks = {1: []}
    backward_tracks = {1: []}

    for feature in seq_record.features:
        feature_added = False
        current_track = 1

        while current_track <= total_tracks and not feature_added:

            current_strand = feature.location.strand

            if current_strand == 1:
                dict_track = forward_tracks
            elif current_strand == -1:
                dict_track = backward_tracks

            overlap = False

            for loc in dict_track[current_track]:

                if feature.location.start + overlap_amount in loc or feature.location.end - overlap_amount in loc:
                    overlap = True

            if overlap:

                # If we've reached the highest track, make a new track
                if current_track == total_tracks:
                    if max_tracks is not None and current_track == max_tracks:

                        # We've reached the highest track we want to try and add, so just add it to the bottom track (there will be overlap)
                        gd_feature_set1.add_feature(feature, label=True, name=feature.type,
                                                    color=region_colours[feature.type], label_position='middle')
                        dict_track[1].append(feature.location)
                        feature_added = True
                    else:
                        # Add to next highest track
                        current_track += 1
                        total_tracks += 1
                        if current_track not in dict_track.keys():
                            dict_track[current_track] = []
                        if current_track not in forward_tracks.keys():
                            forward_tracks[current_track] = []
                        if current_track not in backward_tracks.keys():
                            backward_tracks[current_track] = []
                        exec("gd_track" + str(
                            current_track) + "= gd_diagram.new_track(0, name=name + ' Track " + str(
                            current_track) + "', greytrack=True, start=0, end=len(seq_record))")
                        exec("gd_feature_set" + str(current_track) + " = gd_track" + str(current_track) + ".new_set()")

                else:
                    # Check on higher track
                    current_track += 1

            else:
                # Add to this track
                exec("gd_feature_set" + str(
                    current_track) + ".add_feature(feature, label=True, name=feature.type, color=region_colours[feature.type], label_position='middle')")
                if current_track in dict_track:
                    dict_track[current_track].append(feature.location)
                else:
                    dict_track[current_track] = [feature.location]
                feature_added = True

    for track in range(1, total_tracks):
        exec("gd_track" + str(track) + ".add_set(gd_feature_set" + str(track) + ")")
    gd_diagram.draw(format="linear", pagesize="A2", fragments=10, start=start, end=len(seq_record))
    gd_diagram.write(output_path, "PNG")
    print("Genome Diagram has been added to file " + output_path)
    fh = open(output_path, 'rb')
    return fh