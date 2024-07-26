import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
import numpy as np
import copy
import itertools
import matplotlib.pyplot as plt
from TrackReconstruction_functions import *



def RunReco(data, part):

    # There seems to be a duplicate row sometimes
    data = data.drop_duplicates()

    # display(data)
    data = data[['x', 'y', 'z',"energy"]]

    # shuffle the data to ensure we dont use g4 ordering
    data = data.sample(frac=1).reset_index(drop=True)

    # then sort it based on the x,y,z
    data = data.sort_values(by=['x', "y", "z"]).reset_index(drop=True)

    # Print energy sum
    # print("Energy Sum: ", data.energy.sum())

    # Calculate the distance matrix
    dist_matrix = distance_matrix(data[['x', 'y', 'z']], data[['x', 'y', 'z']])

    # Initialize connections counter, keeps track of number of connections to each index
    connection_count = np.zeros(len(data), dtype=int)

    # This is a dict, format is
    # index : [connected node 1, connected node 2,...]
    connected_nodes = {}

    # Tunable parameters
    init_dist_thresh = 15 # max distance for initial connections [mm]
    incr_dist_thresh = [2,4,6,8,10,12,14,16,18,20] # Second stage, look for closest nodes, then slowly increase threshold [mm]
    dist_threshold = 15 # Third distance threshold for conecting end nodes [mm]


    connections = []

    # -----------------------------------
    # Find the node that is closest to the vertex
    # Calculate the Euclidean distance from the origin for each row
    data['distance'] = np.sqrt(data['x']**2 + data['y']**2 + data['z']**2)
    vertex_index = data['distance'].idxmin()
    # data = data.drop(columns=['distance'])
    # print("Vertex Index is:" , vertex_index)

    # Make a connection to the two closest nodes
    closest_idx = np.argsort(dist_matrix[vertex_index])[1:3]
    UpdateConnections(vertex_index, closest_idx[0], connected_nodes, connections, connection_count)
    UpdateConnections(vertex_index, closest_idx[1], connected_nodes, connections, connection_count)

    # ------------------------------------
    # Find closest nodes and create connections

    for i in range(len(data)):
        # Find the index of the closest node (excluding itself)
        # closest_idx = np.argpartition(dist_matrix[i], 1)[1]
        closest_idx = np.argsort(dist_matrix[i])[1]
        
        # Check if the connection already exists 
        if closest_idx not in connected_nodes.get(i, []) and i not in connected_nodes.get(closest_idx, []):

            # Check the proposed node has 0 or 1 connection
            if (connection_count[closest_idx] <= 1 and connection_count[i] <= 1 and dist_matrix[i][closest_idx] < init_dist_thresh):
                
                cycle  = Testcycle(i, closest_idx ,connected_nodes, connections, connection_count)
                
                # Add connection between node i and closest_idx if it doesnt form a cycle
                if (not cycle):
                    UpdateConnections(i, closest_idx, connected_nodes, connections, connection_count)

    # Get indices where the value is 1
    single_nodes = np.where(connection_count == 1)[0]

    # Incrementally loop over distance steps looking for connections
    # starting from a small step size helps lock onto the nearest nodes
    for dist in incr_dist_thresh:

        # Connect single nodes to the next closest single node
        for i in single_nodes:
            
            # Connections get updated, so this ensures we dont make a connection to a newly formed connection
            if connection_count[i] == 1:
                
                # Find the index of the closest node with one connection (excluding itself)
                sorted_indices = np.argsort(dist_matrix[i])[1:]
                
                for closest_idx in sorted_indices[:dist]:

                    # Check if the index is not itelf and the connection count of the closest index is 1
                    if closest_idx != i and connection_count[closest_idx] <= 1 and connection_count[i] <= 1 and closest_idx not in connected_nodes.get(i, []) and i not in connected_nodes.get(closest_idx, []): 
                        
                        if dist_matrix[i][closest_idx] < dist:

                            cycle  = Testcycle(i, closest_idx ,connected_nodes, connections, connection_count)
                            
                            if not cycle:
                                UpdateConnections(i, closest_idx, connected_nodes, connections, connection_count)
                                break

    colormap = plt.cm.get_cmap('Dark2')
    color_cycle = itertools.cycle(colormap.colors)

    # Get indices where the value is 1
    single_nodes = np.where(connection_count == 1)[0]

    Tracks = []

    for i,node in enumerate(single_nodes):
        # Check that the track hasnt already been added
        if (check_start_end_exists(node,Tracks)):
            continue

        # Get the track path
        path = GetNodePath(connected_nodes, node, connected_nodes[node][0])

        total_length, total_energy = GetTrackLengthEnergy(path,data)
        color = next(color_cycle)

        Track = {"id":i, "start":path[0], "end":path[-1], "nodes":path, "length":total_length, "energy":total_energy,"label":"track","c":color}
        Tracks.append(Track)


    # for t in Tracks:
    #     print(t)

    # print(GetMeanNodeDist(Tracks, data))

    dist_threshold = 4*GetMeanNodeDist(Tracks, data)

    # ------------------------------------------------------
    # Here we break the track containing the vertex ID in two 
    for t in Tracks:

        # Found the track with the vertex
        if (vertex_index in t["nodes"]):
            # Get the length either side of track
            trk1_path = GetNodePath(connected_nodes, vertex_index, connected_nodes[vertex_index][0])[0:]
            trk2_path = GetNodePath(connected_nodes, vertex_index, connected_nodes[vertex_index][1])[0:]
            # print("vertex:",vertex_index)
            # print("Path1:",trk1_path)
            # print("Path2:",trk2_path)
            CreateVertexandSplit(vertex_index, t["id"], trk1_path, trk2_path, Tracks, data)
            break

    # print("Printing Tracks")
    # for t in Tracks:
    #     print(t)

    # ------------------------------------------------------

    UpdatedTracks = copy.deepcopy(Tracks)

    for idx, Track in enumerate(Tracks):
        curr_track = Track["id"]
        # print(curr_track)
        # print("Num new tracks:", len(UpdatedTracks))

        if (Track["label"] == "vertex"):
            # print("Skipping Vertex...")
            continue
        
        start_node = Track["start"]
        end_node   = Track["end"]

        # dont run this if we only got one track!
        if (len(Tracks) == 1):
            break

        # Get the indexes of closest nodes to start and end
        dist_ind_start = np.argsort(dist_matrix[start_node])[1:]
        dist_ind_end   = np.argsort(dist_matrix[end_node])[1:]

        # Filter nodes that are in the current track
        dist_ind_start = [x for x in dist_ind_start if x not in Track["nodes"]]
        dist_ind_end   = [x for x in dist_ind_end if x not in Track["nodes"]]

        # Distances of the end point to the closest track
        dist_start = dist_matrix[start_node][dist_ind_start[0]]
        dist_end   = dist_matrix[end_node][dist_ind_end[0]]

        # apply threshold
        if (dist_start > dist_threshold and dist_end > dist_threshold):
            # print("Failed distance requirements")
            continue

        # Initialize
        closest_idx = 0
        end_conn_node = 0
        con_point = "start"
        curr_track_path = Track["nodes"]

        # Choose the smallest index
        if dist_start < dist_end:
            closest_idx = dist_ind_start[0]
            end_conn_node = start_node
            
        else:
            closest_idx = dist_ind_end[0]
            end_conn_node = end_node
            con_point = "end"

        # Skip if we are trying to reconnect to a vertex
        if (end_conn_node == vertex_index):
            # print("Trying to connect to vertex, skipping...")
            continue


        # Get the track ID where the connecting node is located
        con_track      = GetTrackwithNode(closest_idx, Tracks)
        con_track_dict = GetTrackDictwithNode(closest_idx, Tracks)

        # if node-node then merge nodes and update track in Tracks
        if (closest_idx == con_track_dict["start"] or closest_idx == con_track_dict["end"]):
            
            newpath = join_tracks(curr_track_path,con_track_dict["nodes"])
            UpdateAndMergeTrack(curr_track, con_track, newpath, UpdatedTracks, data)
            UpdateConnections(closest_idx, end_conn_node, connected_nodes, connections, connection_count)
            print("node-node connection")
            continue

        # Check if the proposed connection will form a cycle
        cycle  = Testcycle(end_conn_node, closest_idx ,connected_nodes, connections, connection_count)
        
        if not cycle:

            if (con_point =="start"):
                curr_track_path.insert(0,closest_idx)
            else:
                curr_track_path.append(closest_idx)

            Track["nodes"] = curr_track_path
            UpdateConnections(closest_idx, end_conn_node, connected_nodes, connections, connection_count)
        else:
            break

        # Combine the tracks
        AddConnectedTracksnoDelta(curr_track, con_track, Track["nodes"], UpdatedTracks, data)



    # Add in any nodes without connections to the tracks as gammas and re-label other tracks as gammas
    AddConnectionlessNodes(connection_count, UpdatedTracks, data)
    FixTrackEnergies(UpdatedTracks,vertex_index, data)


    e_sum = 0
    track_id_list = []
    for t in UpdatedTracks:
        # print(t)
        if t["id"] not in track_id_list:
            track_id_list.append(t["id"])

            
        e_sum+=t["energy"]

    # print("Tot Energy: ",e_sum)
    # print("Total Tracks:", len(UpdatedTracks))
    # print("Total Unique:", len(track_id_list))


    # Reconstruction level quantities

    # Add gamma energy to the lowest energy track

    track1_energy = 0
    track2_energy = 0

    track1 = 0
    track2 = 0

    for t in UpdatedTracks:

        if (t["label"] == "Track1"):
            track1_energy = track1_energy + t["energy"]
            track1=track1+1

        if (t["label"] == "Track2"):
            track2_energy = track2_energy + t["energy"]
            track2=track2+1

    # If there is no Track1 or Track2 then return failed reco value
    if (track1 == 0 or track2 == 0):
        print("Error!! No track1 or track 2 in final reconstruction...")
        return -1, -1, -1, -1, -1, connected_nodes, UpdatedTracks

    # ID the track with the lowest label
    highest_track = "Track2"
    if (track1_energy > track2_energy):
        highest_track = "Track1"


    e_gammas = 0
    for t in UpdatedTracks:

        if (t["label"] != "Track1" and t["label"] != "Track2"):
            # print(t["label"])
            e_gammas = e_gammas + t["energy"]

    print("Gamma Energy:",e_gammas)

    if (highest_track == "Track1"):
        track1_energy = track1_energy + e_gammas
    else:
        track2_energy = track2_energy + e_gammas

    print("Track 1 Energy:", track1_energy)
    print("Track 2 Energy:", track2_energy)
    print("Tot Energy:", track1_energy + track2_energy)
    Reco_T1 = ReturnLargest(track1_energy,track2_energy)


    part = part[part.primary == 1]
    Gen_T1 = max(part.kin_energy.values) # generator T1


    ### ------------- ------------- ------------- -------------
    # Filter particles with particle_id 1 and 2
    particle_1 = part[part.particle_id == 1].copy()
    particle_2 = part[part.particle_id == 2].copy()

    # Merge the DataFrames on 'event' to pair particles from the same event
    merged_particles = pd.merge(particle_1, particle_2, on="event_id", suffixes=('_1', '_2'))

    # Apply the function to each row in the merged DataFrame
    merged_particles['angle'] = merged_particles.apply(calculate_angle_parts, axis=1)
    Gen_cos_theta =  merged_particles['angle'].iloc[0]
    # ------------- ------------- ------------- -------------

    data['Track1'] = 0
    data['Track2'] = 0

    track1_indices = []
    track2_indices = []

    for t in UpdatedTracks:

        if t["label"] == "Track1":
            track1_indices = track1_indices + t["nodes"]

        if t["label"] == "Track2":
            track2_indices = track2_indices + t["nodes"]

    data.loc[track1_indices, 'Track1'] = 1
    data.loc[track2_indices, 'Track2'] = 1


    # Given vertex position
    vertex = data[ (data.Track1 == 1) & (data.Track2 == 1)]
    vertex = np.array([vertex.iloc[0].x,vertex.iloc[0].y,vertex.iloc[0].z])

    Track1 = data[ (data.Track1 == 1) & (data.Track2 != 1)]
    Track1 = Track1.reindex(track1_indices)


    Track2 = data[ (data.Track2 == 1) & (data.Track1 != 1)]
    Track2 = Track2.reindex(track2_indices)

    Track1_node1 = np.array(Track1.iloc[1][0:3])
    Track2_node1 = np.array(Track2.iloc[1][0:3])

    direction_vector1 = Track1_node1 - vertex
    direction_vector2 = Track2_node1 - vertex

    # # Compute cosine of the angle between the vectors
    Reco_cos_theta = cosine_angle(direction_vector1, direction_vector2)

   

    return Gen_T1, Gen_cos_theta, Reco_T1, Reco_cos_theta, e_gammas, connected_nodes, UpdatedTracks


# hits = pd.read_hdf('../../NEXT_nudobe/files/data/Leptoquark_SM_nexus.h5',"MC/hits")
hits = pd.read_hdf('../../NEXT_nudobe/files/data/Leptoquark_SM_1mm_smear.h5',"MC/hits")
parts = pd.read_hdf('../../NEXT_nudobe/files/data/Leptoquark_SM_1mm_smear.h5',"MC/particles")

events = [1834, 1835, 1836, 1837, 1838, 1839, 1840]

counter = 0

for event_num in parts.event_id.unique():

    if (counter > 30):
        break

    hit = hits[hits.event_id == event_num]
    part = parts[parts.event_id == event_num]

    # print(hit)

    Gen_T1, Gen_cos_theta, Reco_T1, Reco_cos_theta, e_gammas, connected_nodes, UpdatedTracks = RunReco(hit, part)

    print("Event: ",event_num)
    print("Gen  T1:",Gen_T1)
    print("Reco T1:", Reco_T1)

    print("Gen  Cos Theta:",Gen_cos_theta)
    print("Reco Cos Theta:", Reco_cos_theta)
    print("")

    counter = counter+1

