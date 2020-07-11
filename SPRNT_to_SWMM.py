import random
import numpy as np
import os.path
import re

"""
This script constitutes the write step of the read-process-write workflow for converting SPRNT data into SWMM inputfile
format.

The database result from the processing script is read into variables that will be written into a text file configured 
to mimic the layout of a SWMM .inp file
"""
# Globals
sprnt_file = r"Small_Network_for_SWMM_Test.spt"
swmm_example = r"Lavaca_Sample.inp"
swmm_final = r"Small_Network_Complete.inp"
path_to_external_files = 'C:/Users/edt489/Box Sync/Current Projects/SWMM_SPRNT/Conversion_Scripts/LavacaRB/External_Files'
root_raw = '5707298_630030057'
root_clip = root_raw.split('_')[-1]
root_bc = 0.0

max_segments_per_comid_digits = 4
cross_section_transformation_tolerance = 0.001
max_inflows = 500

f_sprnt = open(sprnt_file, "r")
sprnt_contents = f_sprnt.readlines()


def read_sprnt_for_counters():
    """
    This function inspects the contents of the SPRNT input file and creates arrays for each of the data types that are
    sized correctly.  This allows us to use np.arrays rather than 2D lists. The second dimension of the arrays is
    hard-coded to reflect that each data type in SPRNT contains a different number of parameters.

    In the San Antonio - Guadalupe River basin the SPRNT input file is 3.8 millions lines.
    """
    global total_count, null_count, node_count, segment_count, lateral_count, junction_count, qsource_count, \
        nodes_name, nodes_geo, segments, lateralsources, junctions, qsources, boundaryconditions, sprnt_contents

    total_count = len(sprnt_contents)
    null_count, node_count, segment_count, lateral_count, junction_count, qsource_count = 0, 0, 0, 0, 0, 0

    for ii in range(total_count):
        if sprnt_contents[ii].find("#") != -1:
            null_count = null_count + 1
        elif sprnt_contents[ii].find('node') != -1:
            node_count = node_count + 1
        elif sprnt_contents[ii].find("segment") != -1:
            segment_count = segment_count + 1
        elif sprnt_contents[ii].find("lateralsource") != -1:
            lateral_count = lateral_count + 1
        elif sprnt_contents[ii].find("junction") != -1:
            junction_count = junction_count + 1
        elif sprnt_contents[ii].find("qsource") != -1:
            qsource_count = qsource_count + 1

    # The sizes of these arrays is network dependent, may need to write some functions to generalize them
    nodes_geo = np.empty((node_count, 2), dtype=object)
    nodes_name = np.empty((node_count, 6), dtype=object)
    segments = np.empty((segment_count, 3), dtype=object)
    lateralsources = np.empty((lateral_count, 1491), dtype=object)
    junctions = np.empty((junction_count, 5), dtype=object)
    qsources = np.empty((qsource_count, 1491), dtype=object)
    # The boundary condition is hard-coded with the value obtained from the SPRNT input file
    boundaryconditions = (root_raw, root_bc)

    return


read_sprnt_for_counters()


def ay_hec2_transform():
    global station, elev
    station = [0.0]
    elev = []
    for yy in Y:
        elev.append(yy)
    elev.append(0.0)
    for yy in range(1, len(Y)):
        elev.append(Y[-1 - yy])
    for nn in range(1, len(W)):
        offset = (W[nn - 1] - W[nn]) / 2
        station.append(station[nn - 1] + offset)
    reversed_W = list(reversed(W))
    station.append(station[-1] + W[-1])
    for ww in range(1, len(reversed_W)):
        offset = (reversed_W[ww] - reversed_W[ww - 1]) / 2
        station.append(station[-1] + offset)

    for stat in range(1, len(station)):
        if station[stat] < station[stat - 1]:
            print("Station needs to be monotonically increasing")
            print(station[stat], station[stat - 1])
            quit()

    return elev, station


def transform_check(elev, station, A_max, P_max):
    """This function only needs to be called when a new Network is being analyzed."""
    global minimum_neg, node_counter
    station_cap = [station[0], station[-1]]
    elev_cap = [elev[0], elev[-1]]
    area_channel = np.trapz(elev, station)
    area_cap = np.trapz(elev_cap, station_cap)
    area_difference = abs((A_max - area_cap + area_channel)) / A_max
    if abs(area_difference) > cross_section_transformation_tolerance:
        print("The overall area has not been conserved by ", area_difference, ' percent')
        quit()

    wetted_perimeter = 0
    for ii in range(len(station) - 1):
        x_diff = abs(station[ii + 1] - station[ii])
        y_diff = abs(elev[ii + 1] - elev[ii])
        inc_length = np.sqrt(x_diff ** 2 + y_diff ** 2)
        wetted_perimeter = wetted_perimeter + inc_length

    perimeter_diff = abs((wetted_perimeter - P_max)) / P_max
    if abs(perimeter_diff) > cross_section_transformation_tolerance:
        print("The wetted perimeter has not been conserved by ", perimeter_diff, ' percent')
        quit()

    return


def ay_main():
    global total_count, A, P, Y, W
    node_counter = 0

    for ii in range(total_count):
        if sprnt_contents[ii].find("intrinsic") != -1:
            A = []
            P = []
            Y = []
            W = []
            jj = ii + 1
            while sprnt_contents[jj].find("end") == -1:
                line = sprnt_contents[jj].split()
                for elem in range(len(line)):
                    line[elem] = line[elem].split('=')[1]
                # print(line)
                A.insert(0, float('{:.3f}'.format(float(line[0]))))
                P.insert(0, float('{:.3f}'.format(float(line[1]))))
                Y.insert(0, float('{:.3f}'.format(float(line[2]))))
                W.insert(0, float('{:.3f}'.format(float(line[3]))))
                jj = jj + 1

            for ww in range(1, len(W)):
                if W[ww - 1] < W[ww]:
                    W[ww] = W[ww - 1]

            nodes_geo[node_counter] = ay_hec2_transform()
            transform_check(nodes_geo[node_counter, 0], nodes_geo[node_counter, 1], A[0], P[0])
            node_counter = node_counter + 1

    return


ay_main()


def read_sprnt_contents():
    """
    This function populates the various data type arrays with the values parsed from the SPRNT input file database.

    The function loops through each line in the database; when it finds the keyword that indicates a certain
    data-type, the corresponding values are added to the np.array
    """

    global nodes_geo, nodes_name, segments, lateralsources, junctions, qsources, sprnt_contents, root_clip, \
        node_elevation, upstream, length, id, node_elevation

    # node parameters
    id = 0
    sR = 1
    n = 2
    zR = 3
    hR = 4
    # bw and slope are trapezoidal approximation parameters
    # bw = 5
    # slope = 6
    node_elevation = 5

    # segment parameters
    upstresam = 0
    downstream = 1
    length = 2

    # lateralsource parameters
    source = 0
    timeunit = 1

    # junction parameters
    downstream_node = 0
    upstream1 = 1
    coeff1 = 2
    upstream2 = 3
    coeff2 = 4

    # qsource parameters
    source = 0
    timeunit = 1

    node_number, segment_number, lateralsource_number, junction_number, qsource_number = 0, 0, 0, 0, 0
    for line in range(len(sprnt_contents) - 1):
        # print(contents[line])

        # If the line is commented out, ignore it
        if sprnt_contents[line].find("#") != -1:
            continue

        # When you find a node definition - do this
        # Need to replace the nodes population with the AY-transform
        elif sprnt_contents[line].find("node") != -1:
            print("Found Node")
            spline = re.split(' id=| sR=| n=| zR=| hR=|\n', sprnt_contents[line])
            nodes_name[node_number, id] = spline[1]
            nodes_name[node_number, sR] = spline[2]
            nodes_name[node_number, n] = spline[3]
            nodes_name[node_number, zR] = spline[4]
            nodes_name[node_number, hR] = spline[5]

            # nextspline = re.split('bottomwidth=|slope=|end|\n', sprnt_contents[line+1])
            # nodes_name[node_number, bw] = nextspline[1]
            # nodes_name[node_number, slope] = nextspline[2]
            # #print(nextspline)
            node_number = node_number + 1

        # When a segment definition is found - do this
        elif sprnt_contents[line].find("segment") != -1:
            print("Found Segment")
            spline = re.split(' up=| down=| length=| ', sprnt_contents[line])
            segments[segment_number] = spline[2:5]
            segment_number = segment_number + 1

        # When a lateralsource is defined - add the time series
        elif sprnt_contents[line].find("lateralsource") != -1:
            print("Found Lateral Source")
            spline = re.split("=|\n", sprnt_contents[line + 1])
            lateralsources[lateralsource_number, source] = spline[1]
            nextspline = re.split("=|\n", sprnt_contents[line + 3])
            lateralsources[lateralsource_number, timeunit] = nextspline[1]
            temp_line = line + 4
            time_number = 2
            # Each time series has the source, timeunit, time_number series.  So the time_number series starts at 2
            while sprnt_contents[temp_line].find('t=') != -1:
                timespline = re.split("=|\n| ", sprnt_contents[temp_line])
                # print(timespline)
                lateralsources[lateralsource_number, time_number] = float(timespline[-2])
                temp_line = temp_line + 1
                time_number = time_number + 1
            lateralsource_number = lateralsource_number + 1

        # When a junction is found - add the connecting nodes and coefficients
        elif sprnt_contents[line].find("junction") != -1:
            print("Found Junction")
            spline = re.split("down=| up1=| coeff1=| up2=| coeff2=|\n", sprnt_contents[line + 1])
            spline = spline[1:-1]
            try:
                junctions[junction_number] = spline
            except ValueError:
                junctions[junction_number, downstream_node:coeff1 + 1] = spline
            junction_number = junction_number + 1

        # When a qsource is found - add the time series
        elif sprnt_contents[line].find("qsource") != -1:
            print("Found Qsource")
            spline = re.split("=|\n", sprnt_contents[line + 1])
            qsources[qsource_number, source] = spline[1]
            nextspline = re.split("=|\n", sprnt_contents[line + 3])
            qsources[qsource_number, timeunit] = nextspline[1]
            temp_line = line + 4
            time_number = 2
            # Similarly to lateral sources, the qsource row has source, timeunit, time_number series data
            while sprnt_contents[temp_line].find('t=') != -1:
                timespline = re.split("=|\n| ", sprnt_contents[temp_line])
                qsources[qsource_number, time_number] = float(timespline[-2])
                temp_line = temp_line + 1
                time_number = time_number + 1

            qsource_number = qsource_number + 1

    # This code segment was written to allow Eric to test the upstream boundary condition inflows
    # qsource_names = []
    # for ii in qsources[:, source]:
    #     qsource_names.append(ii)
    # print(qsource_names)
    # with open('qsources.csv', 'w') as result_file:
    #     for ii in qsource_names:
    #         result_file.write(ii)
    #         result_file.write('\n')

    return


read_sprnt_contents()


def read_SWMM_contents(inputfilename=swmm_example):
    """
    The general format of a SWMM input file is adopted from a sample file which has been independently populated
    with the Outfall node
    """
    global swmm_contents, count

    with open(inputfilename, 'r') as swmmput:
        swmm_contents = swmmput.readlines()

    count = len(swmm_contents)
    return


read_SWMM_contents()


def create_inputfile(trialfile):
    """
    This function combines the processed data from the np.arrays (organized into the data types from SPRNT) and inserts
    the data line by line into the SWMM inputfile template

    When a header is found in the swmm_contents, the lines beneath that header are populated with the appropriate
    information from the np.arrays. All of the tokens for which no SPRNT corollary exists are be populated with a
    default value.

    **Where ambiguous, values for conduit/xsection dimensions are defined by the values at the upstream SPRNT node.
    """
    global swmm_contents, nodes_name, nodes_geo, segments, lateralsources, junctions, qsources, comID, count,\
        lateralseries_names

    """
    This first set of for-loops fix the issue that the same "node" in SPRNT is assigned two different names depending on
    which comID it appears.  Removing the comID prefix makes it so that the same location in space only has one name.
    
    This also has the effect of eliminating the need for the Junctions definitions from SPRNT, as the Junctions simply
    acted to connect the same location in space with each identifier.
    
    A similar ID trimming must be done wherever nodeID's exist, namely, Segments and Qsources.
    """
    for node in range(len(nodes_name)):
        if nodes_name[node, 0] not in segments[:, 0]:
            print(nodes_name[node, 0])
            nodes_name[node, 0] = '-998877'

    for segment in range(len(segments)):
        from_node_id = segments[segment, 0].split('_')[-1]
        to_node_id = segments[segment, 1].split('_')[-1]
        if len(from_node_id) > max_segments_per_comid_digits:
            segments[segment, 0] = from_node_id
        if len(to_node_id) > max_segments_per_comid_digits:
            segments[segment, 1] = to_node_id

    for node in range(len(nodes_name)):

        node_id = nodes_name[node, 0].split('_')[-1]
        # This if statement says that if the node_id is very long, it probably corresponds to a node connecting segments
        # which have different names depending on which segment they're on.
        if len(node_id) > max_segments_per_comid_digits:
            nodes_name[node, 0] = node_id

    node_id_list = []

    for source in range(len(qsources)):
        q_id = qsources[source, 0].split('_')[-1]
        if len(q_id) > max_segments_per_comid_digits:
            qsources[source, 0] = q_id

    for source in range(len(lateralsources)):
        l_id = lateralsources[source, 0].split('_')[-1]
        if len(l_id) > max_segments_per_comid_digits:
            lateralsources[source, 0] = l_id


    # SWMM Junctions = SPRNT Nodes
    junction_index = swmm_contents.index('[JUNCTIONS]\n')
    # Iterates backwards so they show up in a familiar order in the .inp file
    for node in nodes_name[::-1]:
        # This logical checks the node name, to ensure each node is only added once (i.e. the comID endpoints)
        if (node[0] not in node_id_list) and (node[0] not in root_clip) and (node[0] != '-998877'):
            # The elevation is converted from centimeters to meters
            try:
                templine = [node[0], '{:.3f}'.format(float(node[3])/100), '0', '0', '0', '0']
            except TypeError:
                print(node)
                continue
            splitline = "      ".join(templine) + "\n"
            swmm_contents.insert(junction_index+3, splitline)
        node_id_list.append(node[0])

        """
        COORDINATES
        SWMM requires coordinates for each node bc it functions through a GUI.  SPRNT, as an api, doesn't have the same
        requirement, so coordinate data does not exist for the Texas river basins
        
        -for each node, randomly select values between 0-1000 for x and y
        
        **The index of each keyword must be updated before insertion to reflect the size change due to the previous
        section's insertions
        """
        coordinates_index = swmm_contents.index('[COORDINATES]\n')
        x_coor = random.randrange(0, 1000)
        y_coor = random.randrange(0, 1000)
        templine = [node[0], str(x_coor), str(y_coor)]
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(coordinates_index + 3, splitline)

    """
    CONDUITS
    -The name of the conduit can be from an arbitrary counter
    -The from_node, to_node, and length are from the segments array
    -The roughness will be assigned by the roughness in the nodes array with the same id as the segment's upstream node
    """
    conduit_name = 1
    conduits_index = swmm_contents.index('[CONDUITS]\n')

    for segment in segments[::-1]:
        from_node = segment[0]
        to_node = segment[1]
        length = segment[2]
        # This variable designates the row of the nodes np.array where the from_node can be found
        node_loc = np.where(nodes_name[:, 0] == from_node)
        mannings_n = nodes_name[node_loc[0][0], 2]
        templine = [str(conduit_name), from_node, to_node, length, str(mannings_n), '0', '0', '', '']
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(conduits_index + 3, splitline)

        """    
        XSECTIONS
        -The name of the XSECTION will be the same as the arbitrary counter for conduits.  Meaning they're populated at
        the same time
        
        For Lavaca, the shape is the transformed HEC2 irregular, and the Tsect name needs to be unique to the conduit
        - Tsect name will be the upstream_node
        # -The shape is TRAPEZOIDAL
        # -Geom1 is the max depth, also arbitrarily defined (probably 100)
        # -Geom2 is the bottom width, which is defined at the upstream node for each segment
        # -Geom3 is the side slope, also defined at the upstream node for each segment
        """
        xsections_index = swmm_contents.index('[XSECTIONS]\n')
        #bottom_width = nodes[mannings_n_loc[0][0], 5]
        #side_slope = nodes[mannings_n_loc[0][0], 6]
        templine = [str(conduit_name), 'IRREGULAR', from_node]
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(xsections_index + 3, splitline)

        conduit_name = conduit_name + 1

    """
    INFLOWS
    -The node name will be the node where the qsource/lateralsource is located
    -Constituent is 'FLOW'
    -Time Series is variable bc there are multiple time series
    -Type 'FLOW'
    -Mfactor, Sfactor = 1.0
    """
    inflows_index = swmm_contents.index('[INFLOWS]\n')

    qseries_names = []
    for source in qsources:
        name = source[0]
        # In Lavaca basin, about half have a constant inflow of 1m3/s, so these can be combined into one timeseries
        if all(source[2:] == 1.0):
            templine = [name, 'FLOW', 'T1', 'FLOW', '1.0', '1.0', '', '']
        else:
            templine = [name, 'Flow', name, 'FLOW', '1.0', '1.0', '', '']
            qseries_names.append(name)
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(inflows_index + 3, splitline)

    lateralseries_names = []
    for source in lateralsources:
        name = source[0]
        # In Lavaca basin, about half have a constant inflow of 1m3/s, so these can be combined into one timeseries
        if all(source[2:] == 1.0):
            templine = [name, 'FLOW', 'T1', 'FLOW', '1.0',  '', '']
            # print(name)
        else:
            templine = [name, 'FLOW', name, 'FLOW', '1.0', '', '']
            lateralseries_names.append(name)
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(inflows_index + 3, splitline)

    """
    [TIMESERIES]
    - need to include flow timeseries from lateral sources and qsources
    - we have checked that there is no overlap
    - each timeseries is named after the node at which it occurs
    - on a given line: the name, (time since simulation began, flow value) pairs
    - The lateral sources make the size of the input file explode
    - We also include the T1 series which describes a subset of qsource and lateralsource data
    """
    timeseries_index = swmm_contents.index('[TIMESERIES]\n')

    # First I should check if the lateralsources and qsources have any overlap

    for name in lateralsources[:, 0]:
        if name in qsources[:, 0]:
            print("There are multiple inflows to node " + name)

    #This is where I need to put a filter for the max_inflows
    for name in lateralseries_names:
        # file = str(name) + '.dat'
        # completeName = os.path.join(path_to_external_files, file)
        # templine = [name, 'FILE', file]
        # splitline = "      ".join(templine) + "\n"
        # swmm_contents.insert(timeseries_index + 3, splitline)

        templine = [name]
        lateral_data_loc = np.where(lateralsources[:, 0] == name)
        # For each time, the range is hardcoded as the number of hours of the simulation
        for time in range(1488):
            templine.append(str(time))
            templine.append(str(lateralsources[lateral_data_loc[0][0], time + 2]))
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(timeseries_index + 3, splitline)

    for name in qseries_names:
        templine = [name]
        qsource_data_loc = np.where(qsources[:, 0] == name)
        # For each time, the range is hardcoded as the number of hours of the simulation
        for time in range(1488):
            templine.append(str(time))
            if qsources[qsource_data_loc[0][0], time + 2] == None:
                templine.append('0.0')
            else:
                templine.append(str(qsources[qsource_data_loc[0][0], time + 2]))
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(timeseries_index + 3, splitline)

    templine = ['T1']
    for time in range(1488):
        templine.append(str(time))
        templine.append('1.0')
    splitline = "      ".join(templine) + "\n"
    swmm_contents.insert(timeseries_index + 3, splitline)

    """
    [TRANSECTS] - this section is going to be a bitch as well
    
    """
    transect_index = swmm_contents.index('[TRANSECTS]\n')

    error_count = 0
    for segment in segments[::-1]:
        from_node = segment[0]
        node_loc = np.where(nodes_name[:, 0] == from_node)
        mannings_n = nodes_name[node_loc[0][0], 2]
        templine_NC = ['NC', mannings_n, mannings_n, mannings_n]
        station = nodes_geo[node_loc[0][0], 1]
        elevation = nodes_geo[node_loc[0][0], 0]
        stat_past = 0
        for stat in station:
            if stat >= stat_past:
                stat_past = stat
            else:
                error_count = error_count + 1
                break
        templine_X1 = ['X1', from_node, str(len(station)), str(station[0]), str(station[-1]), '0.0', '0.0', '0.0', '0.0', '0.0']
        templine_GR = ['GR']
        for index in range(len(station)):
            templine_GR.append(str(elevation[index]))
            templine_GR.append(str(station[index]))
        splitline_NC = "      ".join(templine_NC) + "\n"
        splitline_X1 = "      ".join(templine_X1) + "\n"
        splitline_GR = "      ".join(templine_GR) + "\n"
        swmm_contents.insert(transect_index + 3, ';\n')
        swmm_contents.insert(transect_index + 3, splitline_GR)
        swmm_contents.insert(transect_index + 3, splitline_X1)
        swmm_contents.insert(transect_index + 3, splitline_NC)
    print(error_count)

    # Once the SWMM .inp template has been populated, open a new file and write the contents to it
    with open(trialfile, 'w') as newfile:

        for i in range(len(swmm_contents)):
            newfile.write(swmm_contents[i])

    return


# The function argument specifies the name of the SWMM input file that is created
create_inputfile(swmm_final)

def create_externalfiles(extfile, toFile):
    file = str(extfile) + ".dat"
    completeName = os.path.join(path_to_external_files, file)
    with open(completeName, "w") as file1:
        for i in range(len(toFile)):
            file1.write(toFile[i])
    return

for name in lateralseries_names:
    toFile = []
    lateral_data_loc = np.where(lateralsources[:, 0] == name)
    # For each time, the range is hardcoded as the number of hours of the simulation
    for time in range(1488):
        templine = []
        templine.append(str(time))
        templine.append(str(lateralsources[lateral_data_loc[0][0], time + 2]))
        splitline = "      ".join(templine) + "\n"
        toFile.append(splitline)
    create_externalfiles(name, toFile)

