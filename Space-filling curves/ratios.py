import math
from itertools import combinations
import copy
import time
import pickle
import csv
import multiprocessing as mp
import argparse

parser = argparse.ArgumentParser(description="Maximum depth level")
parser.add_argument('max_depth', type=int)

args = parser.parse_args()

pattern_dict = {}   # number > pattern > number > pattern
id_count = 1
filename = './data/ratios.csv'


def get_curves(prev_curve_list):
    global id_count
    curves = []
    for curve in prev_curve_list:
        # r = no. of refinements to make
        for r in range(1, len(curve.refinement_left_list) + len(curve.refinement_right_list) + 1):
            # find all refinement location combinations for n choose r
            refinements, duplicates = get_refinements(
                curve.refinement_left_list, curve.refinement_right_list, r, curve.duplication)

            # parallel routine
            pool = mp.Pool(mp.cpu_count())
            # create a curve for each refinement
            curves.extend([pool.apply(create_curves, args=(id_count + i, curve.production_list, curve.depths_list,
                                                           load_obj(curve.curve_id),
                                                           refinements[i], curve.midpoint, duplicates[i])) for i in
                           range(0, len(refinements))])
            pool.close()
            id_count += len(refinements)
            # non-parallel routine:
            # for i in range(0, len(refinements)):
            #    curves.append(create_curves(id_count + i, curve.production_list, curve.depths_list,
            #                                                 load_obj(curve.curve_id),
            #                                                refinements[i], curve.midpoint, duplicates[i]))
            # id_count += len(refinements)

            print(id_count)

    return curves


def create_curves(curve_id, production_list, depths_list, cell_neighbours, refinement_list, midpoint, duplicates):
    return Curve(curve_id, production_list, depths_list, cell_neighbours, refinement_list, midpoint, duplicates)


def get_refinements(refinement_left_list, refinement_right_list, no_of_refinements, is_prev_symmetric):
    refinement_list = []
    duplicates_list = []

    if is_prev_symmetric == 1:  # if prev quad-tree is not symmetric, then new quad-tree cannot be symmetric
        left_refinement_max = math.floor(no_of_refinements / 2)
        for left_refinements in range(0, left_refinement_max+1):    # exclude symmetry
            for combo_r in combinations(refinement_right_list, no_of_refinements - left_refinements):
                if left_refinements > 0:
                    # if left refinements = right refinements
                    if left_refinements == no_of_refinements - left_refinements:
                        for combo_l in combinations(refinement_left_list, left_refinements):
                            symmetric = True
                            for i in range(0, left_refinements):
                                left = refinement_left_list.index(combo_l[i])
                                right = refinement_right_list.index(combo_r[len(combo_r) - i - 1])
                                # if left and right locs dont mirror then no refinement
                                if not len(refinement_right_list) - right - 1 == left:
                                    symmetric = False
                                    break
                            refinement_list.append(list(combo_r) + list(combo_l))
                            if symmetric:
                                duplicates_list.append(1)
                                break       # all left combinations complete for given right combination.
                                # carrying on would need to duplicates
                            else:
                                duplicates_list.append(2)
                    else:
                        for combo_l in combinations(refinement_left_list, left_refinements):
                            refinement_list.append(list(combo_r) + list(combo_l))
                            duplicates_list.append(2)

                else:
                    refinement_list.append(list(combo_r))
                    duplicates_list.append(2)
    else:
        for combo in combinations(refinement_right_list+refinement_left_list, no_of_refinements):
            refinement_list.append(list(combo))
            duplicates_list.append(2)

    return [refinement_list, duplicates_list]        # return all new bitstreams found


def load_obj(curve_id):
    with open('obj/' + str(curve_id) + '.pkl', 'rb') as f:
        return pickle.load(f)


# boundaries([left,up,right,down])
# depth 1 = Curve({0: [None, 1, 3, None], 1: [None, None, 2, 3], 2: [1, None, None, 3], 3: [1, 2, None, None]})
class Curve:
    def __init__(self, curve_id, prev_production_list, prev_depths_list, prev_cell_neighbours, refinement,
                 prev_midpoint, duplicates):
        self.curve_id = curve_id
        self.duplication = duplicates
        self.production_list, self.depths_list = self.apply_production_rules(
            copy.deepcopy(prev_production_list), copy.deepcopy(prev_depths_list), refinement)
        cell_neighbours, self.refinement_left_list, self.refinement_right_list, self.midpoint = \
            self.get_cell_neighbours(prev_cell_neighbours, refinement, prev_midpoint)
        self.quadtree_averages = self.calculate_average(cell_neighbours)
        self.save_obj(cell_neighbours, curve_id)
        self.cell_count = len(cell_neighbours)
        self.max_depth = max(self.depths_list)
        with open(filename, "a", newline='') as csv_file:
            my_file = csv.writer(csv_file, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            my_file.writerow([self.max_depth, str(self.quadtree_averages[0]), str(self.cell_count)])
            if self.duplication == 2:
                my_file.writerow([self.max_depth, str(self.quadtree_averages[0]), str(self.cell_count)])

    def save_obj(self, obj, curve_id):
        with open('obj/' + str(curve_id) + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    def apply_production_rules(self, production_list, depths_list, refinement_list):
        refinement_letters_list = []
        for i in range(0, len(refinement_list)):
            if production_list[refinement_list[i]] == 'h':
                refinement_letters_list.append('h')
                production_list[refinement_list[i]:refinement_list[i] + 1] = ['a', 'h', 'h', 'b']
            elif production_list[refinement_list[i]] == 'a':
                refinement_letters_list.append('a')
                production_list[refinement_list[i]:refinement_list[i] + 1] = ['h', 'a', 'a', 'c']
            elif production_list[refinement_list[i]] == 'b':
                refinement_letters_list.append('b')
                production_list[refinement_list[i]:refinement_list[i] + 1] = ['c', 'b', 'b', 'h']
            else:
                refinement_letters_list.append('c')
                production_list[refinement_list[i]:refinement_list[i] + 1] = ['b', 'c', 'c', 'a']
            depths_list[refinement_list[i]:refinement_list[i] + 1] = [depths_list[refinement_list[i]] + 1 for _ in
                                                                      range(4)]
        return production_list, depths_list

    def get_cell_neighbours(self, prev_cell_neighbours, refinement_loc, prev_midpoint):
        midpoint = copy.deepcopy(prev_midpoint)
        cell_neighbours = copy.deepcopy(prev_cell_neighbours)
        shifted_r_loc = copy.deepcopy(refinement_loc)
        replaced_neighbours = {}
        shift = 0
        count = 0
        for r in range(len(refinement_loc) - 1, -1, -1):
            shifted_r_loc[r] += shift  # shift refinement_list to new positions
            shift += 3
            if refinement_loc[r] < midpoint:  # update midpoint location
                count += 3
        midpoint += count
        for neighbour in cell_neighbours.values():
            for i in range(0, 4):  # shift values of cell_neighbours where necessary
                if neighbour[i]:
                    for j in range(0, len(neighbour[i])):
                        shift = 0
                        for r in reversed(refinement_loc):
                            if neighbour[i][j] > r:
                                shift += 3
                            else:
                                break
                        neighbour[i][j] += shift
        for i in range(len(cell_neighbours), len(cell_neighbours) + 3 * len(refinement_loc)):
            cell_neighbours[i] = [[], [], [], []]
        r = 0
        neighbour = len(prev_cell_neighbours) - 1
        for i in range(len(cell_neighbours) - 1, -1, -1):
            if shifted_r_loc[r] <= i < shifted_r_loc[r] + 4:
                if shifted_r_loc[r] == i:
                    replaced_neighbours[refinement_loc[r]] = copy.deepcopy(cell_neighbours[neighbour])
                    neighbour -= 1
                    if r < len(refinement_loc) - 1:
                        r += 1
                cell_neighbours[i] = [[], [], [], []]

            else:
                for d in range(0, 4):
                    if cell_neighbours[neighbour][d]:
                        n_list = []
                        for edge in range(len(cell_neighbours[neighbour][d])):
                            if cell_neighbours[neighbour][d][edge] not in shifted_r_loc:
                                n_list.append(cell_neighbours[neighbour][d][edge])
                        cell_neighbours[i][d] = n_list
                neighbour -= 1

        future_r_left_list = []
        future_r_right_list = []
        shift = 0
        for r in reversed(refinement_loc):
            top_right, top_left, bottom_left, bottom_right = self.get_shape(r + shift)
            # add the two newly refined cell neighbours for each of the four newly refined cells.
            cell_neighbours[r + shift + top_left][2].append(r + shift + top_right)
            cell_neighbours[r + shift + top_left][3].append(r + shift + bottom_left)
            cell_neighbours[r + shift + top_right][0].append(r + shift + top_left)
            cell_neighbours[r + shift + top_right][3].append(r + shift + bottom_right)
            cell_neighbours[r + shift + bottom_left][1].append(r + shift + top_left)
            cell_neighbours[r + shift + bottom_left][2].append(r + shift + bottom_right)
            cell_neighbours[r + shift + bottom_right][0].append(r + shift + bottom_left)
            cell_neighbours[r + shift + bottom_right][1].append(r + shift + top_right)

            for i in range(0, 4):
                if prev_cell_neighbours[r][i]:
                    if i == 0:      # left
                        pos1 = top_left
                        pos2 = bottom_left
                    elif i == 1:    # up
                        pos1 = top_left
                        pos2 = top_right
                    elif i == 2:    # right
                        pos1 = top_right
                        pos2 = bottom_right
                    else:           # down
                        pos1 = bottom_left
                        pos2 = bottom_right
                    if prev_cell_neighbours[r][i][0] not in refinement_loc:
                        self.update_neighbours(pos1, pos2, i, r, shift, cell_neighbours, replaced_neighbours)
                    else:
                        if r > prev_cell_neighbours[r][i][0]:
                            r2 = refinement_loc.index(prev_cell_neighbours[r][i][0])
                            top_right2, top_left2, bottom_left2, bottom_right2 = self.get_shape(shifted_r_loc[r2])
                            top_right_cell = prev_cell_neighbours[r][i][0] + (
                                    len(refinement_loc) - r2 - 1) * 3 + top_right2
                            top_left_cell = prev_cell_neighbours[r][i][0] + (
                                    len(refinement_loc) - r2 - 1) * 3 + top_left2
                            bottom_left_cell = prev_cell_neighbours[r][i][0] + (
                                    len(refinement_loc) - r2 - 1) * 3 + bottom_left2
                            bottom_right_cell = prev_cell_neighbours[r][i][0] + (
                                    len(refinement_loc) - r2 - 1) * 3 + bottom_right2
                            if i == 0:
                                pos3 = top_right_cell
                                pos4 = bottom_right_cell
                            elif i == 1:
                                pos3 = bottom_left_cell
                                pos4 = bottom_right_cell
                            elif i == 2:
                                pos3 = top_left_cell
                                pos4 = bottom_left_cell
                            else:
                                pos3 = top_left_cell
                                pos4 = top_right_cell
                            if i < 2:
                                j = i + 2
                            else:
                                j = i - 2
                            self.update_neighbours_v2(pos1, pos2, pos3, pos4, i, j, r, shift, cell_neighbours)

            refinable = [True, True, True, True]  # [left, top, right, bottom]
            for i in range(0, 4):
                if prev_cell_neighbours[r][i]:
                    if len(prev_cell_neighbours[r][i]) == 1:    # neighbouring cell must also be refined
                        if not prev_cell_neighbours[r][i][0] in refinement_loc:
                            refinable[i] = False
                    else:                                       # if one edge has more than one neighbour
                        refinable[i] = False
            future_refinement_list = []
            if refinable[0] and refinable[1]:  # left and top
                future_refinement_list.append(r + shift + top_left)
            if refinable[0] and refinable[3]:  # left and bottom
                future_refinement_list.append(r + shift + bottom_left)
            if refinable[2] and refinable[1]:  # right and top
                future_refinement_list.append(r + shift + top_right)
            if refinable[2] and refinable[3]:  # right and bottom
                future_refinement_list.append(r + shift + bottom_right)
            if future_refinement_list:
                max_future_r = max(future_refinement_list)
                if max_future_r < midpoint:
                    future_r_left_list.extend(future_refinement_list)
                else:
                    future_r_right_list.extend(future_refinement_list)

            shift += 3
        future_r_left_list.sort(reverse=True)
        future_r_right_list.sort(reverse=True)
        return cell_neighbours, future_r_left_list, future_r_right_list, midpoint

    def update_neighbours(self, pos1, pos2, i, r, shift, cell_neighbours, replaced_neighbours):
        # refined cells with non-refined cells
        cell_neighbours[r + shift + pos1][i].append(replaced_neighbours[r][i][0])
        cell_neighbours[r + shift + pos2][i].append(replaced_neighbours[r][i][0])
        cell_neighbours[replaced_neighbours[r][i][0]][2].extend([r + shift + pos1, r + shift + pos2])

    def update_neighbours_v2(self, pos1, pos2, pos3, pos4, i, j, r, shift, cell_neighbours):
        # refined cells with refined cells
        cell_neighbours[r + shift + pos1][i].append(pos3)
        cell_neighbours[pos3][j].append(r + shift + pos1)
        cell_neighbours[r + shift + pos2][i].append(pos4)
        cell_neighbours[pos4][j].append(r + shift + pos2)

    def get_shape(self, r):
        if self.production_list[r + 1] == 'h':
            bottom_left = 0
            top_left = 1
            top_right = 2
            bottom_right = 3
        elif self.production_list[r + 1] == 'a':
            bottom_left = 0
            top_left = 3
            top_right = 2
            bottom_right = 1
        elif self.production_list[r + 1] == 'b':
            bottom_left = 2
            top_left = 1
            top_right = 0
            bottom_right = 3
        else:
            bottom_left = 2
            top_left = 3
            top_right = 0
            bottom_right = 1
        return top_right, top_left, bottom_left, bottom_right

    def calculate_average(self, cell_neighbours):
        cell_count = len(cell_neighbours)
        ratio_list = [4 for _ in range(0, cell_count)]  # all surface:volumes where n = 1 is  4
        for s in range(0, cell_count):  # start partition at s
            surfaces = 4                # setup as if n = 1
            depths = [self.depths_list[s]]
            min_depth = self.depths_list[s]
            mapped_pattern = ['h']         # all partitions start with h
            for n in range(2, cell_count - (s - 1)):  # number of cells in partition
                # append depth of next cell in path upon incrementing
                depths.append(self.depths_list[s + n - 1])
                # if first non-terminal starts with 'a','b' or 'c' map to 'h'
                if self.production_list[s] == 'a':
                    mapped_pattern.append(self.map_to_h(s + n - 1, 'a', 'h', 'c', ))  # c = 'b'
                elif self.production_list[s] == 'b':
                    mapped_pattern.append(self.map_to_h(s + n - 1, 'b', 'c', 'h', ))  # c = 'a'
                elif self.production_list[s] == 'c':
                    mapped_pattern.append(self.map_to_h(s + n - 1, 'c', 'b', 'a'))  # c = 'h'
                else:  # h = h, a = a, b = b, c = c, d = d
                    mapped_pattern.append(self.production_list[s + n - 1])
                min_depth = min(min_depth, self.depths_list[s + n - 1])
                if min_depth > 1:  # upscale
                    shifted_depths = [d + 1 - min_depth for d in depths]
                    surfaces = self.find_pattern(shifted_depths, s, s + n, surfaces, mapped_pattern, cell_neighbours)
                else:
                    surfaces = self.find_pattern(depths, s, s + n, surfaces, mapped_pattern, cell_neighbours)
                ratio_list.append(surfaces / n)
        return [sum(ratio_list) / len(ratio_list), len(ratio_list)]

    def map_to_h(self, i, h, a, b):
        if self.production_list[i] == h:
            return 'h'
        elif self.production_list[i] == a:
            return 'a'
        elif self.production_list[i] == b:
            return 'b'
        else:  # c
            return 'c'

    # pattern_dict = {depth -> pattern -> depth -> pattern ->...}
    def find_pattern(self, depths, start, end, old_surfaces, pattern, cell_neighbours):
        global pattern_dict
        path = pattern_dict
        for i in range(0, len(depths)):  # find or create empty item in pattern_dict for the specific production list
            if depths[i] in path:
                path = path[depths[i]]
            else:  # nest dictionaries if pattern not in dictionary
                path[depths[i]] = {}
                path = path[depths[i]]
            if pattern[i] in path:
                path = path[pattern[i]]
            else:
                path[pattern[i]] = {}
                path = path[pattern[i]]

        if '-' not in path:  # if pattern not in dictionary, find surfaces brute force and
            surfaces = self.subpartition_average(start, end, old_surfaces, cell_neighbours)
            path['-'] = surfaces  # add pattern to dictionary
        return path['-']

    def subpartition_average(self, start, end, surfaces, cell_neighbours):
        # use surface length from previous partition
        for i in range(start, end - 1):
            for boundary in cell_neighbours[i]:  # check neighbours of each cell in partition
                if end - 1 in boundary:  # if newly appended cell has a boundary with another cell in the partition
                    remove = True
                    for j in range(0, len(boundary)):
                        # if any of the neighbours along the boundary are not in the partition
                        if not (start <= boundary[j] < end):
                            remove = False
                            break
                    if remove:  # if all neighbours along the boundary are in the partition, lower surface count
                        surfaces -= 1
        for boundary in cell_neighbours[end - 1]:  # for each boundary that has at least one neighbour
            if boundary:
                for i in range(0, len(boundary)):
                    if not (start <= boundary[i] < end):
                        surfaces += 1  # add one to the surface count
                        break
            else:
                surfaces += 1
        return surfaces


if __name__ == "__main__":
    open(filename, 'w').close()     # clear file
    max_depth = args.max_depth
    start_time = time.time()
    temp2_time = time.time()

    curve_list = [Curve(0, ['h'], [0], {0: [[], [0], [], []]}, [0], 0, 1)]
    curve_list[0].refinement_left_list = [1, 0]  # correct variables in the initial curve
    curve_list[0].refinement_right_list = [3, 2]
    curve_list[0].midpoint = 2
    for depth in range(1, max_depth+1):  # run through depth levels
        if not depth == 1:
            curve_list = get_curves(curve_list)  # list of all the curves in current depth level
        minimum = min(item.quadtree_averages[0] for item in curve_list)
        print("min: " + str(minimum))       # value of the curve with the minimum s/v ratio
        maximum = max(item.quadtree_averages[0] for item in curve_list)
        print("max: " + str(maximum))       # value of the curve with the maximum s/v ratio
        average = sum([c.quadtree_averages[0] * c.quadtree_averages[1] * c.duplication for c in curve_list]
                      ) / sum([c.quadtree_averages[1] * c.duplication for c in curve_list])
        print("average: " + str(average))   # value of average s/v for all curves in the given depth level
        print("--- %s seconds ---" % (time.time() - temp2_time))
        temp_time2 = time.time()
