import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import argparse

parser = argparse.ArgumentParser(description="Bitstream")
parser.add_argument('bitstream', type=str)
args = parser.parse_args()


def h(depth):
    global bitstream, extra_depth
    if len(bitstream) > 0:
        if bitstream[0] == 1:
            del bitstream[0]

            a(depth + 1)
            draw_line(depth + extra_depth,'up')
            extra_depth = 0

            h(depth + 1)
            draw_line(depth + extra_depth, 'right')
            extra_depth = 0

            h(depth + 1)
            draw_line(depth + extra_depth, 'down')
            extra_depth = 0

            b(depth + 1)

            if not bitstream:
                draw_line(depth + extra_depth, 'end')
            extra_depth += 1

        else:
            del bitstream[0]


def a(depth):
    global bitstream, extra_depth
    if len(bitstream) > 0:
        if bitstream[0] == 1:
            del bitstream[0]

            h(depth + 1)
            draw_line(depth + extra_depth, 'right')
            extra_depth = 0

            a(depth + 1)
            draw_line(depth + extra_depth, 'up')
            extra_depth = 0

            a(depth + 1)
            draw_line(depth + extra_depth, 'left')
            extra_depth = 0

            c(depth + 1)
            extra_depth += 1
        else:
            del bitstream[0]


def b(depth):
    global bitstream, extra_depth
    if len(bitstream) > 0:
        if bitstream[0] == 1:
            del bitstream[0]

            c(depth + 1)
            draw_line(depth + extra_depth, 'left')
            extra_depth = 0

            b(depth + 1)
            draw_line(depth + extra_depth, 'down')
            extra_depth = 0

            b(depth + 1)
            draw_line(depth + extra_depth, 'right')
            extra_depth = 0

            h(depth + 1)
            extra_depth += 1
        else:
            del bitstream[0]


def c(depth):
    global bitstream, extra_depth
    if len(bitstream) > 0:
        if bitstream[0] == 1:
            del bitstream[0]

            b(depth + 1)
            draw_line(depth + extra_depth, 'down')
            extra_depth = 0

            c(depth + 1)
            draw_line(depth + extra_depth, 'left')
            extra_depth = 0

            c(depth + 1)
            draw_line(depth + extra_depth, 'up')
            extra_depth = 0

            a(depth + 1)
            extra_depth += 1

        else:
            del bitstream[0]


def draw_line(depth, direction):
    global oldx, oldy, first, old_depth, prev_direction, processor, count, cells, remainder, line_count
    length = math.pow(2, depth)
    new_depth = grid_size / length
    depth_change = (new_depth - old_depth) / 2

    if prev_direction == 'start':
        # first = False
        newx = grid_size/(length*2)
        newy = grid_size/(length*2)
    elif prev_direction == 'right':
        newx = oldx + new_depth - depth_change
        newy = oldy + depth_change
        plt.plot([oldx, newx], [oldy, newy], 'r', linewidth=4)
    elif prev_direction == 'up':
        newx = oldx + depth_change
        newy = oldy + new_depth - depth_change
        plt.plot([oldx, newx], [oldy, newy], 'r', linewidth=4)
    elif prev_direction == 'left':
        newx = oldx - new_depth + depth_change
        newy = oldy - depth_change
        plt.plot([oldx, newx], [oldy, newy], 'r', linewidth=4)
    elif prev_direction == 'down':
        newx = oldx - depth_change
        newy = oldy - new_depth + depth_change
        plt.plot([oldx, newx], [oldy, newy], 'r', linewidth=4)
    else:
        newx = oldx  # should be unreachable
        newy = oldy
    p = patches.Rectangle((newx - (grid_size/(length * 2)), newy - (grid_size/(length * 2))), grid_size / length,
                          grid_size / length, fill=False, edgecolor='k', linewidth=2)
    plt.gca().add_patch(p)
    count += 1

    # plt.gca().text(newx, newy, processor[line_count], horizontalalignment='center', verticalalignment='center',
                   # fontsize=32, color='black')
    line_count += 1
    oldx = newx
    oldy = newy
    old_depth = new_depth
    prev_direction = direction


if __name__ == "__main__":
    grid_size = 64
    plt.xticks(np.arange(0, grid_size + 1, 8))  # x labels
    plt.yticks(np.arange(0, grid_size + 1, 8))  # y labels
    plt.axis([0, grid_size, 0, grid_size])

    count = 0
    oldx = 0
    oldy = 0
    extra_depth = 0
    prev_direction = 'start'
    old_depth = 0
    line_count = 0

    str_bitstream = args.bitstream
    cells = math.floor(str_bitstream.count('0')/3)
    remainder = str_bitstream.count('0') % 3
    print(str_bitstream.count('0'), " + ", remainder)
    bitstream = list(map(int, str_bitstream))   # 110000000 101000000 1100001000000 1100000010000
    h(1)
    plt.plot([64, 64], [0, 64], 'k', linewidth=4)
    plt.plot([0, 0], [0, 64], 'k', linewidth=4)
    plt.plot([0, 64], [0, 0], 'k', linewidth=4)
    plt.plot([0, 64], [64, 64], 'k', linewidth=4)

    plt.axis('off')

    plt.show()
