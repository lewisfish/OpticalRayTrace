from typing import List

import matplotlib.pyplot as plt  # type: ignore[import]
import numpy as np  # type: ignore[import]


def get_rays(file: str) -> List[List[List[float]]]:
    """Read rays from file.

    Args:
        file: str
            file from which rays are read. Format of file is
                x11, y11, z11
                x12, y12, z12
                .   .   .
                x1n, y1n, z1n



                x21, y21, z21
                x22, y22, z22
                .   .   .
                x2n, y2n, z2n

    Returns:
        List[List[List[float]]]: Description
    """
    rays = []
    with open(file, "r") as f:
        lines = f.readlines()
        skip_counter = 0
        tmp = []
        for line in lines:
            length = len(line)
            if length > 3:
                skip_counter = 0
                line = line.lstrip().rstrip()
                line = line.replace("    ", " ")
                line = line.replace("   ", " ")
                line = line.replace("  ", " ")
                data = line.split(" ")
                tmp.append(list(map(float, data)))
            else:
                if skip_counter == 0:
                    if len(tmp) != 0:
                        rays.append(tmp)
                        tmp = []
                skip_counter += 1
                if skip_counter == 3:
                    skip_counter = 0
                tmp = []
                continue
    return rays


def plot_rays(rays: List[List[List[float]]]) -> None:
    """Plot rays in 3D line plot with each ray the same colour

    Args:
        rays (List[List[List[float]]]): Description
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for ray in rays:
        ray = np.array(ray).T
        ax.plot(*ray, color="#1f77b4")
        ax.scatter(*ray, color="#ff7f0e")
    plt.show()


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-f", "--file", type=str,
                        help="Name of file to be plotted.")

    args = parser.parse_args()
    file = args.file

    rays = get_rays(file)
    plot_rays(rays)
