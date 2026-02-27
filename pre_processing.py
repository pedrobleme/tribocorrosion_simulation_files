import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import weibull_min
import time
from collections import defaultdict

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "lines.linewidth": 2.0,
    "axes.grid": True,
    "grid.linestyle": "--",
    "grid.alpha": 0.6
})


def read_inp_file(filename):
    """
    Reads an Abaqus .inp file and extracts nodal coordinates and element connectivity.

    Parameters
    ----------
    filename : str
        Path to the Abaqus input file.

    Returns
    -------
    nodes : ndarray
        Array containing node IDs and their spatial coordinates.
    ele : ndarray
        Array containing element IDs and their associated node connectivity.
    """
    with open(filename, 'r') as fid:
        lines = [line.strip() for line in fid if line.strip()]

    # Identify relevant sections in the input file
    node_start = next(i for i, line in enumerate(lines) if line.startswith('*Node'))
    elem_start = next(i for i, line in enumerate(lines) if line.startswith('*Element'))
    nset_start = next(i for i, line in enumerate(lines) if line.startswith('*Nset'))

    # Read nodal data
    nodes = []
    for line in lines[node_start + 1:elem_start]:
        parts = [p.strip() for p in line.split(',') if p.strip()]
        if len(parts) >= 4:
            nodes.append([
                float(parts[0]),
                float(parts[1]),
                float(parts[2]),
                float(parts[3])
            ])

    # Read element connectivity (8-node elements)
    ele = []
    for line in lines[elem_start + 1:nset_start]:
        parts = [p.strip() for p in line.split(',') if p.strip()]
        if len(parts) >= 9:
            ele.append([int(parts[0])] + [int(p) for p in parts[1:9]])

    return np.array(nodes), np.array(ele)


def calculate_centroids(ele, nodes):
    """
    Computes element centroids based on nodal coordinates.
    A fallback strategy consistent with a legacy MATLAB implementation
    is applied if the centroid falls outside the element bounds.

    Parameters
    ----------
    ele : ndarray
        Element connectivity array.
    nodes : ndarray
        Nodal coordinates array.

    Returns
    -------
    coordx, coordy, coordz : ndarray
        Centroid coordinates of each element.
    """
    node_dict = {int(node[0]): node[1:4] for node in nodes}
    centroids = np.zeros((len(ele), 3))

    for i, elem in enumerate(ele):
        node_coords = np.array([node_dict[nid] for nid in elem[1:9]])

        # Primary centroid calculation using arithmetic mean
        centroid = np.mean(node_coords, axis=0)

        # Check if centroid lies within the element bounding box
        if not (
            np.all(centroid >= np.min(node_coords, axis=0)) and
            np.all(centroid <= np.max(node_coords, axis=0))
        ):
            # Fallback method
            dx = node_coords[0, 0] - node_coords[:, 0]
            dy = node_coords[0, 1] - node_coords[:, 1]
            dz = node_coords[0, 2] - node_coords[:, 2]

            centroid = [
                node_coords[0, 0] + (np.max(dx) if np.max(dx) > 0 else np.min(dx)) / 2,
                node_coords[0, 1] + (np.max(dy) if np.max(dy) > 0 else np.min(dy)) / 2,
                node_coords[0, 2] + (np.max(dz) if np.max(dz) > 0 else np.min(dz)) / 2
            ]

        centroids[i] = centroid

    return centroids[:, 0], centroids[:, 1], centroids[:, 2]


def find_neighbors(ele):
    """
    Identifies neighboring elements based on shared nodes.
    Two elements are considered neighbors if they share at least four nodes.

    Parameters
    ----------
    ele : ndarray
        Element connectivity array.

    Returns
    -------
    neighbors : list of lists
        List containing the neighboring element IDs for each element.
    """
    node_to_elems = defaultdict(list)

    for i, elem in enumerate(ele):
        for node in elem[1:9]:
            node_to_elems[node].append(i)

    neighbors = [[] for _ in range(len(ele))]

    for i, elem in enumerate(ele):
        neighbor_counts = defaultdict(int)
        for node in elem[1:9]:
            for neighbor_elem in node_to_elems[node]:
                if neighbor_elem != i:
                    neighbor_counts[neighbor_elem] += 1

        valid_neighbors = [
            ele[k, 0] for k, count in neighbor_counts.items() if count >= 4
        ]

        neighbors[i] = sorted(valid_neighbors)[:6]

    return neighbors


def read_volume_file(filename):
    """
    Reads element volume data from an Abaqus output file.

    Parameters
    ----------
    filename : str
        Path to the volume file.

    Returns
    -------
    volumes : ndarray
        Element volume values.
    """
    volumes = []
    with open(filename, 'r') as fid:
        start_reading = False
        for line in fid:
            if 'Element' in line and 'EVOL' in line:
                start_reading = True
                next(fid)
                next(fid)
                continue
            if start_reading:
                parts = line.strip().split()
                if len(parts) == 2 and parts[0].isdigit():
                    volumes.append(float(parts[1]))

    return np.array(volumes)


def main():
    """
    Main execution routine:
    - Reads mesh data
    - Computes element centroids
    - Identifies neighboring elements
    - Assigns stochastic pitting corrosion values
    - Produces a Weibull distribution validation plot
    - Writes processed data to output file
    """
    start_time = time.time()

    print("Loading input data...")
    nodes, ele = read_inp_file('finite_element_model.inp') # Abaqus input file

    print("Computing element centroids...")
    coordx, coordy, coordz = calculate_centroids(ele, nodes)

    print("Identifying neighboring elements...")
    neighbors = find_neighbors(ele)

    # Connectivity matrix
    max_neighbors = 6
    con_list = np.zeros((len(ele), max_neighbors + 1))
    con_list[:, 0] = ele[:, 0]

    for i, neigh in enumerate(neighbors):
        con_list[i, 1:len(neigh) + 1] = neigh

    # Surface identification
    surface_flag = np.array([1 if len(neigh) < 6 else 0 for neigh in neighbors])
    con_list_final = np.column_stack((con_list, surface_flag))

    print("Reading element volumes...")
    volume_values = read_volume_file('element_volumes.txt') # INPUT - Element volume file

    print("Generating stochastic pitting corrosion values...")
    shape_param = 4.0
    scale_param = 10.0
    pitting_values = weibull_min.rvs(
        shape_param, scale=scale_param, size=len(ele)
    )

    # --- Weibull validation plot  ---
    x = np.linspace(0.0, np.max(pitting_values) * 1.1, 300)
    weibull_pdf = weibull_min.pdf(x, shape_param, scale=scale_param)

    plt.figure(figsize=(7, 5))
    plt.hist(
        pitting_values,
        bins=30,
        density=True,
        alpha=0.6,
        label="Generated pitting values"
    )
    plt.plot(
        x,
        weibull_pdf,
        linewidth=2.0,
        label="Theoretical Weibull PDF"
    )
    plt.xlabel("Pitting corrosion parameter")
    plt.ylabel("Probability density")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    # ------------------------------------------

    max_corrosion = np.max(pitting_values)
    normalized_pitting = (
        pitting_values * con_list_final[:, 7] / max_corrosion
    )

    final_data = np.column_stack((
        con_list_final[:, :8].astype(int),
        normalized_pitting,
        volume_values,
        coordx, coordy, coordz,
        pitting_values
    ))

    print("Saving output file...")
    np.savetxt(
        'element_properties.txt', # Output file name
        final_data,
        delimiter=',',
        fmt=['%d'] * 8 + ['%.5f'] * 6
    )

    print(f"Execution completed in {time.time() - start_time:.2f} seconds.")


if __name__ == "__main__":
    main()
