import os
import sys
import numpy as np
import multiprocessing
import argparse
import cooler
import h5py


def getMatrix(path):
    mcool = h5py.File(path, 'r')
    cool = cooler.Cooler(mcool)
    cool_mat = cool.matrix(balance=False)
    return cool_mat


def load_ref_mat(path, genome):
    if genome == 'mm10':
        chr_arr = list(range(1, 20)) + ['X']
    elif genome == 'hg19':
        chr_arr = list(range(1, 23)) + ['X']
    else:
        print("Genome:", genome, "is not supported.")
        sys.exit(1)

    chr_all = []
    for k in chr_arr:
        chr_all.append("chr" + str(k))

    matrix_map = {}

    mat = getMatrix(path)

    for chr in chr_all:
        chr_mat = mat.fetch(chr)
        matrix_map[chr] = np.zeros(chr_mat.shape)

    return matrix_map


def write_mat(k, mat_out, out_f):
    indices = np.nonzero(mat_out)

    for i, j in zip(indices[0], indices[1]):
        out_f.write('{}\t{}\t{}\t{}\n'.format(k, i, j, mat_out[i, j]))


def normalize(arr, max_bin_dist):
    for i in range(1, max_bin_dist + 1):
        diag_entry = np.diagonal(arr, offset=i)
        diag_entry_copy = np.copy(diag_entry)
        diag_entry_stat = np.copy(diag_entry).tolist()

        diag_entry_stat.sort()
        diag_entry_stat.reverse()
        trim_index = round(0.01 * len(diag_entry_stat)) - 1
        remaining = diag_entry_stat[(trim_index):]

        # thresh = np.percentile(diag_entry_copy, 99)
        # diag_entry_copy[diag_entry_copy > thresh] = 0

        # print(diag_entry_copy)

        mean = np.mean(remaining)
        std_dev = np.std(remaining)

        if std_dev < 1e-6:
            np.fill_diagonal(arr[:, i:], 0)
        else:
            z_scores = [(x - mean) / std_dev for x in diag_entry_copy]
            if np.isnan(z_scores).any():
                print(z_scores)
            np.fill_diagonal(arr[:, i:], z_scores)

    return arr


def main(type, prefix, in_suffix, out_suffix, cell_file, genome, res, max_dist):
    outdir = prefix + type + '-' + out_suffix

    n_file = cell_file

    os.makedirs(outdir, exist_ok=True)

    max_bin_dist = int(max_dist / res)

    cells = []
    with open(n_file) as f_f:
        for line_f in f_f:
            cells.append(line_f.rstrip())
            # process_cell(prefix, line_f, outdir)

    params = [(prefix, in_suffix, type, cell, outdir, genome, max_bin_dist) for cell in cells]
    with multiprocessing.Pool(4) as pool:
        pool.starmap(process_cell, params)


def process_cell(prefix, in_suffix, type, line_f, outdir, genome, max_bin_dist):
    print(prefix, in_suffix, type, line_f, outdir, genome)
    f_path = prefix + type + '-' + in_suffix + '/' + line_f.rstrip()
    matrix = load_ref_mat(prefix + type + '/' + line_f.rstrip(), genome)

    with open(f_path) as f:
        for line in f:
            split = line.rstrip().split()

            chr = split[0]
            bin1 = int(split[1])
            bin2 = int(split[2])
            val = float(split[3])

            matrix[chr][bin1][bin2] = val

    out_f = open(outdir + '/' + line_f.rstrip(), 'w')

    for k, v in matrix.items():
        # print(k)
        mat_out = normalize(v, max_bin_dist)

        write_mat(k, mat_out, out_f)

    out_f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', action='store', required=True)
    parser.add_argument('--prefix', action='store', required=True)
    parser.add_argument('--res', action='store', type=int, required=False, default=100000)
    parser.add_argument('--in-suffix', action='store', required=True)
    parser.add_argument('--out-suffix', action='store', required=True)
    parser.add_argument('--cell-file', action='store', required=True)
    parser.add_argument('--genome', action='store', required=True)
    parser.add_argument('--max-dist', action='store', type=float, default=2e6, required=False)

    args = parser.parse_args()

    main(args.type, args.prefix + '/', args.in_suffix, args.out_suffix, args.cell_file, args.genome, args.res,
         args.max_dist)

