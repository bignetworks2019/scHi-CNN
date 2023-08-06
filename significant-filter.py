import sys
import numpy as np
from statsmodels.stats.multitest import multipletests
import argparse
import cooler
import h5py


def getMatrix(path):
    mcool = h5py.File(path, 'r')
    cool = cooler.Cooler(mcool)
    cool_mat = cool.matrix(balance=False)
    return cool_mat


def load_ref_mat(path, chr_all):
    matrix_map = {}

    mat = getMatrix(path)

    for chr in chr_all:
        chr_mat = mat.fetch(chr)
        matrix_map[chr] = np.zeros(chr_mat.shape)

    return matrix_map


def check_z_thresh(val_arr):
    n_greater = np.sum(val_arr > 1.96)
    prop_greater = n_greater / len(val_arr)

    if prop_greater > 0.1:
        return True
    else:
        return False


def main(type, genome, prefix, in_suffix, center_dist, cell_file, fdr, tstat, res, max_dist):
    n_file = cell_file

    num_cells = 0
    with open(n_file, 'r') as file:
        for line in file:
            num_cells += 1
            ref_cell = line.rstrip()

    matrix_chr_arr = {}

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

    matrix = load_ref_mat(prefix + type + '/' + ref_cell, chr_all)

    for chr in chr_all:
        matrix_chr_arr[chr] = np.zeros((num_cells, matrix[chr].shape[0], matrix[chr].shape[1]))

    count = 0

    with open(n_file) as f_f:
        for line_f in f_f:
            f_path = prefix + type + '-' + in_suffix + '/' + line_f.rstrip()
            matrix = load_ref_mat(prefix + type + '/' + line_f.rstrip(), chr_all)

            with open(f_path) as f:
                for line in f:
                    split = line.rstrip().split()

                    chr = split[0]
                    bin1 = int(split[1])
                    bin2 = int(split[2])
                    val = float(split[3])

                    matrix[chr][bin1][bin2] = val

            for chr in chr_all:
                matrix_chr_arr[chr][count] = matrix[chr]

            count += 1

    outfile = open(prefix + type + '-significant.txt', 'w')

    max_bin_dist = int(max_dist / res)

    for chr in chr_all:

        sum_array_chr = np.sum(matrix_chr_arr[chr], axis=0)

        distance_map = {}
        pval_map = {}
        tstat_map = {}

        for i in range(center_dist, sum_array_chr.shape[0] - center_dist):
            for j in range(i + center_dist, sum_array_chr.shape[1] - center_dist):

                if (j - i) > max_bin_dist:
                    break

                val_arr = matrix_chr_arr[chr][:, i, j]

                if np.mean(val_arr) < 0 or not check_z_thresh(val_arr):
                    continue

                def surrounding_mean(matrix, row, col, distance=2):
                    start_row, end_row = row - distance, row + distance + 1
                    start_col, end_col = col - distance, col + distance + 1

                    region = matrix[:, start_row:end_row, start_col:end_col]

                    mask = np.full((region.shape[1], region.shape[2]), True)
                    mask[1:mask.shape[0] - 1, 1:mask.shape[1] - 1] = False

                    avg_neigh = np.zeros(matrix.shape[0])
                    count = 0
                    for l in range(0, mask.shape[0]):
                        for m in range(0, mask.shape[1]):
                            if mask[l][m]:
                                # print(l, m)
                                avg_neigh = avg_neigh + region[:, l, m]
                                count += 1

                    # print(count)
                    avg_neigh = avg_neigh / count

                    return avg_neigh

                avg_neigh = surrounding_mean(matrix_chr_arr[chr], i, j, distance=center_dist)

                t_stat, pval = paired_ttest(val_arr, avg_neigh)

                dist = j - i

                if dist not in distance_map:
                    distance_map[dist] = []
                    pval_map[dist] = []
                    tstat_map[dist] = []

                distance_map[dist].append((chr, i, j))
                pval_map[dist].append(pval)
                tstat_map[dist].append(t_stat)

        count_chr_thresh = 0
        count_all = 0
        for d, values in pval_map.items():
            q_vals = multipletests(values, method='fdr_bh')[1]
            bin_list = distance_map[d]
            tstat_li = tstat_map[d]
            pval_li = pval_map[d]

            for z in range(0, len(bin_list)):
                count_all += 1
                if q_vals[z] <= fdr and tstat_li[z] >= tstat:
                    count_chr_thresh += 1

                    chr = bin_list[z][0]
                    bin1_x = bin_list[z][1] * res
                    bin2_x = bin_list[z][2] * res

                    outfile.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr, bin1_x, bin1_x + res, chr, bin2_x,
                                                                     bin2_x + res,
                                                                     tstat_li[z], pval_li[z], q_vals[z]))


                    # outfile.write("{}\t{}\t{}\n".format(chr, bin1_x, bin2_x))
                                                                      
                                                                     

        print(chr, sum_array_chr.shape, count_all, count_chr_thresh)

        # print(len(q_vals), count_chr_thresh)

        # outfile.flush()

    outfile.close()


def paired_ttest(val_arr, param):
    from scipy.stats import ttest_rel

    t_stat, p_val = ttest_rel(val_arr, param)

    return t_stat, p_val


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', action='store', required=True)
    parser.add_argument('--genome', action='store', required=True)
    parser.add_argument('--res', action='store', type=int, required=False, default=100000)
    parser.add_argument('--prefix', action='store', required=True)
    parser.add_argument('--in-suffix', action='store', required=True)
    parser.add_argument('--center-dist', action='store', type=int, default=2, required=False)
    parser.add_argument('--max-dist', action='store', type=float, default=2e6, required=False)
    parser.add_argument('--cell-file', action='store', required=True)
    parser.add_argument('--fdr', action='store', type=float, default=0.1, required=False)
    parser.add_argument('--tstat', action='store', type=float, default=3.0, required=False)
    args = parser.parse_args()
    main(args.type, args.genome, args.prefix + '/', args.in_suffix, args.center_dist, args.cell_file, args.fdr,
         args.tstat, args.res, args.max_dist)

