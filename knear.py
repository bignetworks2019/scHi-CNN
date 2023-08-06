import sys

import numpy as np
import os
import argparse
import traceback
import heapq
from scipy.sparse import csr_matrix
from mpi4py import MPI
import gc
import cooler
import h5py


def getMatrix(path):
    mcool = h5py.File(path, 'r')
    cool = cooler.Cooler(mcool)
    cool_mat = cool.matrix(balance=False)
    return cool_mat


def load_chr_mat(cells, prefix, type, chr):
    matrix_map = {}
    for cell in cells:
        print(cell)
        mat = getMatrix(prefix + type + '/' + cell)

        chr_mat = mat.fetch(chr)
        np.fill_diagonal(chr_mat, 0)
        matrix_map[cell] = csr_matrix(chr_mat)
        # del chr_mat

    return matrix_map


def load_all_mat(cells, prefix, type, chr_all):
    matrix_map = {}
    for cell in cells:
        print(cell)
        mat = getMatrix(prefix + type + '/' + cell)
        matrix_map[cell] = {}

        for chr in chr_all:
            chr_mat = mat.fetch(chr)
            np.fill_diagonal(chr_mat, 0)
            matrix_map[cell][chr] = csr_matrix(chr_mat)
            # del chr_mat

    return matrix_map


def main(type, prefix, resol, out_suffix, cutoff, rank, size, comm, topk, pad_local, cell_file_path, genome,
         non_zero_flag):
    # cell_file = prefix + type + '-orig.txt'
    # cell_file = prefix + type + '.txt'
    if cell_file_path:
        cell_file = cell_file_path
    else:
        cell_file = prefix + type + '.txt'
    out_dir = prefix + '/' + type + '-' + out_suffix + '/'

    if genome == 'mm10':
        chr_arr = list(range(1, 20)) + ['X']
    elif genome == 'hg19':
        chr_arr = list(range(1, 23)) + ['X']
    else:
        print("Genome:", genome, "is not supported.")
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)

    chr_all = []
    for k in chr_arr:
        chr_all.append("chr" + str(k))
    # chr_all = list(range(5, 8)) + ['X']

    cells = []
    if rank == 0:
        with open(cell_file) as f:
            for line in f:
                cells.append(line.rstrip())
    else:
        cells = None

    cells = comm.bcast(cells, root=0)

    print('start')

    files_per_process, remainder = divmod(len(cells), size)

    start_idx = rank * files_per_process + min(rank, remainder)
    if rank < remainder:
        end_idx = start_idx + files_per_process + 1
    else:
        end_idx = start_idx + files_per_process

    for i in range(start_idx, end_idx):
        process_chr(chr_all, cells[i], resol, cutoff, topk, out_dir, pad_local, cells, prefix, type, non_zero_flag)


def process_chr(chr_all, cell, resol, cutoff, topk, out_dir, pad_local, cells, prefix, type, non_zero_flag):
    try:
        out = open(out_dir + cell, 'w')
        out_val = open(out_dir + cell + '.vals', 'w')
        imputed = 0
        not_imputed = 0
        for chr in chr_all:
            matrix_map = load_chr_mat(cells, prefix, type, chr)
            print(cell, chr)
            mat = matrix_map[cell]

            gc.collect()

            new_mat = mat.toarray().astype(np.float32).copy()

            for i in range(0, new_mat.shape[0]):
                for j in range(i + 1, new_mat.shape[1]):
                    val = mat[i, j]
                    if val != 0:
                        new_mat[i, j] = val
                    else:
                        dist = (abs(i - j)) * resol

                        if dist > cutoff:
                            new_mat[i, j] = val
                        else:
                            # continue
                            new_mat[i, j], vals_out = get_local_val(i, j, mat, matrix_map, topk, pad_local, cell,
                                                                    non_zero_flag)
                            if len(vals_out) > 0:
                                out_val.write(
                                    "{}\t{}\t{}\t{}\t{}\n".format(chr, cell, i, j, '\t'.join(map(str, vals_out))))
                                not_imputed += 1
                            else:
                                imputed += 1

            out_val.flush()

            indices = np.nonzero(np.triu(new_mat))

            for i, j in zip(indices[0], indices[1]):
                out.write('{}\t{}\t{}\t{}\n'.format(chr, i, j, new_mat[i, j]))

            out.flush()

        out.close()
        out_val.close()

        print(cell, imputed, not_imputed)

    except Exception as e:
        print(e)
        traceback.print_exc()


def get_spear_csr(mat1, mat2):
    arr1 = mat1.flatten()
    arr2 = mat2.flatten()

    coefficient = np.corrcoef(arr1, arr2)[0, 1]
    # coefficient = 0
    if np.isnan(coefficient):
        coefficient = -1

    return coefficient


def check_submat_nonzero_elements(mat, non_zero_flag):
    non_zero_count = mat.count_nonzero()

    if non_zero_flag:
        if non_zero_count >= 1:
            # print(mat.toarray())
            return True
        else:
            return False

    if non_zero_count >= (0.05 * mat.shape[0] * mat.shape[1]):
        # print(mat.toarray())
        return True
    else:
        return False


def get_local_val(i, j, mat, matrix_map, topk, pad, cell_mat, non_zero_flag):
    n_i_1 = i - pad
    n_i_2 = i + pad + 1

    n_j_1 = j - pad
    n_j_2 = j + pad + 1

    if n_i_1 < 0 or n_j_1 < 0 or n_i_2 > mat.shape[0] or n_j_2 > mat.shape[1]:
        return 0, []

    sub_mat = mat[n_i_1:n_i_2, n_j_1:n_j_2]

    if not check_submat_nonzero_elements(sub_mat, non_zero_flag):
        return 0, []

    val_arr = []

    for cell, chr_mats in matrix_map.items():
        if cell == cell_mat:
            continue
        sub_mat_comp = chr_mats[n_i_1:n_i_2, n_j_1:n_j_2]

        coef = get_spear_csr(sub_mat.toarray(), sub_mat_comp.toarray())

        if coef <= 0:
            continue

        val_arr.append((coef, chr_mats[i, j]))

    if len(val_arr) < 4:
        return 0, []

    sorted_val = heapq.nlargest(4, val_arr)

    vals_out = []

    sum = 0
    sum_coef = 0
    for i in range(0, topk):
        sum += sorted_val[i][1] * sorted_val[i][0]
        sum_coef += sorted_val[i][0]
        vals_out.append(sorted_val[i][0])

    avg = sum / sum_coef
    return avg, vals_out


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--resol', type=int, required=True)
    parser.add_argument('--cutoff', type=int, required=True)
    parser.add_argument('--pad-local', type=int, default=5, required=False)
    parser.add_argument('--topk', type=int, default=4, required=False)
    parser.add_argument('--type', action='store', required=True)
    parser.add_argument('--prefix', action='store', required=True)
    parser.add_argument('--out-suffix', action='store', required=True)
    parser.add_argument('--cell-file', action='store', required=False)
    parser.add_argument('--genome', action='store', required=True)
    parser.add_argument('--non-zero', action='store_true')
    return parser


if __name__ == '__main__':
    parser = get_args()
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    main(args.type, args.prefix + '/', args.resol, args.out_suffix, args.cutoff, rank, size, comm, args.topk,
         args.pad_local, args.cell_file, args.genome, args.non_zero)

