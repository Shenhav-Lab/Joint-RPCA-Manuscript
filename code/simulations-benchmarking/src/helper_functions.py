import numpy as np
import pandas as pd
from scipy.stats import norm
from skbio.stats.composition import (closure, alr, alr_inv)
from biom import Table
from numpy.random import (poisson, lognormal, randint, normal)

def blocks(ncols, nrows, nblocks, overlap=0, minval=0, sigma=2, maxval=1.0):
    """
    Generate block diagonal with Gaussian distributed values within blocks.
    Parameters
    ----------
    ncol : int
        Number of columns
    nrows : int
        Number of rows
    nblocks : int
        Number of blocks, mucst be greater than one
    overlap : int
        The Number of overlapping columns (Default = 0)
    minval : int
        The min value output of the table (Default = 0)
    maxval : int
        The max value output of the table (Default = 1)
    Returns
    -------
    np.array
        Table with a block diagonal where the rows represent samples
        and the columns represent features.  The values within the blocks
        are gaussian distributed between 0 and 1.
    Note
    ----
    The number of blocks specified by `nblocks` needs to be greater than 1.
    """

    if nblocks <= 1:
        raise ValueError('`nblocks` needs to be greater than 1.')
    mat = np.zeros((nrows, ncols))
    gradient = np.linspace(0, 10, nrows)
    mu = np.linspace(0, 10, ncols)
    xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
          for i in range(len(mu))]
    mat = np.vstack(xs).T
    block_cols = ncols // nblocks #(nblocks * 2)
    block_rows = nrows // nblocks 

    for b in range(nblocks - 1):

        gradient = np.linspace(5, 5, block_rows)  # samples (bock_rows)
        # features (block_cols+overlap)
        mu = np.linspace(0, 10, block_cols + overlap)
        xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
              for i in range(len(mu))]

        B = np.vstack(xs).T * maxval
        lower_row = block_rows * b
        upper_row = min(block_rows * (b + 1), nrows)
        lower_col = block_cols * b
        upper_col = min(block_cols * (b + 1), ncols)

        if b == 0:
            mat[lower_row:upper_row,
                lower_col:int(upper_col + overlap)] = B
        else:
            ov_tmp = int(overlap / 2)
            if (B.shape) == (mat[lower_row:upper_row,
                                 int(lower_col - ov_tmp):
                                 int(upper_col + ov_tmp + 1)].shape):
                mat[lower_row:upper_row, int(
                    lower_col - ov_tmp):int(upper_col + ov_tmp + 1)] = B
            elif (B.shape) == (mat[lower_row:upper_row,
                                   int(lower_col - ov_tmp):
                                   int(upper_col + ov_tmp)].shape):
                mat[lower_row:upper_row, int(
                    lower_col - ov_tmp):int(upper_col + ov_tmp)] = B
            elif (B.shape) == (mat[lower_row:upper_row,
                                   int(lower_col - ov_tmp):
                                   int(upper_col + ov_tmp - 1)].shape):
                mat[lower_row:upper_row, int(
                    lower_col - ov_tmp):int(upper_col + ov_tmp - 1)] = B

    upper_col = int(upper_col - overlap) 
    # Make last block fill in the remainder
    gradient = np.linspace(5, 5, nrows - upper_row)
    mu = np.linspace(0, 10, ncols - upper_col)
    xs = [norm.pdf(gradient, loc=mu[i], scale=sigma)
          for i in range(len(mu))]
    B = np.vstack(xs).T * maxval

    mat[upper_row:, upper_col:] = B

    return mat

def input_matrix_validation(mat, depths):

    if np.any(depths <= 0):
        raise ValueError("Read depth cannot have values "
                         "less than or equal to zero")
    if depths.ndim != 2:
        raise ValueError("Read depth can only have two dimensions")
    if depths.shape[0] != mat.shape[0]:
        raise ValueError("Number of est. read depth does not match number of "
                         "samples in the input matrix")
    # check matrix and ensure
    # data is proportions
    mat = closure(mat)

    return mat

def output_matrix_validation(sim):

    # ensure no negative counts
    sim[sim < 0.0] = 0.0
    # remove zero sums and return a mask (if needed)
    zero_sum_mask_rows = sim.sum(1) > 0
    sim = sim[zero_sum_mask_rows]
    zero_sum_mask_columns = sim.sum(0) > 0
    sim = sim[:, zero_sum_mask_columns]

    return sim, zero_sum_mask_rows, zero_sum_mask_columns

def poisson_lognormal(mat, depths, kappa=1):
    """
    Simulate from counts, probabilities, or
    proportions of input matrix with a
    Poisson Log-Normal distribution.

    Parameters
    ----------
    mat: array_like
        matrix of strictly positive counts
        or probabilities/proportions.
        columns = features (components)
        rows = samples (compositions)
    depth : array_like
        Read depth of the simulation
        for each sample (row).
    kappa: float
        Over-dispersion parameter.
        Default is 1.

    Returns
    -------
    array_like, np.int
       A matrix of counts simulated from
       the input mat by the distribution.
    list, bool
        Mask of rows that summed to zero
    list, bool
        Mask of columns that summed to zero

    Raises
    ------
    ValueError
       Raises an error if any depths are equal
       or less than zero.
    ValueError
       Raises an error if any depths does not have
       exactly 2 dimensions.
    ValueError
       Raises an error if any depths shape does not match the
       input matrix.
    ValueError
       Raises an error if any values are negative.
    ValueError
       Raises an error if the matrix has more than 2 dimension.
    ValueError
       Raises an error if there is a row that has all zeros.

    """

    # check matrix and ensure
    # data is proportions
    mat = input_matrix_validation(mat, depths)
    # simulate from proportions
    mu = depths * mat
    sim = np.vstack([poisson(lognormal(np.log(mu[i, :]), kappa))
                     for i in range(mat.shape[0])])

    return output_matrix_validation(sim)

def add_noise(mat,
              pseudocount=1,
              percent_normal=0.1,
              percent_random=0.1,
              random_count=1,
              add_missing_at_random=False,
              percent_missing=0.1):
    """
    This function transforms count data into
    the simplex with the ALR. This gives the
    data an approximate normal distibution,
    allowing for the addition of normal
    and randomly distributed noise. The data
    is transformed into proportions with the
    inverse alr and then missing values can be
    added, either to match the input or at random.

    Parameters
    ----------
    mat: array_like
        matrix of strictly positive counts
        or probabilities/proportions.
        columns = features (components)
        rows = samples (compositions)
    pseudocount: float
        Pseudocount to add for ALR.
        Default is 1.
    percent_normal: float
        Percent of data to add homoscedastic
        noise. Default is 0.1 (i.e. 10%)
    percent_random: float
        Percent of data to add heteroscedastic
        noise. Default is 0.1 (i.e. 10%)
    random_count: float
        Intensity of random data added.
        Default is 1 (ten would be large).
    add_missing_at_random: bool
        If missing values should match the input
        mat. If True percent_missing is ignored.
    percent_missing: float
        Percent of data to add missing (zero)
        values. Default is 0.1 (i.e. 10%)

    Returns
    -------
    array_like, np.float
       A matrix of noisy proportions.

    Raises
    ------
    ValueError
       Raises an error if any values are negative.
    ValueError
       Raises an error if the matrix has more than 2 dimension.
    ValueError
       Raises an error if there is a row that has all zeros.
    """

    # transform mat into ALR space for
    # adding normal dist. noise
    mat_noise = alr(mat + pseudocount)

    # add homo-scedastic noise
    err = percent_normal * np.ones_like(mat_noise)
    mat_noise = normal(mat_noise, err)

    # add hetero-scedastic noise
    err = percent_random * np.ones_like(mat_noise)
    n_entries = int(percent_random * np.count_nonzero(mat_noise))
    i = randint(0, err.shape[0], n_entries)
    j = randint(0, err.shape[1], n_entries)
    err[i, j] = random_count
    mat_noise = normal(mat_noise, err)

    # transform back
    mat_noise = alr_inv(mat_noise)

    # finally add sparsity
    # Note: there will be no zeros after
    #       using the pseudocount
    if add_missing_at_random:
        n_entries = int(percent_missing * np.count_nonzero(mat_noise))
        i = randint(0, mat_noise.shape[0], n_entries)
        j = randint(0, mat_noise.shape[1], n_entries)
        mat_noise[i, j] = 0
    else:
        i, j = np.nonzero(mat == 0)
        mat_noise[i, j] = 0

    return mat_noise

def simple_blocks(n_samp,
                  n_feat,
                  omic_id='omic1',
                  n_blocks=3,
                  kappa=0.5,
                  sigma=2,
                  overlap=0,
                  pseudocount=1,
                  percent_normal=0.5,
                  percent_random=0.5,
                  random_count=8,
                  add_missing_at_random=True,
                  percent_missing=0.2):

    if n_samp % n_blocks != 0:
        raise ValueError('n_samp must be a multiple of n_blocks')
    block_sim = blocks(n_feat, n_samp, n_blocks,
                       sigma=sigma, overlap=overlap,
                       minval=0, maxval=10000)
    block_sim = block_sim[::-1, ::-1]
    #add omic in front of sample and feature ids
    samp_ids = ['s%i' % i for i in range(n_samp)]
    feat_ids = ['%s_f%i' % (omic_id, i) for i in range(n_feat)]
    block_sim_df = pd.DataFrame(block_sim, samp_ids, feat_ids)
    block_sim_bt = Table(block_sim_df.values.T, block_sim_df.columns, block_sim_df.index)
    table = Table(block_sim_df.values, block_sim_df.index, block_sim_df.columns)
    mat = table.matrix_data.toarray().T
    depths = table.sum('sample').reshape(table.shape[1], -1) * 100
    #add noise
    mat = add_noise(mat, pseudocount, percent_normal,
                    percent_random, random_count,
                    add_missing_at_random, percent_missing)
    #PLN subsampling
    sim_res = poisson_lognormal(mat, depths, kappa=kappa)
    simulation_table = Table(sim_res[0],
                             table.ids()[sim_res[1]],
                             table.ids("observation")[sim_res[2]])
    chunck = n_samp // n_blocks
    mf_sim = pd.DataFrame(['g%i' % (i) for i in range(n_blocks)
                           for j in range(chunck)], samp_ids, ['groups'])

    
    return block_sim_bt, simulation_table, mf_sim