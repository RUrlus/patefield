/* patefield.hpp -- Patefield generators.
 * Copyright 2022 R. Urlus
 */

#ifndef INCLUDE_PATEFIELD_PATEFIELD_HPP_
#define INCLUDE_PATEFIELD_PATEFIELD_HPP_

#if defined(PATEFIELD_HAS_OPENMP_SUPPORT)
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <utility>
#include <stdexcept>
#include <type_traits>

#include <patefield/commons.hpp>
#include <patefield/rcont.hpp>

namespace patefield {
namespace details {

template<typename T, isInt<T> = true>
inline int64_t check_inputs(const T n_row, const T n_col, const T* n_row_sums, const T* n_col_sums) {
    if (n_row <= 1) {
        throw InputError("patefield: number of rows is less than 2.\n");
    }

    if (n_col <= 1) {
        throw InputError("patefield: number of columns is less than 2.\n");
    }

    for (T i = 0; i < n_row; i++) {
        if (n_row_sums[i] <= 0) {
            throw InputError("patefield: an entry in the row sum vector is not positive.\n");
        }
    }

    for (T j = 0; j < n_col; j++) {
        if (n_col_sums[j] <= 0) {
            throw InputError("patefield: an entry in the column sum vector is not positive.\n");
        }
    }

    int64_t n_total = std::accumulate(n_col_sums, n_col_sums + n_col, 0);
    if (n_total != std::accumulate(n_row_sums, n_row_sums + n_row, 0)) {
        throw InputError("patefield: the row and column sum vectors do not have the same sum.\n");
    }
    return n_total;
}  // check_inputs

template <typename T, isInt<T> = true>
double* create_factorial_table(const T n_total) {
    auto table = reinterpret_cast<double*>(std::malloc(static_cast<size_t>(n_total + 1) * sizeof(double)));
    //  Calculate log-factorials.
    double x = 0.0;
    table[0] = 0.0;
    for (T i = 1; i <= n_total; i++) {
        x += std::log(static_cast<double>(i));
        table[i] = x;
    }
    return table;
}  // create_factorial_table

/* Generate a random two-way contingency table with given sums.
 *
 *  It is possible to specify row and column sum vectors which
 *  correspond to no table at all.
 *
 * Parameters
 * ----------
 * n_row : number of rows in the table, must be >= 2;
 * n_col : number of columns in the table, must be >= 2;
 * n_row_sums : the row sums, must be > 0;
 * n_col_sums : the column sums, must be > 0;
 * n_total : optional, default = 0, the sum of the column or row sums
 * seed : optional, a seed for the random number generator,
 *        default = 0, i.e. drawn by random_device.
 * factorial_table : optional, pointer to table containing the log factorials,
 *                   can be created using `patefield::create_factorial_table`
 *                   Otherwise the factorial_table will be computed.
 * result : optional, pointer to the memory where the table will be stored,
 *          must be of size [n_row * n_col] and is expected to be C-contiguous
 *
 *  Returns
 *  -------
 *  result : pointer to the memory where the table has been stored
 */
template<typename T, isInt<T> = true>
inline T* generate_contingency_table(
    const T n_row,
    const T n_col,
    const T* const n_row_sums,
    const T* const n_col_sums,
    int64_t n_total = 0,
    uint64_t seed = 0,
    double* factorial_table = nullptr,
    T* result = nullptr
) {
    if (n_total == 0) {
        n_total = details::check_inputs<T>(n_row, n_col, n_row_sums, n_col_sums);
    }
    if (!factorial_table) {
        factorial_table = details::create_factorial_table<T>(n_total);
    }
    if (!result) {
        result = reinterpret_cast<T*>(std::malloc(n_row * n_col * sizeof(T)));
        if (!result) throw std::bad_alloc();
    }
    pcg64_dxsm rng;
    if (seed == 0) {
        pcg_seed_seq seed_source;
        rng.seed(seed_source);
    } else {
        rng.seed(seed);
    }

    rcont2<T>(
        n_row, n_col, n_total, n_row_sums, n_col_sums, result, factorial_table, rng
    );
    return result;
} // generate_contingency_table

/* Generate a random two-way contingency table with given sums.
 *
 *  It is possible to specify row and column sum vectors which
 *  correspond to no table at all.
 *
 * Parameters
 * ----------
 * n_row : number of rows in the table, must be >= 2;
 * n_col : number of columns in the table, must be >= 2;
 * n_row_sums : the row sums, must be > 0;
 * n_col_sums : the column sums, must be > 0;
 * n_total : optional, default = 0, the sum of the column or row sums
 * seed : optional, a seed for the random number generator,
 *        default = 0, i.e. drawn by random_device.
 * factorial_table : optional, pointer to table containing the log factorials,
 *                   can be created using `patefield::create_factorial_table`
 *                   Otherwise the factorial_table will be computed.
 * result : optional, pointer to the memory where the table will be stored,
 *          must be of size [n_row * n_col] and is expected to be C-contiguous
 *
 *  Returns
 *  -------
 *  result : pointer to the memory where the table has been stored
 */
template<typename T, isInt<T> = true>
inline T* generate_contingency_tables(
    const size_t n_tables,
    const T n_row,
    const T n_col,
    const T* n_row_sums,
    const T* n_col_sums,
    int64_t n_total,
    const size_t n_threads = 1,
    const uint64_t seed = 0,
    double* factorial_table = nullptr,
    T* result = nullptr
) {
    if (n_total == 0) {
        n_total = details::check_inputs<T>(n_row, n_col, n_row_sums, n_col_sums);
    }
    if (!factorial_table) {
        factorial_table = details::create_factorial_table<T>(n_total);
    }
    const size_t block_size = static_cast<size_t>(n_row) * static_cast<size_t>(n_col);
    if (!result) {
        result = reinterpret_cast<T*>(std::malloc(block_size * sizeof(T) * n_tables));
        if (!result) throw std::bad_alloc();
    }

    pcg64_dxsm global_rng;
    if (seed == 0) {
        pcg_seed_seq seed_source;
        global_rng.seed(seed_source);
    } else {
        global_rng.seed(seed);
    }

    #pragma omp parallel num_threads(n_threads) shared(global_rng, n_row, n_col, n_row_sums, n_col_sums, factorial_table, result)
    {
        pcg64_dxsm rng = global_rng;
        rng.set_stream(omp_get_thread_num() + 1);

        #pragma omp for
        for (size_t i = 0; i < n_tables; i++) {
            rcont2<T>(
                n_row, n_col, n_total, n_row_sums, n_col_sums, result + (i * block_size), factorial_table, rng
            );
        }
    }  // pragma parallel

    return result;
} // generate_contingency_tables

}  // namespace details

/* Generate a random two-way contingency table with given sums.
 *
 *  It is possible to specify row and column sum vectors which
 *  correspond to no table at all.
 *
 * Parameters
 * ----------
 * n_row : number of rows in the table, must be >= 2;
 * n_col : number of columns in the table, must be >= 2;
 * n_row_sums : the row sums, must be > 0;
 * n_col_sums : the column sums, must be > 0;
 * n_total : optional, default = 0, the sum of the column or row sums
 * seed : optional, a seed for the random number generator,
 *        default = 0, i.e. drawn by random_device.
 * factorial_table : optional, pointer to table containing the log factorials,
 *                   can be created using `patefield::create_factorial_table`
 *                   Otherwise the factorial_table will be computed.
 * result : optional, pointer to the memory where the table will be stored,
 *          must be of size [n_row * n_col] and is expected to be C-contiguous
 *
 *  Returns
 *  -------
 *  result : pointer to the memory where the table has been stored
 */
int* generate_contingency_table(
    const int n_row,
    const int n_col,
    const int* const n_row_sums,
    const int* const n_col_sums,
    int64_t n_total,
    uint64_t seed = 0,
    double* factorial_table = nullptr,
    int* result = nullptr
);

/* Generate a random two-way contingency table with given sums.
 *
 *  It is possible to specify row and column sum vectors which
 *  correspond to no table at all.
 *
 * Parameters
 * ----------
 * n_row : number of rows in the table, must be >= 2;
 * n_col : number of columns in the table, must be >= 2;
 * n_row_sums : the row sums, must be > 0;
 * n_col_sums : the column sums, must be > 0;
 * n_total : optional, default = 0, the sum of the column or row sums
 * seed : optional, a seed for the random number generator,
 *        default = 0, i.e. drawn by random_device.
 * factorial_table : optional, pointer to table containing the log factorials,
 *                   can be created using `patefield::create_factorial_table`
 *                   Otherwise the factorial_table will be computed.
 * result : optional, pointer to the memory where the table will be stored,
 *          must be of size [n_row * n_col] and is expected to be C-contiguous
 *
 *  Returns
 *  -------
 *  result : pointer to the memory where the table has been stored
 */
int64_t* generate_contingency_table(
    const int64_t n_row,
    const int64_t n_col,
    const int64_t* const n_row_sums,
    const int64_t* const n_col_sums,
    int64_t n_total,
    uint64_t seed = 0,
    double* factorial_table = nullptr,
    int64_t* result = nullptr
);

/* Generate a random two-way contingency table with given sums.
 *
 *  It is possible to specify row and column sum vectors which
 *  correspond to no table at all.
 *
 * Parameters
 * ----------
 * n_row : number of rows in the table, must be >= 2;
 * n_col : number of columns in the table, must be >= 2;
 * n_row_sums : the row sums, must be > 0;
 * n_col_sums : the column sums, must be > 0;
 * n_total : optional, default = 0, the sum of the column or row sums
 * seed : optional, a seed for the random number generator,
 *        default = 0, i.e. drawn by random_device.
 * factorial_table : optional, pointer to table containing the log factorials,
 *                   can be created using `patefield::create_factorial_table`
 *                   Otherwise the factorial_table will be computed.
 * result : optional, pointer to the memory where the table will be stored,
 *          must be of size [n_row * n_col] and is expected to be C-contiguous
 *
 *  Returns
 *  -------
 *  result : pointer to the memory where the table has been stored
 */
int* generate_contingency_tables(
    const size_t n_tables,
    const int n_row,
    const int n_col,
    const int* n_row_sums,
    const int* n_col_sums,
    int64_t n_total,
    const size_t n_threads = 1,
    const uint64_t seed = 0,
    double* factorial_table = nullptr,
    int* result = nullptr
);

/* Generate a random two-way contingency table with given sums.
 *
 *  It is possible to specify row and column sum vectors which
 *  correspond to no table at all.
 *
 * Parameters
 * ----------
 * n_row : number of rows in the table, must be >= 2;
 * n_col : number of columns in the table, must be >= 2;
 * n_row_sums : the row sums, must be > 0;
 * n_col_sums : the column sums, must be > 0;
 * n_total : optional, default = 0, the sum of the column or row sums
 * seed : optional, a seed for the random number generator,
 *        default = 0, i.e. drawn by random_device.
 * factorial_table : optional, pointer to table containing the log factorials,
 *                   can be created using `patefield::create_factorial_table`
 *                   Otherwise the factorial_table will be computed.
 * result : optional, pointer to the memory where the table will be stored,
 *          must be of size [n_row * n_col] and is expected to be C-contiguous
 *
 *  Returns
 *  -------
 *  result : pointer to the memory where the table has been stored
 */
int64_t* generate_contingency_tables(
    const size_t n_tables,
    const int64_t n_row,
    const int64_t n_col,
    const int64_t* n_row_sums,
    const int64_t* n_col_sums,
    int64_t n_total,
    const size_t n_threads = 1,
    const uint64_t seed = 0,
    double* factorial_table = nullptr,
    int64_t* result = nullptr
);

}  // namespace patefield
#endif  // INCLUDE_PATEFIELD_PATEFIELD_HPP_
