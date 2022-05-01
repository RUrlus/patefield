/* patefield.cpp -- Public API for Patefield generators.
 * Copyright 2022 R. Urlus
 */
#include <patefield/patefield.hpp>

namespace patefield {

double* create_factorial_table(const int n_total) {
    return details::create_factorial_table(n_total);
}

double* create_factorial_table(const int64_t n_total) {
    return details::create_factorial_table(n_total);
}

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
    uint64_t seed,
    double* factorial_table,
    int* result
) {
    return details::generate_contingency_table<int>(
        n_row,
        n_col,
        n_row_sums,
        n_col_sums,
        n_total,
        seed,
        factorial_table,
        result
    );
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
int64_t* generate_contingency_table(
    const int64_t n_row,
    const int64_t n_col,
    const int64_t* const n_row_sums,
    const int64_t* const n_col_sums,
    int64_t n_total,
    uint64_t seed,
    double* factorial_table,
    int64_t* result
) {
    return details::generate_contingency_table<int64_t>(
        n_row,
        n_col,
        n_row_sums,
        n_col_sums,
        n_total,
        seed,
        factorial_table,
        result
    );
} // generate_contingency_table

/* Generate `n_tables` random two-way contingency tables with given sums.
 *
 *  It is possible to specify row and column sum vectors which
 *  correspond to no table at all.
 *
 * Parameters
 * ----------
 * n_tables : number of tables to generate
 * n_row : number of rows in the table, must be >= 2;
 * n_col : number of columns in the table, must be >= 2;
 * n_row_sums : the row sums, must be > 0;
 * n_col_sums : the column sums, must be > 0;
 * n_total : optional, default = 0, the sum of the column or row sums
 * n_threads : optional, default = 1, the number of threads used to generate
 *             the tables
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
    const size_t n_threads,
    const uint64_t seed,
    double* factorial_table,
    int* result
) {
    return details::generate_contingency_tables<int>(
        n_tables,
        n_row,
        n_col,
        n_row_sums,
        n_col_sums,
        n_total,
        n_threads,
        seed,
        factorial_table,
        result
    );
}

/* Generate `n_tables` random two-way contingency tables with given sums.
 *
 *  It is possible to specify row and column sum vectors which
 *  correspond to no table at all.
 *
 * Parameters
 * ----------
 * n_tables : number of tables to generate
 * n_row : number of rows in the table, must be >= 2;
 * n_col : number of columns in the table, must be >= 2;
 * n_row_sums : the row sums, must be > 0;
 * n_col_sums : the column sums, must be > 0;
 * n_total : optional, default = 0, the sum of the column or row sums
 * n_threads : optional, default = 1, the number of threads used to generate
 *             the tables
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
    const size_t n_threads,
    const uint64_t seed,
    double* factorial_table,
    int64_t* result
) {
    return details::generate_contingency_tables<int64_t>(
        n_tables,
        n_row,
        n_col,
        n_row_sums,
        n_col_sums,
        n_total,
        n_threads,
        seed,
        factorial_table,
        result
    );
}

}  // namespace patefield
