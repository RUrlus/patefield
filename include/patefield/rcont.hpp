/* patefield.hpp -- header only implementation of Patefield algorithm
 * Copyright (C) 2021 Ralph Urlus
 *
 *   Original FORTRAN77 version by WM Patefield.
 *   C++ version by John Burkardt from:
 *   https://people.sc.fsu.edu/~jburkardt/cpp_src/asa159/asa159.html
 *   Modified by Ralph Urlus
 *
 * Reference:
 *
 *   WM Patefield,
 *   Algorithm AS 159:
 *   An Efficient Method of Generating RXC Tables with
 *   Given Row and Column Totals,
 *   Applied Statistics,
 *   Volume 30, Number 1, 1981, pages 91-97.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INCLUDE_PATEFIELD_RCONT_HPP_
#define INCLUDE_PATEFIELD_RCONT_HPP_
#include <cmath>
#include <limits>
#include <random>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <patefield/commons.hpp>

namespace patefield {
namespace details {

/* rcont2 constructs a random two-way contingency table with given sums.
 *
 *   It is possible to specify row and column sum vectors which
 *   correspond to no table at all.  As far as I can see, this routine does
 *   not detect such a case.
 *
 * Parameters
 * ----------
 * n_row : number of rows in the table, must be >= 2;
 * n_col : number of columns in the table, must be >= 2;
 * n_total : the sum of the column or row sums
 * n_row_sums : the row sums, must be > 0;
 * n_col_sums : the column sums, must be > 0;
 * result : pointer to the memory where the table will be stored,
 *          must be of size [n_row * n_col] and is expected to be C-contiguous;
 * factorial_table : pointer to table containing the log factorials,
 *                   can be created using `patefield::create_factorial_table`
 *  rng : instantiated random number generator of PCG RNG family
 */
template<typename T, isInt<T> = true>
void rcont2(
    const T n_row,
    const T n_col,
    const T n_total,
    const T* n_row_sums,
    const T* n_col_sums,
    T* result,
    double* factorial_table,
    pcg64_dxsm& rng
) {
    bool done1;
    bool done2;
    bool lsm;
    bool lsp;
    T i;
    T ia;
    T iap;
    T ib;
    T ic;
    T id;
    T idp;
    T ie;
    T igp;
    T ihp;
    T ii;
    T iip;
    T j;
    T jc;
    T l;
    T m;
    T nll;
    T nlm;
    T nlmp;
    T n_row_sumsl;
    double r;
    double sumprb;
    double x;
    double y;

    // the distribution should be a uniform over the open set (0, 1)
    // uniform_real_distribution is [a, b)
    std::uniform_real_distribution<double> uni_dist(0.0 + std::numeric_limits<double>::epsilon(), 1.0);

    //  Construct a random matrix.
    std::unique_ptr<int64_t[]> jwork_arr(new int64_t[n_col]);
    int64_t* jwork = jwork_arr.get();

    for (i = 0; i < n_col - 1; i++) {
        jwork[i] = n_col_sums[i];
    }

    jc = n_total;

    for (l = 0; l < n_row - 1; l++) {
        n_row_sumsl = n_row_sums[l];
        ia = n_row_sumsl;
        ic = jc;
        jc -= n_row_sumsl;

        for (m = 0; m < n_col - 1; m++) {
            id = jwork[m];
            ie = ic;
            ic = ic - id;
            ib = ie - ia;
            ii = ib - id;

            //  Test for zero entries in matrix.
            if (ie == 0) {
                ia = 0;
                for (j = m; j < n_col; j++) {
                    result[l + j * n_row] = 0;
                }
                break;
            }

            //  Generate a pseudo-random number.
            r = uni_dist(rng);

            //  Compute the conditional expected value of MATRIX(L,M).
            done1 = false;

            for (;;) {
                nlm = static_cast<int>(static_cast<double>(ia * id) / static_cast<double>(ie) + 0.5);
                iap = ia + 1;
                idp = id + 1;
                igp = idp - nlm;
                ihp = iap - nlm;
                nlmp = nlm + 1;
                iip = ii + nlmp;
                x = std::exp(
                    factorial_table[iap - 1] + factorial_table[ib] + factorial_table[ic] + factorial_table[idp - 1]
                    - factorial_table[ie] - factorial_table[nlmp - 1] - factorial_table[igp - 1]
                    - factorial_table[ihp - 1] - factorial_table[iip - 1]);

                if (r <= x) {
                    break;
                }

                sumprb = x;
                y = x;
                nll = nlm;
                lsp = false;
                lsm = false;

                //  Increment entry in row L, column M.
                while (!lsp) {
                    j = (id - nlm) * (ia - nlm);

                    if (j == 0) {
                        lsp = true;
                    } else {
                        nlm += 1;
                        x = x * static_cast<double>(j) / static_cast<double>(nlm * (ii + nlm));
                        sumprb += x;

                        if (r <= sumprb) {
                            done1 = true;
                            break;
                        }
                    }

                    done2 = false;

                    //  Decrement the entry in row L, column M.
                    while (!lsm) {
                        j = nll * (ii + nll);

                        if (j == 0) {
                            lsm = true;
                            break;
                        }

                        nll -= 1;
                        y = y * static_cast<double>(j) / static_cast<double>((id - nll) * (ia - nll));
                        sumprb += y;

                        if (r <= sumprb) {
                            nlm = nll;
                            done2 = true;
                            break;
                        }

                        if (!lsp) {
                            break;
                        }
                    }

                    if (done2) {
                        break;
                    }
                }

                if (done1) {
                    break;
                }

                if (done2) {
                    break;
                }

                r = sumprb * uni_dist(rng);
            }

            result[l + m * n_row] = nlm;
            ia -= nlm;
            jwork[m] = jwork[m] - nlm;
        }
        result[l + (n_col - 1) * n_row] = ia;
    }
    //  Compute the last row.
    for (j = 0; j < n_col - 1; j++) {
        result[n_row - 1 + j * n_row] = jwork[j];
    }
    result[n_row - 1 + (n_col - 1) * n_row] = ib - result[n_row - 1 + (n_col - 2) * n_row];

    return;
}  // rcont2

}  // namespace details
}  // namespace patefield
#endif  // INCLUDE_PATEFIELD_RCONT_HPP_
