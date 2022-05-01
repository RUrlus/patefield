/* patefield/commons.hpp -- Utility functions and definitions
 * Copyright (C) 2022 Ralph Urlus
 */
#ifndef INCLUDE_PATEFIELD_COMMONS_HPP_
#define INCLUDE_PATEFIELD_COMMONS_HPP_

#include <random>  // random_device

#include <pcg_random.hpp>
#include <pcg_extras.hpp>

namespace patefield {

template<typename T>
using isInt = std::enable_if_t<std::is_integral<T>::value, bool>;

class InputError : public std::runtime_error {
 public:
    explicit InputError(const char* what) : std::runtime_error(what) {}
};

typedef pcg_engines::setseq_dxsm_128_64 pcg64_dxsm;
typedef pcg_extras::seed_seq_from<std::random_device> pcg_seed_seq;

}  // namespace patefield
#endif  // INCLUDE_PATEFIELD_COMMONS_HPP_
