#ifndef _NUMPY_CPP_CONVERSION_H_
#define _NUMPY_CPP_CONVERSION_H_

#include <pybind11/numpy.h>

namespace py = pybind11;

template <typename T>
std::vector<T> numpy_to_vector_1d(const py::array_t<T, py::array::c_style | py::array::forcecast> &array);

template <typename T>
std::vector<std::vector<T>> numpy_to_vector_2d(const py::array_t<T, py::array::c_style | py::array::forcecast> &array);

template <typename T, uint32_t D>
std::vector<std::array<T, D>> numpy_to_vector_of_arrays(const py::array_t<T, py::array::c_style | py::array::forcecast> &array);

template <typename T>
py::array_t<T, py::array::c_style | py::array::forcecast> vector_to_numpy_1d(const std::vector<T> &vector);

template <typename T>
py::array_t<T, py::array::c_style | py::array::forcecast> vector_to_numpy_2d(const std::vector<std::vector<T>> &vector_of_vectors);

template <typename T>
py::array_t<T, py::array::c_style | py::array::forcecast> vector_of_pairs_to_numpy(const std::vector<std::pair<T, T>> &vector_of_pairs);

#include "numpy_cpp_conversion.inl"

#endif
