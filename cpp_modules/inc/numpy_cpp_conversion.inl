#ifndef _NUMPY_CPP_CONVERSION_INL_
#define _NUMPY_CPP_CONVERSION_INL_

#include <vector>
#include <pybind11/numpy.h>

namespace py = pybind11;

template <typename T>
std::vector<T> numpy_to_vector_1d(const py::array_t<T, py::array::c_style | py::array::forcecast> &array)
{
    py::buffer_info buffer_info{array.request()};
    assert(buffer_info.ndim == 1U);
    return std::vector<T>(array.data(), array.data() + array.size());
}

template <typename T>
std::vector<std::vector<T>> numpy_to_vector_2d(const py::array_t<T, py::array::c_style | py::array::forcecast> &array)
{
    py::buffer_info buffer_info{array.request()};
    assert(buffer_info.ndim == 2U);
    const auto data_pointer = static_cast<T *>(buffer_info.ptr);

    size_t num_of_rows = buffer_info.shape[0];
    size_t num_of_columns = buffer_info.shape[1];

    std::vector<std::vector<T>> result;
    for (auto i = 0U; i < num_of_rows; ++i)
    {
        result.emplace_back(data_pointer + i * num_of_columns, data_pointer + (i + 1) * num_of_columns);
    }

    return result;
}

template <typename T, uint32_t N>
std::vector<std::array<T, N>> numpy_to_vector_of_arrays(const py::array_t<T, py::array::c_style | py::array::forcecast> &array)
{
    const auto buffer_info{array.request()};
    assert(buffer_info.ndim == 2U);
    const auto data_pointer{static_cast<T *>(buffer_info.ptr)};

    const auto num_of_rows{buffer_info.shape[0]};
    const auto num_of_columns{N};

    std::vector<std::array<T, N>> result{};
    result.reserve(num_of_rows);
    for (auto row{0U}; row < num_of_rows; ++row)
    {
        std::array<T, N> next_row{};
        for (auto column{0U}; column < num_of_columns; ++column)
        {
            next_row[column] = data_pointer[row * num_of_columns + column];
        }
        result.emplace_back(next_row);
    }

    return result;
}

template <typename T>
py::array_t<T, py::array::c_style | py::array::forcecast> vector_to_numpy_1d(const std::vector<T> &data)
{
    const auto elements{data.size()};
    const std::vector<ssize_t> shape{static_cast<ssize_t>(elements)};

    py::array_t<T, py::array::c_style | py::array::forcecast> result(shape);

    const auto resultPtr{result.mutable_data()};

    // Copy the data from the vector of vectors to the NumPy array
    for (auto i{0U}; i < elements; ++i)
    {
        resultPtr[i] = data[i];
    }

    return result;
}

template <typename T>
py::array_t<T, py::array::c_style | py::array::forcecast> vector_to_numpy_2d(const std::vector<std::vector<T>> &data)
{
    const auto rows{data.size()};
    const auto cols{(rows > 0) ? data[0].size() : 0};
    const std::vector<ssize_t> shape{rows, cols};

    py::array_t<T, py::array::c_style | py::array::forcecast> result(shape);

    const auto resultPtr{result.mutable_data()};

    // Copy the data from the vector of vectors to the NumPy array
    for (auto i{0U}; i < rows; ++i)
    {
        for (auto j{0U}; j < cols; ++j)
        {
            resultPtr[i * cols + j] = data[i][j];
        }
    }

    return result;
}

template <typename T>
py::array_t<T, py::array::c_style | py::array::forcecast> vector_of_pairs_to_numpy(const std::vector<std::pair<T, T>> &vector_of_pairs)
{
    const auto size{vector_of_pairs.size()};
    const std::vector<uint64_t> shape{size, 2U};
    py::array_t<T, py::array::c_style | py::array::forcecast> numpy_array(shape);

    const auto numpy_data{numpy_array.mutable_data()};

    memcpy(numpy_data, vector_of_pairs.data(), size * sizeof(std::pair<T, T>));
    return numpy_array;
}

inline double profile_function(const double argument, const double alpha)
{
    return argument <= alpha ? 1. / (2. * alpha) : 0.;
}

#endif
