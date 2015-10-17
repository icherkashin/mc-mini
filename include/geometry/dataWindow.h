#pragma once

#include <iostream>
#include <cassert>

#include <Eigen/Dense>

namespace e = Eigen;

/** @brief A data array wrapper class
 *
 *  The DataWindow class was designed as a simple and lightweight wrapper
 *  around the raw scalar data arrays which performs bounds-checking on
 *  access.
 */
template<typename T>
class DataWindow {
  public:
    /** @brief Constructs a DataWindow around a contiguous region in memory
     *
     *  Constructs a DataWindow around a contiguous block of memory of size
     *  **columns** \* **rows**, pointed to by **basePtr**.
     */
    DataWindow(T *basePtr,
               unsigned int columns,
               unsigned int rows) :
        __basePtr(basePtr),
        __cols(columns),
        __rows(rows) {};

    /** @brief Subscripts the memory region wrapped by DataWindow, with bounds
     *         -checking
     *
     *  Accesses the scalar value at the position (**col**, **row**) by
     *  calculating the offset from the base data pointer.
     */
    T& operator()(unsigned int _col, unsigned int _row) {
      // Ensure we haven't gone out-of-bounds on memory. There may be cases
      // where we actually want to do that, but we can remove the assertion
      // if that actually happens.
      assert(_row * __cols + _col < __cols * __rows);

      // Column-major memory layout per Eigen.
      return __basePtr[_row * __cols + _col];
    }

    /** @brief Displays the data array wrapped by DataWindow */
    const std::string displayMatrix() {
      /*  TODO: Actually return a string version of the array rather than
       *  outputting it here.
       */
      std::cout << e::Map<e::Matrix<T, e::Dynamic, e::Dynamic, e::RowMajor> >(__basePtr, __rows, __cols).colwise().reverse();

      return "";
    }

  private:
    T *const           __basePtr;
    const unsigned int __cols;
    const unsigned int __rows;
};
