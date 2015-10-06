#pragma once

#include <iostream>
#include <cassert>

#include <Eigen/Dense>

namespace e = Eigen;

template<typename T>
class DataWindow {
  public:
    DataWindow(T *basePtr,
               unsigned int columns,
               unsigned int rows) :
        __basePtr(basePtr),
        __cols(columns),
        __rows(rows) {};
		
    T& operator()(unsigned int _col, unsigned int _row) {
      // Ensure we haven't gone out-of-bounds on memory. There may be cases
      // where we actually want to do that, but we can remove the assertion
      // if that actually happens.
      assert(_row * __cols + _col < __cols * __rows)
      
      // Column-major memory layout per Eigen.  
      return __basePtr[_row * __cols + _col];
    }
    
    void displayMatrix() {
      std::cout << e::Map<e::Matrix<T, e::Dynamic, e::Dynamic, e::RowMajor> >(__basePtr, __rows, __cols).colwise().reverse();
    }

  private:
    T *const           __basePtr;
    const unsigned int __cols;
    const unsigned int __rows;
};
