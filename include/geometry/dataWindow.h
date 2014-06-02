#pragma once

#include <string>
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

template<class T>
class DataWindow {
  public:
    DataWindow (T* _basePtr, unsigned int _columns, unsigned int _rows);
    T& operator() (unsigned int _col, unsigned int _row);
    const string displayMatrix();
  private:
    T*           __basePtr;
    unsigned int __cols;
    unsigned int __rows;
};

// Column-major data window, as per both Eigen and Visit.
template<class T>
DataWindow<T>::DataWindow (T* _basePtr, unsigned int _columns, unsigned int _rows) :
    __basePtr (_basePtr),
    __cols (_columns),
    __rows (_rows) {};

template<class T>
T& DataWindow<T>::operator() (unsigned int _col, unsigned int _row) {
  return __basePtr[_row * __cols + _col];
}

template<class T>
const string DataWindow<T>::displayMatrix() {
  std::cout << Eigen::Map<Matrix<T, Dynamic, Dynamic, RowMajor> >(__basePtr, __rows, __cols).colwise().reverse();
  return "";
}
