
#include <iostream>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cmath>


long double const delta = 1e-9;

#pragma once

class TMatrix {
private:
    size_t rows_;
    size_t cols_;
    long double** data_;
    
    

public:
    explicit TMatrix(size_t = 0, long double = 0.);
    TMatrix(size_t , size_t , long double = 0.);
    TMatrix(TMatrix const &);
    ~TMatrix();
    long double* operator [] (size_t i);
    long double const* operator [] (size_t i) const;
    TMatrix operator + (TMatrix const &) const;
    TMatrix operator - (TMatrix const &) const;
    TMatrix operator * (TMatrix const &) const;
    TMatrix operator * (int) const;
    TMatrix& operator = (TMatrix const &);
    bool operator == (TMatrix const&) const;
    friend std::ostream& operator << (std::ostream&, TMatrix const& );
    friend std::istream& operator >> (std::istream&, TMatrix&);
    std::tuple<TMatrix, TMatrix, std::vector<std::pair<size_t, size_t>>> LUdecomposition() const;
    long double Determinant() const;
    TMatrix Inverse() const;
    size_t Size() const;
    size_t Get_Rows() const;
    size_t Get_Cols() const;
    long double Norm() const;
    std::pair<size_t, size_t> Change_Without_Zero(size_t i);
    void Swap_Rows(size_t i, size_t j);
    


};



