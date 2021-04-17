#include "matrix.hpp"
#include "calculation.hpp"



TMatrix::TMatrix(size_t n, long double elem) : TMatrix(n, n, elem) {}



TMatrix::TMatrix(size_t i, size_t j, long double elem) {
    rows_ = i;
    cols_ = j;
    data_ = new long double* [rows_];
    for (size_t c = 0;  c != rows_; c++)
        data_[c] = new long double [cols_];

    for (size_t i = 0; i != rows_; i++) {
        for (size_t j = 0; j != cols_; j++)
            data_[i][j] = elem;
    }
}

TMatrix::TMatrix(TMatrix const & other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    data_ = new long double* [rows_];
    for (size_t i = 0; i != rows_; i++)
        data_[i] = new long double [cols_];
    
    for (size_t i = 0; i != rows_; i++) {
        for (size_t j = 0; j != cols_; j++)
            data_[i][j] = other.data_[i][j];
    }
}


TMatrix::~TMatrix() {
    for (size_t i = 0; i != rows_; i++)
        delete [] data_[i];
    
    delete [] data_;
}


long double* TMatrix::operator [] (size_t i) {
    return data_[i];
}

long double const* TMatrix::operator [] (size_t i) const {
    return data_[i];
}


TMatrix TMatrix::operator+ (TMatrix const & other) const {
    assert(rows_ == other.rows_ and cols_ == other.cols_);
    TMatrix tmp = *this;
    for (size_t i = 0; i != rows_; i++) {
        for (size_t j = 0; j != cols_; j++)
            tmp.data_[i][j] += other.data_[i][j];
    }
    return tmp;
}


TMatrix TMatrix::operator- (TMatrix const & other) const{
    assert(rows_ == other.rows_ and cols_ == other.cols_);
    TMatrix tmp = *this;
    for (size_t i = 0; i != rows_; i++) {
        for (size_t j = 0; j != cols_; j++)
            tmp.data_[i][j] -= other.data_[i][j];
    }
    return tmp;
}

TMatrix TMatrix::operator* (TMatrix const & other) const{
    assert(cols_ == other.rows_);
    size_t l = rows_;
    size_t n = other.cols_;
    size_t m = cols_;
    TMatrix res(l, n);

    for (size_t i = 0; i != l; i++) {
        for (size_t j = 0; j != n; j++) {
            long double tmp = 0.;
            for (size_t k = 0; k != m; k++)
                tmp += data_[i][k] * other.data_[k][j];
            res.data_[i][j] = tmp;
        }
    }
    return res;
}


TMatrix TMatrix::operator* (int c) const {
    TMatrix res = *this;
    for (size_t i = 0; i != rows_; i++) {
        for (size_t j = 0; i != cols_; j++)
            res.data_[i][j] *= c;
    }
    return res;
}


TMatrix& TMatrix::operator= (TMatrix const& other) {
    this->~TMatrix();

    rows_ = other.rows_;
    cols_ = other.cols_;

    data_ = new long double* [rows_];
    for (size_t i = 0; i != rows_; i++)
        data_[i] = new long double [cols_];

    for (size_t i = 0; i != rows_; i++) {
        for (size_t j = 0; j != cols_; j++)
            data_[i][j] = other.data_[i][j];
    }

    return *this;




}

std::ostream& operator << (std::ostream& out, TMatrix const& m) {
    out << std::setprecision(5) << std::fixed;
    for (size_t i = 0; i != m.rows_; i++){
        for (size_t j = 0; j != m.cols_; j++) {
            out << m.data_[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}




std::istream& operator >> (std::istream& in, TMatrix & m) {
    for (size_t i = 0; i != m.rows_; i++) {
        for (size_t j = 0; j != m.cols_; j++)
            in >> m.data_[i][j];
    }
    return in;
}




bool TMatrix::operator== (TMatrix const& other) const {
    if (rows_ != other.rows_ or cols_ != other.cols_)
        return false;
    
    for (size_t i = 0; i != rows_; i++) {
        for (size_t j = 0; j != cols_; j++) {
            if (std::abs(data_[i][j] - other.data_[i][j]) > delta)
                return false;
        }
    }
    return true;
}

std::tuple<TMatrix, TMatrix, std::vector<std::pair<size_t, size_t>>> TMatrix::LUdecomposition() const {
    assert(rows_ == cols_);
    size_t n = rows_;
    TMatrix U = *this;
    TMatrix L(n);
    std::vector<std::pair<size_t, size_t>> P;
    for (size_t i = 0; i != n; i++)
        L.data_[i][i] = 1;
    
    for (size_t i = 0; i != n - 1; i++) {
        if (U.data_[i][i] == 0) {
            auto p = U.Change_Without_Zero(i);
            P.push_back(p);
        }
        for (size_t j = i + 1; j != n; j++) {
            long double t = U.data_[j][i] / U.data_[i][i];
            for (size_t k = i; k != n; k++) {
                U.data_[j][k] -= t * U.data_[i][k];
            }
            L.data_[j][i] = t;
        }
    }

    return  std::make_tuple(L, U, P);
}

long double TMatrix::Determinant() const {
    auto t = this->LUdecomposition();
    TMatrix U = std::get<1>(t);
    int p = std::get<2>(t).size();
    long double det = 1;
    for (size_t i = 0; i != cols_; i++)
        det *= U.data_[i][i];
    return std::pow(-1, p) * det;
}


TMatrix TMatrix::Inverse() const {
    auto [L, U, P] = this->LUdecomposition();
    return LU_Inverse_Matrix(L, U, P);
}

size_t TMatrix::Size() const {
    assert(rows_ == cols_);
    return rows_;
}

size_t TMatrix::Get_Rows() const {
    return rows_;
}

size_t TMatrix::Get_Cols() const {
    return cols_;
}


long double  TMatrix::Norm() const {
    long double norm = -1.;
    
    for (size_t i = 0; i != rows_; i++) {
        long double sum = 0.;
        for (size_t j = 0; j != cols_; j++)
            sum += std::abs(data_[i][j]);
        norm = sum > norm ? sum : norm;
    }
    return norm;
}


std::pair<size_t, size_t> TMatrix::Change_Without_Zero(size_t i) {
    size_t j = 0;
    for (size_t j = i + 1; j < Size(); j++)
        if (data_[j][i] != 0.) {
            Swap_Rows(i ,j);
            break;
        }
    return std::make_pair(i, j);    
}

void TMatrix::Swap_Rows(size_t i, size_t j) {
    long double * tmp = new long double [cols_];
    for (size_t k = 0; k != cols_; k++)
        tmp[k] = data_[i][k];
    for (size_t k = 0; k != cols_; k++)
        data_[i][k] = data_[j][k];
    for (size_t k = 0; k != cols_; k++)
        data_[j][k] = tmp[k];
    delete [] tmp;
}