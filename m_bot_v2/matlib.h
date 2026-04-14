#pragma once

#include <vector>
#include <iostream>
#include <cmath>




struct csr {
    std::vector<double> values;
    std::vector<int> col_index;
    std::vector<int> row_ptr;
};

struct Matrix {
    std::vector<double> data;
    int nrows, ncols;

    //double& operator()(int i, int j) { return data[i * ncols + j]; }
    //const double& operator()(int i, int j) const { return data[i * ncols + j]; }
    double& operator()(int i, int j) {
        if (i < 0 || j < 0 || i >= nrows || j >= ncols) {
            std::cout << "INDEX ERROR: i=" << i << " j=" << j
                << " nrows=" << nrows << " ncols=" << ncols << std::endl;
            throw std::out_of_range("bad index");
        }

        size_t idx = static_cast<size_t>(i) * ncols + j;

        if (idx >= data.size()) {
            std::cout << "DATA ERROR: idx=" << idx
                << " size=" << data.size()
                << " i=" << i << " j=" << j
                << " ncols=" << ncols << std::endl;
            throw std::out_of_range("bad index");
        }

        return data[idx];
    }
    const double& operator()(int i, int j) const  {
        if (i < 0 || j < 0 || i >= nrows || j >= ncols) {
            std::cout << "INDEX ERROR: i=" << i << " j=" << j
                << " nrows=" << nrows << " ncols=" << ncols << std::endl;
            throw std::out_of_range("bad index");
        }

        size_t idx = static_cast<size_t>(i) * ncols + j;

        if (idx >= data.size()) {
            std::cout << "DATA ERROR: idx=" << idx
                << " size=" << data.size()
                << " i=" << i << " j=" << j
                << " ncols=" << ncols << std::endl;
            throw std::out_of_range("bad index");
        }

        return data[idx];
    }
    size_t size() const { return data.size(); }  // gibt nrows*ncols zurück
    size_t alloc_size() const { return nrows*ncols; }  // gibt nrows*ncols zurück
    int row_size() const { return nrows; }
    int col_size() const { return ncols; }
    std::vector<double> getVec() const { return data; }
};

//-------cpp functions------


Matrix mat_multiplication(const Matrix& A, const Matrix& B) {
    unsigned int rows = A.row_size();
    unsigned int cols = B.col_size();
    if (A.col_size() != B.row_size()) {
        std::cout << "Matrizies doesnt match!\n";
    }
    Matrix result;
    result.ncols = cols;
    result.nrows = rows;
    std::vector<double> data(result.alloc_size(), 0);
    result.data = data;
    // get vecs
    for (int i = 0; i < rows; i++) {//getting indices for result
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < A.col_size(); k++) { //sum loop
                result(i, j) += A(i, k) * B(k, j);    
            }
        }
    }
    return result;
}

void transpose(Matrix& mat) {
    int rows = mat.row_size();
    int cols = mat.col_size();
    Matrix newMat;
    newMat.nrows = cols; //vertauscht wegen transponiert!!!
    newMat.ncols = rows;
    std::vector<double> data(newMat.alloc_size(), 0);
    newMat.data = data;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            newMat(j, i) = mat(i, j);
        }
    }
    mat = newMat;
}

std::vector<double> mat_vec_multiplication(const matlib::Matrix& A, const std::vector<double>& b) {
    if (A.col_size() != b.size()) {
        std::cout << "Mat and Vec doesnt match!\n";
    }
    unsigned int cols = A.col_size();
    unsigned int rows = A.row_size();
    std::vector<double> result(rows, 0);
    for (int i = 0; i < rows; i++) { //result vec

        for (int k = 0; k < cols; k++) { //calculation
            result[i] += A(i, k) * b[k];
        }
    }
    return result;
}

double vec_multiplication(const std::vector<double>& v1, const std::vector<double>& v2) {

    if (v1.size() != v2.size()) {
        throw std::runtime_error("ERROR: Vectors doesnt have compatible dimension");
    }
    int n = v1.size();
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += v1[i] * v2[i];
    }
    return result;

}


std::vector<double> vec_scale(const std::vector<double>& v1, const double& scale) {
    int n = v1.size();
    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = v1[i] * scale;
    }
    return result;
}

std::vector<double> vec_substraction(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw std::runtime_error("ERROR: Vectors doesnt have compatible dimension");
    }
    int n = v1.size();
    std::vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

std::vector<double> vec_addition(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw std::runtime_error("ERROR: Vectors doesnt have compatible dimension");
    }
    int n = v1.size();
    std::vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

std::vector<double> csr_vec_multiplication(const csr& A, const std::vector<double>& v) {
    int n = A.row_ptr.size() - 1;
    std::vector<double> result(n, 0.0);
    if (A.row_ptr.size() - 1 != v.size()) {
        throw std::runtime_error("ERROR: Matrix and Vector doesnt have compatible dimension");
    }
    for (int i = 0; i < n; i++) {
        for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; j++) {
            result[i] += A.values[j] * v[A.col_index[j]];
        }
    }
    return result;

}


csr sparse_csr(const Matrix& A) { //[values, col_index, row_ptr]

    //error handling
    if (A.row_size() != A.col_size()) {
        throw std::runtime_error("ERROR: Matrix is not quadratic");
    }
    int n = A.row_size();
    csr result;
    result.row_ptr.push_back(0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A(i, j) != 0) { //element val is not 0
                result.values.push_back(A(i, j));
                result.col_index.push_back(j);
            }
        }
        result.row_ptr.push_back(result.values.size());
    }


    return result;

}



bool AUTO_linsolve_CG(const Matrix A, const std::vector<double> b, std::vector<double>& x, const long double convergence = 0.0000000000001, const bool printProgress = true) {

    //-------python-->cpp-------
    int n = b.size();

    std::vector<double> r_new(n, 0.0);
    std::vector<double> r_old(n, 0.0);
    std::vector<double> r_glob(n, 0.0);
    std::vector<double> p_old(n, 0.0);
    std::vector<double> p_new(n, 0.0);
    std::vector<double> local(n, 0.0);
    double beta_glob = 0;
    double alpha = 0;
    double convergenceReached = 0;
    //progress

    const int iteration_steps_update = 200;  //1 percent steps
    int progress_counter = 1;
    int iter = 0;
    //double percentage = 0.0;

    //init r_old
    std::vector<double> Ax = mat_vec_multiplication(A, x);
    r_old = vec_substraction(b, Ax);
    p_old = r_old;

    double tempRes1 = 1;
    double tempRes2 = 0;
    double tempRes3 = 1;


    while (tempRes3 / tempRes1 > convergence * convergence) {
        if (iter > n + 500) { //prevent endless loops
            break;
        }

        tempRes1 = 0;
        tempRes2 = 0;
        tempRes3 = 0;

        // local 
        for (int k = 0; k < n; k++) {
            local[k] = 0;
            for (int j = 0; j < n; j++) {
                local[k] += A.getVec()[k * n + j] * p_old[j];
            }
        }

        // alpha berechnen
        for (int k = 0; k < n; k++) {
            tempRes1 += r_old[k] * r_old[k];
            tempRes2 += p_old[k] * local[k];
        }

        if (std::abs(tempRes2) < 1e-18) {
            throw std::runtime_error("Division by zero in CG");
        }

        alpha = tempRes1 / tempRes2;

        // x aktualisieren & r_new berechnen
        for (int k = 0; k < n; k++) {
            x[k] += p_old[k] * alpha;
            r_new[k] = r_old[k] - local[k] * alpha;
            tempRes3 += r_new[k] * r_new[k];
        }

        // Konvergenz prüfen
        if (tempRes3 < convergence * convergence) {
            break;
        }

        // beta berechnen
        beta_glob = tempRes3 / tempRes1;

        // p_new berechnen
        for (int k = 0; k < n; k++) {
            p_new[k] = r_new[k] + p_old[k] * beta_glob;
        }

        // Vektoren tauschen
        std::swap(r_old, r_new);
        std::swap(p_old, p_new);

        //progres
        iter++;
        if (printProgress) {

            if (progress_counter >= iteration_steps_update) {
                std::cout << "\rIterations done: " << iter;
                progress_counter = 0;
            }
            else {
                progress_counter++;
            }

        }
    }
    if (tempRes3 >= convergence * convergence) {
        //py::print("WARNING: convergence criterion has not been met");
        return false;
    }
    else {
        return true;
    }
}

