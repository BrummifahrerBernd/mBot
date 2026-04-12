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

    double& operator()(int i, int j) { return data[i * ncols + j]; }
    const double& operator()(int i, int j) const { return data[i * ncols + j]; }
    size_t size() const { return data.size(); }  // gibt nrows*ncols zurück
    int row_size() const { return nrows; }
    int col_size() const { return ncols; }
    std::vector<double> getVec() const { return data; }
};

//-------cpp functions------


Matrix symm_mat_multiplication(const Matrix& A, const Matrix& B) {
    int n = A.col_size();
    Matrix result;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result(i, j) = 0.0;
            for (int k = 0; k < 3; k++) {
                result(i, j) += A(i, j) * B(i, j);
            }
        }
    }
    return result;
}

void transpose(Matrix& mat) {
    int n = mat.row_size();
    Matrix newMat;
    newMat.nrows = n;
    newMat.ncols = n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            newMat(j, i) = mat(i, j);
        }
    }
    mat = newMat;
}


std::vector<double> mat_vec_multiplication(const std::vector<double>& A, const std::vector<double>& b) {
    float n = b.size();
    std::vector<double> result(b.size(), 0);
    for (int k = 0; k < b.size(); k++) {

        for (int j = 0; j < b.size(); j++) {
            result[k] += A[k * n + j] * b[j];
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



bool AUTO_linsolve_CG(const Matrix A, const std::vector<double> b, const long double convergence = 0.0000000000001, const bool printProgress = true, std::vector<double>& x) {

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
    std::vector<double> Ax = mat_vec_multiplication(A.getVec(), x);
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


