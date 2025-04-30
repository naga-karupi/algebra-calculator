#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using Matrix = vector<vector<int>>;

Matrix add(const Matrix& A, const Matrix& B) {
    int n = A.size(), m = A[0].size();
    Matrix C(n, vector<int>(m));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

Matrix subtract(const Matrix& A, const Matrix& B) {
    int n = A.size(), m = A[0].size();
    Matrix C(n, vector<int>(m));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

int nextPowerOfTwo(int n) {
    return pow(2, ceil(log2(n)));
}

Matrix resizeMatrix(const Matrix& A, int newSize) {
    Matrix B(newSize, vector<int>(newSize, 0));
    for (int i = 0; i < A.size(); ++i)
        for (int j = 0; j < A[0].size(); ++j)
            B[i][j] = A[i][j];
    return B;
}

Matrix strassen(const Matrix& A, const Matrix& B) {
    int n = A.size();
    if (n == 1)
        return {{A[0][0] * B[0][0]}};

    int k = n / 2;

    Matrix A11(k, vector<int>(k)), A12(k, vector<int>(k)), A21(k, vector<int>(k)), A22(k, vector<int>(k));
    Matrix B11(k, vector<int>(k)), B12(k, vector<int>(k)), B21(k, vector<int>(k)), B22(k, vector<int>(k));

    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];
            A22[i][j] = A[i + k][j + k];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + k];
            B21[i][j] = B[i + k][j];
            B22[i][j] = B[i + k][j + k];
        }

    Matrix M1 = strassen(add(A11, A22), add(B11, B22));
    Matrix M2 = strassen(add(A21, A22), B11);
    Matrix M3 = strassen(A11, subtract(B12, B22));
    Matrix M4 = strassen(A22, subtract(B21, B11));
    Matrix M5 = strassen(add(A11, A12), B22);
    Matrix M6 = strassen(subtract(A21, A11), add(B11, B12));
    Matrix M7 = strassen(subtract(A12, A22), add(B21, B22));

    Matrix C11 = add(subtract(add(M1, M4), M5), M7);
    Matrix C12 = add(M3, M5);
    Matrix C21 = add(M2, M4);
    Matrix C22 = add(subtract(add(M1, M3), M2), M6);

    Matrix C(n, vector<int>(n));
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + k] = C12[i][j];
            C[i + k][j] = C21[i][j];
            C[i + k][j + k] = C22[i][j];
        }

    return C;
}

Matrix multiply(const Matrix& A, const Matrix& B) {
    int n = A.size(), m = B[0].size(), k = B.size();
    int s = nextPowerOfTwo(std::max(n, std::max(m, k)));

    Matrix A_pad = resizeMatrix(A, s);
    Matrix B_pad = resizeMatrix(B, s);
    Matrix C_pad = strassen(A_pad, B_pad);

    // トリミングして元のサイズに戻す
    Matrix C(n, vector<int>(m));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            C[i][j] = C_pad[i][j];

    return C;
}

void printMatrix(const Matrix& mat) {
    for (const auto& row : mat) {
        for (int val : row)
            cout << val << " ";
        cout << endl;
    }
}

int main() {
    Matrix A = {
        {12, 20, 16, 13, 5,	11,	10,	20,	8, 4},
        {4,	11,	14,	10,	6,	2,	4,	18,	9,	15}, 
        {2,	20,	1,	17,	14,	17,	19,	7,	10,	9}, 
        {7,	6,	13,	2,	20,	5,	11,	18,	2,	11},
        {19, 19, 8,	17,	12,	18,	10,	9,	3,	17},
        {8,	3,	13,	10,	21,	17,	3,	18,	14,	3},
        {6,	13,	4,	2,	15,	7,	2,	12,	5,	2},
        {17, 15, 7,	5,	4,	6,	15,	6,	6,	20}, 
        {2,	3, 4, 15, 10, 9, 2,	13,	8,	15},
        {15,	15,	15,	16,	8,	4,	9,	7,	14,	3},
        {17,	12,	8,	6,	13,	18,	6,	14,	13,	9},
        {2,	20,	10,	8,	17,	14,	18,	7,	6,	2},
        {10,	6,	18,	11,	2,	20,	16,	21,	12,	11},
        {11,	18,	3,	17,	5,	9,	14,	11,	11,	6},
        {4,	7,	21,	2,	21,	5,	16,	13,	4,	15},
        {8,	14,	20,	19,	7,	19,	7,	14,	20,	2},
        {20,	11,	9,	8,	11,	6,	2,	10,	15,	19},
        {8,	8,	18,	20,	17,	6,	14,	9,	16,	17},
        {12,	19,	15,	17,	9,	5,	3,	5,	3,	9},
        {5,	6,	2,	14,	20,	9,	3,	20,	15,	5},
        {2,	5,	5,	20,	5,	11,	10,	16,	21,	16},
        {5,	14,	4,	19,	19,	11,	10,	5,	1,	7}

    };

    Matrix B = {
        {7,	1,	8,	18,	5,	19,	5,	8},
        {1,	17,	7,	7,	14,	14,	15,	14},
        {14,	16,	17,	20,	4,	21,	5,	2},
        {3,	3,	4,	6,	17,	21,	16,	2},
        {21,	12,	9,	15,	2,	17,	4,	16},
        {4,	8,	12,	15,	20,	17,	20,	3},
        {9,	4,	14,	9,	12,	19,	12,	5},
        {1,	5,	19,	6,	13,	15,	15,	14},
        {12,	16,	19,	15,	13,	16,	4,	19},
        {6,	21,	13,	18,	16,	10,	20,	17}
    };

    Matrix C = multiply(A, B);

    cout << "Result:" << endl;
    printMatrix(C);

    return 0;
}
