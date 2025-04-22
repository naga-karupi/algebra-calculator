#include <iostream>
#include <Eigen/Dense>

template<int N>
std::pair<Eigen::Matrix<double, N, N>, Eigen::Matrix<double, N, N>>
qr_householder(const Eigen::Matrix<double, N, N>& A) {
    Eigen::Matrix<double, N, N> Q = Eigen::Matrix<double, N, N>::Identity();
    Eigen::Matrix<double, N, N> R = A;

    for(int i = 0; i < N; ++i) {
        // 部分ベクトル x = R(i:N-1, i)
        Eigen::VectorXd x = R.block(i, i, N-i, 1);

        // α = -sign(x0) * ||x||
        double sigma = x.norm();
        if (sigma == 0.0) continue;
        double alpha = (x(0) >= 0 ? -sigma : sigma);

        // 反射ベクトル v を作成
        Eigen::VectorXd v = x;
        v(0) -= alpha;
        double beta = v.squaredNorm();
        if (beta == 0.0) continue;

        // 部分反射行列 H_sub = I - 2/β * v v^T
        Eigen::MatrixXd H_sub = Eigen::MatrixXd::Identity(N-i, N-i)
                              - (2.0 / beta) * (v * v.transpose());

        // 全体反射行列 H を単位行列に埋め込む
        Eigen::Matrix<double, N, N> H = Eigen::Matrix<double, N, N>::Identity();
        H.block(i, i, N-i, N-i) = H_sub;

        // R, Q を更新
        R = H * R;
        Q = Q * H.transpose();
    }

    return {Q, R};
}

template<int N>
std::pair<Eigen::Matrix<double, N, N>, Eigen::Matrix<double, N, N>>
    calclate_eigenvalues_and_eigenvector(const Eigen::Matrix<double, N, N>& A)
{

    Eigen::Matrix<double, N, N> A_k = A;
    Eigen::Matrix<double, N, N> pre_A_k = Eigen::Matrix<double, N, N>::Zero();
    Eigen::Matrix<double, N, N> X_k = Eigen::Matrix<double, N, N>::Identity();
    Eigen::Matrix<double, N, N> ret_eigenvalues = Eigen::Matrix<double, N, N>::Identity();
    Eigen::Matrix<double, N, N> ret_eigenvectors = Eigen::Matrix<double, N, N>::Identity();

    // QR法による固有値計算
    while(true)
    {
        auto [Q, R] = qr_householder<N>(A_k);
        A_k = R * Q;
        X_k = X_k * Q;
        pre_A_k = A_k;


        // 収束判定
        bool convergence_flag = true;
        for(int i = 0; i < N; i++)
        {
            if(A_k(i, i) - pre_A_k(i, i) > 1e-8)
            {
                convergence_flag = false;
                break;
            }
        }

        if(convergence_flag)
        {
            ret_eigenvalues = A_k;
            ret_eigenvectors = X_k;
            break;
        }
    }

    return {ret_eigenvalues, ret_eigenvectors};
}

int main() {
    constexpr int N = 4;
    Eigen::Matrix<double, N, N> A;
    A << 1, 2, 3, 4,
         4, 5, 6, 6,
         7, 8,10, 3,
         2, 4, 3, 1;

    auto [eigenvalues, eigenvector] 
        = calclate_eigenvalues_and_eigenvector<N>(A);
    std::cout << "Eigenvalues:\n" << eigenvalues << "\n";
    std::cout << "Eigenvectors:\n" << eigenvector << "\n";

    return 0;
}
