#include <iostream>
#include <complex>

#include <random>

#define USE_EIGEN 0
#define USE_MY_MATRIX 1

#if USE_EIGEN
#include <Eigen/Dense>
using c_double = std::complex<double>;
constexpr int X = -1;


template<int N>
Eigen::Matrix<c_double, N, N> hessenberg(const Eigen::Matrix<c_double, N, N>& A) 
{
    Eigen::Matrix<c_double, N, N> Q = Eigen::Matrix<c_double, N, N>::Identity();
    Eigen::Matrix<c_double, N, N> R = A;
    Eigen::Matrix<c_double, N, N> A_star = A;

    for(int i = 0; i < N-2; ++i) 
    {
        Eigen::Matrix<c_double, N, N> H = create_householder<N>(A_star, i);
        if (H.isZero()) continue;
        // R, Q を更新
        A_star = H * A_star * H.adjoint();
        R = H * R;
        Q = Q * H.adjoint();
    }

    return A_star;
}

template<int N>
Eigen::Matrix<c_double, N, N> 
  create_householder(const Eigen::Matrix<c_double, N, N>& A, size_t n) 
  {
    // 部分ベクトル x = R(i:N-1, i)
    Eigen::Vector<c_double, X> x = A.block(n, n, N-n, 1);

    // 実ノルムは double で取る
    double sigma = x.norm();  
    if (sigma == 0.0) return Eigen::Matrix<c_double, N, N>::Zero();

    // x(0) が 0 だときは位相不定になるのでガード
    c_double sign = std::abs(x(0)) == 0.0 
                    ? c_double(1.0, 0.0) 
                    : x(0) / std::abs(x(0));

    // 複素 α を作成
    c_double alpha = - sign * sigma;

    // 反射ベクトル v を作成
    Eigen::Vector<c_double, X> v = x;
    v(0) -= alpha;
    c_double beta = v.squaredNorm();
    if (beta == 0.0) return Eigen::Matrix<c_double, N, N>::Zero();

    // 部分反射行列 H_sub = I - 2/β * v v^T
    Eigen::Matrix<c_double, X, X> H_sub = Eigen::MatrixXcd::Identity(N-n, N-n)
                            - (2.0 / beta) * (v * v.adjoint());

    // 全体反射行列 H を単位行列に埋め込む
    Eigen::Matrix<c_double, N, N> H = Eigen::Matrix<c_double, N, N>::Identity();
    H.block(n, n, N-n, N-n) = H_sub;

    return H;
}

template<int N>
std::pair<Eigen::Matrix<c_double, N, N>, Eigen::Matrix<c_double, N, N>>
qr_householder(const Eigen::Matrix<c_double, N, N>& A) 
{
    Eigen::Matrix<c_double, N, N> Q = Eigen::Matrix<c_double, N, N>::Identity();
    Eigen::Matrix<c_double, N, N> R = A;

    for(int i = 0; i < N; ++i) 
    {

        Eigen::Matrix<c_double, N, N> H = create_householder<N>(R, i);
        if (H.isZero()) continue;
        // R, Q を更新
        R = H * R;
        Q = Q * H.adjoint();
    }

    return {Q, R};
}



template<int N>
std::pair<Eigen::Matrix<c_double, N, N>, Eigen::Matrix<c_double, N, N>>
    calclate_eigenvalues_and_eigenvector(const Eigen::Matrix<c_double, N, N>& A)
{

    Eigen::Matrix<c_double, N, N> A_k = A;
    Eigen::Matrix<c_double, N, N> pre_A_k = Eigen::Matrix<c_double, N, N>::Zero();
    Eigen::Matrix<c_double, N, N> X_k = Eigen::Matrix<c_double, N, N>::Identity();
    Eigen::Matrix<c_double, N, N> ret_eigenvalues = Eigen::Matrix<c_double, N, N>::Identity();
    Eigen::Matrix<c_double, N, N> ret_eigenvectors = Eigen::Matrix<c_double, N, N>::Identity();
    uint64_t iter = 0;

    // QR法による固有値計算
    while(true)
    {
        auto [Q, R] = qr_householder<N>(A_k);
        A_k = R * Q;
        X_k = X_k * Q;

        // 収束判定
        bool convergence_flag = true;
        for(int i = 0; i < N; i++)
        {
            if(std::fabs(A_k(i, i) - pre_A_k(i, i)) > 1e-8)
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

        pre_A_k = A_k;
        iter++;
    }
    std::cout << "iter: " << iter << "\n";
    return {ret_eigenvalues, ret_eigenvectors};
}

std::pair<c_double, c_double> calculate_2x2_eigenvalue(const Eigen::Matrix<c_double, 2, 2>& A)
{
    c_double a = A(0, 0);
    c_double b = A(1, 0);
    c_double c = A(0, 1);
    c_double d = A(1, 1);

    c_double aa = 1.0;
    c_double bb = -(a + d);
    c_double cc = a * d - b * c;

    c_double t = bb * bb - 4.0 * aa * cc;
    c_double u = 2.0 * aa;

    c_double p = (-bb + std::sqrt(t)) / u;
    c_double q = (-bb - std::sqrt(t)) / u;

    return {p, q};
}

template<int N>
c_double wilkinson_shift(const Eigen::Matrix<c_double, N, N>& A)
{
    c_double a = A(N-1, N-1);
    auto [p, q] = calculate_2x2_eigenvalue(A.block(N-2, N-2, 2, 2));

    return std::abs(p - a) < std::abs(q - a) ? p : q;
}

template<int N>
std::pair<Eigen::Matrix<c_double, N, N>, Eigen::Matrix<c_double, N, N>>
    calclate_eigenvalues_and_eigenvector_hessenberg(const Eigen::Matrix<c_double, N, N>& A)
{
    Eigen::Matrix<c_double, N, N> A_k = A;
    Eigen::Matrix<c_double, N, N> pre_A_k = Eigen::Matrix<c_double, N, N>::Zero();
    Eigen::Matrix<c_double, N, N> X_k = Eigen::Matrix<c_double, N, N>::Identity();
    Eigen::Matrix<c_double, N, N> ret_eigenvalues = Eigen::Matrix<c_double, N, N>::Identity();
    Eigen::Matrix<c_double, N, N> ret_eigenvectors = Eigen::Matrix<c_double, N, N>::Identity();
    uint64_t iter = 0;

    // hessenberg
    Eigen::Matrix<c_double, N, N> Hs = hessenberg<N>(A_k);
    A_k = Hs;

    // QR法による固有値計算
    while(true)
    {
        // ウィルキンソンシフト
        c_double mu = wilkinson_shift<N>(A_k);

        // シフト行列を作成
        Eigen::Matrix<c_double, N, N> As_k = A_k - mu * Eigen::Matrix<c_double, N, N>::Identity();
        
        auto [Q, R] = qr_householder<N>(As_k);
        A_k = R * Q + mu * Eigen::Matrix<c_double, N, N>::Identity();
        X_k = X_k * Q;

        // 収束判定
        bool is_convergence = true;
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < i; j++)
            {
                if(is_convergence && i == N-1 && j == i-1)
                {
                    goto LOOP_END;
                }
                if(std::fabs(A_k(i, j)) > 1e-8)
                {
                    is_convergence = false;
                    break;
                }
            }

            if(not is_convergence) break;
        }
        
        // 右下に2x2の小行列だけが残ったとき
        if(false)
        {
        LOOP_END:
            auto [p, q] = calculate_2x2_eigenvalue(A_k.block(N-2, N-2, 2, 2));
            A_k(N-1, N-1) = p;
            A_k(N-2, N-2) = q;
            A_k(N-1, N-2) = 0.0;
            A_k(N-2, N-1) = 0.0;
        }


        if(is_convergence)
        {
            ret_eigenvalues = A_k;
            ret_eigenvectors = X_k;
            break;
        }

        pre_A_k = A_k;
        iter++;
    }
    
    std::cout << "iter: " << iter << "\n";
    return {ret_eigenvalues, ret_eigenvectors};
}

int main() 
{
    constexpr int N = 25;

    Eigen::Matrix<c_double, N, N> A;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    for(int i = 0; i < N; ++i) 
    {
        for(int j = 0; j < N; ++j) 
        {
            A(i, j) = dis(gen);
        }
    }
    Eigen::Matrix<c_double, N, N> B = A + A.transpose(); // 対称行列にする
          
    auto [eigenvalues, eigenvector] 
        = calclate_eigenvalues_and_eigenvector_hessenberg<N>(B);
    std::cout << "Eigenvalues:\n" << eigenvalues << "\n";
    std::cout << "Eigenvectors:\n" << eigenvector << "\n";

    return 0;
}

#endif // USE_EIGEN

#if USE_MY_MATRIX
#include "matrix.hpp"
using c_double = std::complex<double>;
constexpr int X = 0;


template<int N>
my_mt::Matrix<c_double, N, N> hessenberg(const my_mt::Matrix<c_double, N, N>& A) 
{
    my_mt::Matrix<c_double, N, N> Q = my_mt::Matrix<c_double, N, N>::Identity();
    my_mt::Matrix<c_double, N, N> R = A;
    my_mt::Matrix<c_double, N, N> A_star = A;

    for(int i = 0; i < N-2; ++i) 
    {
        std::cout << "a\n";
        my_mt::Matrix<c_double, N, N> H = create_householder<N>(A_star, i);
        std::cout << "b\n";
        if (H.isZero()) continue;
        std::cout << "c\n";
        // R, Q を更新
        A_star = H * A_star * H.adjoint();
        R = H * R;
        Q = Q * H.adjoint();
    }

    return A_star;
}

template<int N>
my_mt::Matrix<c_double, N, N> 
  create_householder(const my_mt::Matrix<c_double, N, N>& A, size_t n) 
  {
    // 部分ベクトル x = R(i:N-1, i)
    my_mt::VectorA<c_double, X> x = A.block(n, n, N-n, 1);

    // 実ノルムは double で取る
    auto sigma = x.norm();  
    if (std::fabs(sigma) == 0.0) return my_mt::Matrix<c_double, N, N>::Zero();


    // x(0) が 0 だときは位相不定になるのでガード
    c_double sign = std::abs(x(0)) == 0.0 
                    ? c_double(1.0, 0.0) 
                    : x(0) / std::abs(x(0));
 
    // 複素 α を作成
    c_double alpha = - sign * sigma;

    // 反射ベクトル v を作成
    my_mt::VectorA<c_double, X> v = x;
    v(0) -= alpha;
    c_double beta = v.squaredNorm();
    if (beta == 0.0) return my_mt::Matrix<c_double, N, N>::Zero();

    std::cout << "as\n";
    //--------------------------------koko
    // 部分反射行列 H_sub = I - 2/β * v v^T
    //my_mt::Matrix<c_double, X, X> 
    auto H_sub = //my_mt::Matrix<c_double>::Identity(N-n, N-n)
                             /*- (2.0 / beta) */ (v * v.adjoint());
    //--------------------------------koko
    std::cout << "aa\n";
    // 全体反射行列 H を単位行列に埋め込む
    my_mt::Matrix<c_double, N, N> H = my_mt::Matrix<c_double, N, N>::Identity();
    H.block(n, n, N-n, N-n) = H_sub;

    return H;
}

template<int N>
std::pair<my_mt::Matrix<c_double, N, N>, my_mt::Matrix<c_double, N, N>>
qr_householder(const my_mt::Matrix<c_double, N, N>& A) 
{
    my_mt::Matrix<c_double, N, N> Q = my_mt::Matrix<c_double, N, N>::Identity();
    my_mt::Matrix<c_double, N, N> R = A;

    for(int i = 0; i < N; ++i) 
    {

        my_mt::Matrix<c_double, N, N> H = create_householder<N>(R, i);
        if (H.isZero()) continue;
        // R, Q を更新
        R = H * R;
        Q = Q * H.adjoint();
    }

    return {Q, R};
}



template<int N>
std::pair<my_mt::Matrix<c_double, N, N>, my_mt::Matrix<c_double, N, N>>
    calclate_eigenvalues_and_eigenvector(const my_mt::Matrix<c_double, N, N>& A)
{

    my_mt::Matrix<c_double, N, N> A_k = A;
    my_mt::Matrix<c_double, N, N> pre_A_k = my_mt::Matrix<c_double, N, N>::Zero();
    my_mt::Matrix<c_double, N, N> X_k = my_mt::Matrix<c_double, N, N>::Identity();
    my_mt::Matrix<c_double, N, N> ret_eigenvalues = my_mt::Matrix<c_double, N, N>::Identity();
    my_mt::Matrix<c_double, N, N> ret_eigenvectors = my_mt::Matrix<c_double, N, N>::Identity();
    uint64_t iter = 0;

    // QR法による固有値計算
    while(true)
    {
        auto [Q, R] = qr_householder<N>(A_k);
        A_k = R * Q;
        X_k = X_k * Q;

        // 収束判定
        bool convergence_flag = true;
        for(int i = 0; i < N; i++)
        {
            if(std::fabs(A_k(i, i) - pre_A_k(i, i)) > 1e-8)
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

        pre_A_k = A_k;
        iter++;
    }
    std::cout << "iter: " << iter << "\n";
    return {ret_eigenvalues, ret_eigenvectors};
}

template<size_t N, size_t M>
std::pair<c_double, c_double> calculate_2x2_eigenvalue(const my_mt::Matrix<c_double, N, M> A)
{
    if(A.rows() != 2 || A.cols() != 2)
    {
        throw std::invalid_argument("Matrix must be 2x2.");
    }

    c_double a = A.get(0, 0);
    c_double b = A.get(1, 0);
    c_double c = A.get(0, 1);
    c_double d = A.get(1, 1);

    c_double aa = 1.0;
    c_double bb = -(a + d);
    c_double cc = a * d - b * c;

    c_double t = bb * bb - 4.0 * aa * cc;
    c_double u = 2.0 * aa;

    c_double p = (-bb + std::sqrt(t)) / u;
    c_double q = (-bb - std::sqrt(t)) / u;

    return {p, q};
}

template<int N>
c_double wilkinson_shift(const my_mt::Matrix<c_double, N, N>& A)
{
    c_double a = A.get(N-1, N-1);
    auto [p, q] = calculate_2x2_eigenvalue(A.block(N-2, N-2, 2, 2));

    return std::abs(p - a) < std::abs(q - a) ? p : q;
}

template<int N>
std::pair<my_mt::Matrix<c_double, N, N>, my_mt::Matrix<c_double, N, N>>
    calclate_eigenvalues_and_eigenvector_hessenberg(const my_mt::Matrix<c_double, N, N>& A)
{
    my_mt::Matrix<c_double, N, N> A_k = A;
    my_mt::Matrix<c_double, N, N> pre_A_k = my_mt::Matrix<c_double, N, N>::Zero();
    my_mt::Matrix<c_double, N, N> X_k = my_mt::Matrix<c_double, N, N>::Identity();
    my_mt::Matrix<c_double, N, N> ret_eigenvalues = my_mt::Matrix<c_double, N, N>::Identity();
    my_mt::Matrix<c_double, N, N> ret_eigenvectors = my_mt::Matrix<c_double, N, N>::Identity();
    uint64_t iter = 0;

    std::cout << "a\n";
    // hessenberg
    my_mt::Matrix<c_double, N, N> Hs = hessenberg<N>(A_k);
    A_k = Hs;


    // QR法による固有値計算
    while(true)
    {
        std::cout << "b\n";
        // ウィルキンソンシフト
        c_double mu = wilkinson_shift<N>(A_k);
        std::cout << "c\n";
        // シフト行列を作成
        my_mt::Matrix<c_double, N, N> As_k = A_k - mu * my_mt::Matrix<c_double, N, N>::Identity();
        
        std::cout << "d\n";
        auto [Q, R] = qr_householder<N>(As_k);
        std::cout << "e\n";
        A_k = R * Q + mu * my_mt::Matrix<c_double, N, N>::Identity();
        std::cout << "f\n";
        X_k = X_k * Q;
        std::cout << "g\n";

        // 収束判定
        bool is_convergence = true;
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < i; j++)
            {
                if(is_convergence && i == N-1 && j == i-1)
                {
                    goto LOOP_END;
                }
                if(std::fabs(A_k(i, j)) > 1e-8)
                {
                    is_convergence = false;
                    break;
                }
            }

            if(not is_convergence) break;
        }

        // 右下に2x2の小行列だけが残ったとき
        if(false)
        {
        LOOP_END:
            auto [p, q] = calculate_2x2_eigenvalue(A_k.block(N-2, N-2, 2, 2));
            A_k(N-1, N-1) = p;
            A_k(N-2, N-2) = q;
            A_k(N-1, N-2) = 0.0;
            A_k(N-2, N-1) = 0.0;
        }

        if(is_convergence)
        {
            ret_eigenvalues = A_k;
            ret_eigenvectors = X_k;
            break;
        }
        pre_A_k = A_k;
        iter++;
    }
    
    std::cout << "iter: " << iter << "\n";
    return {ret_eigenvalues, ret_eigenvectors};
}

int main() 
{
    constexpr int N = 25;

    my_mt::Matrix<c_double, N, N> A;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    for(int i = 0; i < N; ++i) 
    {
        for(int j = 0; j < N; ++j) 
        {
            A(i, j) = dis(gen);
        }
    }
    my_mt::Matrix<c_double, N, N> B = A + A.transpose(); // 対称行列にする
          
    auto [eigenvalues, eigenvector] 
        = calclate_eigenvalues_and_eigenvector_hessenberg<N>(B);
    // std::cout << "Eigenvalues:\n" << eigenvalues << "\n";
    // std::cout << "Eigenvectors:\n" << eigenvector << "\n";

    return 0;
}
#endif // MY_MATRIX
