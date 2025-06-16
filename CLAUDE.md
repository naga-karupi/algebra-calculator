# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System & Commands

This is a C++ project using CMake with C++20 features. The project depends on Eigen3 and GoogleTest.

**Build commands:**
```bash
mkdir -p build && cd build
cmake ..
make
```

**Run executables:**
```bash
./determinant        # Matrix determinant calculator
./eigenvalue         # Eigenvalue computation
./multiplication     # Matrix multiplication demo
```

**Run tests:**
```bash
./test_matrix        # Run all GoogleTest unit tests
```

## Code Architecture

### Core Matrix Library (`my_mt` namespace)
- **Template-based Matrix class** (`inc/matrix.hpp`): Supports both compile-time fixed-size matrices (`Matrix<TYPE, ROW, COL>`) and runtime dynamic matrices (`Matrix<TYPE, 0, 0>`)
- **Dual storage strategy**: Fixed-size matrices use stack arrays, dynamic matrices use `std::vector<std::vector<TYPE>>`
- **Implementation files**: Core logic in `inc/detail/matrix.ipp`, helper types in `inc/helper_class.hpp`

### Linear Algebra Features
- **LU Decomposition** (`inc/LU_decomposition.hpp`): Template-based LU factorization for square matrices
- **Matrix operations**: Addition, subtraction, multiplication, scalar operations, determinant, trace, norm
- **Special matrices**: Identity, Zero, transpose, adjoint (for complex matrices)
- **Matrix properties**: Symmetric, zero, identity checks

### Applications
- **Determinant calculation** (`src/det.cpp`): Uses both custom implementation and Eigen3
- **Eigenvalue computation** (`src/eigen-value.cpp`): Eigenvalue/eigenvector calculations
- **Matrix multiplication** (`src/multi-mat.cpp`): Matrix multiplication demonstrations

### Dependencies
- **Eigen3**: External library for reference implementations and performance comparison
- **GoogleTest**: Unit testing framework for comprehensive matrix operation testing
- **C++20**: Uses modern C++ features, compiled with `-Wall -Wextra -pedantic -O2`

### Testing Strategy
All matrix operations are thoroughly tested in `test/test.cpp` including:
- Constructor variations, arithmetic operations, special matrix properties
- Numerical precision testing (using `std::fabs` for floating-point comparisons)
- Complex number support (adjoint operations)