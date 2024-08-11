# Numerical Methods in MATLAB: Mini-Projects

## Table of Contents
1. [Overview](#overview)  
    1.1 [Mathematical Backgrounds](#mathematical-backgrounds)  
    1.2 [Complexity Analysis](#complexity-analysis)  
2. [Mini-Project Descriptions](#mini-project-descriptions)  
    2.1 [Eigenvalues Project](#eigenvalues-project)  
    2.2 [Gaussian Methods Project](#gaussian-methods-project)  
    2.3 [Iterative Methods Project](#iterative-methods-project)  
    2.4 [LU Decomposition Project](#lu-decomposition-project)  
    2.5 [QR Decomposition Project](#qr-decomposition-project)  
    2.6 [Nonlinear Equations Project](#nonlinear-equations-project)  
3. [Contact Information](#contact-information)

---

## Overview
This repository contains a series of mini-projects focused on fundamental numerical methods implemented in MATLAB. These projects cover key areas in numerical linear algebra and root-finding techniques, each designed to provide practical insights into mathematical theories and their computational applications.

### 1.1. Mathematical Backgrounds
Each mini-project is grounded in mathematical theory, using methods widely applied in various scientific and engineering contexts. Below is a brief overview of the mathematical background for each project:

#### Eigenvalues Project
- **Eigenvalues and Eigenvectors:**  
  For a square matrix \( A \), an eigenvalue \( \lambda \) and a corresponding eigenvector \( v \) are defined by the equation:  
  \[
  Av = \lambda v
  \]
  The eigenvalue \( \lambda \) is a scalar that scales the eigenvector \( v \), which is non-zero.
  
- **Power Method:**  
  Used to find the largest eigenvalue (in magnitude) of a matrix. The method iteratively multiplies a vector by the matrix, normalizing it each time:  
  \[
  v_{k+1} = \frac{Av_k}{\|Av_k\|}
  \]
  The sequence \( v_k \) converges to the eigenvector associated with the largest eigenvalue.
  
- **Inverse Power Method:**  
  Used to find the smallest eigenvalue (in magnitude) by applying the power method to the inverse of the matrix:  
  \[
  v_{k+1} = \frac{A^{-1}v_k}{\|A^{-1}v_k\|}
  \]
  This method can also be modified to find eigenvalues near a specific value using a shift \( \sigma \), applied as \( (A - \sigma I)^{-1} \).
  
- **Deflation:**  
  After finding an eigenvalue \( \lambda_1 \) and eigenvector \( v_1 \), deflation reduces the matrix to simplify finding the next eigenvalue. The deflated matrix \( B \) is:  
  \[
  B = A - \lambda_1 v_1 v_1^T
  \]
  This process is repeated to find subsequent eigenvalues.

#### Gaussian Methods Project
- **Gaussian Elimination:**  
  A method to solve a system of linear equations \( Ax = b \) by transforming the system into an upper triangular form, then performing back substitution:
  \[
  A \rightarrow U
  \]
  where \( U \) is an upper triangular matrix.

- **Partial Pivoting:**  
  In Gaussian elimination, partial pivoting involves swapping rows to place the largest available element (in absolute value) in the pivot position to reduce numerical errors.

- **Scaled Partial Pivoting:**  
  Before performing partial pivoting, matrix rows are scaled by their largest absolute element, further reducing the chances of numerical instability.

- **Total Pivoting:**  
  Total pivoting involves swapping both rows and columns to ensure the best pivot element is used, minimizing errors and ensuring stability.

#### Iterative Methods Project
- **Iterative Methods for Linear Systems:**  
  Used to solve systems of linear equations \( Ax = b \) where direct methods (like Gaussian elimination) are computationally expensive. These methods iteratively improve an initial guess \( x_0 \).

- **Jacobi Method:**  
  Updates each component of the solution vector independently:
  \[
  x_i^{(k+1)} = \frac{1}{a_{ii}}\left( b_i - \sum_{j \neq i} a_{ij} x_j^{(k)} \right)
  \]
  Convergence depends on the matrix being diagonally dominant or symmetric positive definite.

- **Gauss-Seidel Method:**  
  An improvement over Jacobi, using the latest values of the variables as soon as they are available:
  \[
  x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)} \right)
  \]

- **Successive Over-Relaxation (SOR):**  
  Accelerates the Gauss-Seidel method by introducing a relaxation factor \( \omega \):
  \[
  x_i^{(k+1)} = (1 - \omega) x_i^{(k)} + \omega \frac{1}{a_{ii}} \left( b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)} \right)
  \]
  Choosing an optimal \( \omega \) can significantly speed up convergence.

#### LU Decomposition Project
- **LU Decomposition:**  
  Factors a matrix \( A \) into a lower triangular matrix \( L \) and an upper triangular matrix \( U \):
  \[
  A = LU
  \]
  Solving \( Ax = b \) involves solving two simpler systems \( Ly = b \) and \( Ux = y \).

- **Cholesky Decomposition:**  
  For symmetric positive definite matrices, \( A \) can be decomposed as:
  \[
  A = LL^T
  \]
  This decomposition is more efficient and numerically stable than general LU decomposition.

- **Crout and Doolittle Methods:**  
  - **Crout:** Constructs \( L \) with ones on the diagonal and \( U \) with non-unit diagonal elements.  
  - **Doolittle:** Constructs \( L \) with unit diagonal elements and \( U \) with the computed values.

#### QR Decomposition Project
- **QR Decomposition:**  
  Factors a matrix \( A \) into an orthogonal matrix \( Q \) and an upper triangular matrix \( R \):
  \[
  A = QR
  \]
  \( Q \) has orthogonal columns (i.e., \( Q^TQ = I \)), and \( R \) is upper triangular.

- **Gram-Schmidt Process:**  
  Orthogonalizes a set of vectors to form the columns of \( Q \):
  \[
  q_i = a_i - \sum_{j=1}^{i-1} (a_i \cdot q_j)q_j
  \]
  This process is susceptible to numerical instability, leading to loss of orthogonality.

- **Givens Rotation:**  
  Introduces zeros below the diagonal in a matrix, used to form \( R \) by applying a series of rotations.

- **Householder Transformation:**  
  A more numerically stable method compared to Gram-Schmidt, using reflections to introduce zeros below the diagonal:
  \[
  H = I - 2\frac{vv^T}{v^Tv}
  \]
  \( H \) is used to reflect a vector onto a coordinate axis.

#### Nonlinear Equations Project
- **Nonlinear Equations:**  
  Finding the roots of nonlinear equations \( f(x) = 0 \) is a common problem in various fields. Closed-form solutions often don't exist, necessitating iterative methods.

- **Bisection Method:**  
  A bracketing method that iteratively narrows down the interval \([a, b]\) where the root lies:
  \[
  c = \frac{a + b}{2}
  \]
  If \( f(c) \) changes sign between \( a \) and \( c \), then the root is in \([a, c]\); otherwise, it's in \([c, b]\).

- **Newton-Raphson Method:**  
  An open method that uses the derivative to find successive approximations to the root:
  \[
  x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}
  \]
  It converges quickly if the initial guess is close to the true root but may fail if \( f'(x) \) is small or the initial guess is poor.

- **Secant Method:**  
  An approximation to Newton-Raphson that doesn't require the derivative:
  \[
  x_{k+1} = x_k - f(x_k)\frac{x_k - x_{k-1}}{f(x_k) - f(x_{k-1})}
  \]
  It uses secant lines to approximate the derivative.

### 1.2. Complexity Analysis
Each mini-project has varying computational complexities depending on the method used and the size of the problem. Below is a brief complexity analysis:

- **Eigenvalues Project:**
  - **Power Method:** \( O(n^2) \) per iteration, with total complexity depending on the number of iterations for convergence.
  - **Inverse Power Method:** Similar to the Power Method, but includes additional complexity due to the matrix inversion.
  
- **Gaussian Methods Project:**
  - **Gaussian Elimination:** \( O(n^3) \).

- **Iterative Methods Project:**
  - **Jacobi and Gauss-Seidel:** \( O(n^2) \) per iteration, with the total complexity dependent on the number of iterations required for convergence.

- **LU Decomposition Project:**
  - **LU Decomposition:** \( O(n^3) \), with optimizations available for special matrix types (e.g., symmetric positive definite matrices in Cholesky decomposition).

- **QR Decomposition Project:**
  - **QR Decomposition:** \( O(n^3) \), with different methods (e.g., Gram-Schmidt vs. Householder) having slightly varying constants.

- **Nonlinear Equations Project:**
  - **Bisection Method:** \( O(\log(\epsilon)) \), where \( \epsilon \) is the desired accuracy.
  - **Newton-Raphson Method:** \( O(n) \) per iteration, with the total complexity depending on the number of iterations.

## Mini-Project Descriptions

### 2.1. Eigenvalues Project
Explore methods for computing eigenvalues and eigenvectors, focusing on the Power Method, Inverse Power Method, and Deflation techniques. The project involves implementing these methods and analyzing their convergence properties.

### 2.2. Gaussian Methods Project
This project covers Gaussian elimination and its variants: partial, scaled partial, and total pivoting. The emphasis is on understanding how pivoting improves numerical stability and prevents issues like division by small numbers.

### 2.3. Iterative Methods Project
Investigate iterative methods for solving linear systems, including the Jacobi, Gauss-Seidel, and Successive Over-Relaxation (SOR) methods. The project includes implementing these methods and comparing their convergence rates for different types of matrices.

### 2.4. LU Decomposition Project
Implement LU decomposition using different methods: Cholesky (for symmetric positive definite matrices), Crout, and Doolittle. The project involves solving linear systems using LU decomposition and comparing it to other methods like Gaussian elimination.

### 2.5. QR Decomposition Project
Explore QR decomposition through different approaches: Gram-Schmidt, Givens Rotations, and Householder transformations. The project includes applications of QR decomposition in solving linear systems and finding eigenvalues.

### 2.6. Nonlinear Equations Project
Focus on root-finding methods for nonlinear equations, including the Bisection Method, Newton-Raphson Method, and Secant Method. The project involves implementing these methods and analyzing their convergence and stability.

## Contact Information
For any questions or further information, please contact [Your Name] at [Your Email Address].
