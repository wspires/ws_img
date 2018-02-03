/**
	@file   Linear_Algebra.cpp
	@author Wade Spires
	@date   2005/08/03
	@brief  Solve linear system Ax=b using LU Factorization by Gaussian
		elimination with partial pivoting.

	Solve linear system Ax=b using LU Factorization by Gaussian
	elimination with partial pivoting.
	First, A is factored into lower- and upper-triangular matrices L and U,
	respectively. L and U overwrite A in memory, and L's diagonal is not
	explicitly stored as it consists entirely of 1's.  To solve the system
	LUx=b, we let y=Ux and solve for y in the system Ly=b by
	forward-substitution. We then solve for x in Ux=y by
	backward-substitution.  The forward-substitution routine is not general
	purpose as it does not divide by the diagonal element since L's diagonal
	holds only 1's. In both substitution procedures, row sums are computed
	and then subtracted from the corresponding element in the right-hand
	side, so the procedures do not perform a traditional Gaussian reduction.
	Finally, the pivot vector holds the location of each row in A to reduce
	the possible division by zero that could result in dividing by the
	multipliers when performing LU Factorization.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Linear_Algebra.hpp"

using namespace ws_img;

typedef Lin_Alg::value_type value_type;
typedef Lin_Alg::size_type  size_type;

const double    Lin_Alg::DEFAULT_TOL  = .0000001;
const unsigned  Lin_Alg::DEFAULT_ITER = 100;
const double    Lin_Alg::EPS          = DBL_EPSILON;

using std::cout; using std::cerr; using std::endl;
using std::vector;
using std::domain_error;
using std::plus;
using std::minus;
using std::multiplies;
using std::transform;

/**
	Compute determinant of matrix.
	@param[in] A Matrix
	@return determinant Determinant of matrix A
 */
value_type
Lin_Alg::det( const Matrix& A )
{
	value_type determinant = 0;

	// put A in upper-triangular form 
	Matrix U = A;
	// perform Gaussian elimination

	// product of diagonal entries is determinant
	for( unsigned i = 0; i != U.row_size(); ++i )
	{
		determinant *= U(i, i);
	}

	return( determinant );
}

/**
	Find inverse of 2-by-2 matrix directly.
	@param[in] A Matrix
	@retval A_inv Inverse of A
 */
Matrix
Lin_Alg::inverse_2_by_2( const Matrix& A )
{
	// verify matrix is 2 x 2
	if( A.row_size() != 2 || A.col_size() != 2 )
	{
		throw std::domain_error( "Matrix must be 2x2 to compute inverse" );
	}

	// compute determinant: det(A) = ad - bc
	double determinant = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);

	/*
		Given matrix	[ a b ],
							[ c d ]
		change the following elements for its inverse:
		a <-> d, b -> -b, and c -> -c.
	 */
	Matrix A_inv = A;
	A_inv(0, 0) =  A(1, 1);
	A_inv(0, 1) = -A(0, 1);
	A_inv(1, 0) = -A(1, 0);
	A_inv(1, 1) =  A(0, 0);

	// divide each element by det(A)
	A_inv /= determinant;

	return( A_inv );
}

/**
	Find inverse of 3-by-3 matrix A directly.
	@param[in] A Matrix to invert
	@retval A_inv Inverse of A
 */
Matrix
Lin_Alg::inverse_3_by_3( const Matrix& A )
{
	if( A.sub_row_size() != 3 || A.sub_col_size() != 3 )
	{
		err_quit( "Matrix dimensions (%u x %u) must be 3-by-3",
				A.sub_row_size(), A.sub_col_size() );
	}

	// easier to deal with short names
	const Matrix::value_type a = A(0,0);
	const Matrix::value_type b = A(0,1);
	const Matrix::value_type c = A(0,2);
	const Matrix::value_type d = A(1,0);
	const Matrix::value_type e = A(1,1);
	const Matrix::value_type f = A(1,2);
	const Matrix::value_type g = A(2,0);
	const Matrix::value_type h = A(2,1);
	const Matrix::value_type i = A(2,2);

	// compute determinant
	double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
	if( det == 0 )
	{
		err_quit( "inverse: Determinant is 0\n" );
	}

	// compute inverse
	Matrix A_inv( 3, 3 );
	A_inv(0,0) = e * i - f * h;
	A_inv(0,1) = c * h - b * i;
	A_inv(0,2) = b * f - c * e;
	A_inv(1,0) = f * g - d * i;
	A_inv(1,1) = a * i - c * g;
	A_inv(1,2) = c * d - a * f;
	A_inv(2,0) = d * h - e * g;
	A_inv(2,1) = b * g - a * h;
	A_inv(2,2) = a * e - b * d;

	A_inv /= det;

	return( A_inv );
}

/**
	Solve linear system Ax = b by LU factorization--solve LU*x = b where A = LU.
	@param[in,out] A Square matrix
	@param[in,out] b Right-hand side solution vector
 */ 
void
Lin_Alg::solve_system( Matrix& A, Vector& b, double tolerance )
{
	// verify that matrix is square
	assert( A.row_size() == A.col_size() );

	Pivot pivot;

	// initialize pivot
	for( unsigned i = 0; i < A.row_size(); ++i )
	{
		pivot.push_back( i );
	}

	double d;
	LU_factorization( A, pivot, &d, 1.0E-20 );

	LU_back_sub( A, pivot, b, tolerance );
}

/**
	Given a matrix a [1..n][1..n], this routine replaces it by the LU
	decomposition of a rowwise permutation of itself. a and n are input. a is
	output, arranged with L and U in the same matrix; indx[1..n] is an
	output vector that records the row permutation affected by the partial
	pivoting; d is output as +/-1 depending on whether the number of row
	interchanges was even or odd, respectively.
	This routine is used in combination with LU_back_sub()
	to solve linear equations or invert a matrix.
 */
void
Lin_Alg::LU_factorization( Matrix& A, Pivot& pivot, double *d,
	double small_value )
{
	unsigned n = A.row_size();

	// implicit scaling of each row
	Vector scale_vec( n );

	*d = 1.0;  // no row interchanges yet

	// loop over rows to get the implicit scaling information
	for( unsigned i = 0; i != n; ++i )
	{
		double big = 0.0;
		for( unsigned j = 0; j != n; ++j )
		{
			double val = fabs( A(i,j) ); 
			if( val > big)
			{
				big = val;
			}
		}
		if( big == 0.0 )
		{
			err_quit( "Singular matrix\n" );
		}

		// nonzero largest element: save the scaling
		scale_vec(i) = 1.0 / big;
	}

	// loop over columns of Crout's method
	for( unsigned j = 0; j != n; ++j )
	{
		// sum is form of a triangular matrix except for i=j
		for( unsigned i = 0; i != j; ++i )
		{
			double sum = A(i,j);     
			for( unsigned k = 0; k != i; ++k )
			{
				sum -= A(i,k) * A(k,j);
			}
			A(i,j) = sum;
		}

		// initialize for the search for largest pivot element
		// set default value for imax
		double big = 0.0;
		unsigned pivot_pos = 0;

		for( unsigned i = j; i != n; ++i )
		{
			// i=j of previous sum and i=j+1...N of the rest of the sum
			double sum = A(i,j);
			for( unsigned k = 0; k != j; ++k )
			{
				sum -= A(i,k) * A(k,j);
			}
			A(i,j) = sum;

			// if the figure for the pivot is better than the best so far
			double val = scale_vec(i) * fabs(sum);
			if( val >= big )
			{
				big = val;
				pivot_pos = i;
			}
		}

		// check to see if rows should be interchanged
		if( j != pivot_pos )
		{
			for( unsigned k = 0; k != n; ++k )
			{
				// interchange rows
				double val = A(pivot_pos,k);
				A(pivot_pos,k) = A(j,k);
				A(j,k) = val;
			}

			// change the parity of d and interchange the scale factor
			*d = -(*d);
			scale_vec(pivot_pos) = scale_vec(j);
		}

		// If the pivot element is zero the matrix is singular (at least
		// to the precision of the algorithm). For some applications on
		// singular matrices, it is desiderable to substitute TINY for zero
		pivot[j] = pivot_pos;
		if( A(j,j) == 0.0 )
		{
			A(j,j) = small_value;
		}

		if( j != n )
		{
			// divide by pivot element
			double val = 1.0 / A(j,j);
			for( unsigned i = (j + 1) ; i != n; ++i )
			{
				A(i,j) *= val;
			}
		}
	}  // go back for the next column in the reduction
}
    
/**
	Solves the set of n linear equations A.X = B.
	Here a[1..n][1..n] is input, not as the matrix A
	but rather as its LU decomposition, determined
	by the routine LU_decomp(). indx[1..n] is input as
	the permutation vector returned by LUdcmp().
	b[1..n] is input as the right hand side vector B,
	and returns with the solution vector X. a, n, and
	indx are not modified by this routine and can be
	left in place for successive calls with different
	right-hand sides b. This routine takes into account
	the possibility that b will begin with many zero
	elements, so it it efficient for use in matrix
	inversion.
 */
void
Lin_Alg::LU_back_sub( const Matrix& A, const Pivot& pivot, Vector& b,
	double tolerance )
{
	unsigned n = A.row_size();

	// initially don't perform sum (as long as b_i is 0)
	unsigned idx = 0;

	// forward substitution: traverse columns
	// use pivot to move b's elements in correct positions
	for( unsigned i = 0; i != n; ++i )
	{
		// When idx is set to a positive value, it will
		// become the index of the first nonvanishing
		// element of b. We now do the forward substitution.
		// The only new wrinkle is to unscramble the
		// permutation as we go.

		// use pivot to get actual b[i]
		unsigned pivot_pos = pivot[i];
		double sum = b(pivot_pos);
		b(pivot_pos) = b(i);

		if( idx != 0 )
		{
			// perform forward sub. at row i: 
			assert( idx != 0 );
			for( unsigned j = (idx - 1); j != i; ++j )
			{
				sum -= A(i,j) * b(j);
			}
		}

		// if a nonzero element was encountered, we will have to do the sums
		// in the loop above
		else if( sum != 0.0 )
		{
			idx = i;
		}

		b(i) = sum;  // finish swapping elements
	}

	// perform the actual back-substitution
	for( int i = (n - 1); i >= 0; --i )
	{
		double sum = b(i);

		for( unsigned j = i + 1; j < n; ++j )
		{
			sum -= A(i,j) * b(j);
		}

		// Store a component of the solution vector x
		b(i) = sum / A(i,i);

		// set to 0 if very small (less than our tolerance)
		if( fabs( b(i) ) < tolerance )
		{
			// Verify small numbers
			b(i) = 0.0;
		}
	}
}

/**
	Solve linear system Ax = b by LU factorization--solve LU*x = b where A = LU.
	@param[in,out] A Square matrix
	@param[in] b Right-hand side solution vector
	@param[out] x Solution vector
 */ 
void
Lin_Alg::solve_system( Matrix& A, const Vector& b, Vector& x )
{
	// verify that matrix is square
	assert( A.row_size() == A.col_size() );

	Pivot pivot( A.row_size() );

	// initialize pivot
	for( unsigned i = 0; i < A.row_size(); ++i )
	{
		pivot[i] = i;
	}

	LU_factorization( A, pivot );

	// forward substitution to find solution to L*y=b where y = U*x
	Vector  y( A.row_size(), 0 );
	//Matrix  y( A.row_size(), 1, 0 );
	LU_forw_sub( A, b, y, pivot );

	// backward substitution to find solution to U*x=y
	Vector  x_nonpiv( A.row_size(), 0 );
	//Matrix  x( A.row_size(), 1, 0 );
	int err_element = LU_back_sub( A, y, x_nonpiv, pivot );

	if( err_element >= 0 )
	{
		err_quit( "Matrix is singular (diagonal entry %d = 0)\n",
			err_element );
		//throw std::domain_error( "Singular matrix" );
	}

	// resize vector if too small
	if( x.size() < pivot.size() )
	{
		x.resize( pivot.size() );
	}

	// create solution x by storing pivoted elements
	for( unsigned i = 0; i != pivot.size(); ++i )
	{
		x(i) = x_nonpiv( pivot[i] );
	}
}

/**
	Solve linear system Ax = b by LU factorization--solve LU*x = b where A = LU.
	@param[in] A Square matrix
	@param[in] input_b Right-hand side solution vector
	@return solution Solution vector to system
 */ 
Vector
Lin_Alg::solve_system( const Matrix& A, const Vector& b )
{
/*
	Matrix A_copy = A;
	Vector b_copy = b;
	Vector x = Vector( b.size() );
	Lin_Alg::solve_system( Matrix& A, Vector& b, Vector& x )
	return( x );
 */

	// verify that matrix is square
	assert( A.row_size() == A.col_size() );

	Matrix A_copy = A;
	Vector b_copy = b;
	Pivot pivot;

	// initialize pivot
	for( unsigned i = 0; i < A_copy.row_size(); ++i )
	{
		pivot.push_back( i );
	}

	LU_factorization( A_copy, pivot );

	// forward substitution to find solution to L*y=b where y = U*x
	Vector  y( A_copy.row_size(), 0 );
	//Matrix  y( A_copy.row_size(), 1, 0 );
	LU_forw_sub( A_copy, b_copy, y, pivot );

	// backward substitution to find solution to U*x=y
	Vector  x( A_copy.row_size(), 0 );
	//Matrix  x( A_copy.row_size(), 1, 0 );
	int err_element = LU_back_sub( A_copy, y, x, pivot );

	if( err_element >= 0 )
	{
		throw std::domain_error( "Singular matrix" );
	}

	// permute final solution
	Vector solution( pivot.size() );
	for( unsigned i = 0; i != pivot.size(); ++i )
	{
		solution(i) =  x( pivot[i] );
	}

	return( solution );
}

/**
	Factorize matrix A into L and U matrices.
	Notes: L and U write over A,
	       A's diagonal is replaced by U's
	       L's diagonal is not stored since it is all 1s

	@param[in,out] A Square matrix to factorize
 */ 
void Lin_Alg::LU_factorization( Matrix& A, Pivot& pivot )
{
	// verify that matrix is square
	assert( A.row_size() == A.col_size() );

	// loop over columns
	for( unsigned k = 0; k != A.col_size(); k++ )
	{
		// search for pivot in current column
		unsigned p = k;
		for( unsigned i = k + 1; i < A.row_size(); i++ )
		{
			if( fabs( A(pivot[i], k) ) > fabs( A(pivot[p], k) ) )
			{
				p = i;
			}
		}

		// interchange rows, if necessary
		if( p != k )
		{
			int tmp  = pivot[k];
			pivot[k] = pivot[p];
			pivot[p] = tmp;	
		}

		// skip zero column
		if( A(pivot[k], k) == 0 )
		{
			continue;
		}

		// compute multiplier for current column
		for( unsigned i = k + 1; i < A.row_size(); i++ )
		{
			double multiplier = A(pivot[i], k) / A(pivot[k], k);

			// set L element as multiplier
			A(pivot[i], k) = multiplier;

			// apply transformation to remaining matrix
			for( unsigned j = k + 1; j < A.col_size(); j++ )
			{
				A(pivot[i], j) = A(pivot[i], j) - multiplier * A( pivot[k], j);
			}
		}
	}
}

/**
	Perform forward substitution to obtain solution y to Ly=b.
	Note: Not a general forward-substitution routine--matrix must have unit
		diagonal.
 */ 
void Lin_Alg::LU_forw_sub( const Matrix& A, const Vector& b, Vector& y, 
	const Pivot& pivot )
{
	// traverse columns
	for( unsigned i = 0; i != A.row_size(); i++ )
	{
		double row_sum = 0.0;
		for( unsigned j = 0; j != i; j++ )
		{
			// calculate row sum
			row_sum += A(pivot[i], j) * y( pivot[j] );
		}

		// calculate i_th solution
		y( pivot[i] ) = b( pivot[i] ) - row_sum;
	}
}

/**
	Perform backward substitution to obtain solution x to Ux=y.
 */ 
int Lin_Alg::LU_back_sub( const Matrix& A, const Vector& y, Vector& x, 
	const Pivot& pivot )
{
	for( int j = int(A.row_size()) - 1; j >= 0; j-- )
	{
		// if diagonal entry is zero, return position
		if( A(pivot[j], j) == 0 )
		{
			return( pivot[j] );
		}

		// initialize row sum
		double row_sum = 0.0;
		for( unsigned k = j + 1; k < A.row_size(); k++ )
		{
			// calculate row sum
			row_sum += A(pivot[j], k) * x( pivot[k] );
		}

		// calculate j_th solution
		x( pivot[j] ) = ( y( pivot[j] ) - row_sum ) / A(pivot[j], j);
	}
	
	return( -1 );
}

/**
	Jacobit iterative method to solve Ax=b.

	@param[in] tolerance Threshold in which to stop computation when error
		reduces to less than
	@param[in] max_iterations Maximum number of iterations to perform in
		algorithm
 */
Vector Lin_Alg::jacobi( Matrix& A, Vector& b, unsigned max_iterations,
	double tolerance, int* number_of_iterations )
{
	// current and previous solutions
	Vector  x( A.row_size() );
	Vector  x_old( A.row_size() );

	for( unsigned k = 0; k != max_iterations; k++ )
	{
		for( unsigned i = 0; i != A.row_size(); i++ )
		{
			double sum = 0;
			for( unsigned j = 0; j != A.row_size(); j++ )
			{
				if( j == i )
				{
					continue;
				}
				sum += A(i, j) * x_old(j);
			}

			x(i) = ( -sum + b( i ) ) / A(i, i);
		}

		Vector difference = x - x_old;
		if( difference.inf_norm() < tolerance )
		{
			if( number_of_iterations != NULL )
			{
				*number_of_iterations = k + 1;
			}
			return( x );
		}

		x_old = x;
	}

	err_warn( "Maximum iterations of %u was exceeded with tolerance of %lf\n",
		max_iterations, tolerance );

	return( x );
}

/**
	Gauss-Seidel iterative method to solve Ax=b.
 */
Vector 
Lin_Alg::gauss_seidel( Matrix& A, Vector& b, unsigned max_iterations, 
	double tolerance, int* number_of_iterations )
{
	// current and previous solutions
	Vector  x( A.row_size() );
	Vector  x_old( A.row_size() );

	for( unsigned k = 0; k != max_iterations; k++ )
	{
		for( unsigned i = 0; i != A.row_size(); i++ )
		{
			double sum_1 = 0;
			for( unsigned j = 0; j != i; j++ )
			{
				sum_1 += A(i, j) * x(j);
			}

			double sum_2 = 0;
			for( unsigned j = i + 1; j != A.row_size(); j++ )
			{
				sum_2 += A(i, j) * x_old(j);
			}

			x(i) = ( -sum_1 - sum_2 + b( i ) ) / A(i, i);
		}

		Vector difference = x - x_old;
		if( difference.inf_norm() < tolerance )
		{
			if( number_of_iterations != NULL )
			{
				*number_of_iterations = k + 1;
			}
			return( x );
		}

		x_old = x;
	}

	err_warn( "Maximum iterations of %u was exceeded with tolerance of %lf\n",
		max_iterations, tolerance );

	return( x );
}

/**
	SOR iterative method to solve Ax=b.
 */
Vector 
Lin_Alg::SOR( Matrix& A, Vector& b, double omega,
	unsigned max_iterations, double tolerance, int* number_of_iterations )
{
	// current and previous solutions
	Vector  x( A.row_size() );
	Vector  x_old( A.row_size() );

	for( unsigned k = 0; k != max_iterations; k++ )
	{
		for( unsigned i = 0; i != A.row_size(); i++ )
		{
			double sum_1 = 0;
			for( unsigned j = 0; j != i; j++ )
			{
				sum_1 += A(i, j) * x(j);
			}

			double sum_2 = 0;
			for( unsigned j = i + 1; j != A.row_size(); j++ )
			{
				sum_2 += A(i, j) * x_old(j);
			}

			x(i) = ( 1 - omega ) * x_old(i) + 
				( omega * ( -sum_1 - sum_2 + b(i) ) ) / A(i, i);
		}

		Vector difference = x - x_old;
		if( difference.inf_norm() < tolerance )
		{
			if( number_of_iterations != NULL )
			{
				*number_of_iterations = k + 1;
			}
			return( x );
		}

		x_old = x;
	}

	err_warn( "Maximum iterations of %u was exceeded with tolerance of %lf\n",
		max_iterations, tolerance );

	return( x );
}

/**
	Singular Value Decomposition (SVD).

	A general, rectangular M-by-N matrix A has a singular value decomposition
	(SVD) into the product of an M-by-N orthogonal matrix U, an N-by-N diagonal
	matrix of singular values S, and the transpose of an N-by-N orthogonal square
	matrix V: A = U S V^T 

	The diagonal elements of the singular value matrix S are stored in the
	vector S. The singular values are non-negative and form a non-increasing
	sequence from S_1 to S_N. The matrix V contains the elements of V in
	untransposed form. To form the product U S V^T, it is necessary to take the
	transpose of V.

	@param[in,out] A Input matrix (transformed to orthogonal matrix U)
	@param[out] S Vector of singular values
	@param[out] V Square orthogonal matrix
	@param[in] method Algorithm for computing SVD (default is Golub_Reinsch)
 */
void
Lin_Alg::SVD( Matrix& A, Vector& S, Matrix& V, const SVD_Method& method )
{
#ifdef HAVE_GSL
	// number of rows and columns in a subregion of A
	size_type M = A.sub_row_size();
	size_type N = A.sub_col_size();

	// allocate output matrices
	S.resize( N, false );
	V.resize( N, N, false );

	// allocate U, V, and S
	// Note the following:
	// 1. convert A to a form that GSL can use: 'matrix views' use the actual
	// data being pointed to, so we should not call 'delete[] data' until the
	// view is no longer used
	// 2. S is a _vector_ of singular values representing a diagonal
	gsl_matrix_view m = gsl_matrix_view_array( A.data(), M, N );
	gsl_matrix* U_gsl = &m.matrix;
	gsl_vector* S_gsl = gsl_vector_alloc( N );
	gsl_matrix* V_gsl = gsl_matrix_alloc( N, N );

	// workspace for Golub-Reinsch method (unused for Jacobi method)
	gsl_vector* work = 0;
	gsl_matrix* X    = 0;

	// perform singular value decomposition using the given algorithm
	switch( method )
	{
		// computes singular values to higher relative accuracy than Golub-Reinsch
		case Jacobi:
			gsl_linalg_SV_decomp_jacobi( U_gsl, V_gsl, S_gsl );
			break;

		case Golub_Reinsch:
			work = gsl_vector_alloc( N );
			gsl_linalg_SV_decomp( U_gsl, V_gsl, S_gsl, work );
			gsl_vector_free( work );
			break;

		// modified to be faster for M >> N
		case Golub_Reinsch_mod:
			work = gsl_vector_alloc( N );
			X    = gsl_matrix_alloc( N, N );
			gsl_linalg_SV_decomp_mod( U_gsl, X, V_gsl, S_gsl, work );
			gsl_vector_free( work );
			gsl_matrix_free( X );
			break;

		default:
			err_quit( "SVD: Unknown SVD algorithm '%u'\n", method );
			break;
	}

	// store GSL matrices and vector in my matrix/vector objects
	for( size_type i = V.row_begin(); i != V.row_end(); ++i )
	{
		for( size_type j = V.col_begin(); j != V.col_end(); ++j )
		{
			V(i,j) = gsl_matrix_get( V_gsl, i, j );
		}
		S(i) = gsl_vector_get( S_gsl, i );
	}

	// deallocate U, S, and V
	gsl_vector_free( S_gsl );
	gsl_matrix_free( V_gsl );

#else
	err_warn( "Function Lin_Alg::SVD() not implemented:"
			" GNU Scientific Laboratory (GSL) is not installed.\n" );
	return;
#endif // HAVE_GSL
}

/**
	Calculate pseudo-inverse of matrix A.

	Compute the pseudo-inverse using SVD since
		A = U * W * V^t
		=> pinv = (U * W * V^t)^-1 = (V^t)^-1 * W^-1 * U^-1
		= V * W^-1 * U^t
	This matrix has the property that
		A = A * P * A

	@param[in,out] A Input matrix (transformed to orthogonal matrix U)
	@param[in] A Input matrix
	@retval A_inv Pseudo-inverse of A
 */
Matrix
Lin_Alg::pseudo_inverse( const Matrix& A, const SVD_Method& method )
{
#ifdef HAVE_GSL
	// number of rows and columns in a subregion of A
	size_type M = A.sub_row_size();
	size_type N = A.sub_col_size();

	// allocate U, V, and S
	// Note the following:
	// 1. convert A to a form that GSL can use: 'matrix views' use the actual
	// data being pointed to, so we should not call 'delete data' until view
	// is no longer used
	// 2. S is _vector_ of singular values => S = diag(W)
	double* data = A.to_array();
	gsl_matrix_view m = gsl_matrix_view_array( data, M, N );
	gsl_matrix* U = &m.matrix;
	gsl_matrix* V = gsl_matrix_alloc( N, N );
	gsl_vector* S = gsl_vector_alloc( N );

	// workspace for Golub-Reinsch method (unused for Jacobi method)
	gsl_vector* work = 0;
	gsl_matrix* X    = 0;

	// perform singular value decomposition using the given algorithm
	switch( method )
	{
		// computes singular values to higher relative accuracy than Golub-Reinsch
		case Jacobi:
			gsl_linalg_SV_decomp_jacobi( U, V, S );
			break;

		case Golub_Reinsch:
			work = gsl_vector_alloc( N );
			gsl_linalg_SV_decomp( U, V, S, work );
			gsl_vector_free( work );
			break;

		// modified to be faster for M >> N
		case Golub_Reinsch_mod:
			work = gsl_vector_alloc( N );
			X    = gsl_matrix_alloc( N, N );
			gsl_linalg_SV_decomp_mod( U, X, V, S, work );
			gsl_vector_free( work );
			gsl_matrix_free( X );
			break;

		default:
			err_quit( "SVD: Unknown SVD algorithm '%u'\n", method );
			break;
	}

	// compute V * W^-1
	for( unsigned j = 0; j != V->size2; ++j )
	{
		// set W_inv = 1 / s, but check for a singular matrix (values close to 0)
		if( fabs( gsl_vector_get(S,j) >= 0.01 ) )
		{
			double s_inverse = (1.0 / gsl_vector_get(S, j));
			for( unsigned i = 0; i != V->size1; ++i )
			{
				gsl_matrix_set(V, i, j, s_inverse * gsl_matrix_get(V, i, j) );
			}
		}
		else // fabs( gsl_vector_get(S,i) < 0.01 )
		{
			// S(j) is approximately 0, so V(i,j) * S(j) = 0 for all i
			for( unsigned i = 0; i != V->size1; ++i )
			{
				gsl_matrix_set(V, i, j, 0 );
			}
		}
	}

	// A_inv = (V * W^-1) * U^T
	Matrix A_inv( N, M );
	for( unsigned i = 0; i != N; ++i )
	{
		for( unsigned j = 0; j != M; ++j )
		{
			A_inv(i,j) = 0;
			for( unsigned k = 0; k != N; ++k )
			{
				// Note: we use U(j,k) instead of U(k,j) since we need U^T
				A_inv(i,j) += gsl_matrix_get(V, i, k) * gsl_matrix_get(U, j, k);
			}
		}
	}

	// deallocate U, V, and S
	delete [] data;
	gsl_vector_free( S );
	gsl_matrix_free( V );

	return( A_inv );

#else
	err_quit( "Function Lin_Alg::pseudo_inverse() not implemented:"
			" GNU Scientific Laboratory (GSL) is not installed.\n" );

	// not reached
	Matrix A_inv;
	return( A_inv );
#endif // HAVE_GSL
}

/**
	Solve linear system Ax = b by SVD.

	@param[in] A Square matrix
	@param[in] b Right-hand side solution vector
	@param[out] x Solution vector
 */ 
void
Lin_Alg::solve_system_SVD( const Matrix& A, const Vector& b, Vector& x,
	const SVD_Method& method )
{
#ifdef HAVE_GSL
	// number of rows and columns in a subregion of A
	size_type M = A.sub_row_size();
	size_type N = A.sub_col_size();

	if( b.size() != M )
	{
		err_quit( "solve_system_SVD: Invalid size for b (%u)\n", b.size() );
	}

	// allocate U, V, and S
	// Note the following:
	// 1. convert A to a form that GSL can use: 'matrix views' use the actual
	// data being pointed to, so we should not call 'delete[] data' until the
	// view is no longer used
	// 2. S is a _vector_ of singular values representing a diagonal
	double* A_data = A.to_array();
	gsl_matrix_view A_view = gsl_matrix_view_array( A_data, M, N );
	gsl_matrix* U = &A_view.matrix;
	gsl_vector* S = gsl_vector_alloc( N );
	gsl_matrix* V = gsl_matrix_alloc( N, N );

	// workspace for Golub-Reinsch method (unused for Jacobi method)
	gsl_vector* work = 0;

	// perform singular value decomposition using the given algorithm
	switch( method )
	{
		// computes singular values to higher relative accuracy than Golub-Reinsch
		case Jacobi:
			gsl_linalg_SV_decomp_jacobi( U, V, S );
			break;

		case Golub_Reinsch:
			work = gsl_vector_alloc( N );
			gsl_linalg_SV_decomp( U, V, S, work );
			gsl_vector_free( work );
			break;

		default:
			err_quit( "SVD: Unknown SVD algorithm '%u'\n", method );
			break;
	}

	// copy b into a GSL vector
	double* b_data = b.to_array();
	gsl_vector_view b_view = gsl_vector_view_array( b_data, M );
	gsl_vector* b_gsl = &b_view.vector;

	// create GSL solution vector
	gsl_vector* x_gsl = gsl_vector_alloc( N );

	// solve system using SVD
	gsl_linalg_SV_solve( U, V, S, b_gsl, x_gsl );

	// copy solution into output vector x
	if( x.size() != N )
	{
		x.resize( N, false );
	}
	unsigned ix = 0;
	for( size_type i = x.vec_begin(); i != x.vec_end(); ++i, ++ix )
	{
		x(i) = gsl_vector_get(x_gsl, ix);
	}

	// deallocate GSL matrices and vectors
	delete[] A_data;
	gsl_vector_free( S );
	gsl_matrix_free( V );
	delete[] b_data;
	gsl_vector_free( x_gsl );
#else
	err_warn( "Function Lin_Alg::solve_system_SVD() not implemented:"
			" GNU Scientific Laboratory (GSL) is not installed.\n" );
	return;
#endif // HAVE_GSL
}

/**
	Singular Value Decomposition (SVD).

	A general, rectangular M-by-N matrix A has a singular value decomposition
	(SVD) into the product of an M-by-N orthogonal matrix U, an N-by-N diagonal
	matrix of singular values S, and the transpose of an N-by-N orthogonal square
	matrix V: A = U S V^T 

	The diagonal elements of the singular value matrix S are stored in the
	vector S. The singular values are non-negative and form a non-increasing
	sequence from S_1 to S_N. The matrix V contains the elements of V in
	untransposed form, so to form the product U S V^T, it is necessary to take
	the transpose of V first. The columns of U are the eigenvectors of AA',
	while the columns of V are the eigenvectors of A'A.

	This code follows the GSL function gsl_linalg_SV_decomp_jacobi() by G.
	Jungman, which was itself based on code by Arthur Kosowsky at Rutgers
	University (kosowsky@physics.rutgers.edu).
	Support was added in this implementation for matrix regions, so the SVD of a
	matrix sub-region can be computed instead of for the entire matrix.
	Additionally, the number of non-zero singular values is returned as well.

	This method also implements two versions of computing the 2-norm of a vector:
		|| v || = sqrt( v_1^2, v_2^2, ..., v_n^2 ).
	Version 1 simply computes this directly, while version 2 follows the Fortran
	code in the CBLAS routine dnrm2. Version 1 is faster but less numerically
	accurate than version 2, which is slower than version 1.  
	For example, suppose we have the matrix A(i,j) = i^2 + j^2.
	Version 1 takes 14 seconds while version 2 takes 19 seconds if A is
	500-by-500. The differences in accuracy were observed in the matrix V
	for the columns that corresponded to singular values that are 0 only. Since
	these columns should probably be avoided anyway, it may be preferable to use
	the faster one. Version 1 is enabled by leaving the line '#define
	FAST_NORM' uncommented, while version 2 can be used by commenting this line
	in the function and recompiling.

	The SVD algorithm comes from J.C. Nash, Compact Numerical Methods for
	Computers (New York: Wiley and Sons, 1979), chapter 3.
	See also Algorithm 4.1 in James Demmel, Kresimir Veselic, "Jacobi's Method
	is more accurate than QR", Lapack Working Note 15 (LAWN15), October 1989.
	Available from netlib.

	Another relevant paper is P.P.M. De Rijk, "A One-Sided Jacobi Algorithm for
	computing the singular value decomposition on a vector computer", SIAM
	Journal of Scientific and Statistical Computing, Vol 10, No 2, pp 359-371,
	March 1989.

	@param[in,out] A Matrix to compute SVD where the matrix U replaces the
		contents of A
	@param[out] S Vector of singular values in non-increasing order
	@param[out] V Orthogonal matrix V such that A = U S V^T
	@param[in] tol Tolerance for singular values: singular values within the
		tolerance are set to 0 along with the corresponding columns 
	@retval num_nonzero_values Number of singular values that are not zero, which
		also corresponds to the number of non-zero columns of A
 */
unsigned
Lin_Alg::SVD_jacobi( Matrix& A, Vector& S, Matrix& V, value_type tol )
{
	const size_type M = A.sub_row_size();
	const size_type N = A.sub_col_size();

	if( A.sub_row_size() < A.sub_col_size() )
	{
		// TODO Implement this case
		err_quit( "SVD of MxN matrix with M < N is not implemented\n" );
	}

	value_type tolerance = 10 * M * tol; // readjust tolerance
	// value_type tolerance = tol;          // no readjustment

	// length of vector S must match the second dimension of matrix A
	S.resize(N, false);

	// set V to the identity matrix (V must be square and match the second
	// dimension of A)
	V = Matrix::identity( N );

	// store the column error estimates in S for use during the
	// orthogonalization
	Vector::size_type Sj = S.vec_begin();
	for( size_type j = A.col_begin(); j != A.col_end(); ++j, ++Sj )
	{
		A.push_region( A.row_begin(), j, A.sub_row_size(), 1 );
		S(Sj) = EPS * A.two_norm();
		A.pop_region();
	}

	// initialize rotation counter and sweep counter
	int      rot_counter   = 1;
	unsigned current_sweep = 0;

	// do a minimum of 12 sweeps through the matrix
	unsigned num_sweeps = std::max<unsigned>( 5 * N, 12 );

	// orthogonalize A by plane rotations
	while( rot_counter > 0 && current_sweep <= num_sweeps )
	{
		// initialize rotation counter
		rot_counter = N * (N - 1) / 2;

		// A's columns are in the range [A.col_begin, A.col_end), but the
		// columns of V are in the range [0,N), so we use j and k for A's columns
		// and jj and kk for V's columns
		for( size_type j = A.col_begin(), jj = 0; jj != (N - 1); ++j, ++jj )
		{
			for( size_type k = j + 1, kk = jj + 1; kk != N; ++k, ++kk )
			{
				value_type p = 0;         // dot product of columns j and k
				value_type a = 0, b = 0;  // 2-norm of columns j and k

// uncomment to use faster, less accurate 2-norm; comment to use slower, more
// accurate 2-norm using the method in cblas's dnrm2
#define FAST_NORM

#ifdef FAST_NORM
				for( size_type i = A.row_begin(); i != A.row_end(); ++i )
				{
					// compute dot product of columns j and k
					p += A(i,j) * A(i,k);

					// compute two norm of columns j and k
					a += A(i,j) * A(i,j);
					b += A(i,k) * A(i,k);
				}
				a = sqrt(a);  // finish computing norm
				b = sqrt(b);
#else
				value_type scale_a = 0, scale_b = 0;
				value_type ssq_a   = 1, ssq_b   = 1;
				for( size_type i = A.row_begin(); i != A.row_end(); ++i )
				{
					// compute dot product of columns j and k
					p += A(i,j) * A(i,k);

					// compute two norm of columns j and k using cblas dnrm2 algorithm
					if( A(i,j) != 0 )
					{
						double abs_ij = fabs( A(i,j) );
						if( scale_a < abs_ij )
						{
							value_type frac = (scale_a / abs_ij);
							ssq_a = 1 + ssq_a * (frac * frac);
							scale_a = abs_ij;
						}
						else
						{
							value_type frac = (abs_ij / scale_a);
							ssq_a += (frac * frac);
						}
					}

					if( A(i,k) != 0 )
					{
						double abs_ik = fabs( A(i,k) );
						if( scale_b < abs_ik )
						{
							value_type frac = (scale_b / abs_ik);
							ssq_b = 1 + ssq_b * (frac * frac);
							scale_b = abs_ik;
						}
						else
						{
							value_type frac = (abs_ik / scale_b);
							ssq_b += (frac * frac);
						}
					}

				}
				a = scale_a * sqrt( ssq_a );  // finish computing norm
				b = scale_b * sqrt( ssq_b );
#endif

				p *= 2.0;     // from equation 9a in literature:  p = 2 * x.y

				value_type q = (a * a) - (b * b);
				value_type v = hypot(p, q);

				// test whether columns j and k are orthogonal, or dominant errors
				value_type abserr_a = S(j);
				value_type abserr_b = S(k);

				bool sorted = (a >= b);
				bool orthog = (fabs(p) <= tolerance * (a * b));
				bool noisya = (a < abserr_a);
				bool noisyb = (b < abserr_b);

				if( sorted && (orthog || noisya || noisyb) )
				{
					--rot_counter;
					continue;
				}

				// calculate rotation angles
				value_type cosine = 0, sine = 0;
				if( v == 0 || !sorted )
				{
					cosine = 0.0;
					sine   = 1.0;
				}
				else
				{
					cosine = sqrt( (v + q) / (2.0 * v) );
					sine   = p / (2.0 * v * cosine);
				}

				// apply rotation to A
				for( size_type i = A.row_begin(); i != A.row_end(); ++i )
				{
					// save original values to use in computing the rotation
					const value_type A_ij = A(i,j);
					const value_type A_ik = A(i,k);

					A(i,j) = (A_ij * cosine) + (A_ik * sine);
					A(i,k) = (-A_ij * sine)  + (A_ik * cosine);
				}

				// apply rotation to singular values
				S(j) = (abserr_a * fabs(cosine) + abserr_b * fabs(sine));
				S(k) = (abserr_a * fabs(sine)   + abserr_b * fabs(cosine));

				// apply rotation to V
				for( size_type i = V.row_begin(); i != V.row_end(); ++i )
				{
					// save original values to use in computing the rotation
					const value_type V_ij = V(i,jj);
					const value_type V_ik = V(i,kk);

					V(i,jj) = (V_ij * cosine) + (V_ik * sine);
					V(i,kk) = (-V_ij * sine)  + (V_ik * cosine);
				}
			}
		}

		++current_sweep;  // sweep completed
	}

	// orthogonalization is complete, so compute singular values and apply
	// normalization

	value_type prev_norm = -1.0;  // initialize to some non-zero norm

	unsigned num_nonzero_values = 0;

	Sj  = S.vec_begin();
	size_type col_begin = A.col_begin();
	for( size_type j = A.col_begin(); j != A.col_end(); ++j, ++Sj )
	{
		// get column j of A
		A.push_region( A.row_begin(), j, A.sub_row_size(), 1 );

		// compute norm of column
		value_type norm = A.two_norm();

		// determine if the singular value is zero according to the
		// criteria used in the main loop above (i.e., compare the value
		// with the norm of the previous column)
		// note: once a singular value is found to be too close to 0, then all
		// subsequent singular values are set to 0 as well since the values are
		// ordered in S in non-increasing order
		if( norm == 0 || prev_norm == 0 
				|| (j > col_begin && norm <= (tolerance * prev_norm)) )
		{
			// mark singular value as 0 (matrix A is singular) 
			S(Sj) = 0;

			// annihilate the column
			for( size_type i = A.row_begin(); i != A.row_end(); ++i )
			{
				A(i,j) = 0;
			}

			// ensure that subsequent singular values are also set to 0
			prev_norm = 0;
		}
		else
		{
			// save singular value (matrix A is non-singular so far)
			S(Sj) = norm;
			++num_nonzero_values;

			// normalize the column
			A /= norm;

			// save previous norm
			prev_norm = norm;
		}

		A.pop_region();
	}

	// reached sweep limit
	if( rot_counter > 0 )
	{
		err_quit( "Jacobi iterations did not reach the desired tolerance\n" );
	}

	return( num_nonzero_values );
}

#if false

#ifdef HAVE_GSL

/*
	Factorise a general M x N matrix A into,

	  A = U D V^T

	where U is a column-orthogonal M x N matrix (U^T U = I), 
	D is a diagonal N x N matrix, 
	and V is an N x N orthogonal matrix (V^T V = V V^T = I)

	U is stored in the original matrix A, which has the same size

	V is stored as a separate matrix (not V^T). You must take the
	transpose to form the product above.

	The diagonal matrix D is stored in the vector S,  D_ii = S_i
 */
int
gsl_linalg_SV_decomp( Matrix& A, Vector& S, Matrix& V )
{
	size_t a, b, i, j;

	const size_type M = A.row_size();
	const size_type N = A.col_size();
	const size_type K = min<size_type>(M, N);

	if( A.sub_row_size() < A.sub_col_size() )
	{
		// TODO Implement this case
		err_quit( "SVD of MxN matrix with M < N is not implemented\n" );
	}

	V.resize(N, N, false);
	S.resize(N, false);
	work.resize(N, false);

	/* Handle the case of N = 1 (SVD of a column vector) */

	if( N == 1 )
	{
		value_type norm = A.two_norm();

		S(0) = norm;
		V(0) = 1;

		if( norm != 0.0 )
		{
			A /= norm;
			return( 1 );
		}
		return( 0 );
	}

	work.set_region( 0, K - 1 );
	gsl_vector_view f = gsl_vector_subvector (work, 0, K - 1);

	/* bidiagonalize matrix A, unpack A into U S V */

	gsl_linalg_bidiag_decomp (A, S, &f.vector);
	gsl_linalg_bidiag_unpack2 (A, S, &f.vector, V);

	/* apply reduction steps to B=(S,Sd) */

	chop_small_elements (S, &f.vector);

	/* Progressively reduce the matrix until it is diagonal */

	b = N - 1;

	while (b > 0)
	{
		double fbm1 = gsl_vector_get (&f.vector, b - 1);

		if (fbm1 == 0.0 || gsl_isnan (fbm1))
		{
			b--;
			continue;
		}

		/* Find the largest unreduced block (a,b) starting from b
			and working backwards */

		a = b - 1;

		while (a > 0)
		{
			double fam1 = gsl_vector_get (&f.vector, a - 1);

			if (fam1 == 0.0 || gsl_isnan (fam1))
			{
				break;
			}

			a--;
		}

		{
			const size_t n_block = b - a + 1;
			gsl_vector_view S_block = gsl_vector_subvector (S, a, n_block);
			gsl_vector_view f_block = gsl_vector_subvector (&f.vector, a, n_block - 1);

			gsl_matrix_view U_block =
				gsl_matrix_submatrix (A, 0, a, A->size1, n_block);
			gsl_matrix_view V_block =
				gsl_matrix_submatrix (V, 0, a, V->size1, n_block);

			qrstep (&S_block.vector, &f_block.vector, &U_block.matrix, &V_block.matrix);

			/* remove any small off-diagonal elements */

			chop_small_elements (&S_block.vector, &f_block.vector);
		}
	}
	/* Make singular values positive by reflections if necessary */

	for (j = 0; j < K; j++)
	{
		double Sj = gsl_vector_get (S, j);

		if (Sj < 0.0)
		{
			for (i = 0; i < N; i++)
			{
				double Vij = gsl_matrix_get (V, i, j);
				gsl_matrix_set (V, i, j, -Vij);
			}

			gsl_vector_set (S, j, -Sj);
		}
	}

	/* Sort singular values into decreasing order */

	for (i = 0; i < K; i++)
	{
		double S_max = gsl_vector_get (S, i);
		size_t i_max = i;

		for (j = i + 1; j < K; j++)
		{
			double Sj = gsl_vector_get (S, j);

			if (Sj > S_max)
			{
				S_max = Sj;
				i_max = j;
			}
		}

		if (i_max != i)
		{
			/* swap eigenvalues */
			gsl_vector_swap_elements (S, i, i_max);

			/* swap eigenvectors */
			gsl_matrix_swap_columns (A, i, i_max);
			gsl_matrix_swap_columns (V, i, i_max);
		}
	}

	return GSL_SUCCESS;
}

int 
gsl_linalg_bidiag_decomp (gsl_matrix * A, gsl_vector * tau_U, gsl_vector * tau_V)  
{
  if (A->size1 < A->size2)
    {
      GSL_ERROR ("bidiagonal decomposition requires M>=N", GSL_EBADLEN);
    }
  else if (tau_U->size  != A->size2)
    {
      GSL_ERROR ("size of tau_U must be N", GSL_EBADLEN);
    }
  else if (tau_V->size + 1 != A->size2)
    {
      GSL_ERROR ("size of tau_V must be (N - 1)", GSL_EBADLEN);
    }
  else
    {
      const size_t M = A->size1;
      const size_t N = A->size2;
      size_t i;
  
      for (i = 0 ; i < N; i++)
        {
          /* Apply Householder transformation to current column */
          
          {
            gsl_vector_view c = gsl_matrix_column (A, i);
            gsl_vector_view v = gsl_vector_subvector (&c.vector, i, M - i);
            double tau_i = gsl_linalg_householder_transform (&v.vector);
            
            /* Apply the transformation to the remaining columns */
            
            if (i + 1 < N)
              {
                gsl_matrix_view m = 
                  gsl_matrix_submatrix (A, i, i + 1, M - i, N - (i + 1));
                gsl_linalg_householder_hm (tau_i, &v.vector, &m.matrix);
              }

            gsl_vector_set (tau_U, i, tau_i);            

          }

          /* Apply Householder transformation to current row */
          
          if (i + 1 < N)
            {
              gsl_vector_view r = gsl_matrix_row (A, i);
              gsl_vector_view v = gsl_vector_subvector (&r.vector, i + 1, N - (i + 1));
              double tau_i = gsl_linalg_householder_transform (&v.vector);
              
              /* Apply the transformation to the remaining rows */
              
              if (i + 1 < M)
                {
                  gsl_matrix_view m = 
                    gsl_matrix_submatrix (A, i+1, i+1, M - (i+1), N - (i+1));
                  gsl_linalg_householder_mh (tau_i, &v.vector, &m.matrix);
                }

              gsl_vector_set (tau_V, i, tau_i);
            }
        }
    }
        
  return GSL_SUCCESS;
}

/* Form the orthogonal matrices U, V, diagonal d and superdiagonal sd
   from the packed bidiagonal matrix A */

int
gsl_linalg_bidiag_unpack (const gsl_matrix * A, 
                          const gsl_vector * tau_U, 
                          gsl_matrix * U, 
                          const gsl_vector * tau_V,
                          gsl_matrix * V,
                          gsl_vector * diag, 
                          gsl_vector * superdiag)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  const size_t K = GSL_MIN(M, N);

  if (M < N)
    {
      GSL_ERROR ("matrix A must have M >= N", GSL_EBADLEN);
    }
  else if (tau_U->size != K)
    {
      GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    }
  else if (tau_V->size + 1 != K)
    {
      GSL_ERROR ("size of tau must be MIN(M,N) - 1", GSL_EBADLEN);
    }
  else if (U->size1 != M || U->size2 != N)
    {
      GSL_ERROR ("size of U must be M x N", GSL_EBADLEN);
    }
  else if (V->size1 != N || V->size2 != N)
    {
      GSL_ERROR ("size of V must be N x N", GSL_EBADLEN);
    }
  else if (diag->size != K)
    {
      GSL_ERROR ("size of diagonal must match size of A", GSL_EBADLEN);
    }
  else if (superdiag->size + 1 != K)
    {
      GSL_ERROR ("size of subdiagonal must be (diagonal size - 1)", GSL_EBADLEN);
    }
  else
    {
      size_t i, j;

      /* Copy diagonal into diag */

      for (i = 0; i < N; i++)
        {
          double Aii = gsl_matrix_get (A, i, i);
          gsl_vector_set (diag, i, Aii);
        }

      /* Copy superdiagonal into superdiag */

      for (i = 0; i < N - 1; i++)
        {
          double Aij = gsl_matrix_get (A, i, i+1);
          gsl_vector_set (superdiag, i, Aij);
        }

      /* Initialize V to the identity */

      gsl_matrix_set_identity (V);

      for (i = N - 1; i > 0 && i--;)
        {
          /* Householder row transformation to accumulate V */
          gsl_vector_const_view r = gsl_matrix_const_row (A, i);
          gsl_vector_const_view h = 
            gsl_vector_const_subvector (&r.vector, i + 1, N - (i+1));
          
          double ti = gsl_vector_get (tau_V, i);
          
          gsl_matrix_view m = 
            gsl_matrix_submatrix (V, i + 1, i + 1, N-(i+1), N-(i+1));
          
          gsl_linalg_householder_hm (ti, &h.vector, &m.matrix);
        }

      /* Initialize U to the identity */

      gsl_matrix_set_identity (U);

      for (j = N; j > 0 && j--;)
        {
          /* Householder column transformation to accumulate U */
          gsl_vector_const_view c = gsl_matrix_const_column (A, j);
          gsl_vector_const_view h = gsl_vector_const_subvector (&c.vector, j, M - j);
          double tj = gsl_vector_get (tau_U, j);
          
          gsl_matrix_view m = 
            gsl_matrix_submatrix (U, j, j, M-j, N-j);
          
          gsl_linalg_householder_hm (tj, &h.vector, &m.matrix);
        }

      return GSL_SUCCESS;
    }
}

double
gsl_linalg_householder_transform (gsl_vector * v)
{
  /* replace v[0:n-1] with a householder vector (v[0:n-1]) and
     coefficient tau that annihilate v[1:n-1] */

  const size_t n = v->size ;

  if (n == 1)
    {
      return 0.0; /* tau = 0 */
    }
  else
    { 
      double alpha, beta, tau ;
      
      gsl_vector_view x = gsl_vector_subvector (v, 1, n - 1) ; 
      
      double xnorm = gsl_blas_dnrm2 (&x.vector);
      
      if (xnorm == 0) 
        {
          return 0.0; /* tau = 0 */
        }
      
      alpha = gsl_vector_get (v, 0) ;
      beta = - (alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, xnorm) ;
      tau = (beta - alpha) / beta ;
      
      gsl_blas_dscal (1.0 / (alpha - beta), &x.vector);
      gsl_vector_set (v, 0, beta) ;
      
      return tau;
    }
}

int
gsl_linalg_householder_hm (double tau, const gsl_vector * v, gsl_matrix * A)
{
  /* applies a householder transformation v,tau to matrix m */

  size_t i, j;

  if (tau == 0.0)
    {
      return GSL_SUCCESS;
    }

#ifdef USE_BLAS
  {
    gsl_vector_const_view v1 = gsl_vector_const_subvector (v, 1, v->size - 1);
    gsl_matrix_view A1 = gsl_matrix_submatrix (A, 1, 0, A->size1 - 1, A->size2);
    
    for (j = 0; j < A->size2; j++)
      {
        double wj = 0.0;
        gsl_vector_view A1j = gsl_matrix_column(&A1.matrix, j);
        gsl_blas_ddot (&A1j.vector, &v1.vector, &wj);
        wj += gsl_matrix_get(A,0,j);

        {
          double A0j = gsl_matrix_get (A, 0, j);
          gsl_matrix_set (A, 0, j, A0j - tau *  wj);
        }

        gsl_blas_daxpy (-tau * wj, &v1.vector, &A1j.vector);
      }
  }
#else
  for (j = 0; j < A->size2; j++)
    {
      /* Compute wj = Akj vk */

      double wj = gsl_matrix_get(A,0,j);  

      for (i = 1; i < A->size1; i++)  /* note, computed for v(0) = 1 above */
        {
          wj += gsl_matrix_get(A,i,j) * gsl_vector_get(v,i);
        }

      /* Aij = Aij - tau vi wj */

      /* i = 0 */
      {
        double A0j = gsl_matrix_get (A, 0, j);
        gsl_matrix_set (A, 0, j, A0j - tau *  wj);
      }

      /* i = 1 .. M-1 */

      for (i = 1; i < A->size1; i++)
        {
          double Aij = gsl_matrix_get (A, i, j);
          double vi = gsl_vector_get (v, i);
          gsl_matrix_set (A, i, j, Aij - tau * vi * wj);
        }
    }
#endif

  return GSL_SUCCESS;
}

int
gsl_linalg_householder_mh (double tau, const gsl_vector * v, gsl_matrix * A)
{
  /* applies a householder transformation v,tau to matrix m from the
     right hand side in order to zero out rows */

  size_t i, j;

  if (tau == 0)
    return GSL_SUCCESS;

  /* A = A - tau w v' */

#ifdef USE_BLAS
  {
    gsl_vector_const_view v1 = gsl_vector_const_subvector (v, 1, v->size - 1);
    gsl_matrix_view A1 = gsl_matrix_submatrix (A, 0, 1, A->size1, A->size2-1);

    for (i = 0; i < A->size1; i++)
      {
        double wi = 0.0;
        gsl_vector_view A1i = gsl_matrix_row(&A1.matrix, i);
        gsl_blas_ddot (&A1i.vector, &v1.vector, &wi);
        wi += gsl_matrix_get(A,i,0);  
        
        {
          double Ai0 = gsl_matrix_get (A, i, 0);
          gsl_matrix_set (A, i, 0, Ai0 - tau *  wi);
        }
        
        gsl_blas_daxpy(-tau * wi, &v1.vector, &A1i.vector);
      }
  }
#else
  for (i = 0; i < A->size1; i++)
    {
      double wi = gsl_matrix_get(A,i,0);  

      for (j = 1; j < A->size2; j++)  /* note, computed for v(0) = 1 above */
        {
          wi += gsl_matrix_get(A,i,j) * gsl_vector_get(v,j);
        }
      
      /* j = 0 */
      
      {
        double Ai0 = gsl_matrix_get (A, i, 0);
        gsl_matrix_set (A, i, 0, Ai0 - tau *  wi);
      }

      /* j = 1 .. N-1 */
      
      for (j = 1; j < A->size2; j++) 
        {
          double vj = gsl_vector_get (v, j);
          double Aij = gsl_matrix_get (A, i, j);
          gsl_matrix_set (A, i, j, Aij - tau * wi * vj);
        }
    }
#endif

  return GSL_SUCCESS;
}

#endif // HAVE_GSL

#endif

/**
	Print contents of pivot.
	@param[in] pivot Vector of pivots
 */ 
void
Lin_Alg::print_pivot( const Pivot& pivot )
{
   // print nothing if empty
   if( pivot.size() == 0 )
   {
      return;
   }
                                                                                
   // see 10 digits of number
   cout << std::setprecision( 10 );
                                                                                
   // print initial element
   Pivot::const_iterator iter = pivot.begin();
   cout << "[ " << *iter;
                                                                                
   // print each element in vector
   for( iter++;
      iter != pivot.end();
      iter++ )
   {
      cout << ", " << *iter;
   }
   cout << " ]" << endl;
}

/**
	Approximate the dominant eigenvalue and an associated eigenvector for the
	matrix A, given a nonzero vector x.
	Note: Eigenvalue is the returned value, and the eigenvector overwrites x.

	@param[in,out] x Output eigenvector
	@param[in] tolerance Threshold in which to stop computation when error
		reduces to less than
	@param[in] max_iterations Maximum number of iterations to perform in
		algorithm
	@return eigenvalue Eigenvalue of matrix with the largest magnitude
 */
value_type
Lin_Alg::power_method( Matrix& A, Matrix& x, value_type tolerance,
	unsigned max_iterations )
{
	value_type eigenvalue = 0;
	size_type p;

	inf_norm( x, p );
	x = x / x( p, 0 );

	for( size_type k = 0; k != max_iterations; k++ )
	{
		Matrix y = A * x;
		eigenvalue = y(p,0);

		// find the smallest integer p with |y_p| = ||y||_inf
		inf_norm( y, p );

		// if || y ||_2 = 0, then quit
		if( y(p,0) == 0 )
		{
		/*
			std::cerr << "A has eigenvalue 0, select new vector x and restart" 
				<< std::endl;
		 */
			return( eigenvalue );
		}

		// error = || x - ( y / || y ||_2 ) ||_2
		value_type error = ( x - y / y(p,0) ).inf_norm();

		// x = y / || y ||_2
		x = y / y(p,0);

		if( error < tolerance )
		{
			return( eigenvalue );
		}
	}

	err_warn( "Maximum iterations of %u was exceeded with tolerance of %lf\n",
		max_iterations, tolerance );

	return( eigenvalue );
}

/**
// TODO: Does not work yet! (worked on it on 8/15/05)
	Approximate the dominant eigenvalue and an associated eigenvector for the
	symmetric matrix A, given a nonzero vector x.
	Note: Eigenvalue is the returned value, and the eigenvector overwrites x.

	@param[in,out] x Output eigenvector
	@param[in] tolerance Threshold in which to stop computation when error
		reduces to less than
	@param[in] max_iterations Maximum number of iterations to perform in
		algorithm
 */
value_type
Lin_Alg::inverse_power_method( Matrix& A, Vector& x,
	value_type tolerance, unsigned max_iterations )
{
	value_type eigenvalue = 0;
	size_type p = 0;

	// q = (x^t)Ax / (x^t)x
	Matrix numerator = (x * (A * x));
	Matrix denominator = x * x;

	if( denominator(0, 0) == 0 )
	{
		return( eigenvalue );
	}
value_type q = numerator(0,0) / denominator(0,0);

	// find smallest p s.t. p = |x|_inf
	value_type norm = x.inf_norm();
	for( unsigned i = 0; i != x.size(); ++i )
	{
		if( x(i) == norm )
		{
			p = i;
			break;
		}
	}
	x = x / x(p);

	for( unsigned k = 0; k != max_iterations; k++ )
	{
		// Solve (A - qI)y = x.
		Matrix A_qI = A - (q * Matrix::identity( A.row_size() ));
		Vector y;
		Lin_Alg::solve_system( A_qI, x, y );

		//Matrix y( x.row_size(), 1 );

		// if system does not have a unique solution, then stop
		if( y(p) == 0 )
		{
			err_warn( "System does not have a unique solution\n" );
			return( eigenvalue );
		}

		eigenvalue = y(p);

		// find smallest p s.t. p = |y|_inf
		value_type norm = y.inf_norm();
		for( unsigned i = 0; i != y.size(); ++i )
		{
			if( y(i) == norm )
			{
				p = i;
				break;
			}
		}

		// error = || x - ( y / || y ||_2 ) ||_2
		value_type error = ( x - y / y(p) ).inf_norm();

		// x = y / || y ||_2
		x = y / y(p);

		if( error < tolerance )
		{
			return( eigenvalue );
		}
	}

	err_warn( "Maximum iterations of %u was exceeded with tolerance of %lf\n",
		max_iterations, tolerance );

	return( eigenvalue );
}

/**
	Compute infinity norm of matrix and save row that it was found in for use
	in the power method.
	@param[in] A
	@param[out] row_position Row with largest infinity norm
	@return max_value
 */
value_type
Lin_Alg::inf_norm( const Matrix& A, size_type& row_position )
{
	value_type max_value = 0;
	row_position = 0;
	for( size_type i = 0; i != A.row_size(); i++ )
	{
		value_type row_sum = 0;
		for( size_type j = 0; j != A.col_size(); j++ )
		{
			row_sum += fabs( A(i, j) );
		}

		// compute new max and new row position
		if(  row_sum > max_value )
		{
			max_value = row_sum;
			row_position = i;
		}
	}

	// return location of maximum element
	return( max_value );
}

/**
	Approximate the dominant eigenvalue and an associated eigenvector for the
	symmetric matrix A, given a nonzero vector x.
	Note: Eigenvalue is the returned value, and the eigenvector overwrites x.

	@param[in,out] x Output eigenvector
	@param[in] tolerance Threshold in which to stop computation when error
		reduces to less than
	@param[in] max_iterations Maximum number of iterations to perform in
		algorithm
 */
value_type
Lin_Alg::symmetric_power_method( Matrix& A, Matrix& x,
	value_type tolerance, unsigned max_iterations )
{
	value_type eigenvalue = 0;

	Matrix x_transpose = x;
	x_transpose.transpose();

	for( unsigned k = 0; k != max_iterations; k++ )
	{
		Matrix y = A * x;
		eigenvalue = (x_transpose * y)(0,0);

		// if || y ||_2 = 0, then quit
		if( y.two_norm() == 0 )
		{
			err_warn( "A has eigenvalue 0, select new vector x and restart\n" );
			return( eigenvalue );
		}

		// error = || x - ( y / || y ||_2 ) ||_2
		value_type error = ( x - y / y.two_norm() ).two_norm();

		// x = y / || y ||_2
		x = y / y.two_norm();

		if( error < tolerance )
		{
			return( eigenvalue );
		}
	}

	err_warn( "Maximum iterations of %u was exceeded with tolerance of %lf\n",
		max_iterations, tolerance );

	return( eigenvalue );
}

/**
	Obtain principal-component transformation matrix and sample mean of
	the vector population contained in the rows of X, a matrix of size K-by-n
	where K is the number of vectors and n is their dimensionality. Q, with
	values in the range [0, n], is the number of eigenvectors used in
	constructing the principal-components transformation matrix.

	@param[in] Q Number of components to use
	@param[in,out] X Matrix representing rows of data--this matrix is altered
		upon returning such that the mean is subtracted from each column
	@param[out] A Q-by-n principal components transformation matrix whose rows
		are the Q eigenvectors of Cx corresponding to the Q largest eigenvalues.
	@param[out] mean n-by-1 mean vector of the population in X
 */
void
Lin_Alg::prin_comp_matrix( unsigned Q, Matrix& X, Matrix& A, Vector& mean )
{
	if( Q > X.sub_col_size() )
	{
		err_quit( "Matrix::prin_comp_matrix: The number of principal components"
				" requested (%u) exceeds the number of data items (%u)\n",
				Q, X.sub_col_size() );
	}

	// obtain mean vector and covariance matrix of the vectors in X
	Matrix cov;
	Lin_Alg::covariance( X, cov, mean );

	// obtain the eigenvectors and corresponding eigenvalues of cov. The
	// eigenvectors are the columns of n-by-n matrix V. D is an n-by-n diagonal
	// matrix whose elements along the main diagonal are the eigenvalues
	// corresponding to the eigenvectors in V, so that X*V = D*V.
	// Note: The eigenvalues are sorted in decreasing order with the
	// eigenvectors rearranged to match.
	Matrix V;
	Vector D;
	Lin_Alg::eigen( cov, D, V );

	// form the q rows of A from first q columns of V
	V.set_region( V.row_begin(), V.col_begin(), V.sub_row_size(), Q );
	A = V.transpose();
}

/**
	Computes the covariance matrix 'cov' and the mean vector 'mean' of a vector
	population organized as the rows of this matrix. 'cov' is of size N-by-N and
	'mean' is of size N-by-1, where N is the dimension of the vectors (i.e., the
	number of columns of X).  

	@param[in,out] M Matrix representing rows of data--this matrix is altered
		upon returning such that the mean is subtracted from each column
	@param[out] cov Covariance matrix
	@param[out] mean Mean/average of each column of this matrix
 */
void
Lin_Alg::covariance( Matrix& M, Matrix& cov, Vector& mean )
{
	// compute an unbiased estimate of mean for each column
	mean.resize( M.sub_col_size(), false );
	for( size_type j = M.col_begin(), k = mean.vec_begin();
			j != M.col_end(); ++j, ++k )
	{
		unsigned   item_num = 0;   // current item number
		value_type unused_var = 0; // unused, but we must pass a variance pointer

		for( size_type i = M.row_begin(); i != M.row_end(); ++i )
		{
			ws_tools::stat( M(i,j), item_num++, &mean(k), &unused_var );
		}
	}

	// for each column, subtract the column's mean from each component of column
	for( size_type j = M.col_begin(), k = mean.vec_begin();
			j != M.col_end(); ++j, ++k )
	{
		for( size_type i = M.row_begin(); i != M.row_end(); ++i )
		{
			M(i,j) -= mean(k);
		}
	}

	// compute an unbiased estimate of 'cov': X' * X where X = *this
	// Note that the product is X'*X since the population vectors are the rows of
	// X, and this equivalent to the following inefficient code:
	// cov = (transpose() * (*this)) / (sub_row_size() - 1);
	// which does not make use of the symmetry of X' * X and creates unnecessary
	// temporary objects.

	// same as above but avoids computations by computing multiplication
	// directly and makes use of the symmetry of the covariance matrix
	cov.resize( M.sub_col_size(), M.sub_col_size(), false );

	// for each row and column
	for( size_type j = M.col_begin(); j != M.col_end(); ++j )
	{
		for( size_type i = j; i != M.col_end(); ++i )
		{
			// compute dot product of a row in X' with a column in X, but using
			// only X (i.e., without explicitly computing X') 
			value_type sum = 0;
			for( size_type k = M.row_begin(); k != M.row_end(); ++k )
			{
				sum += M(k, j) * M(k, i);
			}

			// store value (i,j) (note that the covariance matrix is symmetric)
			cov(i, j) = cov(j, i) = sum;
		}
	}

	// finish computation of the covariance matrix (the -1 in the division
	// provides an unbiased estimation) 
	cov /= (M.sub_row_size() - 1);
}

/**
	Project data Y onto basis vectors using principal-component analysis.

	@param[in] A Q-by-n principal components transformation matrix whose rows
		are the Q eigenvectors of Cx corresponding to the Q largest eigenvalues.
	@param[in] mean n-by-1 mean vector of the population in X
	@param[in,out] Y Data to project
 */
void
Lin_Alg::prin_comp_project( const Matrix& A, const Vector& mean, Vector& Y )
{
	if( Y.sub_size() != mean.sub_size() )
	{
		err_quit( "Matrix::prin_comp_project: The dimension of the input vector"
				" (%u) does not match the dimension of the mean vector (%u)\n",
				Y.sub_size(), mean.sub_size() );
	}
	if( Y.sub_size() != A.sub_col_size() )
	{
		err_quit( "Matrix::prin_comp_project: The dimension (%u) of the input"
				" vector does not match the column dimension (%u) of the"
				" transformation matrix\n", Y.sub_size(), mean.sub_size() );
	}

	// project data onto basis vectors
	Y = A * Matrix(Y - mean);
}

/**
	Project data Y onto basis vectors using principal-component analysis.

	@param[in] A Q-by-n principal components transformation matrix whose rows
		are the Q eigenvectors of Cx corresponding to the Q largest eigenvalues.
	@param[in] mean n-by-1 mean vector of the population in X
	@param[in,out] Y Data to project
 */
void
Lin_Alg::prin_comp_project( const Matrix& A, const Vector& mean, Matrix& Y )
{
	if( Y.sub_col_size() != mean.sub_size() )
	{
		err_quit( "Matrix::prin_comp_project: The dimension of the input vector"
				" (%u) does not match the dimension of the mean vector (%u)\n",
				Y.sub_col_size(), mean.sub_size() );
	}
	if( Y.sub_col_size() != A.sub_col_size() )
	{
		err_quit( "Matrix::prin_comp_project: The dimension (%u) of the input"
				" vector does not match the column dimension (%u) of the"
				" transformation matrix\n", Y.sub_col_size(), mean.sub_size() );
	}

	Matrix Mean( A.sub_row_size(), mean.sub_size() );
	Vector::size_type k = mean.vec_begin();
	for( size_type j = Mean.col_begin(); j != Mean.col_end(); ++j )
	{
		for( size_type i = Mean.row_begin(); i != Mean.row_end(); ++i )
		{
			Mean(i,j) = mean(k);
		}
		++k;
	}

	// project data onto basis vectors
	Y = A * (Y.transpose() - Mean);
}

/**
	Obtain principal-component vectors and related quantities.
	P = princomp(M, Q) computes the principal-component vectors of the vector
	population contained in the rows of M, a matrix of size K-by-n where K is
	the number of vectors and n is their dimensionality. Q, with values in the
	range [0, n], is the number of eigenvectors used in constructing the
	principal-components transformation matrix. P is a structure with the
	following fields:

	P.Y   K-by-Q matrix whose columns are the principal-component vectors.

	P.A   Q-by-n principal components transformation matrix whose rows are the Q
			eigenvectors of Cx corresponding to the Q largest eigenvalues.

	P.X   K-by-n matrix whose rows are the vectors reconstructed from the
			principal-component vectors. P.X and P.Y are identical if Q = n.

	P.ems The mean square error incurred in using only the Q eigenvectors
			corresponding to the largest eigenvalues. P.ems is 0 if Q = n.

	P.Cx  The n-by-n covariance matrix of the population in X.

	P.mx  The n-by-1 mean vector of the population in X.

	P.Cy  The Q-by-Q covariance matrix of the population in Y. The main diagonal
			contains the eigenvalues (in descending order) corresponding to the Q
			eigenvectors.
 */
void
Lin_Alg::prin_comp( unsigned Q, Matrix& M )
{
	if( Q > M.sub_col_size() )
	{
		err_quit( "Matrix::prin_comp_matrix: The number of principal components"
				" requested (%u) exceeds the number of data items (%u)\n",
				Q, M.sub_col_size() );
	}

	// [K, n] = size(X);
	// X = double(X);

	// obtain mean vector and covariance matrix of the vectors in X
	// [P.Cx, P.mx] = covmatrix(X);
	// P.mx = P.mx';
	Matrix Cx;
	Vector mx;
	covariance( M, Cx, mx );

	// obtain the eigenvectors and corresponding eigenvalues of Cx. The
	// eigenvectors are the columns of n-by-n matrix V. D is an n-by-n diagonal
	// matrix whose elements along the main diagonal are the eigenvalues
	// corresponding to the eigenvectors in V, so that X*V = D*V.
	// Note: The eigenvalues are sorted in decreasing order with the
	// eigenvectors rearranged to match.
	// [V, D] = eig( P.Cx );
	// d = diag(D);
	// [d, idx] = sort(D);
	// d = flipud(d);
	// idx = flipud(idx);
	// D = diag(d);
	// V = V(:, idx);
	Matrix V;
	Vector D;
	Lin_Alg::eigen( Cx, D, V );

	fprintf( stderr, "\n" );
	V.write( "", 5 );
	fprintf( stderr, "\n" );
	D.write( "", 5 );
	fprintf( stderr, "\n" );

	// form the q rows of A from first q columns of V
	// P.A = V(:, 1:q)';
	V.set_region( V.row_begin(), V.col_begin(), V.sub_row_size(), Q );
	Matrix A = V.transpose();

	// compute principal component vectors (note that the mean is already
	// subtracted from *this due to the call to covariance()
	// Mx = repmat(P.mx, K, 1);
	// P.Y = P.A * (X - Mx)';
	Matrix Y = A * M.transpose();

	// obtain reconstructed matrix
	// P.X = (P.A' * P.Y)' + Mx;
	Matrix Mx( M.sub_row_size(), mx.sub_size() );
	Vector::size_type k = mx.vec_begin();
	for( size_type j = Mx.col_begin(); j != Mx.col_end(); ++j )
	{
		for( size_type i = Mx.row_begin(); i != Mx.row_end(); ++i )
		{
			Mx(i,j) = mx(k);
		}
		++k;
	}
	Matrix X = (A.transpose() * Y).transpose() + Mx;
	X.write( "", 4 );

	// convert P.Y to K-by-q array and P.mx to n-by-1 vector
	// P.Y = P.Y';
	// P.mx = P.mx';
	Y = Y.transpose();

	// mean square error is given by the sum of all the eigenvalues minus the
	// sum of the q largest eigenvalues
	// d = diag(D);
	// P.ems = sum( d( q + 1: end ) );
	double ems = 0;
	D.set_region( Q, D.sub_size() - Q );
	for( size_type i = D.vec_begin(); i != D.vec_end(); ++i )
	{
		ems += D(i);
	}

	// covariance matrix of the Y's
	// P.Cy = P.A * P.Cx * P.A';
	Matrix Cy = A * Cx * A.transpose();
}

/**
	Find eigenvalues of symmetric matrix A using GSL library.
	@param[in] A Symmetric matrix
	@param[out] eig_val Eigenvalues of A that are sorted in descending order
	@param[out] eig_vec Eigenvectors of A ordered such that each column i of 
		eig_vec corresponds to the i_th eigenvalue in eig_val, which is sorted
 */
void
Lin_Alg::eigen( const Matrix& A, Vector& eig_val, Matrix& eig_vec )
{
#ifdef HAVE_GSL
	if( A.sub_row_size() != A.sub_col_size() )
	{
		err_warn( "eigen: Invalid matrix dimension: must be symmetric\n" );
		return;
	}

	// convert A to a form that GSL can use: 'matrix views' use the actual data
	// being pointed to, so we should not call 'delete data' until view is no
	// longer used
	double* data = A.to_array();
	gsl_matrix_view m = gsl_matrix_view_array( data, A.sub_row_size(),
			A.sub_col_size());

	// allocate space for eigenvalues, eigenvectors, and a workspace
	gsl_vector* e_val = gsl_vector_alloc( A.sub_row_size() );
	gsl_matrix* e_vec = gsl_matrix_alloc( A.sub_row_size(), A.sub_row_size() );
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc( A.row_size() );

	// find eigenvalues and eigenvectors of A
	gsl_eigen_symmv( &m.matrix, e_val, e_vec, w );

	// deallocate matrix view and workspace since they are no longer used
	delete[] data;
	gsl_eigen_symmv_free( w );

	// sort eigenvalues and eigenvectors by absolute value in descending order
	gsl_eigen_symmv_sort( e_val, e_vec, GSL_EIGEN_SORT_ABS_DESC );

	// save eigenvalues
	eig_val.resize( A.sub_row_size(), false );
	unsigned k = 0;
	for( Vector::size_type i = eig_val.vec_begin(); i != eig_val.vec_end(); ++i )
	{
		eig_val(i) = gsl_vector_get( e_val, k );
		++k;
	}
	gsl_vector_free( e_val );

	// save eigenvectors
	eig_vec.resize( A.sub_row_size(), A.sub_col_size(), false );
	k = 0;
	for( Matrix::size_type i = eig_vec.row_begin(); i != eig_vec.row_end(); ++i )
	{
		unsigned l = 0;
		for( Matrix::size_type j = eig_vec.col_begin();
				j != eig_vec.col_end(); ++j )
		{
			eig_vec(i, j) = gsl_matrix_get( e_vec, k, l );
			++l;
		}
		++k;
	}
	gsl_matrix_free( e_vec );

#else
	err_warn( "Function Lin_Alg::eigen() not implemented:"
			" GNU Scientific Laboratory (GSL) is not installed.\n" );
	return;
#endif // HAVE_GSL
}

/**
	Find eigenvalues of symmetric matrix A.
	@param[in,out] A Symmetric matrix
 */
Vector Lin_Alg::eig_values2( const Matrix& A )
{
	Matrix A_copy = A;
	Lin_Alg::householder( A_copy );
	Vector eigenvalues = Lin_Alg::QR( A_copy );
	return( eigenvalues );
}

/**
	Convert symmetric matrix A into a similar symmetric tridiagonal matrix.
	@param[in,out] A Symmetric matrix
 */
void Lin_Alg::householder( Matrix& A )
{
	if( A.row_size() != A.col_size() )
	{
		err_quit( "Invalid matrix dimension--must be symmetric\n" );
	}
	if( A.row_size() < 2 )
	{
		err_quit( "Invalid matrix dimension--must be at least 2x2\n" );
	}

	unsigned n = A.row_size();

	// Step 1:
	for( unsigned k = 0; k != n - 2; ++k )
	{
		// Step 2: q = sum from j = k to n of A(j,k)^2
		double q = 0;
		for( unsigned j = k + 1; j != n; ++j )
		{
			q += A(j, k) * A(j, k);
		}

		// Step 3:
		double alpha;
		if( A(k + 1, k) == 0 )
		{
			alpha = -sqrt( q );
		}
		else
		{
			alpha = -sqrt( q ) * ( A(k + 1, k) / fabs( A(k + 1, k) ) );
		}

		// Step 4:
		double RSQ = (alpha * alpha) - (alpha * A(k + 1, k));

		// Step 5:
		Vector v( n );
		v(k) = 0;
		v(k + 1) = A(k + 1, k) - alpha;
		for( unsigned j = k + 2; j != n; j++ )
		{
			v(j) = A(j, k);
		}

		// Step 6:
		Vector u( n );
		double one_over_RSQ = (1 / RSQ);
		for( unsigned j = k; j != n; j++ )
		{
			double sum = 0;
			for( unsigned i = k + 1; i != n; i++ )
			{
				sum += A(j, i) * v(i);
			}
			u(j) = one_over_RSQ * sum;
		}

		// Step 7:
		double PROD = 0;
		for( unsigned i = k + 1; i != n; i++ )
		{
			PROD += v(i) * u(i);
		}

		// Step 8:
		Vector z( n );
		for( unsigned j = k; j != n; j++ )
		{
			z(j) = u(j) - (PROD / (2 * RSQ) ) * v(j);
		}

		// Step 9:
		for( unsigned l = k + 1; l != n - 1; ++l )
		{
			// Step 10:
			for( unsigned j = l + 1; j != n; j++ )
			{
				A(j, l) = A(j, l) - (v(l) * z(l)) - (v(j) - z(l));
				A(l, j) = A(j, l);
			}

			// Step 11:
			A(l, l) = A(l, l) - (2 * v(l) * z(l));
		}

		// Step 12:
		A(n - 1, n - 1) = A(n - 1, n - 1) - (2 * v(n - 1) * z(n - 1));

		// Step 13:
		for( unsigned j = k + 2; j != n; ++j )
		{
			A(k, j) = A(j, k) = 0;
		}

		// Step 14:
		A(k + 1, k) = A(k + 1, k) - (v(k + 1) * z(k));
		A(k, k + 1) = A(k + 1, k);

		// Step 15:
		// the process is complete
	}
}

/**
	Perform QR algorithm to find eigenvalues of symmetric tridiagonal matrix.
	Method comes from Algorithm 9.6 in "Numerical Analysis" on 
		pages 592 - 594.

	Note: Diagonal of matrix is stored in vector a, and one of the off 
		diagonals is stored in b.

	@param[in] tolerance Threshold in which to stop computation when error
		reduces to less than
	@param[in] max_iterations Maximum number of iterations to perform in
		algorithm
 */
Vector 
Lin_Alg::QR( const Matrix& A, value_type tolerance, unsigned max_iterations )
{
	// ensure that matrix is at least square and 3x3
	if( A.row_size() != A.col_size() )
	{
		err_quit( "Invalid matrix dimension--must be symmetric, tridiagonal\n" );
	}
	else if( A.row_size() < 3 )
	{
		err_quit( "Invalid matrix dimension (%ux%u)--must be at at least 3x3\n",
			A.row_size(), A.col_size() );
	}

	// save frequently used size--sizes are actually n + 1, so indexing is okay
	size_type n = A.row_size();

	// eigenvalues and number of eigenvalues currently found
	Vector eigenvalues( n );
	unsigned eig_found = 0;

	// store diagonal and off-diagonal separately--this is less efficient but
	// is simpler to implement and does not overwrite this matrix
	Vector a = initialize_diagonal( A );
	Vector b = initialize_off_diagonal( A );

	// arrays to hold intermediate results
	Vector x( n + 1 );
	Vector y( n + 1 );
	Vector z( n + 1 );
	Vector c( n + 1 );
	Vector s( n + 1 );
	Vector q( n + 1 );

	// initialize accumulated shift
	value_type accumulated_shift = 0;

	for( unsigned k = 0; k != max_iterations; k++ )
	{
		// if |b_n| <= TOL
		if( fabs( b[ n ] ) <= tolerance )
		{
			// save eigenvalue
			eigenvalues[ eig_found++ ] = a[ n ] + accumulated_shift;
			--n;
		}

		// if |b_2| <= TOL
		if( fabs( b[ 2 ] ) <= tolerance )
		{
			// save eigenvalue
			eigenvalues[ eig_found++ ] = a[ 1 ] + accumulated_shift;	
			--n;

			// this is set outside loop since b array is updated fewer times
			a[ 1 ] = a[ 2 ];

			// shift data up one element
			for( size_type j = 2; j <= n; j++ )
			{
				a[ j ] = a[ j + 1 ];
				b[ j ] = b[ j + 1 ];
			}
		}

		if( n == 0 )
		{
			return( eigenvalues );
		}

		if( n == 1 )
		{
			eigenvalues[ eig_found++ ] = a[ 1 ] + accumulated_shift;
			return( eigenvalues );
		}
		// determine if matrix needs to be split
		for( size_type j = 3; j <= n - 1; j++ )
		{
			// if |b[j]| <= TOL, split
			if( fabs( b[ j ] ) <= tolerance )
			{
				/* probably throw an exception here to split matrix */
				return( eigenvalues );
			}
		}

		Vector sub_eigenvalues = compute_2_by_2_eigenvalues( a, b, n );
		value_type mu_1 = sub_eigenvalues[ 0 ];
		value_type mu_2 = sub_eigenvalues[ 1 ];

		// if the entire matrix left is 2x2, then we've computed the last
		// eigenvalues
		if( n == 2 )
		{
			eigenvalues[ eig_found++ ] = mu_1 + accumulated_shift;
			eigenvalues[ eig_found++ ] = mu_2 + accumulated_shift;
			return( eigenvalues );
		}

		/*
			choose s s.t. |s - a_n| = min{ |mu_1 - a_n|, |mu_2 - a_n| }
			i.e., choose eigenvalue closest to last diagonal element 
		 */
		value_type current_shift = 
			( fabs( mu_1 - a[ n ] ) <= fabs( mu_2 - a[ n ] ) ) ?  mu_1 : mu_2;

		accumulated_shift += current_shift;

		// perform shift
		for( size_type j = 1; j <= n; j++ )
		{
			a[ j ] -= current_shift;
		}

		x[ 1 ] = a[ 1 ];
		y[ 1 ] = b[ 2 ];

		for( size_type j = 2; j <= n; j++ )
		{
			// z[ j-1 ] = ( x[ j-1 ]^2 + b[ j ]^2 )^1/2
			z[ j - 1 ] =  sqrt( (x[ j - 1 ] * x[ j - 1 ]) + (b[ j ] * b[ j ]) );

			c[ j ]     =  x[ j - 1 ] / z[ j - 1 ];

			s[ j ]     =  b[ j ] / z[ j - 1 ];

			q[ j - 1 ] =  c[ j ] * y[ j - 1 ] + s[ j ] * a[ j ];

			x[ j ]     = -s[ j ] * y[ j - 1 ] + c[ j ] * a[ j ];

			if( j != n )
			{
				// technically part of algorithm, although it's not used anywhere
				// r[ j - 1 ] = s[ j ] * b[ j + 1 ];

				y[ j ]  =  c[ j ] * b[ j + 1 ];
			}
		}

		/*
			compute A^(k+1)
		 */

		// compute first diagonal element and off-diagonal element
		z[ n ] = x[ n ];
		a[ 1 ] = s[ 2 ] * q[ 1 ] + c[ 2 ] * z[ 1 ];
		b[ 2 ] = s[ 2 ] * z[ 2 ];

		// compute rest of elements
		for( size_type j = 2; j <= n - 1; j++ )
		{
			a[ j ]     = s[ j + 1 ] * q[ j ] + c[ j ] * c[ j + 1 ] * z[ j ];
			b[ j + 1 ] = s[ j + 1 ] * z[ j + 1 ];
		}

		// compute last diagonal element
		a[ n ] = c[ n ] * z[ n ];

		//Matrix A = to_matrix( a, b );
		//A.write();
	}

	err_warn( "Maximum iterations of %u was exceeded with tolerance of %lf\n",
		max_iterations, tolerance );

	return( eigenvalues );
}

/**
	Initialize diagonal elements for use in QR method.	
	@return a Vector
 */
Vector
Lin_Alg::initialize_diagonal( const Matrix& A )
{
	Vector a( A.row_size() + 1 );
	for( size_type i = 1; i != a.size(); i++ )
	{
		a[ i ] = A( i - 1, i - 1 );
	}
	return( a );
}

/**
	Initialize off-diagonal elements for use in QR method.	
	@return b Vector
 */
Vector
Lin_Alg::initialize_off_diagonal( const Matrix& A )
{
	Vector b( A.row_size() + 1 );
	for( size_type i = 2; i != b.size(); i++ )
	{
		b[ i ] = A( i - 2, i - 1 );
	}
	return( b );
}

/**
	Create matrix from vectors a and b where a is the new matrix's diagonal
	and b is the off diagonal on both sides. The resulting matrix is
	tridiagonal.
 */
Matrix 
Lin_Alg::to_matrix( Vector& a, Vector& b )
{
	if( a.size() != b.size() )
	{
		throw std::domain_error( "Dimensions incompatible in to_matrix" );
	}

	// a is one size larger than the matrix
	Matrix result( a.size() - 1 );

	// fill diagonal
	for( size_type i = 1; i != a.size(); i++ )
	{
		result( i - 1, i - 1) = a(i);
	}

	// fill off diagonal
	for( size_type j = 2; j != b.size(); j++ )
	{
		result( j - 2, j - 1) = b(j);
		result( j - 1, j - 2) = b(j);
	}

	return( result );
}

/**
	Compute eigenvalues of 2x2 submatrix in QR method.
	@param a Vector
	@param b Vector
	@param n
	@return eigenvalues Eigenvalues
 */
Vector
Lin_Alg::compute_2_by_2_eigenvalues( Vector& a, Vector& b, size_type n )
{
	Vector eigenvalues( 2 );

	/*
		use quadratic formula to solve det( A - lambda * I) = 0 for eigenvalue
		Note: underscore added to b since we already have a vector named b
	 */
	value_type _b  =  -( a[ n - 1 ] + a[ n ] );
	value_type c   =  a[ n ] * a[ n - 1 ] - ( b[ n ] * b[ n ] );
	value_type d   =  sqrt( _b * _b - 4 * c );

	// compute eigenvalues, mu_1 and mu_2
	value_type mu_1, mu_2;
	if( _b > 0 )
	{
		mu_1 = (-2 * c) / (_b + d);
		mu_2 = -(_b + d) / 2;
	}
	else
	{
		mu_1 = (d - _b) / 2;
		mu_2 = (2 * c) / (d - _b);
	}

	eigenvalues[ 0 ] = mu_1;
	eigenvalues[ 1 ] = mu_2;

	return( eigenvalues );
}
