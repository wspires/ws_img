/**
	@file   Linear_Algebra.hpp
	@author Wade Spires
	@date   2005/08/03
	@brief  Numerical algorithms for matrices.
	
	Solve linear system Ax=b using LU Factorization by Gaussian
		elimination with partial pivoting.

	Solve linear system Ax=b using LU Factorization by Gaussian
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

#ifndef LINEAR_ALG_H
#define LINEAR_ALG_H

#include  <algorithm>
#include  <cfloat>
#include  <cmath>
#include  <functional>
#include  <stdexcept>
#include  <iomanip>
#include  <iostream>
#include  <iterator>
#include  <vector>

#ifdef HAVE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>
#endif // HAVE_GSL

#include "ws_tools.hpp"

#include "Matrix.hpp"
#include "Vector.hpp"

namespace ws_img
{

/**
	@brief Class of static methods for solving linear systems.
 */
class Linear_Algebra
{
	typedef std::vector<unsigned>  Pivot;

   // default tolerance and number of iterations for iterative methods
	static const double    DEFAULT_TOL;
	static const unsigned  DEFAULT_ITER;
	static const double    EPS;

public:

	typedef Matrix::value_type value_type;
	typedef Matrix::size_type  size_type;

	enum SVD_Method { Golub_Reinsch, Golub_Reinsch_mod, Jacobi };

	Linear_Algebra( )
	{ }

	/*
		Static member functions
	 */

	static value_type det( const Matrix& );

	// solving matrices
	static Matrix  inverse_2_by_2( const Matrix& );
	static Matrix  inverse_3_by_3( const Matrix& );
	static void    solve_system( Matrix&, Vector&, double = DEFAULT_TOL );
	static void    solve_system( Matrix&, const Vector&, Vector& );
	static Vector  solve_system( const Matrix&, const Vector& );
	static Vector  jacobi( Matrix&, Vector&, unsigned, double, int* );
	static Vector  gauss_seidel( Matrix&, Vector&, unsigned, double, int* );
	static Vector  SOR( Matrix&, Vector&, double, unsigned, double, int* );

	// singular value decomposition
	static void     SVD( Matrix&, Vector&, Matrix&,
			const SVD_Method& = Golub_Reinsch );
	static Matrix   pseudo_inverse( const Matrix&,
			const SVD_Method& = Golub_Reinsch );
	static void     solve_system_SVD( const Matrix&, const Vector&, Vector&,
			const SVD_Method& = Golub_Reinsch );
	static unsigned SVD_jacobi( Matrix&, Vector&, Matrix&, value_type = EPS );

	// eigensystems
	static value_type power_method( Matrix&, Matrix&,
		value_type = DEFAULT_TOL, unsigned = DEFAULT_ITER );
	// TODO: does not work yet
	static value_type inverse_power_method( Matrix&, Vector&,
		value_type = DEFAULT_TOL, unsigned = DEFAULT_ITER );
	static value_type inf_norm( const Matrix&, size_type& );
	static value_type symmetric_power_method( Matrix&, Matrix&,
		value_type = DEFAULT_TOL, unsigned = DEFAULT_ITER );

	static void prin_comp_matrix( unsigned, Matrix&, Matrix&, Vector& );
	static void covariance( Matrix&, Matrix&, Vector& );
	static void prin_comp_project( const Matrix&, const Vector&, Vector& );
	static void prin_comp_project( const Matrix&, const Vector&, Matrix& );
	static void prin_comp( unsigned, Matrix& );

	static void eigen( const Matrix&, Vector&, Matrix& );
	static Vector eig_values2( const Matrix& );
	static void householder( Matrix& );
	static Vector QR( const Matrix&,
		value_type = DEFAULT_TOL, unsigned = DEFAULT_ITER );

protected:

	// internal methods used by solve_system()
	static void LU_factorization( Matrix&, Pivot&, double*,
		double = 1.0E-20 );
	static void LU_back_sub( const Matrix&, const Pivot&, Vector&,
		double = DEFAULT_TOL );

	static void LU_factorization( Matrix&, Pivot& );
	static void LU_forw_sub( const Matrix&, const Vector&, Vector&,
		const Pivot& );
	static int  LU_back_sub( const Matrix&, const Vector&, Vector&,
		const Pivot& );

	// for debugging issues involving pivot
	static void print_pivot( const Pivot& );

	static Vector initialize_diagonal( const Matrix& );
	static Vector initialize_off_diagonal( const Matrix& );
	static Matrix to_matrix( Vector&, Vector& );
	static Vector compute_2_by_2_eigenvalues( Vector&, Vector&, size_type );
};

typedef Linear_Algebra Lin_Alg;

} // namespace ws_img

#endif // LINEAR_ALG_H
