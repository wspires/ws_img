/**
	@file   main.cpp
	@author Wade Spires
	@date   2005/10/25
	@brief  Test program for using matrix class.
 */

#include <iostream>
#include <string>
#include <vector>

#include <cstdio>
#include <ctime>

#include "Gauss.hpp"
#include "Image.hpp"
#include "Matrix.hpp"
#include "Linear_Algebra.hpp"
#include "Vector.hpp"
#include "err_mesg.h"
#include "ws_tools.hpp"

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

using namespace ws_img;

using std::cerr;  using std::cout;  using std::endl;
using std::string;
using std::vector;

void print_usage( );

void print_test( );
void sort_test( );
void iter_test( );
void vec_test( );
void const_iter_test( );

void comp_matrix( );

void conv_test( );
void cov_test( );
void mult_test( );

void svd_test( );
void pseudo_test( );
void svd_solve_test( );
void jacobi_test( );
void jacobi_test2( );

void read_test( );
void scale_test( );

void log_test( );

void toeplitz_test( );

void range_test( );

void special_matrix_test( );

int gsl_jacobi(gsl_matrix * A, gsl_matrix * Q, gsl_vector * S);

void kurtosis_test( );

void fft_test( );

typedef Matrix::size_type size_type;
typedef Matrix::value_type value_type;

int main( int argc, char** argv )
{
	if( argc != 1 )
	{
		print_usage();
	}

	// print_test();

	sort_test();

	// iter_test();
	// conv_test();
	// cov_test();
	// vec_test();
	// const_iter_test();
	// comp_matrix();
	// mult_test();

	// SVD function tests
	// svd_test();
	// pseudo_test();
	// svd_solve_test();
	// jacobi_test();
	// jacobi_test2( );

	// read_test();
	// scale_test();

	// log_test();

	// toeplitz_test();
	// range_test();
	// special_matrix_test();

	// kurtosis_test();
	// fft_test();

	return( EXIT_SUCCESS );
}

void print_usage( )
{
	cerr << "usage: test_matrix" << endl;
	exit( EXIT_FAILURE );
}

void print_test( )
{
	unsigned size = 5;
	Matrix m = Matrix::random( size, size, -5, 10 );
	m.write();

	for( unsigned i = 0; i != m.row_size(); ++i )
	{
		for( unsigned j = 0; j != m.col_size(); ++j )
		{
			m(i,j) = 10;
		}
	}
	m.write();

	// assign each element to 100
	m = 100;
	m.write();
}

void sort_test( )
{
	unsigned size = 5;
	Matrix m = Matrix::random( size, size, 1000 );
	// Vector m = Vector::random( size, 1000 );
	m.write();
	fprintf( stderr, "\n" );

	m.sort();
	m.write();
	fprintf( stderr, "\n" );

	m.reverse_sort();
	m.write();
}

/**
	Test iterators.
 */
void iter_test( )
{
	Matrix mat( "in.m" );
	mat.write();
	fprintf( stderr, "\n" );

	Matrix_Region region1( 0, 2, 2, 3 );
	Matrix_Region region2( 1, 1, 2, 2 );

	mat.set_region( region1 );
	Matrix mat2 = mat;

	mat2.set_region( region2 );
	mat.write();
	mat2.write();

	// mat2 += mat;
	mat2.reset_region();
	mat2.write();
	fprintf( stderr, "\n" );

	// zero-out region and view the complete matrix to show the hole
	mat.clear();
	mat.write();
	mat.reset_region();
	mat.write();
	fprintf( stderr, "\n" );


	// iterate across each column
	fprintf( stderr, "Column Iterator:\n" );
	for( Matrix::size_type i = mat.col_begin(); i != mat.col_end(); ++i )
	{
		fprintf( stderr, "Column %u / %u\n", i, mat.col_size() - 1 );
		Matrix::col_iterator col_iter = mat.col_begin( i );
		Matrix::col_iterator col_end  = mat.col_end( i );

		fprintf( stderr, "Iterator -: Size: %u Iter diff: %d\n",
				mat.sub_col_size(), int(col_end - col_iter) );

		while( col_iter != col_end )
		{
			fprintf( stderr, "  %lf\n", *col_iter );
			++col_iter;
		}
	}
	fprintf( stderr, "\n" );

	// iterate across each row
	fprintf( stderr, "Row Iterator:\n" );
	for( Matrix::size_type i = mat.row_begin(); i != mat.row_end(); ++i )
	{
		fprintf( stderr, "Row %u / %u\n", i, mat.col_size() - 1 );
		Matrix::row_iterator row_iter = mat.row_begin( i );
		Matrix::row_iterator row_end  = mat.row_end( i );
		while( row_iter != row_end )
		{
			fprintf( stderr, "  %lf\n", *row_iter );
			++row_iter;
		}
	}
	fprintf( stderr, "\n" );

	// iterate across every other row
	fprintf( stderr, "Row Iterator (offset = 2):\n" );
	for( Matrix::size_type i = mat.row_begin(); i != mat.row_end(); ++i )
	{
		fprintf( stderr, "Row %u / %u\n", i, mat.col_size() - 1 );
		Matrix::row_iterator row_iter = mat.row_begin( i );
		Matrix::row_iterator row_end  = mat.row_end( i );
		while( row_iter < row_end )
		{
			fprintf( stderr, "  %lf\n", *row_iter );
			row_iter += 2;
		}
	}
	fprintf( stderr, "\n" );
}

/*
class Vec
{
public:
	typedef int              value_type;
	typedef value_type&      reference;
	typedef const reference  const_reference;
	typedef value_type*      pointer;
	typedef const pointer    const_pointer;

	template <typename V, typename R, typename P>
	struct Nonconst_traits
	{
		typedef typename V  value_type;
		typedef typename R  reference;
		typedef typename P  pointer;

		typedef random_access_iterator_tag iterator_category;
	};

	template <typename V, typename R, typename P>
	struct Const_traits
	{
		typedef typename       V  value_type;
		typedef typename const R  reference;
		typedef typename const P  pointer;

		typedef random_access_iterator_tag iterator_category;
	};

	// typedef value_type* iterator;
	// typedef const value_type* const_iterator;

	template <typename _traits>
	class iterator
	{
	public:
		typedef typename _traits::value_type         value_type;
		typedef typename _traits::const_reference    reference;
		typedef typename _traits::const_pointer      pointer;

		typedef typename _traits::iterator_category  iterator_category;

	private:
		pointer _iter;

	public:

		iterator( )
		: _iter(0)
		{ }

		iterator( Vec::pointer iter )
		: _iter( iter )
		{ }
	};

private:
	pointer data;

public:
	Vec() {}

	typedef _iterator<
		Nonconst_traits< value_type, reference, pointer > >
	iterator;

	typedef _iterator<
		Const_traits< value_type, const_reference, const_pointer> >
	const_iterator;

	iterator begin( )
	{
		return( iterator(data) );
	}

	const_iterator begin( ) const
	{
		return( const_iterator(data) );
	}
};
 */

void const_iter_test( )
{
	// Vec vec;
	// Vec::iterator iter        = vec.begin();
	// Vec::const_iterator citer = vec.begin();

	//Vec::const_iterator iter;

	//Vec::iterator iter = vec.begin();
	//Vec::iterator iter;

	//Matrix::const_col_iterator const_col_iter
		//= mat.col_begin( mat.col_size() - 1 );
/*
	Matrix::const_col_iterator col_end  = mat.col_end( mat.col_size() - 1 );

	vector<double> v2;
	vector<double>::iterator v_iter = v2.begin();
	vector<double>::const_iterator vc_iter = v2.begin();
 */

}

/**
	Compare matrix versions for speed.
	Should have identical function in example/example.cpp
 */
void comp_matrix( )
{
	time_t start_time = time(NULL);

	unsigned size = 5000;
	//Matrix m = Matrix::random( size );
	Matrix m( size );
	Matrix x( size );

	// iterate through matrix
	// for( Matrix::iterator iter = m.begin(); iter != m.end(); ++iter )
	// {
		// *iter = 100;
	// }

	for( unsigned i = 0; i != m.row_size(); ++i )
	{
		for( unsigned j = 0; j != m.col_size(); ++j )
		{
			m(i,j) = 100;
			x(i,j) = 100;
		}
	}

	// for( unsigned i = m.row_begin(); i != m.row_end(); ++i )
	// {
		// for( unsigned j = m.col_begin(); j != m.col_end(); ++j )
		// {
			// m(i,j) = 100;
		// }
	// }

	// m *= m;
	//m.write();

	time_t end_time = time(NULL);
	fprintf( stderr, "Time: %u\n", unsigned(end_time - start_time) );
}

/**
	Test Vector class.
 */
void vec_test( )
{
	/*
	Vector v( "vec.m" );
	v.write();

	v.set_region( 0, 2 );
	Matrix M = v;
	M.write();
	// Vector v2 = M.to_vector();
	// Vector v3 = v.to_vector( M );

	fprintf( stderr, "%lf\n", v[2] );
	 */

	// test convolution
	Vector x( 5 );
	x(0) = -2;
	x(1) =  0;
	x(2) =  1;
	x(3) = -1;
	x(4) =  3;
	Vector h( 4 );
	h(0) =  1;
	h(1) =  2;
	h(2) =  0;
	h(3) = -3;

	// test vector
	// x.set_region(1, 2);
	// h.set_region(1, 3);
	x.write();
	h.write();
	Vector y = x.convolve( h );
	y.write();

	// test matrix
	Matrix X = x;
	X = X.transpose();
	X.set_region(0, 1, 1, 2);
	X.write();
	X.convolve_row(h);
	X.write();

	fprintf( stderr, "\n" );
	h.write();
	h /= 3;
	h.write();
}

/**
	Test convolution.
 */
void conv_test( )
{
	unsigned size = 10;
	Matrix m( size );
	for( unsigned i = m.row_begin(); i != m.row_end(); ++i )
	{
		for( unsigned j = m.col_begin(); j != m.col_end(); ++j )
		{
			m(i,j) = i + j;
		}
	}
	// Matrix m = Matrix::random( size );

	m.write();
	fprintf( stderr, "\n" );

	Vector::size_type f_size = 9;
	Vector h( f_size );
	h(0) = 1;
	h(1) = 0;
	h(2) = -1;
	h(3) = -1;
	h(4) = -1;
	h(5) = .25;
	h(6) = .125;
	h(7) = -.25;
	h(8) = .5;
	// h.write();

	Vector v( size );
	v(0) = 2;
	v(1) = 3;
	v(2) = 4;
	v(3) = 5;
	v(4) = 6;
	v(5) = 7;
	v(6) = 8;
	v(7) = 9;
	v(8) = 10;
	// v(3) = 1;
	// v(4) = 1;

	// Vector y = v.convolve( h );
	// y.write();

	// m = v;
	// m = m.transpose();
	// m.write();
	// m.convolve_row( v );
	m.convolve_row( h );
	// m.convolve( h, h );
	m.write();

	/*
	Matrix f( f_size );
	f(0,0) =  1;
	f(0,1) =  0;
	f(0,2) = -1;
	f(1,0) =  1;
	f(1,1) =  0;
	f(1,2) = -1;
	f(2,0) =  1;
	f(2,1) =  0;
	f(2,2) = -1;
	f.write();

	// m.convolve( h, v );
	m.convolve( f );

	m.write();
	 */
}

/**
	Test covariance.
 */
void cov_test( )
{
	// unsigned size = 5;
	// Matrix m( size, size - 1 );

	Matrix m( 6, 3 );

	unsigned k = 0;
	for( unsigned j = m.col_begin(); j != m.col_end(); ++j )
	{
		for( unsigned i = m.row_begin(); i != m.row_end(); ++i )
		{
			m(i,j) = ++k;
		}
	}
	m(0,0) = .7271;
	m(0,1) = .5466;
	m(0,2) = .5226;
	m(1,0) = .3093;
	m(1,1) = .4449;
	m(1,2) = .8801;
	m(2,0) = .8385;
	m(2,1) = .6946;
	m(2,2) = .1730;
	m(3,0) = .5681;
	m(3,1) = .6213;
	m(3,2) = .9797;
	m(4,0) = .3710;
	m(4,1) = .7904;
	m(4,2) = .2714;
	m(5,0) = .7027;
	m(5,1) = .9568;
	m(5,2) = .2523;
	m.write( "", 4 );

	Matrix cov;
	Vector mean;
	// m.covariance( cov, mean );
	// cov.write("", 4);
	// mean.write();
	// m.write();

	Vector v = Vector::random( 3, 5 );
	v.write( "", 4 );

	// m.prin_comp( 3 );
	Matrix A;
	Lin_Alg::prin_comp_matrix( 3, m, A, mean );
	Lin_Alg::prin_comp_project( A, mean, v );
	v.write( "", 4 );
}

void mult_test( )
{
	unsigned size = 5;
	Matrix m( size );

	for( unsigned i = m.row_begin(); i != m.row_end(); ++i )
	{
		for( unsigned j = m.col_begin(); j != m.col_end(); ++j )
		{
			m(i,j) = i + i * j;
		}
	}
	// Matrix m = Matrix::random( size );

	// double avg, var, std;
	// m.stat( &avg, &var, &std );
	// fprintf( stderr, "%lf %lf %lf\n", avg, var, std );

	m.write();
	// m = m.transpose();
	// m.write();
	// m = m * m;
	// m.write();

	Vector v( size );
	for( unsigned i = v.vec_begin(); i != v.vec_end(); ++i )
	{
		v(i) = i;
	}
	Matrix x( v );

	v.write();
	x.write();

	Vector mv = m * v;
	Matrix mx = m * x;

	mv.write();
	mx.write();
	double norm = mv.two_norm_2();
	fprintf( stderr, "%lf\n", norm );
}

void svd_test( )
{
	unsigned size = 3;
	Matrix A( size );
	for( size_type i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( size_type j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = i + i * j;
			A(i,j) = i * i + j * j;
		}
	}
	A.write();
	fprintf( stderr, "\n" );

	// perform SVD and output results
	Matrix U = A;
	Matrix V;
	Vector S;
	// Lin_Alg::SVD( U, S, V, Lin_Alg::Golub_Reinsch );
	Lin_Alg::SVD( U, S, V, Lin_Alg::Jacobi );
	U.write();
	S.write();
	V.write();
	fprintf( stderr, "\n" );

	// calculate U * S * V^T
	// U * S where S diagonal matrix (stored here as a vector only)
	for( size_type i = U.row_begin(); i != U.row_end(); ++i )
	{
		for( size_type j = U.col_begin(); j != U.col_end(); ++j )
		{
			U(i,j) *= S(j);
		}
	}
	Matrix A_new = U * V.transpose();  // finish computation

	A_new.write();  // write the result

	// compare reconstructed matrix with the original to verify they are equal
	Matrix diff = A - A_new;
	for( size_type i = diff.row_begin(); i != diff.row_end(); ++i )
	{
		for( size_type j = diff.col_begin(); j != diff.col_end(); ++j )
		{
			if( fabs( diff(i,j) ) > .0001 )
			{
				fprintf( stderr, "Not equal at position (%u, %u): "
						"%lf (original) versus %lf (reconstructed)\n",
						i, j, A(i,j), A_new(i,j) );
			}
		}
	}
}

/**
	Test pseudo-inverse.
 */
void pseudo_test( )
{
	unsigned size = 3;
	Matrix A( size );

	for( unsigned i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( unsigned j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = i + i * j;
			A(i,j) = i * i + j * j;
		}
	}
	A.write();

	Matrix P_inv = Lin_Alg::pseudo_inverse( A );
	P_inv.write();

	Matrix A_new = A * P_inv * A;
	A_new.write();

	// compare reconstructed matrix with the original to verify they are equal
	Matrix diff = A - A_new;
	for( size_type i = diff.row_begin(); i != diff.row_end(); ++i )
	{
		for( size_type j = diff.col_begin(); j != diff.col_end(); ++j )
		{
			if( fabs( diff(i,j) ) > .0001 )
			{
				fprintf( stderr, "Not equal at position (%u, %u): "
						"%lf (original) versus %lf (reconstructed)\n",
						i, j, A(i,j), A_new(i,j) );
			}
		}
	}
}

/**
	Solve Ax = b using SVD.
 */
void svd_solve_test( )
{
	unsigned size = 3;

	// create matrix A
	Matrix A( size );
	for( unsigned i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( unsigned j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = i * i + j * j;  // singular (no solution)
			A(i,j) = i + i * j;      // non-singular
		}
	}
	A.write();
	fprintf( stderr, "\n" );

	// create vector b
	Vector b( size );
	for( unsigned i = b.vec_begin(); i != b.vec_end(); ++i )
	{
		b(i) = i;
	}
	b.write();
	fprintf( stderr, "\n" );

	// solve for x
	Vector x;
	Lin_Alg::solve_system_SVD( A, b, x );
	x.write();

	// compare solutions
	Vector solution = A * x;
	solution.write();

	fprintf( stderr, "\n" );

	// compare reconstructed matrix with the original to verify they are equal
	fprintf( stderr, "Comparing solution (may not be exact if A is singular:\n" );
	Vector diff = b - solution;
	for( size_type i = diff.vec_begin(); i != diff.vec_end(); ++i )
	{
		if( fabs( diff(i) ) > .0001 )
		{
			fprintf( stderr, "Not equal at position %u: "
					"%lf (original) versus %lf (reconstructed):\n %lf (difference)\n",
					i, b(i), solution(i), diff(i) );
		}
	}
}

/**
	Test re-implementation of SVD Jacobi algorithm.
 */
void jacobi_test( )
{
	unsigned size = 5;
	Matrix A( size );
	for( size_type i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( size_type j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = i + i * j;
			A(i,j) = i * i + j * j;
		}
	}
	// A.set_region( 1, 1, 3, 3 );
	// A.write();
	fprintf( stderr, "\n" );

	Matrix U = A;
	Matrix V;
	Vector S;

time_t start, end;

	fprintf( stderr, "---- GSL ----\n" );
start = time(NULL);
	Lin_Alg::SVD( U, S, V, Lin_Alg::Jacobi );
	// Lin_Alg::SVD( U, S, V, Lin_Alg::Golub_Reinsch_mod );
end = time(NULL);
	U.write();
	S.write();
	V.write();
	fprintf( stderr, "%u seconds\n", unsigned(end - start) );
	fprintf( stderr, "\n" );

	// perform SVD and output results
	fprintf( stderr, "---- Reimplementation ----\n" );
start = time(NULL);
	// unsigned num_nonzero = Lin_Alg::SVD_jacobi( A, S, V );
end = time(NULL);
	A.write();
	S.write();
	V.write();
	fprintf( stderr, "%u seconds\n", unsigned(end - start) );
	fprintf( stderr, "\n" );

	// calculate U * S * V^T
	// U * S where S diagonal matrix (stored here as a vector only)
	for( size_type i = U.row_begin(); i != U.row_end(); ++i )
	{
		for( size_type j = U.col_begin(); j != U.col_end(); ++j )
		{
			A(i,j) *= S(j);
		}
	}
	Matrix A_new = A * V.transpose();  // finish computation

	// A_new.write();  // write the result
}

/**
	Test copy of GSL SVD Jacobi algorithm.
 */
void jacobi_test2( )
{
	fprintf( stderr, "---- Local SVD Jacobi ----\n" );

	unsigned size = 5;
	Matrix A( size );
	for( size_type i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( size_type j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = i + i * j;
			A(i,j) = i * i + j * j;
		}
	}

	// number of rows and columns in a subregion of A
	size_type M = A.sub_row_size();
	size_type N = A.sub_col_size();

	double* data = A.to_array();
	gsl_matrix_view m = gsl_matrix_view_array( data, M, N );
	gsl_matrix* U_gsl = &m.matrix;
	gsl_vector* S_gsl = gsl_vector_alloc( N );
	gsl_matrix* V_gsl = gsl_matrix_alloc( N, N );

	// use local SVD Jacobi
	gsl_jacobi( U_gsl, V_gsl, S_gsl );

	Matrix U, V;
	Vector S;

	// allocate output matrices
	U = Matrix( M, N );
	S = Vector( N );
	V = Matrix( N );

	// store GSL matrices and vector in my matrix/vector objects
	for( size_type i = U.row_begin(); i != U.row_end(); ++i )
	{
		for( size_type j = U.col_begin(); j != U.col_end(); ++j )
		{
			U(i,j) = gsl_matrix_get( U_gsl, i, j );
		}
	}
	for( size_type i = V.row_begin(); i != V.row_end(); ++i )
	{
		for( size_type j = V.col_begin(); j != V.col_end(); ++j )
		{
			V(i,j) = gsl_matrix_get( V_gsl, i, j );
		}
		S(i) = gsl_vector_get( S_gsl, i );
	}

	// deallocate U, S, and V
	delete[] data;
	gsl_vector_free( S_gsl );
	gsl_matrix_free( V_gsl );

	// output results
	U.write();
	S.write();
	V.write();
	fprintf( stderr, "\n" );
}

/* This is a the jacobi version */
/* Author:  G. Jungman */

/*
 * Algorithm due to J.C. Nash, Compact Numerical Methods for
 * Computers (New York: Wiley and Sons, 1979), chapter 3.
 * See also Algorithm 4.1 in
 * James Demmel, Kresimir Veselic, "Jacobi's Method is more
 * accurate than QR", Lapack Working Note 15 (LAWN15), October 1989.
 * Available from netlib.
 *
 * Based on code by Arthur Kosowsky, Rutgers University
 *                  kosowsky@physics.rutgers.edu  
 *
 * Another relevant paper is, P.P.M. De Rijk, "A One-Sided Jacobi
 * Algorithm for computing the singular value decomposition on a
 * vector computer", SIAM Journal of Scientific and Statistical
 * Computing, Vol 10, No 2, pp 359-371, March 1989.
 * 
 */

int
gsl_jacobi(gsl_matrix * A, gsl_matrix * Q, gsl_vector * S)
{
  if (A->size1 < A->size2)
    {
      /* FIXME: only implemented  M>=N case so far */

      GSL_ERROR ("svd of MxN matrix, M<N, is not implemented", GSL_EUNIMPL);
    }
  else if (Q->size1 != A->size2)
    {
      GSL_ERROR ("square matrix Q must match second dimension of matrix A",
                 GSL_EBADLEN);
    }
  else if (Q->size1 != Q->size2)
    {
      GSL_ERROR ("matrix Q must be square", GSL_ENOTSQR);
    }
  else if (S->size != A->size2)
    {
      GSL_ERROR ("length of vector S must match second dimension of matrix A",
                 GSL_EBADLEN);
    }
  else
    {
      const size_t M = A->size1;
      const size_t N = A->size2;
      size_t i, j, k;

      /* Initialize the rotation counter and the sweep counter. */
      int count = 1;
      int sweep = 0;
      int sweepmax = 5*N;

      double tolerance = 10 * M * GSL_DBL_EPSILON;

      /* Always do at least 12 sweeps. */
      sweepmax = GSL_MAX (sweepmax, 12);

      /* Set Q to the identity matrix. */
      gsl_matrix_set_identity (Q);

      /* Store the column error estimates in S, for use during the
         orthogonalization */

      for (j = 0; j < N; j++)
        {
          gsl_vector_view cj = gsl_matrix_column (A, j);
          double sj = gsl_blas_dnrm2 (&cj.vector);
          gsl_vector_set(S, j, GSL_DBL_EPSILON * sj);
        }
    
      /* Orthogonalize A by plane rotations. */

		int iter = 0;
      while (count > 0 && sweep <= sweepmax)
		{
// fprintf( stderr, "Count: %d\n", count );

			/* Initialize rotation counter. */
			count = N * (N - 1) / 2;

			for (j = 0; j < N - 1; j++)
			{
				for (k = j + 1; k < N; k++)
				{
		// fprintf( stderr, "%u %u: ", j, k );
					double a = 0.0;
					double b = 0.0;
					double p = 0.0;
					double q = 0.0;
					double cosine, sine;
					double v;
					double abserr_a, abserr_b;
					int sorted, orthog, noisya, noisyb;

					gsl_vector_view cj = gsl_matrix_column (A, j);
					gsl_vector_view ck = gsl_matrix_column (A, k);

					gsl_blas_ddot (&cj.vector, &ck.vector, &p);
					p *= 2.0 ;  /* equation 9a:  p = 2 x.y */

					a = gsl_blas_dnrm2 (&cj.vector);
					b = gsl_blas_dnrm2 (&ck.vector);

					q = a * a - b * b;
					v = hypot(p, q);


			// fprintf( stderr, "(p,a,b,q,v): %lf %lf %lf %lf %lf\n", p, a, b, q, v );
					/* test for columns j,k orthogonal, or dominant errors */

					abserr_a = gsl_vector_get(S,j);
					abserr_b = gsl_vector_get(S,k);

			// fprintf( stderr, "\tabs_err: %.16lf %.16lf\n", abserr_a, abserr_b );

					// sorted = (gsl_coerce_double(a) >= gsl_coerce_double(b));
					// orthog = (fabs (p) <= tolerance * gsl_coerce_double(a * b));
					sorted = ((a) >= (b));
					orthog = (fabs (p) <= tolerance * (a * b));
					noisya = (a < abserr_a);
					noisyb = (b < abserr_b);

		if( iter == 35 || iter == 36 )
		{

			// fprintf( stderr, "Iteration: %u\n", iter );
			// fprintf( stderr, "%.16lf %.16lf %lf %lf %lf\n", a, b, p, q, v );

				// fprintf( stderr, "cs: %lf %lf\n", cosine, sine );
			// for( size_type i = 0; i != M; ++i )
			// {
				// fprintf( stderr, "%.16lf ", gsl_vector_get(S,i) );
			// }
// fprintf( stderr, "\t(bool): %d %d %d %d\n", sorted, orthog, noisya, noisyb );
			// fprintf( stderr, "------------\n" );
		}
		// fprintf( stderr, "%u\n", iter );

					if (sorted && (orthog || noisya || noisyb))
					{
						count--;
						continue;
					}

					/* calculate rotation angles */
					if (v == 0 || !sorted)
					{
						cosine = 0.0;
						sine = 1.0;
	// fprintf( stderr, "\t unsorted: (cos,sin): %lf %lf\n", cosine, sine );
					}
					else
					{
						cosine = sqrt((v + q) / (2.0 * v));
						sine   = p / (2.0 * v * cosine);
	// fprintf( stderr, "\t sorted: (cos,sin): %lf %lf\n", cosine, sine );
					}


// fprintf( stderr, "\n" );

					/* apply rotation to A */
				// fprintf( stderr, "\tA:\n" );
					for (i = 0; i < M; i++)
					{
					// fprintf( stderr, "\t\t%lf %lf -> ",
						// gsl_matrix_get (A, i, j), gsl_matrix_get (A, i, k) );

						const double Aik = gsl_matrix_get (A, i, k);
						const double Aij = gsl_matrix_get (A, i, j);
						gsl_matrix_set (A, i, j, Aij * cosine + Aik * sine);
						gsl_matrix_set (A, i, k, -Aij * sine + Aik * cosine);

					// fprintf( stderr, "%lf %lf\n",
						// gsl_matrix_get (A, i, j), gsl_matrix_get (A, i, k) );
					}

					gsl_vector_set(S, j,
							fabs(cosine) * abserr_a + fabs(sine) * abserr_b);
					gsl_vector_set(S, k,
							fabs(sine) * abserr_a + fabs(cosine) * abserr_b);


					/* apply rotation to Q */
				// fprintf( stderr, "\tV:\n" );
					for (i = 0; i < N; i++)
					{
						const double Qij = gsl_matrix_get (Q, i, j);
						const double Qik = gsl_matrix_get (Q, i, k);

					// fprintf( stderr, "\t\t%lf %lf -> ",
						// gsl_matrix_get (Q, i, j), gsl_matrix_get (Q, i, k) );

						gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
						gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);

					// fprintf( stderr, "%lf %lf\n",
						// gsl_matrix_get (Q, i, j), gsl_matrix_get (Q, i, k) );
					}
					// if( sweep == 4 )
					if( iter == 35 || iter == 36 )
					{
		// for( size_type i = 0; i != M; ++i )
		// {
			// for( size_type j = 0; j != N; ++j )
			// {
				// fprintf( stderr, " %lf", gsl_matrix_get( Q, i, j ) );
			// }
			// fprintf( stderr, "\n" );
		// }
		// fprintf( stderr, "-----------\n" );
					}
			++iter;
				}
			}

		// fprintf( stderr, "Sweep %d\n", sweep );


			/* Sweep completed. */
			sweep++;
		}
		/*
		for( size_type i = 0; i != M; ++i )
		{
			for( size_type j = 0; j != N; ++j )
			{
				fprintf( stderr, " %lf", gsl_matrix_get( A, i, j ) );
			}
			fprintf( stderr, "\n" );
		}
		fprintf( stderr, "\n" );
		for( size_type i = 0; i != M; ++i )
		{
			for( size_type j = 0; j != N; ++j )
			{
				fprintf( stderr, " %lf", gsl_matrix_get( Q, i, j ) );
			}
			fprintf( stderr, "\n" );
		}
		return(0);
		 */

      /* 
       * Orthogonalization complete. Compute singular values.
       */

      {
        double prev_norm = -1.0;

        for (j = 0; j < N; j++)
          {
            gsl_vector_view column = gsl_matrix_column (A, j);
            double norm = gsl_blas_dnrm2 (&column.vector);

            /* Determine if singular value is zero, according to the
               criteria used in the main loop above (i.e. comparison
               with norm of previous column). */

            if (norm == 0.0 || prev_norm == 0.0 
                || (j > 0 && norm <= tolerance * prev_norm))
              {
                gsl_vector_set (S, j, 0.0);     /* singular */
                gsl_vector_set_zero (&column.vector);   /* annihilate column */

                prev_norm = 0.0;
              }
            else
              {
                gsl_vector_set (S, j, norm);    /* non-singular */
                gsl_vector_scale (&column.vector, 1.0 / norm);  /* normalize column */

                prev_norm = norm;
              }
          }
      }

      if (count > 0)
        {
          /* reached sweep limit */
          GSL_ERROR ("Jacobi iterations did not reach desired tolerance",
                     GSL_ETOL);
        }

      return GSL_SUCCESS;
    }
}

/**
	Test reading matrix.
	U_T is very large, so our original read method did not work correctly.
 */
void read_test( )
{
	Matrix A( "U_T.m" );
	A.write( "A.m" );
	Matrix B( "in.m" );
	B.write();
}

/**
	Test scaling function.
 */
void scale_test( )
{
	unsigned size = 5;
	Matrix m( size, size );

	for( unsigned i = m.row_begin(); i != m.row_end(); ++i )
	{
		for( unsigned j = m.col_begin(); j != m.col_end(); ++j )
		{
			m(i,j) = i * j;
		}
	}
	m.write();

	// m.scale(2.75);
	m.scale(.5);
	// m.scale(2);
	m.write();
}

/**
	Test Laplacian of Gaussian.
 */
void log_test( )
{
	fprintf( stderr, "LoG Test\n" );
	Image img( "Germany2.pgm" );

	Vector f, g;
	make_log_filter( f, g, 1.5 );
	fprintf( stderr, "Gaussian 2nd derivative:\n" );
	f.write();
	fprintf( stderr, "Gaussian:\n" );
	g.write();

	// equivalent to filter_log( img );
	// img.convolve( f, g ) += Image( img ).convolve( g, f );

	filter_log( img );
	img.write( "log.pgm" );

	// fprintf( stderr, "LoG: %lf\n", laplace_gauss(1,2,1.3) );
}

void toeplitz_test( )
{
	fprintf( stderr, "Toeplitz Test\n" );

	// column
	Vector c( 5 );
	for( unsigned i = 0; i != c.sub_size(); ++i )
	{
		c(i) = i;  // 0, 1, ...
	}

	// row
	Vector r( 3 );
	for( unsigned i = 0; i != r.sub_size(); ++i )
	{
		r(i) = i + 1;  // l, 2, ...
	}

	Matrix A = Matrix::toeplitz( c, r );
	A.write();
}

/**
	A simple struct used in the next function.
 */
struct Range
{
	double begin;  //< First value in the range
	double inc;    //< Amount to increment each successive value
	double end;    //< Last value in the range

	Range( const double b, const double i, const double e )
	: begin(b), inc(i), end(e)
	{ }
};

void range_test( )
{
	fprintf( stdout, "Colon Test\n" );
	fprintf( stdout, "----------\n" );

	// add each test
	vector<Range> ranges;
	ranges.push_back( Range( 1, 1, 10 ) );
	ranges.push_back( Range( 0, 2, 10 ) );
	ranges.push_back( Range( 0, 11, 10 ) );
	ranges.push_back( Range( 1, -.1, .1 ) );
	ranges.push_back( Range( 1, -.2, .1 ) );
	ranges.push_back( Range( 1, .2, 2.1 ) );
	ranges.push_back( Range( 1, 0, 1 ) );
	ranges.push_back( Range( 2, -2, -4 ) );
	ranges.push_back( Range( -1, 1.5, 3 ) );
	ranges.push_back( Range( -2.5, -.5, -4.7 ) );
	ranges.push_back( Range( 1, -.05, .55 ) );
	ranges.push_back( Range( 1, -.05, .5499999999 ) );

	// perform each test
	for( unsigned i = 0; i != ranges.size(); ++i )
	{
		fprintf( stdout, "%lf: %lf: %lf\n",
				ranges[i].begin, ranges[i].inc, ranges[i].end );
		Vector v( ranges[i].begin, ranges[i].inc, ranges[i].end );
		v.write( "", 3 );
		fprintf( stdout, "\n" );
	}
}

void special_matrix_test( )
{
	fprintf( stderr, "Special Matrix Test\n" );
	fprintf( stderr, "-------------------\n" );

	// add each test
	vector<Matrix> matrices;
	matrices.push_back( Matrix::ones( 5, 6 ) );
	matrices.push_back( Matrix::ones( 6 ) );
	matrices.push_back( Matrix::zeros( 5, 6 ) );
	matrices.push_back( Matrix::zeros( 6 ) );
	matrices.push_back( Matrix::identity( 4 ) );
	matrices.push_back( Matrix::identity( 4 ) );
	matrices.push_back( Matrix::random( 4, 4, 100 ) );
	matrices.push_back( Matrix::random( 4, 4, 5 ) );
	matrices.push_back( Matrix::random_normal( 4, 4 ) );
	matrices.push_back( Matrix::random_normal( 3, 4 ) );

	// perform each test
	for( unsigned i = 0; i != matrices.size(); ++i )
	{
		Matrix& A = matrices[i];
		A.write( "", 3 );
		fprintf( stdout, "\n" );
	}
}

void kurtosis_test( )
{
	fprintf( stdout, "Kurtosis Test\n" );
	fprintf( stdout, "-------------\n" );

	Matrix A(5,4);
	A(0,0) =  1.1650;
	A(0,1) =  1.6961;
	A(0,2) = -1.4462;
	A(0,3) = -0.3600;

	A(1,0) =  0.6268;
	A(1,1) =  0.0591;
	A(1,2) = -0.7012;
	A(1,3) = -0.1356;

	A(2,0) =  0.0751;
	A(2,1) =  1.7971;
	A(2,2) =  1.2460;
	A(2,3) = -1.3493;

	A(3,0) =  0.3516;
	A(3,1) =  0.2641;
	A(3,2) = -0.6390;
	A(3,3) = -1.2704;

	A(4,0) = -0.6965;
	A(4,1) =  0.8717;
	A(4,2) =  0.5774;
	A(4,3) =  0.9846;

	A.write();

	Vector kurtosis = A.kurtosis();
	kurtosis.write();

	fprintf( stdout, "\nNormally distributed data\n" );

	// normally distributed data
	Matrix B = Matrix::random_normal( 5, 5 );
	B.write();
	kurtosis = B.kurtosis();
	kurtosis.write();
}

void fft_test( )
{
	fprintf( stdout, "FFT Test\n" );
	fprintf( stdout, "-------------\n" );

	Matrix A(5,4);
	A(0,0) =  1.1650;
	A(0,1) =  1.6961;
	A(0,2) = -1.4462;
	A(0,3) = -0.3600;

	A(1,0) =  0.6268;
	A(1,1) =  0.0591;
	A(1,2) = -0.7012;
	A(1,3) = -0.1356;

	A(2,0) =  0.0751;
	A(2,1) =  1.7971;
	A(2,2) =  1.2460;
	A(2,3) = -1.3493;

	A(3,0) =  0.3516;
	A(3,1) =  0.2641;
	A(3,2) = -0.6390;
	A(3,3) = -1.2704;

	A(4,0) = -0.6965;
	A(4,1) =  0.8717;
	A(4,2) =  0.5774;
	A(4,3) =  0.9846;
	A.write();

	// A.write( "A.m" );

	// apply forward FFT
	Matrix real, imag;
	A.fft( real, imag );
	// real.write();
	// imag.write();

	// display real and imaginary components
	fprintf( stdout, "\n" );
	for( size_type i = real.row_begin(); i != real.row_end(); ++i )
	{
		for( size_type j = real.col_begin(); j != real.col_end(); ++j )
		{
			fprintf( stdout, "%lf %s %lf i\n",
					real(i,j),
					imag(i,j) > 0 ? "+" : "-",
					fabs( imag(i,j) )
				);
		}
	}
	fprintf( stdout, "\n" );

	// apply inverse FFT
	Matrix B;
	B.ifft( real, imag );
	B.write();

	// verify the inverse yielded the original matrix
	Matrix C = A - B;
	const value_type tolerance = .0000001;
	for( size_type i = C.row_begin(); i != C.row_end(); ++i )
	{
		for( size_type j = C.col_begin(); j != C.col_end(); ++j )
		{
			if( C(i,j) > tolerance )
			{
				err_quit( "(%u, %u) = %.50lf\n", C(i,j) );
			}
		}
	}
}
