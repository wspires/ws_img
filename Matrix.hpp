/**
	@file   Matrix.hpp
	@author Wade Spires
	@date   2005/10/28
	@brief  Definition of matrix type.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef MATRIX_HPP
#define MATRIX_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include  <algorithm>
#include  <functional>
#include  <stdexcept>
#include  <iomanip>
#include  <iostream>
#include  <iterator>
#include  <string>
#include  <vector>

// c headers
#include  <cfloat>
#include  <climits>
#include  <cmath>

// This determines how the matrix internally stores its data.
// If HAVE_BOOST_SCOPED_PTR_WS_IMG is defined, then a boost::scoped_array is
// used. If USE_STL is defined (but HAVE_BOOST_SCOPED_PTR_WS_IMG is not), then a
// std::vector is used. Otherwise, a raw pointer is used.
// #define HAVE_BOOST_SCOPED_PTR_WS_IMG
#ifdef HAVE_BOOST_SCOPED_PTR_WS_IMG
	#include  <boost/scoped_array.hpp>
#else
	// whether to use the STL vector class (use raw pointer otherwise)
	// #define USE_STL
#endif // HAVE_BOOST_SCOPED_PTR_WS_IMG

// local headers
#include "Matrix_Iterator.hpp"
#include "Matrix_Region.hpp"
#include "ws_tools.hpp"

#ifndef PI
#define PI M_PI
#endif

namespace ws_img
{

class Vector;  // forward declaration

/**
	@brief Class representing matrix of numeric data.
 */
class Matrix
{

public:  // define basic types

	typedef double             value_type;

	// define the internal array type--either use a boost smart pointer, a STL
	// vector, or a raw pointer
#ifdef HAVE_BOOST_SCOPED_PTR_WS_IMG
	typedef boost::scoped_array<value_type> matrix_type;
#else
	#ifdef USE_STL
	typedef std::vector<value_type>         matrix_type;
	#else
	typedef value_type*                     matrix_type;
	#endif // USE_STL
#endif // HAVE_BOOST_SCOPED_PTR_WS_IMG

	typedef unsigned           size_type;
	typedef value_type&        reference;
	typedef const value_type&  const_reference;
	typedef value_type*        pointer;
	typedef const value_type*  const_pointer;
	typedef ptrdiff_t          difference_type;

	typedef std::vector<Matrix_Region>    Region_List;

	// iterator types are defined below

	typedef int prec_type;  //< Precision type for setting double's precision

	static const value_type MAX_VALUE;   //< Maximum value allowed in matrix

protected:  // types used by multiply()'s private helper functions

	typedef std::vector< int >                 Vector_Int;
	typedef std::vector< std::vector< int > >  Matrix_Int;

public:  // constructors and access functions

	/**
		Construct matrix by reading given matrix file.
		@param[in] matrix_file_name Name of matrix file to read matrix from
	 */
	Matrix( const std::string& );

	/**
		Construct square M-by-M matrix with each element left uninitialized.

		Note that the elements are not initialized to any value, not even 0. Even
		if by inspection this seems to be the case, do not rely upon the values
		being set.

		@param[in] M Number of rows and columns in matrix
	 */
	explicit Matrix( size_type = 1 );

	/**
		Construct M-by-N matrix with each element left uninitialized.

		Note that the elements are not initialized to any value, not even 0. Even
		if by inspection this seems to be the case, do not rely upon the values
		being set.

		@param[in] M Number of rows in matrix
		@param[in] N Number of columns in matrix
	 */
	explicit Matrix( size_type, size_type );

	/**
		Construct M-by-N matrix with each element initialized to the given value.
		@param[in] M Number of rows in matrix
		@param[in] N Number of columns in matrix
		@param[in] default_value Initial value to set each element of matrix
	 */
	explicit Matrix( size_type, size_type, value_type );

	/**
		Copy constructor copies matrix region.

		TODO Check whether better to copy entire matrix and just set the region
		rather than only copy the region

		@param[in] m Matrix to copy
	 */
	Matrix( const Matrix& );

	/**
		Assignment operator copies matrix region.
		@param[in] m Matrix to copy
	 */
	Matrix& operator=( const Matrix& m );

	Matrix& operator=( const value_type );

	Matrix( const Vector& ) ;
	Matrix& operator=( const Vector& ); 
	Matrix( const Vector&, size_type ); 
	Matrix( const Vector&, size_type, size_type ); 

	/**
		Construct M by N matrix using given STL array for values.
		@param[in] matrix Matrix of data represented as C++ vector
	 */
	Matrix( std::vector< std::vector< value_type > >& );

	/**
		Construct M by N matrix using given STL array for values.
		@param[in] matrix Matrix of data represented as C++ vector
		@retval this This matrix
	 */
	Matrix& operator=( std::vector< std::vector< value_type > >& );

	/**
		Construct M by N matrix using given array of values.
		@param[in] matrix_array Matrix of data represented as C-style array
		@param[in] M Number of rows in matrix
		@param[in] N Number of columns in matrix
	 */
	Matrix( value_type**, size_type, size_type );

	/**
	  	Deallocate internal matrix.
	 */
	virtual ~Matrix( );

	/**
		Create a copy of this matrix.

		@retval A Copy of this matrix
	 */
	Matrix copy( ) const;

	/**
		Return number of rows in matrix, ignoring the current matrix region.

		@retval num_rows Number of rows
	 */
	size_type real_row_size( ) const;

	/**
		Return number of rows in matrix.

		The size actually refers to the current submatrix region.

		@retval num_rows Number of rows
	 */
	size_type row_size( ) const;

	/**
		Return number of rows in matrix.
		@retval num_rows Number of rows
	 */
	size_type height( ) const;

	/**
		Return number of columns in matrix, ignoring the current matrix region.

		@retval num_rows Number of columns
	 */
	size_type real_col_size( ) const;

	/**
		Return number of columns in matrix.

		The size actually refers to the current submatrix region.

		@retval num_rows Number of columns
	 */
	size_type col_size( ) const;

	/**
		Return number of columns in matrix.
		@retval num_rows Number of columns
	 */
	size_type width( ) const;

	/**
		Return number of elements in matrix, ignoring the current matrix region.

		@retval num_elements Number of elements
	 */
	size_type real_size( ) const;

	/**
		Return number of elements in matrix.

		The size actually refers to the current submatrix region.

		@retval num_elements Number of elements
	 */
	size_type size( ) const;

	/**
		Access row i and column j in matrix.
		@retval ref Reference to row i in matrix
	 */
	reference at( size_type, size_type );

	/**
		Access row i and column j in matrix.
		@retval ref Reference to row i in matrix
	 */
	const_reference at( size_type, size_type ) const;

	/**
		Access row i and column j in matrix.
		@retval ref Reference to row i in matrix
	 */
	reference operator()( size_type, size_type );

	/**
		Access row i and column j in matrix.
		@retval ref Reference to row i in matrix
	 */
	const_reference operator()( size_type, size_type ) const;

	/**
		Get pointer to internal data.
		@retval _matrix Internal data
	 */
	pointer data( );

	/**
		Get read-only pointer to internal data.
		@retval _matrix Internal data
	 */
	const_pointer data( ) const;

	/**
		Get starting row position of the matrix.
		
		The position actually refers to the current submatrix region, which is the
		first row of the matrix (0) unless explicitly changed.

		@retval row_size Number of rows in region
	 */
	size_type row_begin( ) const;

	/**
		Get starting column position of the matrix.

		The position actually refers to the current submatrix region, which is the
		first column of the matrix (0) unless explicitly changed.

		@retval col_size Number of columns in region
	 */
	size_type col_begin( ) const;

	/**
		Return last row position in region.
		@retval row_end 
	 */
	size_type row_end( ) const;

	/**
		Return last column position in region.
		@retval num_rows Number of rows
	 */
	size_type col_end( ) const;

	/**
		Return number of rows in submatrix marked by region.

		The size actually refers to the current submatrix region, which is the
		entire matrix unless explicitly changed.

		@retval num_rows Number of rows
	 */
	size_type sub_row_size( ) const;

	/**
		Return number of columns in submatrix marked by region.

		The size actually refers to the current submatrix region, which is the
		entire matrix unless explicitly changed.

		@retval num_cols Number of columns
	 */
	size_type sub_col_size( ) const;

	/**
		Get current matrix region.
		@retval region Current matrix region
	 */
	const Matrix_Region& get_region( ) const;

	/**
		Get current matrix region.
		@retval region Current matrix region
	 */
	const Matrix_Region& region( ) const;

	/**
		Set matrix region to the given region.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
	 */
	const Matrix& set_region( size_type, size_type,
		size_type, size_type );

	/**
		Set matrix region to the given region.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
	 */
	const Matrix& region( size_type, size_type,
		size_type, size_type );

	/**
		Set matrix region to the given region.

		The starting positions are set first followed by the size.
		In the case of a region that does not fit inside the matrix exactly, this
		implicitly gives priority to the starting position of the region over the
		size of the region as the size is reduced so that the starting position +
		the size is less than or equal to the entire matrix's size. For example,
		suppose the matrix is 4 x 4 (where [0,0] is the first element) and the
		given region is [1,1] -> [5,5]; that is, the region starts at [1,1] and
		ends at [5,5] giving it a size of 4 in both dimensions. Since the given
		region does not fit inside the matrix, the region is truncated to
		[1,1] -> [4,4]. Also, if the given starting position does not fit inside
		the matrix, then the created region is
		[row_size,col_size] -> [row_size,col_size]; that is, a region with size 0
		that is placed just outside the matrix's bounds is used.

		@param[in] region New matrix region
	 */
	const Matrix& set_region( const Matrix_Region& );

	/**
		Set matrix region to the given region.

		The starting positions are set first followed by the size.
		In the case of a region that does not fit inside the matrix exactly, this
		implicitly gives priority to the starting position of the region over the
		size of the region as the size is reduced so that the starting position +
		the size is less than or equal to the entire matrix's size. For example,
		suppose the matrix is 4 x 4 (where [0,0] is the first element) and the
		given region is [1,1] -> [5,5]; that is, the region starts at [1,1] and
		ends at [5,5] giving it a size of 4 in both dimensions. Since the given
		region does not fit inside the matrix, the region is truncated to
		[1,1] -> [4,4]. Also, if the given starting position does not fit inside
		the matrix, then the created region is
		[row_size,col_size] -> [row_size,col_size]; that is, a region with size 0
		that is placed just outside the matrix's bounds is used.

		@param[in] region New matrix region
	 */
	const Matrix& region( const Matrix_Region& );

	/**
		Save current region.
		This is useful if only the current region needs to be saved before
		changing to multiple regions that do not need saving.
	 */
	void save_region( );

	/**
		Save current region and set new region.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
	 */
	void push_region( size_type, size_type,
		size_type, size_type);

	/**
		Save current region and set new region.
		@param[in] region New matrix region
	 */
	void push_region( const Matrix_Region& );

	/**
		Restore last saved region.
	 */
	void pop_region( );

	/**
		Restore last saved region.
	 */
	void restore_region( );

	/**
		Reset matrix region to the entire matrix.
	 */
	void reset_region( );

	Vector get_col( size_type ) const;
	void set_col( Vector&, size_type );

	/**
	  Determine if file is a matrix or not by its name. Can be passed to
	  ws_tools::dir_traverse() as a file name filter.
	  @param[in] file_name Name of file to check
	 */
	static bool is_matrix( const std::string& );

	/**
	  Determine if file is a ascii matrix or not by its name.
	  @param[in] file_name Name of file to check
	 */
	static bool is_ascii_matlab( const std::string& );

	/**
	  Determine if file is a binary matrix or not by its name.
	  @param[in] file_name Name of file to check
	 */
	static bool is_binary_matlab( const std::string& );

	/**
	  Determine if file is a Lush matrix or not by its name.
	  @param[in] file_name Name of file to check
	 */
	static bool is_lush( const std::string& );

protected:  // section to define iterators

	/**
		@brief Traits for a non-constant matrix, column, or row iterator.
	 */
	struct Nonconst_traits
	{
		typedef Matrix::value_type       value_type;
		typedef Matrix::size_type        size_type;
		typedef Matrix::reference        reference;
		typedef Matrix::pointer          pointer;
		typedef Matrix::difference_type  difference_type;
		typedef Matrix::pointer          iterator;
	};

	/**
		@brief Traits for a constant matrix, column, or row iterator.
	 */
	struct Const_traits
	{
		typedef Matrix::value_type           value_type;
		typedef Matrix::size_type            size_type;
		typedef Matrix::const_reference      reference;
		typedef Matrix::const_pointer        pointer;
		typedef Matrix::difference_type      difference_type;
		typedef Matrix::const_pointer        iterator;
	};

public:  // section to define iterators

	// define iterator type with given traits
	typedef _matrix_iterator<Nonconst_traits>   iterator;
	typedef _matrix_iterator<Const_traits>      const_iterator;

	// define row and column iterator types with given traits
	typedef _row_col_iterator<Nonconst_traits>  row_col_iterator;
	typedef _row_col_iterator<Const_traits>     const_row_col_iterator;
	typedef row_col_iterator                    col_iterator;
	typedef const_row_col_iterator              const_col_iterator;
	typedef row_col_iterator                    row_iterator;
	typedef const_row_col_iterator              const_row_iterator;

	/**
		Get iterator to beginning of matrix region.
		@retval iter Iterator to beginning of matrix region
	 */
	iterator begin( );

	/**
		Get iterator to beginning of matrix region.
		@retval iter Iterator to beginning of matrix region
	 */
	const_iterator const_begin( ) const;

	/**
		Get iterator to end of matrix region.
		@retval iter Iterator to end of matrix region
	 */
	iterator end( );

	/**
		Get iterator to end of matrix region.
		@retval iter Iterator to end of matrix region
	 */
	const_iterator const_end( ) const;

	/**
		Get iterator to start of column j in matrix.

		If the column is not inside the current active region, then an iterator
		to the end of the column is returned.

		@param[in] j Column index
		@retval col_iter Iterator to column
	 */
	col_iterator col_begin( size_type );

	/**
		Get iterator to start of column j in matrix.

		If the column is not inside the current active region, then an iterator
		to the end of the column is returned.

		@param[in] j Column index
		@retval col_iter Iterator to column
	 */
	const_col_iterator const_col_begin( size_type ) const;

	/**
		Get iterator to end of column j in matrix.
		@param[in] j Column index
		@retval col_iter Iterator to column
	 */
	col_iterator col_end( size_type );

	/**
		Get iterator to end of column j in matrix.
		@param[in] j Column index
		@retval col_iter Iterator to column
	 */
	const_col_iterator const_col_end( size_type ) const;

	/**
		Get iterator to start of row i in matrix.

		If the row is not inside the current active region, then an iterator
		to the end of the row is returned.

		@param[in] i Row index
		@retval row_iter Iterator to row
	 */
	row_iterator row_begin( size_type );

	/**
		Get iterator to start of row i in matrix.

		If the row is not inside the current active region, then an iterator
		to the end of the row is returned.

		@param[in] i Row index
		@retval row_iter Iterator to row
	 */
	const_row_iterator const_row_begin( size_type ) const;

	/**
		Get iterator to end of row i in matrix.
		@param[in] i Row index
		@retval row_iter Iterator to row
	 */
	row_iterator row_end( size_type );

	/**
		Get iterator to end of row i in matrix.
		@param[in] i Row index
		@retval row_iter Iterator to row
	 */
	const_row_iterator const_row_end( size_type ) const;

public:

	void read( const std::string& );
	void write( FILE*, prec_type = 5 ) const;
	void write( const std::string& = "", prec_type = 5 ) const;
	void clear( );
	void assign( const value_type );
	Matrix& sort( );
	Matrix& reverse_sort( );

	// unary negation
	Matrix operator-( ) const;

	// binary operators for two matrices
	Matrix& operator+=( const Matrix& );
	Matrix  operator+( const Matrix& ) const;
	Matrix& operator-=( const Matrix& );
	Matrix  operator-( const Matrix& ) const;
	Matrix& operator*=( const Matrix& );
	Matrix  operator*( const Matrix& ) const;

	// binary operators for matrix and value
	Matrix& operator+=( const value_type );
	Matrix  operator+( const value_type ) const;
	Matrix& operator-=( const value_type );
	Matrix  operator-( const value_type ) const;
	Matrix& operator*=( const value_type );
	Matrix  operator*( const value_type ) const;
	Matrix& operator/=( const value_type );
	Matrix  operator/( const value_type ) const;

	/**
		Add matrix by value.
		Some older versions of gcc complain if friend functions are not defined
		in a header file.
		@param[in] rhs Value on right-hand side of +
		@retval matrix Matrix result 
	 */
	friend Matrix operator+( const value_type value, const Matrix& rhs )
	{
		return( rhs + value );
	}

	/**
		Subtract matrix by value.
		Some older versions of gcc complain if friend functions are not defined
		in a header file.
		@param[in] rhs Value on right-hand side of -
		@retval matrix Matrix result 
	 */
	friend Matrix operator-( const value_type value, const Matrix& rhs )
	{
		return( -rhs + value );
	}

	/**
		Multiply matrix by value.
		Some older versions of gcc complain if friend functions are not defined
		in a header file.
		@param[in] rhs Value on right-hand side of *
		@retval matrix Matrix result 
	 */
	friend Matrix operator*( const value_type value, const Matrix& rhs )
	{
		return( rhs * value );
	}

	/**
		Divide matrix by value.
		Some older versions of gcc complain if friend functions are not defined
		in a header file.
		@param[in] rhs Value on right-hand side of /
		@retval matrix Matrix result 
	 */
	friend Matrix operator/( const value_type value, const Matrix& rhs )
	{
		return( rhs / value );
	}

	// Matrix& operator*=( const Vector& );
	// Matrix operator*( const Vector& ) const;

	Vector avg( ) const;
	Vector var( ) const;
	void stat( Vector&, Vector&, Vector& ) const;
	void stat( Vector&, Vector& ) const;
	void stat( double*, double*, double* = 0 ) const;
	Vector kurtosis( ) const;
	Vector sum( ) const;

	value_type dot_product( const Matrix& ) const;

	Matrix transpose( ) const;

	value_type two_norm( bool = true ) const;
	value_type two_norm_2( bool = true ) const;
	value_type inf_norm( );
	Matrix     sin( ) const;
	Matrix     cos( ) const;
	Matrix     tan( ) const;
	void normalize( );
	Matrix     abs( ) const;

	static Matrix identity( size_type );
	static Matrix zeros( size_type );
	static Matrix zeros( size_type, size_type );
	static Matrix ones( size_type );
	static Matrix ones( size_type, size_type );
	static Matrix random( size_type, size_type, value_type = INT_MAX );
	static Matrix random( size_type, size_type, value_type, value_type );
	static Matrix random_normal( size_type, size_type, double = 0, double = 1 );
	static Matrix toeplitz( const Vector&, const Vector& );
	static Matrix toeplitz( const Vector& );

	static Matrix multiply( std::vector< Matrix >& );

	//Vector to_vector( ) const;
	double* to_array( ) const;

	value_type get_min( ) const;
	value_type get_max( ) const;
	void get_min_max( value_type*, value_type* ) const;
	Matrix& convolve_row( const Vector&, Matrix& );
	Matrix& convolve_row( const Vector& );
	Matrix& convolve_col( const Vector&, Matrix& );
	Matrix& convolve_col( const Vector& );
	Matrix& convolve( const Vector&, const Vector& );
	Matrix& convolve( const Matrix& );

	void prin_comp( unsigned );
	void covariance( Matrix&, Vector& );
	Matrix covariance( ) const;

	Matrix reshape( size_type, size_type ) const;

	void fft( Matrix&, Matrix& ) const;
	void ifft( Matrix&, Matrix& );

protected:

	void allocate_matrix( size_type, size_type );
	void write_ascii_matlab( const std::string& = "", prec_type = 1 ) const;
	void write_binary_matlab( const std::string& = "", prec_type = 1 ) const;
	void write_lush( const std::string& = "", prec_type = 1 ) const;

	// Note: The next 4 functions help multiply(vector) function above.
	static Vector_Int initialize_dim( std::vector< Matrix >& );
	static Matrix_Int best_mult_order( std::vector< Matrix >& matrices,
		Vector_Int& dim );
	static std::string make_exp( Matrix_Int&, int, int );
	static Matrix multiply( std::vector< Matrix >&, Matrix_Int&, int, int );

public:  // scale functions
	Matrix& resize( size_type, size_type, bool = false );
	Matrix& scale( double );
	Matrix& scale( double, double );
	Matrix& bicubic_scale( size_type, size_type );
	Matrix& nointerp_scale( double, double );

protected:  // support functions for bicubic_scale()
	void get_scaled_row( double*, double**, int, size_type ) const; 
	void upsample_row( const double*, double*, size_type ) const;
	void downsample_row( const double*, double*, size_type ) const;
	void rotate_pointers( double** ) const;
	double cubic( double, double, double, double, double ) const;
	double cubic_0( double x ) const;
	double cubic_1( double x ) const;
	double cubic_2( double x ) const;
	double cubic_3( double x ) const;

private:

	/*
		Definition of member variables
	 */

	matrix_type    _matrix;        //< Internal matrix of data
	size_type      _row_size;      //< Number of rows in matrix
	size_type      _col_size;      //< Number of columns in matrix

	Matrix_Region  _region;        //< Current working subregion in matrix
	Region_List    _region_stack;  //< List of saved regions in the order saved
};

/*
	Define inlined member functions.
 */

/**
	Create a copy of this matrix.

	@retval A Copy of this matrix
 */
inline Matrix
Matrix::copy( ) const
{
	return( Matrix( *this ) );
}

/**
	Return number of rows in matrix, ignoring the current matrix region.

	@retval num_rows Number of rows
 */
inline Matrix::size_type
Matrix::real_row_size( ) const
{ 
	return( _row_size );	
}

/**
	Return number of rows in matrix.

	The size actually refers to the current submatrix region, which is the
	entire matrix unless explicitly changed.

	@retval num_rows Number of rows
 */
inline Matrix::size_type
Matrix::row_size( ) const
{ 
	return( sub_row_size() );	
	// return( real_row_size() );	
}

/**
	Return number of rows in matrix.
	@retval num_rows Number of rows
 */
inline Matrix::size_type
Matrix::height( ) const
{ 
	return( row_size() );	
}

/**
	Return number of columns in matrix, ignoring the current matrix region.

	@retval num_rows Number of columns
 */
inline Matrix::size_type
Matrix::real_col_size( ) const
{ 
	return( _col_size );	
}

/**
	Return number of columns in matrix.

	The size actually refers to the current submatrix region, which is the
	entire matrix unless explicitly changed.

	@retval num_rows Number of columns
 */
inline Matrix::size_type
Matrix::col_size( ) const
{ 
	return( sub_col_size() );	
	// return( real_col_size() );	
}

/**
	Return number of columns in matrix.
	@retval num_rows Number of columns
 */
inline Matrix::size_type
Matrix::width( ) const
{ 
	return( col_size() );	
}

/**
	Return number of elements in matrix, ignoring the current matrix region.

	@retval num_elements Number of elements
 */
inline Matrix::size_type
Matrix::real_size( ) const
{
	return( real_row_size() * real_col_size() );
}

/**
	Return number of elements in matrix.

	The size actually refers to the current submatrix region.

	@retval num_elements Number of elements
 */
inline Matrix::size_type
Matrix::size( ) const
{
	return( row_size() * col_size() );
}

/**
	Access row i and column j in matrix.
	@retval ref Reference to row i in matrix
 */
inline Matrix::reference
Matrix::at( size_type i, size_type j )
{
	assert( i >= 0 && i < real_row_size() );
	assert( j >= 0 && j < real_col_size() );
	return( _matrix[ (i * real_col_size()) + j ] );
}

/**
	Access row i and column j in matrix.
	@retval ref Reference to row i in matrix
 */
inline Matrix::const_reference
Matrix::at( size_type i, size_type j ) const
{
	assert( i >= 0 && i < real_row_size() );
	assert( j >= 0 && j < real_col_size() );
	return( _matrix[ (i * real_col_size()) + j ] );
}

/**
	Access row i and column j in matrix.
	@retval ref Reference to row i in matrix
 */
inline Matrix::reference
Matrix::operator()( size_type i, size_type j )
{
	return( at(i,j) );
}

/**
	Access row i and column j in matrix.
	@retval ref Reference to row i in matrix
 */
inline Matrix::const_reference
Matrix::operator()( size_type i, size_type j ) const
{
	return( at(i,j) );
}

/**
	Get pointer to internal data.
	@retval _matrix Internal data
 */
inline Matrix::pointer
Matrix::data( )
{
#ifdef HAVE_BOOST_SCOPED_PTR_WS_IMG
	return( _matrix.get() );
#else
#ifdef USE_STL
	return( &_matrix[0] );
#else
	return( _matrix );
#endif // USE_STL
#endif // HAVE_BOOST_SCOPED_PTR_WS_IMG
}

/**
	Get read-only pointer to internal data.
	@retval _matrix Internal data
 */
inline Matrix::const_pointer
Matrix::data( ) const
{
#ifdef HAVE_BOOST_SCOPED_PTR_WS_IMG
	return( _matrix.get() );
#else
#ifdef USE_STL
	return( &_matrix[0] );
#else
	return( _matrix );
#endif // USE_STL
#endif // HAVE_BOOST_SCOPED_PTR_WS_IMG
}

/**
	Get starting row position of the matrix.
	
	The position actually refers to the current submatrix region, which is the
	first row of the matrix (0) unless explicitly changed.

	@retval row_size Number of rows in region
 */
inline Matrix::size_type
Matrix::row_begin( ) const
{
	return( region().row_begin() );
}

/**
	Get starting column position of the matrix.

	The position actually refers to the current submatrix region, which is the
	first column of the matrix (0) unless explicitly changed.

	@retval col_size Number of columns in region
 */
inline Matrix::size_type
Matrix::col_begin( ) const
{
	return( region().col_begin() );
}

/**
	Return last row position in region.
	@retval row_end 
 */
inline Matrix::size_type
Matrix::row_end( ) const
{ 
	return( region().row_end() );
}

/**
	Return last column position in region.
	@retval num_rows Number of rows
 */
inline Matrix::size_type
Matrix::col_end( ) const
{ 
	return( region().col_end() );
}

/**
	Return number of rows in submatrix marked by region.

	The size actually refers to the current submatrix region.

	@retval num_rows Number of rows
 */
inline Matrix::size_type
Matrix::sub_row_size( ) const
{ 
	return( region().row_size() );	
}

/**
	Return number of columns in submatrix marked by region.

	The size actually refers to the current submatrix region.

	@retval num_cols Number of columns
 */
inline Matrix::size_type
Matrix::sub_col_size( ) const
{ 
	return( region().col_size() );	
}

/**
	Get current matrix region.
	@retval region Current matrix region
 */
inline const Matrix_Region&
Matrix::get_region( ) const
{
	return( _region );
}

/**
	Get current matrix region.
	@retval region Current matrix region
 */
inline const Matrix_Region&
Matrix::region( ) const
{
	return( _region );
}

/**
	Set matrix region to the given region.

	@param[in] row_begin Starting row position of the region
	@param[in] col_begin Starting col position of the region
	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region

	@retval this This matrix
 */
inline const Matrix&
Matrix::set_region( size_type row_begin, size_type col_begin,
	size_type row_size, size_type col_size )
{
	return( region( row_begin, col_begin, row_size, col_size ) );
}

/**
	Set matrix region to the given region.

	@param[in] row_begin Starting row position of the region
	@param[in] col_begin Starting col position of the region
	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region

	@retval this This matrix
 */
inline const Matrix&
Matrix::region( size_type row_begin, size_type col_begin,
	size_type row_size, size_type col_size )
{
	const Matrix_Region new_region( row_begin, col_begin,
			row_size, col_size );
	return( region( new_region ) );
}

/**
	Set matrix region to the given region.

	The starting positions are set first followed by the size.
	In the case of a region that does not fit inside the matrix exactly, this
	implicitly gives priority to the starting position of the region over the
	size of the region as the size is reduced so that the starting position +
	the size is less than or equal to the entire matrix's size. For example,
	suppose the matrix is 4 x 4 (where [0,0] is the first element) and the
	given region is [1,1] -> [5,5]; that is, the region starts at [1,1] and
	ends at [5,5] giving it a size of 4 in both dimensions. Since the given
	region does not fit inside the matrix, the region is truncated to
	[1,1] -> [4,4]. Also, if the given starting position does not fit inside
	the matrix, then the created region is
	[row_size,col_size] -> [row_size,col_size]; that is, a region with size 0
	that is placed just outside the matrix's bounds is used.

	@param[in] region New matrix region

	@retval this This matrix
 */
inline const Matrix&
Matrix::set_region( const Matrix_Region& new_region )
{
	return( region( new_region ) );
}

/**
	Set matrix region to the given region.

	The starting positions are set first followed by the size.
	In the case of a region that does not fit inside the matrix exactly, this
	implicitly gives priority to the starting position of the region over the
	size of the region as the size is reduced so that the starting position +
	the size is less than or equal to the entire matrix's size. For example,
	suppose the matrix is 4 x 4 (where [0,0] is the first element) and the
	given region is [1,1] -> [5,5]; that is, the region starts at [1,1] and
	ends at [5,5] giving it a size of 4 in both dimensions. Since the given
	region does not fit inside the matrix, the region is truncated to
	[1,1] -> [4,4]. Also, if the given starting position does not fit inside
	the matrix, then the created region is
	[row_size,col_size] -> [row_size,col_size]; that is, a region with size 0
	that is placed just outside the matrix's bounds is used.

	@param[in] region New matrix region

	@retval this This matrix
 */
inline const Matrix&
Matrix::region( const Matrix_Region& new_region )
{
	// set starting positions 
	_region.set_row_begin(
			std::min<size_type>( new_region.row_begin(), _row_size ) );

	_region.set_col_begin( 
			std::min<size_type>( new_region.col_begin(), _col_size ) );

	// set sizes
	_region.set_row_size(
		std::min<size_type>( new_region.row_size(),
		_row_size - _region.row_begin() )
	);

	_region.set_col_size(
		std::min<size_type>( new_region.col_size(),
		_col_size - _region.col_begin() )
	);
	return( *this );
}

/**
	Save current region.
	This is useful if only the current region needs to be saved before
	changing to multiple regions that do not need saving.
 */
inline void
Matrix::save_region( )
{
	push_region( region() );
}

/**
	Save current region and set new region.

	@param[in] row_begin Starting row position of the region
	@param[in] col_begin Starting col position of the region
	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region
 */
inline void
Matrix::push_region( size_type row_begin, size_type col_begin,
	size_type row_size, size_type col_size )
{
	Matrix_Region region( row_begin, col_begin, row_size, col_size );
	push_region( region );
}

/**
	Save current region and set new region.
	@param[in] region New matrix region
 */
inline void
Matrix::push_region( const Matrix_Region& new_region )
{
	_region_stack.push_back( region() );
	set_region( new_region );
}

/**
	Restore last saved region.
 */
inline void
Matrix::pop_region( )
{
	if( !_region_stack.empty() )
	{
		set_region( _region_stack.back() );
		_region_stack.pop_back();
	}
}

/**
	Restore last saved region.
 */
inline void
Matrix::restore_region( )
{
	pop_region();
}

/**
	Reset matrix region to the entire matrix.
 */
inline void
Matrix::reset_region( )
{
	_region_stack.clear();
	_region.set_row_begin( 0 );
	_region.set_col_begin( 0 );
	_region.set_row_size( _row_size );
	_region.set_col_size( _col_size );
}

/**
	Cubic interpolation function.

	Note: the cubic function no longer clips results. Interpolation using
	Catmull-Rom splines.

	@param[in] x Input value
	@retval value Interpolated value
 */
inline double
Matrix::cubic( double x, double jm1, double j, double jp1, double jp2 ) const
{
	return(
			(
				(
					(
						(
							-jm1 + 3*j - 3*jp1 + jp2
						) * x + (2*jm1 - 5*j + 4*jp1 - jp2)
					) * x + (-jm1 + jp1)
				) * x + (j + j)
			) / 2.0
	);
}

/**
	Cubic interpolation function with fixed parameters such that it is
	equivalent to cubic( x, 1, 0, 0, 0 ).

	@param[in] x Input value
	@retval value Interpolated value
 */
inline double
Matrix::cubic_0( double x ) const
{
	const double x_squared = x * x;
	const double x_cubed   = x * x_squared;

	return( (-x_cubed + (2 * x_squared) - x) / 2.0 );
}

/**
	Cubic interpolation function with fixed parameters such that it is
	equivalent to cubic( x, 0, 1, 0, 0 ).

	@param[in] x Input value
	@retval value Interpolated value
 */
inline double
Matrix::cubic_1( double x ) const
{
	const double x_squared = x * x;
	const double x_cubed   = x * x_squared;

	return( ((3 * x_cubed) - (5 * x_squared) + 2) / 2.0 );
}

/**
	Cubic interpolation function with fixed parameters such that it is
	equivalent to cubic( x, 0, 0, 1, 0 ).

	@param[in] x Input value
	@retval value Interpolated value
 */
inline double
Matrix::cubic_2( double x ) const
{
	const double x_squared = x * x;
	const double x_cubed   = x * x_squared;

	return( ((-3 * x_cubed) + (4 * x_squared) + x) / 2.0 );
}

/**
	Cubic interpolation function with fixed parameters such that it is
	equivalent to cubic( x, 0, 0, 0, 1 ).

	@param[in] x Input value
	@retval value Interpolated value
 */
inline double
Matrix::cubic_3( double x ) const
{
	const double x_squared = x * x;
	const double x_cubed   = x * x_squared;

	return( (x_cubed - x_squared) / 2.0 );
}

/**
	Circularly rotate pointers.

	This is used so we can store 4 rows at a time since we need to interpolate
	down each column, which requires 4 values.  
	Since a given row gets reused 4 times, it is more efficient to just move the
	last 3 row pointers up rather than recopy each row.

	@param[in,out] p Array of pointers to rotate
	@param[in] n Size of p
 */
inline void
// Matrix::rotate_pointers( double** p, unsigned n ) const
Matrix::rotate_pointers( double** p ) const
{
	const unsigned n = 4;

	// shift each pointer down and move the first pointer to the end
	double* tmp = p[0];
	for( unsigned i = 0; i != (n - 1); ++i )
	{
		p[i] = p[i+1];
	}
	p[n - 1] = tmp;
}

} // namespace ws_img

#endif // MATRIX_HPP
