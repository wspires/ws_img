/**
	@file   Matrix.cpp
	@author Wade Spires
	@date   2005/10/28
	@brief  Implementation of matrix type.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Matrix.hpp"
#include "Vector.hpp"

#include  <new>
#include  <cstdlib>
#include  <cmath>

#ifdef HAVE_BOOST_SCOPED_PTR_WS_IMG

// boost headers
#include "boost/random.hpp"
#include <boost/generator_iterator.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>

#endif // HAVE_BOOST_SCOPED_PTR_WS_IMG

#ifdef HAVE_FFTW
	#include <complex>
	#include <fftw3.h>

	using std::complex;
#endif // HAVE_FFTW

using namespace ws_img;

using std::min;
using std::max;
using std::min_element;
using std::max_element;
using std::string;
using std::vector;

typedef Matrix::value_type       value_type;
typedef Matrix::size_type        size_type;
typedef Matrix::reference        reference;
typedef Matrix::const_reference  const_reference;
typedef Matrix::pointer          pointer;
typedef Matrix::const_pointer    const_pointer;
typedef Matrix::prec_type        prec_type;
typedef Matrix::difference_type  difference_type;

typedef Matrix::iterator            iterator;
typedef Matrix::const_iterator      const_iterator;
typedef Matrix::col_iterator        col_iterator;
typedef Matrix::const_col_iterator  const_col_iterator;
typedef Matrix::row_iterator        row_iterator;
typedef Matrix::const_row_iterator  const_row_iterator;

const value_type Matrix::MAX_VALUE = DBL_MAX;

/**
	Construct matrix by reading given matrix file.
	@param[in] matrix_file_name Name of matrix file to read matrix from
 */
Matrix::Matrix( const std::string& matrix_file_name )
: _matrix(0)
{
	read( matrix_file_name );
}

/**
	Construct square M-by-M matrix with each element left uninitialized.

	Note that the elements are not initialized to any value, not even 0. Even
	if by inspection this seems to be the case, do not rely upon the values
	being set.

	@param[in] M Number of rows and columns in matrix
 */
Matrix::Matrix( size_type M )
: _matrix(0), _row_size(M), _col_size(M), _region(0, 0, M, M)
{
	assert( M > 0 );
	allocate_matrix( _row_size, _col_size );
}

/**
	Construct M-by-N matrix with each element left uninitialized.

	Note that the elements are not initialized to any value, not even 0. Even
	if by inspection this seems to be the case, do not rely upon the values
	being set.

	@param[in] M Number of rows in matrix
	@param[in] N Number of columns in matrix
 */
Matrix::Matrix( size_type M, size_type N )
: _matrix(0), _row_size(M), _col_size(N), _region(0, 0, M, N)
{
	assert( M > 0 && N > 0 );
	allocate_matrix( _row_size, _col_size );
}

/**
	Construct M-by-N matrix with each element initialized to the given value.
	@param[in] M Number of rows in matrix
	@param[in] N Number of columns in matrix
	@param[in] default_value Initial value to set each element of matrix
 */
Matrix::Matrix( size_type M, size_type N, value_type default_value )
: _matrix(0), _row_size(M), _col_size(N), _region(0, 0, M, N)
{
	assert( M > 0 && N > 0 );
	allocate_matrix( _row_size, _col_size );
	*this = default_value;
}

/**
	Copy constructor copies matrix region.

	TODO Check whether better to copy entire matrix and just set the region
	rather than only copy the region

	@param[in] m Matrix to copy
 */
Matrix::Matrix( const Matrix& m )
:
#ifdef USE_STL
	_matrix(m.sub_row_size() * m.sub_col_size() ),
#else
	_matrix(0),
#endif // USE_STL
_row_size( m.sub_row_size() ), _col_size( m.sub_col_size() ),
_region(0, 0, m.sub_row_size(), m.sub_col_size() )
{
#ifdef HAVE_BOOST_SCOPED_PTR_WS_IMG
	_matrix.reset( new value_type[ m.sub_row_size() * m.sub_col_size() ] );
	if( _matrix == 0 )
	{
		err_quit( "Unable to construct matrix: malloc failed\n" );
	}
#else
#ifndef USE_STL
	_matrix = new value_type[ m.sub_row_size() * m.sub_col_size() ];
	// _matrix = (value_type*) malloc( _row_size * _col_size
			// * sizeof(value_type) );
	if( _matrix == 0 )
	{
		err_quit( "Unable to construct matrix: malloc failed\n" );
	}
#endif // USE_STL
#endif // HAVE_BOOST_SCOPED_PTR_WS_IMG

	// TODO See if memcpy is faster. We may have to memcpy each row inside a loop
	// in order to only copy a region.
	// for( size_type i = row_begin(), im = m.row_begin();
			// i != row_end(); ++i, ++im )
	// {
			// memcpy( 
					// data() + (i * real_col_size()) + col_begin(),
					// m.data() + (im * m.real_col_size()) + m.col_begin(),
					// sizeof(value_type) * m.sub_col_size()
			// );
	// }

	for( size_type i = row_begin(), im = m.row_begin();
			i != row_end(); ++i, ++im )
	{
		for( size_type j = col_begin(), jm = m.col_begin();
				j != col_end(); ++j, ++jm )
		{
			at(i,j) = m(im,jm);
		}
	}
}

/**
	Assignment operator copies matrix region.
	@param[in] m Matrix to copy
 */
Matrix&
Matrix::operator=( const Matrix& m )
{
	if( &m == this )
	{
		return( *this );
	}
	allocate_matrix( m.sub_row_size(), m.sub_col_size() );

	for( size_type i = row_begin(), im = m.row_begin();
			i != row_end(); ++i, ++im )
	{
		for( size_type j = col_begin(), jm = m.col_begin();
				j != col_end(); ++j, ++jm )
		{
			at(i,j) = m(im,jm);
		}
	}

	return( *this );
}

/**
	Construct M by 1 matrix using given column vector.
	@param[in] v Vector to initialize matrix to
 */
Matrix::Matrix( const Vector& v ) 
: _matrix(0), _row_size(v.sub_size()), _col_size(1),
_region(0, 0, v.sub_size(), 1)
{
	allocate_matrix( v.sub_size(), 1 );
	size_type iv = v.vec_begin();
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			at(i,j) = v(iv++);
		}
	}

	/*
	// iterators to this matrix
	iterator iter     = begin();
	iterator end_iter = end();

	// iterators to vector
	Vector::const_iterator v_iter     = v.const_begin();
	Vector::const_iterator v_end_iter = v.const_end();

	// copy vector
	for( ; iter != end_iter; ++iter, ++v_iter )
	{
		*iter = *v_iter;
	}
	 */
}

/**
	Construct M by 1 matrix using given column vector.
	@param[in] v Vector to initialize matrix to
	@retval this This matrix
 */
Matrix&
Matrix::operator=( const Vector& v ) 
{
	allocate_matrix( v.sub_size(), 1 );

	size_type iv = v.vec_begin();
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			at(i,j) = v(iv++);
		}
	}

	/*
	// iterators to this matrix
	iterator iter     = begin();
	iterator end_iter = end();

	// iterators to vector
	Vector::const_iterator v_iter     = v.const_begin();
	Vector::const_iterator v_end_iter = v.const_end();

	// copy vector
	for( ; iter != end_iter; ++iter, ++v_iter )
	{
		*iter = *v_iter;
	}
	 */

	return( *this );
}

/**
	Construct matrix using the given column vector v such that v is
	copied N times along the columns of this matrix.
	@param[in] v Vector to initialize matrix to
	@param[in] N Number of copies of v
 */
Matrix::Matrix( const Vector& v, size_type N ) 
: _matrix(0), _row_size(v.sub_size()), _col_size(N),
_region(0, 0, v.sub_size(), N)
{
	allocate_matrix( N, v.sub_size() );

	// copy v onto each column of matrix
	size_type iv = v.vec_begin();
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			at(i,j) = v(iv);
		}
		++iv;
	}
}

/**
	Construct matrix consisting of an M-by-N tiling of copies of v.
	@param[in] v Vector to initialize matrix to
	@param[in] M Number of row tilings of v
	@param[in] N Number of column tilings of v
 */
Matrix::Matrix( const Vector& v, size_type M, size_type N ) 
: _matrix(0), _row_size(M * v.sub_size()), _col_size(N),
_region(0, 0, M * v.sub_size(), N)
{
	allocate_matrix( M, N );

	// for each column, copy v along the column
	for( size_type j = col_begin(); j != col_end(); ++j )
	{
		for( size_type i = row_begin(), k = v.vec_begin();
				i != row_end(); ++i, ++k )
		{
			at(i,j) = v(k);

			// return to start of vector to copy another vector
			if( k == v.vec_end() )
			{
				k = v.vec_begin();
			}
		}
	}
}

/**
	Construct M by N matrix using given STL array for values.
	@param[in] matrix Matrix of data represented as C++ vector
 */
Matrix::Matrix( std::vector< std::vector< value_type > >& matrix )
: _matrix(0), _row_size(matrix.size()), _col_size(matrix[0].size()),
	_region(0, 0, matrix.size(), matrix[0].size())
{
	// verify that matrix has positive dimension
	if( matrix.empty() || matrix[0].empty() )
	{
		return;
	}

	allocate_matrix( matrix.size(), matrix[0].size() );

	// copy each element from array into this matrix
	for( size_type i = 0; i != real_row_size(); i++ )
	{
		// verify that each column is the same size
		assert( matrix[i].size() == real_col_size() );

		for( size_type j = 0; j != real_col_size(); j++ )
		{
			at(i, j) = matrix[ i ][ j ];
		}
	}
}

/**
	Construct M by N matrix using given STL array for values.
	@param[in] matrix Matrix of data represented as C++ vector
	@retval this This matrix
 */
Matrix&
Matrix::operator=( std::vector< std::vector< value_type > >& matrix )
{
	// verify that matrix has positive dimension
	if( matrix.empty() || matrix[0].empty() )
	{
		return( *this );
	}

	allocate_matrix( matrix.size(), matrix[0].size() );

	// copy each element from array into this matrix
	for( size_type i = 0; i != real_row_size(); i++ )
	{
		// verify that each column is the same size
		assert( matrix[i].size() == real_col_size() );

		for( size_type j = 0; j != real_col_size(); j++ )
		{
			at(i, j) = matrix[ i ][ j ];
		}
	}

	return( *this );
}

/**
	Construct M by N matrix using given array of values.
	@param[in] matrix_array Matrix of data represented as C-style array
	@param[in] M Number of rows in matrix
	@param[in] N Number of columns in matrix
 */
Matrix::Matrix( value_type** matrix_array, size_type M, size_type N )
: _matrix(0), _row_size(M), _col_size(N), _region(0, 0, M, N)
{ 
	allocate_matrix( _row_size, _col_size );

	// verify that matrix has positive dimension
	if( M == 0 )
	{
		return;
	}
	if( N == 0 )
	{
		return;
	}

	// copy each element from array into this matrix
	for( size_type i = 0; i != real_row_size(); i++ )
	{
		for( size_type j = 0; j != real_col_size(); j++ )
		{
			at(i,j) = matrix_array[ i ][ j ];
		}
	}
}

/**
	Deallocate internal matrix.
 */
Matrix::~Matrix( )
{
#ifndef HAVE_BOOST_SCOPED_PTR_WS_IMG
#ifndef USE_STL
	if( _matrix != 0 )
	{
		delete [] _matrix;
		// free( _matrix );
		_matrix = 0;
	}
#endif // USE_STL
#endif // HAVE_BOOST_SCOPED_PTR_WS_IMG
}

/**
	Allocate space for matrix.

	The previous data is deleted before allocation.
	The values are left uninitialized after allocation.

	@param[in] num_rows Number of rows
	@param[in] num_cols Number of columns
 */
void
Matrix::allocate_matrix( size_type num_rows, size_type num_cols )
{
	// delete previous matrix data
#ifdef HAVE_BOOST_SCOPED_PTR_WS_IMG

	// Boost: nothing to do since deletion is automatic

#else
	#ifdef USE_STL

	// STL: clear vector's data
	_matrix.clear();

	#else

	// Raw pointer: delete array
	if( _matrix != 0 )
	{
		delete [] _matrix;
	}

	#endif // USE_STL
#endif // HAVE_BOOST_SCOPED_PTR_WS_IMG

	// set new dimensions
	_row_size = num_rows;
	_col_size = num_cols;
	reset_region();

	// allocate space for new data
#ifdef HAVE_BOOST_SCOPED_PTR_WS_IMG

	// Boost: reset the smart pointer
	_matrix.reset( new value_type[ real_row_size() * real_col_size() ] );

#else
	#ifdef USE_STL

	// STL: resize the vector
	_matrix.reserve( real_row_size() * real_col_size() );

	#else

	// Raw pointer: allocate new array
	try
	{
		_matrix = new value_type[ real_row_size() * real_col_size() ];
	}
	catch( std::bad_alloc& except )
	{
		err_quit( "Matrix::allocate_matrix: "
				"Unable to construct %u x %u matrix: new failed\n",
				real_row_size(), real_col_size() );
	}

	#endif // USE_STL
#endif // HAVE_BOOST_SCOPED_PTR_WS_IMG
}

/**
  Determine if file is a matrix or not by its name. Can be passed to
  ws_tools::dir_traverse() as a file name filter.
  @param[in] file_name Name of file to check
 */
bool
Matrix::is_matrix( const std::string& file_name )
{
	return( is_ascii_matlab(file_name) || is_binary_matlab(file_name)
		|| is_lush(file_name) );
}

/**
  Determine if file is a ascii matrix or not by its name.
  @param[in] file_name Name of file to check
 */
bool
Matrix::is_ascii_matlab( const std::string& file_name )
{
	std::string ext = ws_tools::get_ext_name( file_name );
	return( ext == "m" );
}

/**
  Determine if file is a binary matrix or not by its name.
  @param[in] file_name Name of file to check
 */
bool
Matrix::is_binary_matlab( const std::string& file_name )
{
	std::string ext = ws_tools::get_ext_name( file_name );
	return( ext == "mat" );
}

/**
  Determine if file is a Lush matrix or not by its name.
  @param[in] file_name Name of file to check
 */
bool
Matrix::is_lush( const std::string& file_name )
{
	std::string ext = ws_tools::get_ext_name( file_name );
	return( ext == "lsh" );
}

/**
	Get iterator to beginning of matrix region.
	@retval iter Iterator to beginning of matrix region
 */
iterator
Matrix::begin( )
{
	// get position in matrix data array that region starts at
	size_type region_start = (row_begin() * real_col_size()) + col_begin();

	iterator iter( &_matrix[0] + region_start,
		col_begin(), col_end(), real_col_size() );

	return( iter );

	// begin at starting point of region
	// return( _matrix.begin() + region_start );

	// iterate through entire matrix, not just region
	// return( _matrix.begin() );  // start of matrix data array
}

/**
	Get iterator to beginning of matrix region.
	@retval iter Iterator to beginning of matrix region
 */
const_iterator
Matrix::const_begin( ) const
{
	// get position in matrix data array that region starts at
	size_type region_start = (row_begin() * real_col_size()) + col_begin();

	const_iterator iter( &_matrix[0] + region_start,
		col_begin(), col_end(), real_col_size() );

	return( iter );

	// begin at starting point of region
	// return( &_matrix[0] + region_start );

	// iterate through entire matrix, not just region
	// return( &_matrix[0] );  // start of matrix data array
}

/**
	Get iterator to end of matrix region.
	@retval iter Iterator to end of matrix region
 */
iterator
Matrix::end( )
{
	// get position in matrix data array that region ends
	size_type region_end = ((row_end()) * real_col_size()) + col_begin();

	iterator iter( &_matrix[0] + region_end,
		col_begin(), col_end(), real_col_size() );

	return( iter );

	// get position in matrix data array that region ends
	// size_type region_end = ((row_end() - 1) * real_col_size()) + (col_end() - 1) ;
	// go one passed the end (we could add it in above, but this is clearer)
	// region_end += 1;

	// iterate through entire matrix, not just region
	// return( &_matrix[0] + region_end );

	// iterate through entire matrix, not just region
	// return( _matrix.end() );  // end of matrix data array
}

/**
	Get iterator to end of matrix region.
	@retval iter Iterator to end of matrix region
 */
const_iterator
Matrix::const_end( ) const
{
	// get position in matrix data array that region ends
	size_type region_end = ((row_end()) * real_col_size()) + col_begin();

	const_iterator iter( &_matrix[0] + region_end,
		col_begin(), col_end(), real_col_size() );

	return( iter );

	// get position in matrix data array that region ends
	// size_type region_end = ((row_end() - 1) * real_col_size()) + (col_end() - 1) ;
	// go one passed the end (we could add it in above, but this is clearer)
	// region_end += 1;

	// iterate through entire matrix, not just region
	// return( &_matrix[0] + region_end );

	// iterate through entire matrix, not just region
	// return( &_matrix[0] );  // end of matrix data array
}

/**
	Get constant iterator to beginning of matrix region.
	@retval iter Iterator to beginning of matrix region
const_iterator
Matrix::begin( ) const
{
	// get position in matrix data array that region starts at
	size_type region_start = (row_begin() * real_col_size()) + col_begin();

	// begin at starting point of region
	return( &_matrix[0] + region_start );

	// iterate through entire matrix, not just region
	// return( &_matrix[0] );  // start of matrix data array
}
 */

/**
	Get constant iterator to end of matrix region.
	@retval iter Iterator to end of matrix region
const_iterator
Matrix::end( ) const
{
	// get position in matrix data array that region ends
	size_type region_end = ((row_end() - 1) * real_col_size()) + (col_end() - 1) ;

	// go one passed the end (we could add it in above, but this is clearer)
	region_end += 1;

	return( &_matrix[0] + region_end );

	// iterate through entire matrix, not just region
	// return( &_matrix[0] );  // end of matrix data array
}
 */

/**
	Get iterator to start of column j in matrix.

	If the column is not inside the current active region, then an iterator
	to the end of the column is returned.

	@param[in] j Column index
	@retval col_iter Iterator to column
 */
col_iterator
Matrix::col_begin( size_type j )
{
	assert( j >= 0 && j < real_col_size() );

	// if outside of region, return iterator to end of column
	if( j < col_begin() || j >= col_end() )
	{
		return( col_end(j) );
	}

	// j_th column + start of row region * number of column positions per row
	col_iterator col_iter( &_matrix[0] + j + (row_begin() * real_col_size()),
		j, real_col_size() );

	// iterate through entire matrix, not just region
	// col_iterator col_iter( begin() + j, j, real_col_size() );

	return( col_iter );
}

/**
	Get iterator to start of column j in matrix.

	If the column is not inside the current active region, then an iterator
	to the end of the column is returned.

	@param[in] j Column index
	@retval col_iter Iterator to column
 */
const_col_iterator
Matrix::const_col_begin( size_type j ) const
{
	assert( j >= 0 && j < real_col_size() );

	// if outside of region, return iterator to end of column
	if( j < col_begin() || j >= col_end() )
	{
		return( const_col_end(j) );
	}

	// j_th column + start of row region * number of column positions per row
	const_col_iterator col_iter(
		&_matrix[0] + j + (row_begin() * real_col_size()), j, real_col_size() );

	// iterate through entire matrix, not just region
	// col_iterator col_iter( begin() + j, j, real_col_size() );

	return( col_iter );
}

/**
	Get iterator to end of column j in matrix.
	@param[in] j Column index
	@retval col_iter Iterator to column
 */
col_iterator
Matrix::col_end( size_type j )
{
	assert( j >= 0 && j < _col_size );

	// if outside of region, return iterator to end of matrix
	if( j < col_begin() || j >= col_end() )
	{
		// return( _matrix.end() );
	}

	// j_th column + end of row region * number of column positions per row
	col_iterator col_iter( &_matrix[0] + j + (row_end() * real_col_size()),
		j, real_col_size() );

	// end of j_th column + next offset position
	// col_iterator col_iter(
		// begin() + j + (real_row_size() * real_col_size()), j, real_col_size() );

	// alternate version calculates iterator position by subtracting from end
	// col_iterator col_iter( _matrix.end() - (real_col_size() - j - 1), j,
		// real_col_size() );

	return( col_iter );
}

/**
	Get iterator to end of column j in matrix.
	@param[in] j Column index
	@retval col_iter Iterator to column
 */
const_col_iterator
Matrix::const_col_end( size_type j ) const
{
	assert( j >= 0 && j < _col_size );

	// if outside of region, return iterator to end of matrix
	if( j < col_begin() || j >= col_end() )
	{
		// return( _matrix.end() );
	}

	// j_th column + end of row region * number of column positions per row
	const_col_iterator col_iter(
		&_matrix[0] + j + (row_end() * real_col_size()), j, real_col_size() );

	// end of j_th column + next offset position
	// col_iterator col_iter(
		// begin() + j + (real_row_size() * real_col_size()), j, real_col_size() );

	// alternate version calculates iterator position by subtracting from end
	// col_iterator col_iter( _matrix.end() - (real_col_size() - j - 1), j,
		// real_col_size() );

	return( col_iter );
}

/**
	Get iterator to start of row i in matrix.

	If the row is not inside the current active region, then an iterator
	to the end of the row is returned.

	@param[in] i Row index
	@retval row_iter Iterator to row
 */
row_iterator
Matrix::row_begin( size_type i )
{
	assert( i >= 0 && i < real_row_size() );

	// if outside of region, return iterator to end of row
	if( i < row_begin() || i >= row_end() )
	{
		return( row_end(i) );
	}

	// i_th row * number of column positions per row + start of column region
	row_iterator row_iter( &_matrix[0] + (i * real_col_size()) + col_begin(),
		i, 1 );

	// beginning of i_th row is i * number of column positions
	// row_iterator row_iter( begin() + (i * real_col_size()), i, 1 );

	return( row_iter );
}

/**
	Get iterator to start of row i in matrix.

	If the row is not inside the current active region, then an iterator
	to the end of the row is returned.

	@param[in] i Row index
	@retval row_iter Iterator to row
 */
const_row_iterator
Matrix::const_row_begin( size_type i ) const
{
	assert( i >= 0 && i < real_row_size() );

	// if outside of region, return iterator to end of row
	if( i < row_begin() || i >= row_end() )
	{
		return( const_row_end(i) );
	}

	// i_th row * number of column positions per row + start of column region
	const_row_iterator row_iter(
		&_matrix[0] + (i * real_col_size()) + col_begin(), i, 1 );

	// beginning of i_th row is i * number of column positions
	// row_iterator row_iter( begin() + (i * real_col_size()), i, 1 );

	return( row_iter );
}

/**
	Get iterator to end of row i in matrix.
	@param[in] i Row index
	@retval row_iter Iterator to row
 */
row_iterator
Matrix::row_end( size_type i )
{
	assert( i >= 0 && i < real_row_size() );

	// if outside of region, return iterator to end of matrix
	if( i < row_begin() || i >= row_end() )
	{
		// return( _matrix.end() );
	}

	// i_th row * number of column positions per row + end of column region
	row_iterator row_iter( &_matrix[0] + (i * real_col_size()) + col_end(),
		i, 1 );

	// end of i_th row is start of next row
	// row_iterator row_iter( begin() + ((i + 1) * real_col_size()), i, 1 );

	return( row_iter );
}

/**
	Get iterator to end of row i in matrix.
	@param[in] i Row index
	@retval row_iter Iterator to row
 */
const_row_iterator
Matrix::const_row_end( size_type i ) const
{
	assert( i >= 0 && i < real_row_size() );

	// if outside of region, return iterator to end of matrix
	if( i < row_begin() || i >= row_end() )
	{
		// return( _matrix.end() );
	}

	// i_th row * number of column positions per row + end of column region
	const_row_iterator row_iter(
		&_matrix[0] + (i * real_col_size()) + col_end(), i, 1 );

	// end of i_th row is start of next row
	// row_iterator row_iter( begin() + ((i + 1) * real_col_size()), i, 1 );

	return( row_iter );
}

/**
	Compare two values for use by qsort.
	@param[in] a Pointer to first value
	@param[in] b Pointer to second value
	@retval compare_result Integer less than, equal to, or greater than zero 
		if the first argument is considered to be respectively less than, equal
		to, or greater than the second
 */
int
qsort_compare( const void* a, const void* b )
{
	const value_type _a = *((value_type*) a);
	const value_type _b = *((value_type*) b);
	if( _a < _b )
	{
		return(-1);
	}
	else if( _a > _b )
	{
		return(1);
	}
	return(0);

	// this is cute, but rounding will cause errors:
	// 0 - .4 = -.4 => (int) -.4 = 0 => 0 == .4
	// return( static_cast<int>( a - b ) );
}

/**
	Sort values in matrix.

	All values are sorted, not just the current region.

	@retval this This matrix after sorting
 */
Matrix&
Matrix::sort( )
{
	qsort( data(), size_t(real_row_size() * real_col_size()), sizeof(value_type),
			qsort_compare );
	return( *this );
}

/**
	Compare two values for use by qsort to reverse the order of elements.
	@param[in] a Pointer to first value
	@param[in] b Pointer to second value
	@retval compare_result Integer less than, equal to, or greater than zero 
		if the first argument is considered to be respectively greater than, equal
		to, or less than the second
 */
int
qsort_reverse_compare( const void* a, const void* b )
{
	return( qsort_compare( b, a ) );

	/*
	const value_type _a = *((value_type*) a);
	const value_type _b = *((value_type*) b);
	if( _a > _b )
	{
		return(-1);
	}
	else if( _a < _b )
	{
		return(1);
	}
	return(0);
	 */
}

/**
	Sort values in matrix in descending order.

	All values are sorted, not just the current region.

	@retval this This matrix after sorting
 */
Matrix&
Matrix::reverse_sort( )
{
	qsort( data(), size_t(real_row_size() * real_col_size()), sizeof(value_type),
			qsort_reverse_compare );
	return( *this );
}

/**
	Get column j in matrix.
	@param[in] j Index of column to retrieve
	@retval column Column j in matrix
 */
Vector
Matrix::get_col( size_type j ) const
{
	Vector column( real_row_size() );
	for( unsigned i = 0; i != real_row_size(); ++i )
	{
		column(i) = at(i, j);
	}
	return( column );
}

/**
	Set column j in matrix with given column vector.
	@param[in] column Column vector
	@param[in] j Index of column to set
	@retval column Column j in matrix
 */
void
Matrix::set_col( Vector& column, size_type j )
{
	if( column.size() != real_row_size() )
	{
		err_warn( "set_col: Column's size (%d) does not match number of rows"
			" in matrix (%d): ignoring request to set new column\n",
			column.size(), real_row_size() );
		return;
	}

	// copy each column entry into j_th column of matrix
	for( unsigned i = 0; i != column.size(); ++i )
	{
		at(i,j) = column(i);
	}
}

/**
	Write contents of matrix to file pointer matrix_file_fp using specified
	precision for each value.
	@param[in] matrix_file_fp Pointer to file to write matrix to
	@param[in] precision Precision to use for floating-point values
 */
void Matrix::write( FILE* matrix_file_fp, prec_type precision ) const
{
	fprintf( matrix_file_fp, "[" );

	if( row_begin() == row_end() )
	{
		fprintf( matrix_file_fp, " ];\n" );
		return;
	}

	for( size_type i = row_begin(); i != row_end() - 1; i++ )
	{
		// write each column of row i
		for( size_type j = col_begin(); j != col_end(); j++ )
		{
			value_type val = at( i, j );
			if( val <= MAX_VALUE )
			{
				fprintf( matrix_file_fp, " %.*lf", precision, val );
			}
			else
			{
				fprintf( matrix_file_fp, " %.*lf", precision, MAX_VALUE );
			}
		}
		fprintf( matrix_file_fp, ";\n" );
	}

	// print last row without a ';' after it
	for( size_type j = col_begin(); j != col_end(); j++ )
	{
		value_type val = at( row_end() - 1, j );
		if( val <= MAX_VALUE )
		{
			fprintf( matrix_file_fp, " %.*lf", precision, val );
		}
		else
		{
			fprintf( matrix_file_fp, " %.*lf", precision, MAX_VALUE );
		}
	}
	fprintf( matrix_file_fp, " ];\n" );
}

/**
	Write contents of matrix to file file_name using specified precision for
	each value.
	@param[in] file_name Name of file to write matrix to
	@param[in] precision Precision to use for floating-point values
 */
void Matrix::write( const std::string& file_name, prec_type precision ) const
{
	// determine file type
	if( is_ascii_matlab(file_name) )
	{
		write_ascii_matlab( file_name, precision );
	}
	else if( is_binary_matlab(file_name) )
	{
		write_binary_matlab( file_name, precision );
	}
	else if( is_lush(file_name) )
	{
		write_lush( file_name, precision );
	}
	else // default to writing as ascii matlab if file type is undetermined
	{
		write_ascii_matlab( file_name, precision );
		// err_warn( "Mode is unsupported: refusing to write matrix\n" );
	}
}

/**
	Write matrix as an ascii MATLAB matrix to file matrix_file_name.
	@param[in] matrix_file_name Name of matrix file to write matrix to
	@param[in] precision Precision to use for floating-point values
 */
void Matrix::write_ascii_matlab( const std::string& matrix_file_name,
	prec_type precision ) const
{
	FILE* matrix_file_fp = stdout;

	if( matrix_file_name != "" )
	{
		if( (matrix_file_fp = fopen( matrix_file_name.c_str(), "wb" )) == NULL )
		{
			err_quit( "Unable to open matrix file '%s'\n", 
					matrix_file_name.c_str() );
		}
	}

	// erase file extension to get matrix name
	std::string matrix_name = matrix_file_name;
	matrix_name = ws_tools::get_base_name( ws_tools::get_file_name(
				matrix_file_name ) );

	if( matrix_name != "" )
	{
		// add "_" so MATLAB doesn't complain about the file and the matrix
		// having the same name
		fprintf( matrix_file_fp, "%s_ = ", matrix_name.c_str() );
	}

	write( matrix_file_fp, precision );

	if( matrix_file_name != "" )
	{
		fclose( matrix_file_fp );
	}
}

/**
	Write matrix as a binary MATLAB matrix to file matrix_file_name.
	@param[in] matrix_file_name Name of matrix file to write matrix to
	@param[in] precision Precision to use for floating-point values
 */
void Matrix::write_binary_matlab( const std::string& matrix_file_name,
	prec_type precision ) const
{
	FILE* matrix_file_fp = stdout;

	if( matrix_file_name != "" )
	{
		if( (matrix_file_fp = fopen( matrix_file_name.c_str(), "wb" )) == NULL )
		{
			err_quit( "Unable to open matrix file '%s'\n", 
					matrix_file_name.c_str() );
		}
	}

	// erase file extension to get matrix name
	std::string matrix_name = matrix_file_name;
	matrix_name =
		ws_tools::get_base_name( ws_tools::get_file_name( matrix_file_name ) );

	if( matrix_name != "" )
	{
		fprintf( matrix_file_fp, "%s = ", matrix_name.c_str() );
	}

	write( matrix_file_fp, precision );

	if( matrix_file_name != "" )
	{
		fclose( matrix_file_fp );
	}
}

/**
	Write matrix as a Lush matrix to file matrix_file_name.
	@param[in] matrix_file_name Name of matrix file to write matrix to
	@param[in] precision Precision to use for floating-point values
 */
void Matrix::write_lush( const std::string & matrix_file_name,
	prec_type precision ) const
{
	FILE* matrix_file_fp = stdout;
	
	if( matrix_file_name != "" )
	{
		if( (matrix_file_fp = fopen( matrix_file_name.c_str(), "wb" )) == NULL )
		{
			err_quit( "Unable to open matrix file '%s'\n", 
					matrix_file_name.c_str() );
		}
	}

	// erase file extension to get matrix name
	std::string matrix_name = matrix_file_name;
	std::string::size_type dot_pos = matrix_name.find_first_of( '.' );
	if( dot_pos != std::string::npos )
	{
		matrix_name.erase( dot_pos );
	}

	fprintf( matrix_file_fp, "(setq %s [\n", matrix_name.c_str() );

	if( row_begin() == row_end() )
	{
		fprintf( matrix_file_fp, "])\n" );
		return;
	}

	for( unsigned i = row_begin(); i != row_end(); i++ )
	{
		fprintf( matrix_file_fp, "[" );

		// write each column of row i
		for( unsigned j = col_begin(); j != col_end(); j++ )
		{
			value_type pixel = at( i, j );
			if( pixel <= MAX_VALUE )
			{
				fprintf( matrix_file_fp, " %.*lf", precision, pixel );
			}
			else
			{
				fprintf( matrix_file_fp, " %.*lf", precision, MAX_VALUE );
			}
		}
		fprintf( matrix_file_fp, " ]\n" );
	}
	fprintf( matrix_file_fp, "])\n" );

	if( matrix_file_name != "" )
	{
		fclose( matrix_file_fp );
	}
}

/**
	Read matrix as a MATLAB matrix from file matrix_file_name.
	@param[in] matrix_file_name Name of matrix file to read matrix from
 */
void
Matrix::read( const string & matrix_file_name )
{
	// open matrix file
	FILE* matrix_file_fp = NULL;
	if( (matrix_file_fp = fopen( matrix_file_name.c_str(), "r" )) == NULL )
	{
		err_quit( "Unable to open matrix file '%s'\n", matrix_file_name.c_str() );
	}

	// variables for parsing the matrix file
	int                c;          // character read from file
	string             value_str;  // single decimal value as a string
	vector<value_type> values;     // sequential list of values in the matrix

	// number of rows and columns in the matrix (updated while the file is read)
	size_type num_rows = 0, num_cols = 0;

	// read name of matrix, e.g., 'A = [' of 'A = [ 1 2; 3 4 ];'
	// Note: anything can be at the front of the matrix as long as '[' is found
	while( (c = fgetc( matrix_file_fp )) != EOF )
	{
		// stop when the start of the matrix is encountered
		if( c == '[' )
		{
			break;
		}
	}

	// check for file errors and that a matrix was found
	if( ferror( matrix_file_fp ) || c != '[' )
	{
		err_quit( "Error while reading file '%s': start of matrix is incorrect\n",
				matrix_file_name.c_str() );
	}

	// process each character at a time (we avoid processing strings at a time
	// since lines containing rows may be very long, which makes it harder to
	// handle buffering) 
	while( (c = fgetc( matrix_file_fp )) != EOF )
	{
		// if c is a character in a decimal number, such as -33.75, then add it
		// to the current value in the string
		if( isdigit(c) || c == '.' || c == '-' || c == '+' )
		{
			value_str += static_cast<string::value_type>( c );
		}
		else // c is not part of a decimal number
		{
			// if reached end of a value, convert it to decimal and add it to the
			// list of all values to be added to the matrix later
			if( value_str != "" )
			{
				values.push_back( ws_tools::string_to_double( value_str ) );
				value_str = "";
			}

			// if at the end of a matrix column
			if( c == ';' )
			{
				++num_rows;  // update the number of rows in the matrix

				// if we were in the first row, then we don't know the number of
				// columns in the matrix, so get this number; otherwise, make sure
				// the number of values in the current row matches the number of
				// columns we are expecting
				if( num_rows == 1 )
				{
					num_cols = values.size();
				}
				else if( values.size() != (num_rows * num_cols) )
				{
					err_quit( "Number of columns (%u) in row %u does not match the"
							" number of columns (%u) in row 1 in matrix file '%s'\n",
							num_cols + (values.size() - (num_rows * num_cols)),
							num_rows, num_cols, matrix_file_name.c_str() );
				}
			}

			// if at the end of matrix, quit reading the matrix file
			else if( c == ']' )
			{
				break;
			}

			// any other character besides whitespace is invalid
			else if( !isspace(c) )
			{
				err_quit( "Invalid character '%c' in matrix file '%s'\n",
						c, matrix_file_name.c_str() );
			}
		}
	}

	// check for file errors and close file (note that we should not have read
	// EOF since we stop reading when ']' is encountered
	if( ferror( matrix_file_fp ) )
	{
		err_quit( "Error while reading file '%s'\n", matrix_file_name.c_str() );
	}
	if( c == EOF )
	{
		err_quit( "Error while reading file '%s': no closing ']' found\n",
				matrix_file_name.c_str() );
	}
	if( matrix_file_name != "" )
	{
		fclose( matrix_file_fp );
	}

	// add one more row since the last row ends with something like "1.2 ];",
	// i.e., no semicolon follows the last decimal value in the last row
	++num_rows;

	// make sure the number of values in the current row matches the number of
	// columns we are expecting
	if( num_rows == 1 )
	{
		num_cols = values.size();
	}
	else if( values.size() != (num_rows * num_cols) )
	{
		err_quit( "Number of columns (%u) in row %u does not match "
				"the number of columns (%u) in row 1 in matrix file '%s'\n",
				num_cols + (values.size() - (num_rows * num_cols)),
				num_rows, num_cols, matrix_file_name.c_str() );
	}

	allocate_matrix( num_rows, num_cols );

	// copy each value into matrix
	vector<value_type>::size_type k = 0;
	for( size_type i = 0; i != num_rows; ++i )
	{
		for( size_type j = 0; j != num_cols; ++j )
		{
			at(i,j) = values[k++];
		}
	}
}

/**
	Set all elements of matrix to 0. The dimensions remain the same.
 */
void
Matrix::clear( )
{
	*this = 0;
}

/**
	Assign all elements of matrix to the given value. The dimensions remain the
	same.  
	@param[in] value Value to assign to each matrix element
 */
Matrix&
Matrix::operator=( const value_type value )
{
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			at(i,j) = value;
		}
	}
	return( *this );
}

/**
	Perform unary negation on matrix.
	@retval neg -A where *this = A
 */
Matrix Matrix::operator-( ) const
{
	Matrix neg( *this );
	iterator iter     = neg.begin();
	iterator end_iter = neg.end();

	// negate each value
	for( ; iter != end_iter; ++iter )
	{
		*iter = -(*iter);
	}

	return( neg );
}

/**
	Add two matrices.
	@param[in] rhs Matrix on right-hand side of +
	@retval this This matrix
 */
Matrix& Matrix::operator+=( const Matrix& rhs )
{
	// verify that region dimensions match
	if( rhs.sub_row_size() != sub_row_size()
		|| rhs.sub_col_size() != sub_col_size() )
	{
		err_quit( "Matrix::operator+: Matrix dimensions are not the same:"
		" (%d x %d) and (%d x %d)\n",
		sub_row_size(), sub_col_size(), rhs.sub_row_size(), rhs.sub_col_size() );
	}

	/*
	iterator iter                 = begin();
	iterator end_iter             = end();
	const_iterator const_rhs_iter = rhs.const_begin();

	// compute sum
	for( ; iter != end_iter; ++iter, ++const_rhs_iter )
	{
		assert( const_rhs_iter != rhs.end() );
		*iter += *const_rhs_iter;
	}
	 */

	// compute sum
	for( size_type i = row_begin(), i_rhs = rhs.row_begin();
			i != row_end();
			++i, ++i_rhs )
	{
		for( size_type j = col_begin(), j_rhs = rhs.col_begin();
				j != col_end();
				++j, ++j_rhs )
		{
			at(i, j) += rhs(i_rhs, j_rhs);
		}
	}

	return( *this );
}

/**
	Add two matrices.
	@param[in] rhs Matrix on right-hand side of +
	@retval matrix Matrix sum 
 */
Matrix Matrix::operator+( const Matrix& rhs ) const
{
	return( Matrix( *this ) += rhs );
}

/**
	Subtract two matrices.
	@param[in] rhs Matrix on right-hand side of -
	@retval this This matrix
 */
Matrix& Matrix::operator-=( const Matrix& rhs )
{
	// verify that region dimensions match
	if( rhs.sub_row_size() != sub_row_size()
		|| rhs.sub_col_size() != sub_col_size() )
	{
		err_quit( "Matrix::operator-: Matrix dimensions are not the same:"
		" (%d x %d) and (%d x %d)\n",
		sub_row_size(), sub_col_size(), rhs.sub_row_size(), rhs.sub_col_size() );
	}

	/*
	iterator iter                 = begin();
	iterator end_iter             = end();
	const_iterator const_rhs_iter = rhs.const_begin();

	// compute difference
	for( ; iter != end_iter; ++iter, ++const_rhs_iter )
	{
		assert( const_rhs_iter != rhs.end() );
		*iter -= *const_rhs_iter;
	}
	 */

	// compute difference
	for( size_type i = row_begin(), i_rhs = rhs.row_begin();
			i != row_end();
			++i, ++i_rhs )
	{
		for( size_type j = col_begin(), j_rhs = rhs.col_begin();
				j != col_end();
				++j, ++j_rhs )
		{
			at(i, j) -= rhs(i_rhs, j_rhs);
		}
	}

	return( *this );
}

/**
	Subtract two matrices.
	@param[in] rhs Matrix on right-hand side of -
	@retval matrix Matrix sum 
 */
Matrix Matrix::operator-( const Matrix& rhs ) const
{
	return( Matrix( *this ) -= rhs );
}

/*
	Multiply two matrices.
	@param[in] rhs Matrix on right-hand side of *
	@retval this This matrix
 */
Matrix& Matrix::operator*=( const Matrix& rhs )
{
	// verify that inner dimensions match for matrix product
	if( sub_col_size() != rhs.sub_row_size() )
	{
		err_quit( "Matrix::operator*: Inner dimensions are not the same: "
				"%ux%u with %ux%u\n", sub_row_size(), sub_col_size(),
				rhs.sub_row_size(), rhs.sub_col_size() );
	}

	Matrix result( sub_row_size(), rhs.sub_col_size() );

	// for each row and column
	for( size_type i = result.row_begin(); i != result.row_end(); ++i )
	{
		for( size_type j = result.col_begin(); j != result.col_end(); ++j )
		{
			// compute sum of products of row entries in this matrix and column
			// entries of rhs
			value_type sum = 0;
			for( size_type k = 0; k != col_end(); k++ )
			{
				sum += at(i, k) * rhs(k, j);
			}
			result(i, j) = sum;
		}
	}

	return( (*this = result) );
}

/**
	Multiply two matrices.
	@param[in] rhs Matrix on right-hand side of *
	@retval matrix Matrix product 
 */
Matrix Matrix::operator*( const Matrix& rhs ) const
{
	return( Matrix( *this ) *= rhs );
}

/**
	Add value to each component of matrix.
	@param[in] rhs Value on right-hand side of +
	@retval this This matrix
 */
Matrix& Matrix::operator+=( const value_type rhs )
{
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			at(i,j) += rhs;
		}
	}

	return( *this );
}

/**
	Add value to each component of matrix.
	@param[in] rhs Value on right-hand side of +
	@retval matrix Matrix result 
 */
Matrix Matrix::operator+( const value_type rhs ) const
{
	return( Matrix( *this ) += rhs );
}

/**
	Subtract value from each component of matrix.
	@param[in] rhs Value on right-hand side of -
	@retval this This matrix
 */
Matrix& Matrix::operator-=( const value_type rhs )
{
	*this += -rhs;
	return( *this );
}

/**
	Subtract value from each component of matrix.
	@param[in] rhs Value on right-hand side of -
	@retval matrix Matrix result 
 */
Matrix Matrix::operator-( const value_type rhs ) const
{
	return( Matrix( *this ) -= rhs );
}

/**
	Multiply matrix by value.
	@param[in] rhs Value on right-hand side of *
	@retval this This matrix
 */
Matrix& Matrix::operator*=( const value_type rhs )
{
	// compute product
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			at(i,j) *= rhs;
		}
	}

	return( *this );
}

/**
	Multiply matrix by value.
	@param[in] rhs Value on right-hand side of *
	@retval matrix Matrix result 
 */
Matrix Matrix::operator*( const value_type rhs ) const
{
	return( Matrix( *this ) *= rhs );
}

/**
	Divide matrix by value.
	@param[in] rhs Value on right-hand side of /
	@retval this This matrix
 */
Matrix& Matrix::operator/=( const value_type rhs )
{
	// division is equivalent to multiplication by reciprocal
	*this *= (1 / rhs);
	return( *this );
}

/**
	Divide matrix by value.
	@param[in] rhs Value on right-hand side of /
	@retval matrix Matrix result 
 */
Matrix Matrix::operator/( const value_type rhs ) const
{
	return( Matrix( *this ) /= rhs );
}

/**
	Transpose matrix--interchange rows and columns.
	@retval transposed_matrix This matrix transposed
 */
Matrix Matrix::transpose( ) const
{
	Matrix transposed_matrix( sub_col_size(), sub_row_size() );

	// copy ith row and jth col. of A into jth row and ith col. of transpose
	for( size_type i = row_begin(), i_t = 0; i != row_end(); ++i, ++i_t )
	{
		for( size_type j = col_begin(), j_t = 0; j != col_end(); ++j, ++j_t )
		{
			transposed_matrix(j_t, i_t) = at(i, j);
		}
	}

	return( transposed_matrix );
}

/**
	Compute average across each column in this matrix.

	Note: To compute the average across each row of a matrix M,
	call M.transpose().avg().

	@retval avg Average value of each column in the matrix
 */
Vector
Matrix::avg( ) const
{
	Vector avg, var; 
	stat( avg, var );
	return( avg );
}

/**
	Compute variance across each column in this matrix.

	Note: To compute the variance across each row of a matrix M,
	call M.transpose().var().

	@retval var Variance of each column in the matrix
 */
Vector
Matrix::var( ) const
{
	Vector avg, var; 
	stat( avg, var );
	return( var );
}

/**
	Compute statistics across each column in this matrix.

	Note: To compute the statistics across each row of a matrix M,
	call M.transpose().stat( avg, var ).

	@param[out] avg Average value of each column in the matrix
	@param[out] var Variance of each column in the matrix
 */
void
Matrix::stat( Vector& avg, Vector& var, Vector& std ) const
{
	stat( avg, var );

	if( std.sub_size() != sub_col_size() )
	{
		// std = Vector( sub_col_size(), 0 );
		std.resize( sub_col_size() );
	}

	for( size_type i = std.vec_begin(); i != std.vec_end(); ++i )
	{
		std(i) = sqrt( var(i) );
	}
}

/**
	Compute statistics across each column in this matrix.

	Note: To compute the statistics across each row of a matrix M,
	call M.transpose().stat( avg, var ).

	@param[out] avg Average value of each column in the matrix
	@param[out] var Variance of each column in the matrix
 */
void
Matrix::stat( Vector& avg, Vector& var ) const
{
	if( avg.sub_size() != sub_col_size() )
	{
		// avg = Vector( sub_col_size(), 0 );
		avg.resize( sub_col_size() );
		avg = 0;
	}
	if( var.sub_size() != sub_col_size() )
	{
		// var = Vector( sub_col_size(), 0 );
		var.resize( sub_col_size() );
		var = 0;
	}

	unsigned k = 0;
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			ws_tools::stat( at(i,j), k, &avg(j), &var(j) );
		}
		++k;
	}
	var /= sub_row_size() - 1;
}

/**
	Compute matrix statistics--average, variance, and standard deviation.

	Algorithm comes from:
	http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

	@param[out] avg Average value of all pixels in matrix
	@param[out] var Variance of all pixels in matrix
	@param[out] std Standard deviation of all pixels in matrix
 */
void
Matrix::stat( double* avg, double* var, double* std ) const
{
	// do nothing if pointers are invalid
	if( avg == 0 || var == 0 )  return;

	double a = *avg = 0;  // use non-pointers to avoid using '*' in code below
	double v = *var = 0;

	// Version 1:
	// algorithm which avoids large numbers in sum_sqr while summing up
	/*
	size_type n = 0;
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			a = ( (a * n) + at( i, j ) ) / (n + 1);
			v = ( (v * n) + (at(i, j) - a) * (at(i, j) - a) ) / (n + 1);
			++n;
		}
	}
	*/

	// Version 2:
	// A similar algorithm, slightly more efficient, which is also numerically
	// stable. It also computes the mean. This algorithm is due to Knuth (The
	// Art of Computer Programming: Semi-Numerical Algorithms), who cites
	// Welford, in Technometrics v. 4, 1962, pp 419-420.  
	size_type n = 0;
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			double delta = at(i,j) - a;
			a += delta / (n + 1);
			v += delta * (at(i,j) - a);
			++n;
		}
	}
	v /= n - 1;

	*avg = a;
	*var = v;
	if( std != 0 )
	{
		*std = sqrt( v );
	}

/*
	// Version 3: alternate (slower) method that computes average and variance
	// separately compute average first
	avg = 0;
	for( unsigned i = 0; i != real_row_size(); i++ )
	{
		for( unsigned j = 0; j != real_col_size(); j++ )
		{
			*avg += at( i, j );
		}
	}
	*avg /= (real_row_size() * real_col_size());

	// compute variance second using average
	var = 0;
	for( unsigned i = 0; i != real_row_size(); i++ )
	{
		for( unsigned j = 0; j != real_col_size(); j++ )
		{
			double val = at( i, j ) - avg;
			*var += (val * val);
		}
	}
	*var /= (real_row_size() * real_col_size()) - 1;

	*std = sqrt( v );
 */
}

/**
	Calculate kurtosis of each column in this matrix.

	The computation follows the formula given here:
		http://en.wikipedia.org/wiki/Kurtosis#Sample_kurtosis

	@retval kurt Kurtosis
 */
Vector
Matrix::kurtosis( ) const
{
	// calculate mean down each column (the variance is unused)
	Vector avg;
	Vector var;
	stat( avg, var );

	// v2 is the variance squared about the mean
	Matrix v2 = *this;
	for( size_type i = v2.row_begin(); i != v2.row_end(); ++i )
	{
		Vector::size_type jj = avg.vec_begin();
		for( size_type j = v2.col_begin(); j != v2.col_end(); ++j, ++jj )
		{
			v2(i,j) -= avg(jj);  // remove mean
			v2(i,j) *= v2(i,j);
		}
	}

	// m4 is 4th moment about the mean (i.e., to the 4th power,
	// so just square v2)
	Matrix m4 = v2;
	for( size_type i = m4.row_begin(); i != m4.row_end(); ++i )
	{
		for( size_type j = m4.col_begin(); j != m4.col_end(); ++j )
		{
			m4(i,j) *= m4(i,j);
		}
	}

	Vector v2_avg;
	Vector v2_var;
	v2.stat( v2_avg, v2_var );

	Vector m4_avg;
	Vector m4_var;
	m4.stat( m4_avg, m4_var );

	/*
fprintf( stderr, "\n" );
v2.write();
fprintf( stderr, "\n" );
m4.write();

fprintf( stderr, "\n" );
v2_avg.write();

fprintf( stderr, "\n" );
m4_avg.write();
	 */

	Vector kurt( sub_col_size() );
	Vector::size_type v2_i = v2_avg.vec_begin();
	Vector::size_type m4_i = m4_avg.vec_begin();
	for( Vector::size_type i = kurt.vec_begin(); i != kurt.vec_end();
			++i, ++v2_i, ++m4_i )
	{
		kurt(i) = m4_avg(m4_i) / (v2_avg(v2_i) * v2_avg(v2_i));
	}

	return( kurt );
}

/**
	Compute sum of each column in this matrix.

	Note: To compute the sum across all elements of a matrix, call M.sum().sum()
	since the Vector class also implements sum().

	@retval col_sum Vector of column sums
 */
Vector
Matrix::sum( ) const
{
	Vector col_sum( sub_col_size() );
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			col_sum(j) += at(i,j);
		}
	}

	return( col_sum );
}

/**
	Compute dot product of two vectors.
	@param[in] rhs Matrix on right-hand side of dot product
	@retval dot_product Dot product of vectors
 */
value_type
Matrix::dot_product( const Matrix& rhs ) const
{
	value_type dot_product = 0;

	// verify matrices are vectors of same size and nonempty--vectors can either
	// be column or row vectors
	if( (sub_row_size() != rhs.sub_row_size()
		|| sub_col_size() != 1 || rhs.sub_col_size() != 1)
	&& (sub_col_size() != rhs.sub_col_size()
		|| sub_row_size() != 1 || rhs.sub_row_size() != 1) )
	{
		return( dot_product );
	}

	if( sub_col_size() == 1 )
	{
		const_col_iterator iter     = const_col_begin( col_begin() );
		const_col_iterator end_iter = const_col_end( col_begin() );
		const_col_iterator rhs_iter = rhs.const_col_begin( rhs.col_begin() );

		// sum the products of the components in each vector
		for( ; iter != end_iter; ++iter, ++rhs_iter )
		{
			dot_product += (*iter) * (*rhs_iter);
		}
	}
	else if( sub_row_size() == 1 )
	{
		const_row_iterator iter     = const_row_begin( row_begin() );
		const_row_iterator end_iter = const_row_end( row_begin() );
		const_row_iterator rhs_iter = rhs.const_row_begin( rhs.row_begin() );

		// sum the products of the components in each vector
		for( ; iter != end_iter; ++iter, ++rhs_iter )
		{
			dot_product += (*iter) * (*rhs_iter);
		}
	}

	return( dot_product );
}

/**
	Compute 2-norm of vector, i.e.,
		|| x || = sqrt( x_1^2, x_2^2, ..., x_n^2 ).
	@param[in] use_fast_method Whether to use the faster, but possibly less
		accurate method (true), or the slower, but possibly more accurate method
		(false), when computing the norm 
	@retval norm 2-norm of vector
 */
value_type Matrix::two_norm( bool use_fast_method ) const
{
	value_type norm = 0;

	// if( use_fast_method )
	{
		// Version 1: directly compute using sum of two_norm_2()
		norm = sqrt( two_norm_2() );
	}
	// else // Version 2: cblas version
	{
		/*
		// verify that matrix is a row or column vector
		if( !(sub_row_size() == 1 || sub_col_size() == 1) )
		{
			return( 0 );
		}

		// Version 2: compute norm using more stable method defined in cblas's
		// dnrm2--slower, but possibly more numerically stable for small values
		value_type scale = 0;
		value_type ssq   = 1;
		for( size_type i = row_begin(); i != row_end(); ++i )
		{
			for( size_type j = col_begin(); j != col_end(); ++j )
			{
				if( at(i,j) != 0 )
				{
					double abs_ij = fabs( at(i,j) );
					if( scale < abs_ij )
					{
						value_type frac = (scale / abs_ij);
						ssq = 1 + ssq * (frac * frac);
						scale = abs_ij;
					}
					else
					{
						value_type frac = (abs_ij / scale);
						ssq += (frac * frac);
					}
				}
			}
		}
		norm = scale * sqrt( ssq );
		 */
	}

	return( norm );
}

/**
	Compute the 2-norm squared of this matrix (taken to be a vector), i.e.,
		|| x ||_2^2 = x_1^2, x_2^2, ..., x_n^2.
	@retval norm 2-norm of vector
 */
value_type Matrix::two_norm_2( bool use_fast_method ) const
{
	value_type norm = 0;

	// Version 1: directly compute norm using literal definition--faster, but
	// possibly numerically unstable for small values
	// if( use_fast_method )
	{
		// verify that matrix is a row or column vector
		if( !(sub_row_size() == 1 || sub_col_size() == 1) )
		{
			return( 0 );
		}

		for( size_type i = row_begin(); i != row_end(); ++i )
		{
			for( size_type j = col_begin(); j != col_end(); ++j )
			{
				norm += at(i,j) * at(i,j);
			}
		}
	}
	// else
	{
		// Version 2: cblas version that also uses cblas version of two_norm()
		// value_type norm_2 = two_norm( false );
		// norm = norm_2 * norm_2;
	}

	return( norm );
}

/**
	Compute infinity norm of matrix, i.e.,
		|| A ||_inf = max{ |a_1,1| + |a_1,2| + ..., 
								 |a_2,1| + |a_2,2| + ...,
								 ..., 
								 |a_n,1| + |a_n,2| + ... }.

	@retval norm Infinity norm of matrix/vector
 */
value_type Matrix::inf_norm( )
{
	value_type max_value = 0;
	for( size_type i = row_begin(); i != row_end(); i++ )
	{
		value_type row_sum = 0;
		for( size_type j = col_begin(); j != col_end(); j++ )
		{
			row_sum += fabs( at(i, j) );
		}

		// compute new max
		max_value = max<value_type>( row_sum, max_value );
	}

	return( max_value );
}

/**
	Compute the sine of each element in this matrix.

	@retval A Result matrix
 */
Matrix
Matrix::sin( ) const
{
	Matrix A( *this );
	for( size_type i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( size_type j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = static_cast<value_type>( std::sin( static_cast<double>(
							at(i,j)
						) ) );
		}
	}
	return( A );
}

/**
	Compute the cosine of each element in this matrix.

	@retval A Result matrix
 */
Matrix
Matrix::cos( ) const
{
	Matrix A( *this );
	for( size_type i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( size_type j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = static_cast<value_type>( std::cos( static_cast<double>(
							at(i,j)
						) ) );
		}
	}
	return( A );
}

/**
	Compute the tangent of each element in this matrix.

	@retval A Result matrix
 */
Matrix
Matrix::tan( ) const
{
	Matrix A( *this );
	for( size_type i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( size_type j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = static_cast<value_type>( std::tan( static_cast<double>(
							at(i,j)
						) ) );
		}
	}
	return( A );
}

/**
	Normalize matrix--make it a unit matrix.
 */
void
Matrix::normalize( )
{
   // divide each element of matrix by its 2-norm
	*this /= two_norm();
}

/**
	Compute absolute value of matrix.
	TODO shouldn't we iterate across A and this with different indices just to
	be safe?
	TODO should we return a new matrix or just compute abs() of this matrix?
 */
Matrix
Matrix::abs( ) const
{
	Matrix A( *this );
	for( size_type i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( size_type j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = fabs( at(i,j) );
		}
	}
	return( A );
}

/**
	Create identity matrix with given dimensions.
	@param[in] dim Dimension of matrix
	@param[in] I Identity matrix
 */
Matrix Matrix::identity( size_type dim )
{
	Matrix I( dim, dim, 0 );
	for( size_type i = I.row_begin(); i != I.row_end(); ++i )
	{
		I(i,i) = 1;
	}
	return( I );
}

/**
	Create NxN matrix initialized to all 0s.
	@param[in] N Dimension of matrix
	@retval A New matrix
 */
Matrix
Matrix::zeros( size_type N )
{
	return( Matrix( N, N, 0 ) );
}

/**
	Create matrix with given dimensions initialized to all 0s.
	@param[in] rows Number of rows in matrix
	@param[in] cols Number of columns in matrix
	@retval A New matrix
 */
Matrix
Matrix::zeros( size_type rows, size_type cols )
{
	return( Matrix( rows, cols, 0 ) );
}

/**
	Create NxN matrix initialized to all 1s.
	@param[in] N Dimension of matrix
	@retval A New matrix
 */
Matrix
Matrix::ones( size_type N )
{
	return( Matrix( N, N, 1 ) );
}

/**
	Create matrix with given dimensions initialized to all 1s.
	@param[in] rows Number of rows in matrix
	@param[in] cols Number of columns in matrix
	@retval A New matrix
 */
Matrix
Matrix::ones( size_type rows, size_type cols )
{
	return( Matrix( rows, cols, 1 ) );
}

/**
	Create random matrix with given dimensions.
	@param[in] rows Number of rows in matrix
	@param[in] cols Number of columns in matrix
	@param[in] max_value Maximum value in range

	@retval A Random matrix
 */
Matrix Matrix::random( size_type rows, size_type cols, value_type max_value )
{
	Matrix A( rows, cols );

	// seed the random number generator and save it since we use the reentrant
	// rand_r() function
	// unsigned seed = time( NULL );
	// srand( seed );
	ws_tools::Uniform_Number un( 0, max_value );

	iterator iter     = A.begin();
	iterator end_iter = A.end();

	// generate random values within the given range for each element
	for( ; iter != end_iter; ++iter )
	{
		// 1. Uses lower order bits, which is not very good
		// *iter = rand_r( &seed ) % int(max_value);

		// 2. Uses higher order bits, which is better
		// *iter = max_value * (rand_r( &seed ) / (RAND_MAX + 1.0));

		// 3. Uses random number object
		*iter = un();
	}
	return( A );
}

/**
	Create random matrix with given dimensions.

	@param[in] rows Number of rows in matrix
	@param[in] cols Number of columns in matrix
	@param[in] min_value Minimum value in range
	@param[in] max_value Maximum value in range

	@retval A Random matrix
 */
Matrix Matrix::random( size_type rows, size_type cols,
		value_type min_value, value_type max_value )
{
	Matrix A( rows, cols );

	iterator iter     = A.begin();
	iterator end_iter = A.end();

	// generate random values within the given range for each element
	ws_tools::Uniform_Number un( min_value, max_value );
	for( ; iter != end_iter; ++iter )
	{
		*iter = un();
	}
	return( A );
}

/**
	Create random matrix with given dimensions whose values are normally
	distributed.

	@param[in] rows Number of rows in matrix
	@param[in] cols Number of columns in matrix
	@param[in] mean Mean of normal distribution (default is 0)
	@param[in] std Standard deviation of normal distribution (default is 1)
	@param[in] A Random matrix
 */
Matrix Matrix::random_normal( size_type rows, size_type cols,
		double mean, double std )
{
	Matrix A( rows, cols );

	ws_tools::Normal_Number normal_number( mean, std );

	// add Gaussian noise to each pixel
	for( size_type i = A.row_begin(); i != A.row_end(); ++i )
	{
		for( size_type j = A.col_begin(); j != A.col_end(); ++j )
		{
			A(i,j) = normal_number();
		}
	}

	return( A );
}

/**
	Construct symmetric toeplitz matrix.

	@param[in] c
	@retval A Toeplitz matrix
 */
Matrix
Matrix::toeplitz( const Vector& c )
{
	return( toeplitz(c, c) );
}

/**
	Construct toeplitz matrix.

	The first element of r is ignored since it sets the main diagonal, but c
	already does.  

	@param[in] c
	@param[in] r
	@retval A Toeplitz matrix
 */
Matrix
Matrix::toeplitz( const Vector& c, const Vector& r )
{
	if( c.sub_size() == 0 )
	{
		err_quit( "Matrix::toeplitz: Column vector c has invalid size %u\n",
				c.sub_size() );
	}
	if( r.sub_size() == 0 )
	{
		err_quit( "Matrix::toeplitz: Row vector r has invalid size %u\n",
				r.sub_size() );
	}

	Matrix A( c.sub_size(), r.sub_size() );

	// set lower triangular region
	// for each diagonal
	for( size_type ci = c.vec_begin(), row_start = A.row_begin();
			ci != c.vec_end();
			++ci, ++row_start )
	{
		// set each diagonal
		for( size_type i = row_start, j = A.col_begin();
				i != A.row_end() && j != A.col_end();
				++i, ++j )
		{
			assert( i < A.real_row_size() && j < A.real_col_size()
					&& ci < c.size() );

			A(i,j) = c(ci);
		}
	}

	// we skip the first element of r below, and the code below assumes we have a
	// second element, which may not be the case
	if( r.sub_size() <= 1 )
	{
		return( A );
	}

	// set upper triangular region
	// for each diagonal--skip first element
	for( size_type ri = (r.vec_begin() + 1), col_start = (A.col_begin() + 1);
			ri != r.vec_end();
			++ri, ++col_start )
	{
		// set each diagonal
		for( size_type i = A.row_begin(), j = col_start;
				i != A.row_end() && j != A.col_end();
				++i, ++j )
		{
			assert( i < A.real_row_size() && j < A.real_col_size()
					&& ci < c.size() );

			A(i,j) = r(ri);
		}
	}

	return( A );
}

/**
	Multiply all the matrices in vector. This is more efficient than simply
	multipling with separate calls to operator* as the dimensions of each matrix
	are examined to minimize the number of operations.

	@param[in] matrices Matrices to multiply
	@retval product Matrix product
 */
Matrix Matrix::multiply( std::vector< Matrix >& matrices )
{
	// matrix p
	Matrix product;

	// compute best order to perform multiplication
	Matrix::Vector_Int dim   = initialize_dim( matrices );
	Matrix::Matrix_Int order = best_mult_order( matrices, dim );

	// print expression that will be evaluated
	// Note: need to start in upper right corner of order matrix to work back
	std::string expression = make_exp( order, 0, matrices.size() - 1 );
	std::cout << expression << std::endl;

	// perform multiplication
	// Note: need to start in upper right corner of order matrix to work back
	return( multiply( matrices, order, 0, matrices.size() - 1 ) );
}

/**
	Initialize array holding dimensions of matrices.
	If n = # of matrices, then dim[i] = # of rows( matrix i ) for i = 0 to n-1
	and dim[n] = # of cols( matrix n-1 )

	@param[in] matrices Matrices to multiply
	@retval dim Matrix initialized with matrices' dimensions
 */
Matrix::Vector_Int
Matrix::initialize_dim( vector< Matrix >& matrices )
{
	// verify that have at least one matrix
	if( matrices.empty() )
	{
		Vector_Int empty; 
		return( empty );
	}

	// storing all row numbers + last column number
	Matrix::Vector_Int dim( matrices.size() + 1 );

	// dim[ i ] = # of rows of matrix i
	dim[ 0 ] = matrices[ 0 ].real_row_size();
	for( unsigned i = 1; i != dim.size() - 1; i++ )
	{
		dim[ i ] = matrices[ i ].real_row_size();

		// verify that # of rows of this matrix equals # of columns of previous
		if( unsigned( dim[ i ] ) != matrices[ i - 1 ].real_col_size() )
		{
			err_quit( "Incompatible dimensions for matrices %u and %u\n",
				i - 1, i );
		}
	}

	// last element is number of columns of last matrix
	dim[ dim.size() - 1 ] = matrices[ dim.size() - 2 ].real_col_size();

	return( dim );
}

/**
	Compute optimal order for multiplication.

	@param[in] matrices Matrices to multiply
	@param[in] dim Matrix initialized with matrices' dimensions
 */
Matrix::Matrix_Int 
Matrix::best_mult_order( std::vector< Matrix >& matrices,
	Matrix::Vector_Int& dim )
{
	// record number of matrices since it will be used heavily (and shorter)
	size_type n = matrices.size();

	if( n == 0 )
	{
		Matrix_Int empty; 
		return( empty );
	}

	Matrix::Matrix_Int order( n, Matrix::Vector_Int( n, 0 ) );

	// place to store counts--initialized to 0
	Matrix::Matrix_Int M( n, Matrix::Vector_Int( n, 0 ) );

	for( int q = n - 2; q >= 0; q-- )  // q must be of type int
	{
		int i = 0;
		int j = n - q - 1;
		for( ; i <= q; i++, j++ )
		{
			/*
			M_i,j = min{ M_i,k + M_k+1,j + d_i * d_k+1 * d_j+1 } for i <= k < j
			 */

			/*
				maximum possible value for first min() comparison
				Note max might need changing depending on M's value type, e.g.,
				long, int, unsigned, etc.
			 */
			//M[ i ][ j ] = ULONG_MAX;
			M[ i ][ j ] = INT_MAX;

			// work along diagonal offsets
			for( int k = i; k < j; k++ )
			{
				// store old minimum for comparison later
				size_type saved_min = M[ i ][ j ];

				// compute number of operations required
				M[ i ][ j ] = min<int>( M[ i ][ j ], 
					M[ i ][ k ] + M[ k + 1 ][ j ]
					+ dim[ i ] * dim[ k + 1 ] * dim[ j + 1 ] );

				// if have new min, save it
				if( M[ i ][ j ] != int(saved_min) )
				{
					order[ i ][ j ] = k;
				}
			}
		}
	}

	std::cout << "Number of multiplications: " << M[ 0 ][ n - 1 ] 
		<< std::endl;

	return( order );
}

/**
	Create expression string to represent best order to multiply matrices.
	(i,k)*(k+1,j)

	@param[in] order
	@param[in] i
	@param[in] j
 */
std::string 
Matrix::make_exp( Matrix::Matrix_Int& order, int i, int j )
{
	std::string ret;

	int k = order[ i ][ j ];

	if( i == j )
	{
		ret.append( ws_tools::int_to_string(i) );
	}
	else if( k == i )
	{
		ret.append( ws_tools::int_to_string(i) );
		ret.append( " * ( " );
		ret.append( make_exp( order, i + 1, j ) );
		ret.append( " )" );
	}
	else if( k + 1 == j )
	{
		ret.append( "( " );
		ret.append( make_exp( order, i, k ) );
		ret.append( " ) * " );
		ret.append( ws_tools::int_to_string(j) );
	}
	else // (i,k)*(k+1,j)
	{
		ret.append( "( " );
		ret.append( make_exp( order, i, k ) );
		ret.append( " ) * ( " );
		ret.append( make_exp( order, k + 1, j ) );
		ret.append( " )" );
	}

	return( ret );
}

/**
	Multiply matrices using most efficient order.
	(i,k)*(k+1,j)
 */
Matrix Matrix::multiply( std::vector< Matrix >& matrices, 
	Matrix::Matrix_Int& order, int i, int j )
{
	int k = order[ i ][ j ];
	if( i == j )
	{
		return( matrices[ i ] );
	}
	else if( k == i )
	{
		Matrix result = multiply( matrices, order, i + 1, j );
		return( matrices[ i ] * result );
	}
	else if( k + 1 == j )
	{
		Matrix result = multiply( matrices, order, i, k );
		return( result * matrices[ j ] );
	}
	else // (i,k)*(k+1,j)
	{
		Matrix result1 = multiply( matrices, order, i, k );
		Matrix result2 = multiply( matrices, order, k + 1, j );

		return( result1 * result2 );
	}
}

/**
	Convert matrix to vector.
Vector Matrix::to_vector( ) const
{
	assert( real_col_size() == 1 || real_row_size() == 1 );

	Vector v;
	if( real_col_size() == 1 )
	{
		v = Vector( real_row_size() );
		for( size_type i = row_begin(); i != row_end(); ++i )
		{
			v[i] = at(i,0);
		}
	}
	else if( real_row_size() == 1 )
	{
		v = Vector( real_col_size() );
		for( size_type j = col_begin(); i != col_end(); ++i )
		{
			v[i] = at(0,i);
		}
	}
	return( v );
}
 */

/**
	Convert matrix to an array (allocated with 'new').
	The returned array must be dellocated by a call to 'delete[]'.
	@retval array
 */
double* Matrix::to_array( ) const
{
	double* array = new double[ unsigned(sub_row_size() * sub_col_size()) ];

	unsigned current_pos = 0;
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			array[current_pos++] = at(i, j);
		}
	}

	return( array );
}

/**
	Find minimum value in matrix.
	@retval max Minimum value in matrix.
 */
value_type
Matrix::get_min( ) const
{
	value_type min, max;
	get_min_max( &min, &max );
	return( min );
}

/**
	Find maximum value in matrix.
	@retval max Maximum value in matrix.
 */
value_type
Matrix::get_max( ) const
{
	value_type min, max;
	get_min_max( &min, &max );
	return( max );
}

/**
	Find minimum and maximum value in matrix.
	@retval max Maximum value in matrix.
 */
void
Matrix::get_min_max( value_type* min, value_type* max ) const
{
	if( min == 0 || max == 0 )
	{
		return;
	}

	*min =  MAX_VALUE;
	*max = -MAX_VALUE;
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			value_type pixel = at(i,j);
			if( pixel < *min )
			{
				*min = pixel;
			}
			if( pixel > *max )
			{
				*max = pixel;
			}
		}
	}
}

/**
	Convolve matrix's rows with given filter.

	The code is very similar to the Vector class's convolve() except only the
	portion of the.

	@param[in] row_filter Filter to convolve along each row with
	@param[out] out Output matrix
	@retval out Output matrix
 */
Matrix&
Matrix::convolve_row( const Vector& filter, Matrix& out )
{
	if( &out == this )
	{
		return( convolve_row( filter ) );
	}
	if( filter.sub_size() > sub_col_size() )
	{
		return( *this );
	}

	out.resize( sub_row_size(), filter.sub_size() + sub_col_size() - 1 );
	out = 0;

	// determine center of filter kernel
	const size_type filter_center = filter.sub_size() / 2;

	// filter along each row
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		// y = filter * h where * is the convolution operator
		// Vector y( filter.sub_size() + sub_col_size() - 1 );

		// first section of convolution: not all of row is used
		for( size_type n = filter_center; n != filter.sub_size() - 1; ++n )
		{
			for( size_type k = filter.vec_begin(); k != (n + 1); ++k )
			{
				out(i, n) += filter(k) * at(i, col_begin() + n - k);
			}
		}

		// middle: all of filter and row are used
		for( size_type n = filter.sub_size() - 1; n != sub_col_size(); ++n )
		{
			for( size_type k = filter.vec_begin(); k != filter.sub_size(); ++k )
			{
				out(i, n) += filter(k) * at(i, col_begin() + n - k);
			}
		}

		// last section of convolution: not all of filter is used
		for( size_type n = sub_col_size(), m = filter.vec_begin();
				m != filter_center; ++n )
		{
			++m;  // shift start of filter one position to the right

			for( size_type k = m; k != filter.sub_size(); ++k )
			{
				out(i, n) += filter(k) * at(i, col_begin() + n - k);
			}
		}
	}

	out.set_region( 0, 0, sub_row_size(), sub_col_size() );

	return( out );
}

/**
	Convolve matrix's rows with given filter.

	The code is very similar to the Vector class's convolve() except only the
	portion of the.

	@param[in] row_filter Filter to convolve along each row with
 */
Matrix&
Matrix::convolve_row( const Vector& filter )
{
	if( filter.sub_size() > sub_col_size() )
	{
		return( *this );
	}

	// determine center of filter kernel
	const size_type filter_center = filter.sub_size() / 2;

	// TODO Make filters work with regions: do we need to add filter(
	// filter.begin() + k) like for Vector?  

	// filter along each row
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		// y = filter * h where * is the convolution operator
		Vector y( filter.sub_size() + sub_col_size() - 1, 0 );

		// first section of convolution: not all of row is used
		for( size_type n = filter_center; n != filter.sub_size() - 1; ++n )
		{
			for( size_type k = filter.vec_begin(); k != (n + 1); ++k )
			// for( size_type k = 0; k != (n + 1); ++k )
			{
				y(n) += filter(k) * at(i, col_begin() + n - k);
				// y(n) += filter(filter.vec_begin() + k) * at(i, col_begin() + n - k);
			}
		}

		// middle: all of filter and row are used
		for( size_type n = filter.sub_size() - 1; n != sub_col_size(); ++n )
		{
			for( size_type k = filter.vec_begin(); k != filter.sub_size(); ++k )
			// for( size_type k = 0; k != filter.sub_size(); ++k )
			{
				y(n) += filter(k) * at(i, col_begin() + n - k);
				// y(n) += filter(filter.vec_begin() + k) * at(i, col_begin() + n - k);
			}
		}

		// last section of convolution: not all of filter is used
		for( size_type n = sub_col_size(), m = filter.vec_begin();
				m != filter_center; ++n )
		{
			++m;  // shift start of filter one position to the right

			for( size_type k = m; k != filter.sub_size(); ++k )
			{
				y(n) += filter(k) * at(i, col_begin() + n - k);
			}
		}

		/*
		// last section of convolution: not all of filter is used
		for( size_type n = sub_col_size(), m = 0; n != y.sub_size(); ++n )
		{
			++m;  // shift start of filter one position to the right

			for( size_type k = m; k != filter.sub_size(); ++k )
			{
				y(n) += filter(filter.vec_begin() + k) * at(i, col_begin() + n - k);
			}
		}
		 */

		// copy convolution result into this row (but only the section that lied
		// within the row originally)
		for( size_type j = col_begin(), j2 = filter_center;
				j != col_end();
				++j, ++j2 )
		{
			at(i,j) = y(j2);
		}
	}

	/*
	// version 2: rely on Vector class's convolve()

	// determine center of filter kernel
	const size_type filter_center = filter.size() / 2;

	// filter along each row
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		// result from convolution along a single row
		Vector row_conv( sub_col_size() );

		for( size_type j = col_begin(), j2 = row_conv.vec_begin();
				j != col_end();
				++j, ++j2 )
		{
			row_conv(j2) = at(i,j);
		}

		row_conv = row_conv.convolve( filter );

		for( size_type j = col_begin(), j2 = row_conv.vec_begin() + filter_center;
				j != col_end();
				++j, ++j2 )
		{
			at(i,j) = row_conv(j2);
		}
	}
	 */
	return( *this );
}

/**
	Convolve matrix's columns with given filter.

	The code is very similar to the Vector class's convolve().

	@param[in] row_filter Filter to convolve along each row with
	@param[out] out Output matrix
	@retval out Output matrix
 */
Matrix&
Matrix::convolve_col( const Vector& filter, Matrix& out )
{
	if( &out == this )
	{
		return( convolve_row( filter ) );
	}
	if( filter.sub_size() > sub_row_size() )
	{
		return( *this );
	}

	out.resize( filter.sub_size() + sub_row_size() - 1, sub_col_size() );
	out = 0;

	// determine center of filter kernel
	const size_type filter_center = filter.sub_size() / 2;

	// filter along each row
	for( size_type j = col_begin(); j != col_end(); ++j )
	{
		// y = filter * h where * is the convolution operator
		// Vector y( filter.sub_size() + sub_row_size() - 1 );

		// first section of convolution: not all of row is used
		for( size_type n = filter_center; n != filter.sub_size() - 1; ++n )
		{
			for( size_type k = filter.vec_begin(); k != (n + 1); ++k )
			{
				out(n,j) += filter(k) * at(row_begin() + n - k, j);
			}
		}

		// middle: all of filter and row are used
		for( size_type n = filter.sub_size() - 1; n != sub_row_size(); ++n )
		{
			for( size_type k = filter.vec_begin(); k != filter.sub_size(); ++k )
			{
				out(n,j) += filter(k) * at(row_begin() + n - k, j);
			}
		}

		// last section of convolution: not all of filter is used
		for( size_type n = sub_row_size(), m = filter.vec_begin();
				m != filter_center; ++n )
		{
			++m;  // shift start of filter one position to the right

			for( size_type k = m; k != filter.sub_size(); ++k )
			{
				out(n,j) += filter(k) * at(row_begin() + n - k, j);
			}
		}
	}

	out.set_region( 0, 0, sub_row_size(), sub_col_size() );

	return( out );
}

/**
	Convolve matrix's columns with given filter.

	The code is very similar to the Vector class's convolve().

	@param[in] filter Filter to convolve along each column with
 */
Matrix&
Matrix::convolve_col( const Vector& filter )
{
	if( filter.sub_size() > sub_row_size() )
	{
		return( *this );
	}

	// determine center of filter kernel
	const size_type filter_center = filter.sub_size() / 2;

	// filter along each row
	for( size_type j = col_begin(); j != col_end(); ++j )
	{
		// y = filter * h where * is the convolution operator
		Vector y( filter.sub_size() + sub_row_size() - 1, 0 );

		// first section of convolution: not all of row is used
		for( size_type n = filter_center; n != filter.sub_size() - 1; ++n )
		{
			for( size_type k = filter.vec_begin(); k != (n + 1); ++k )
			{
				y(n) += filter(k) * at(row_begin() + n - k, j);
			}
		}

		// middle: all of filter and row are used
		for( size_type n = filter.sub_size() - 1; n != sub_row_size(); ++n )
		{
			for( size_type k = filter.vec_begin(); k != filter.sub_size(); ++k )
			{
				y(n) += filter(k) * at(row_begin() + n - k, j);
			}
		}

		// last section of convolution: not all of filter is used
		for( size_type n = sub_row_size(), m = filter.vec_begin();
				m != filter_center; ++n )
		{
			++m;  // shift start of filter one position to the right

			for( size_type k = m; k != filter.sub_size(); ++k )
			{
				y(n) += filter(k) * at(row_begin() + n - k, j);
			}
		}

		// copy convolution result into this row (but only the section that lied
		// within the row originally)
		for( size_type i = row_begin(), i2 = filter_center;
				i != row_end();
				++i, ++i2 )
		{
			at(i,j) = y(i2);
		}
	}

	/*
	// version 2: rely on Vector class's convolve()

	// determine center of filter kernel
	const size_type filter_center = filter.size() / 2;

	// filter along each column
	for( size_type j = col_begin(); j != col_end(); ++j )
	{
		// result from convolution along a single column
		Vector col_conv( sub_row_size() );

		for( size_type i = row_begin(), i2 = row_conv.vec_begin();
				i != row_end();
				++i, ++i2 )
		{
			row_conv(i2) = at(i,j);
		}

		col_conv = col_conv.convolve( filter );

		for( size_type i = row_begin(), i2 = row_conv.vec_begin() + filter_center;
				i != row_end();
				++i, ++i2 )
		{
			at(i,j) = col_conv(i2);
		}
	}
	 */
	return( *this );
}

/**
	Convolve matrix with given separable filters.

	@param[in] row_filter Filter to convolve along each row with
	@param[in] col_filter Filter to convolve along each column with
 */
Matrix&
Matrix::convolve( const Vector& row_filter, const Vector& col_filter )
{
	convolve_row( row_filter );
	convolve_col( col_filter );

	/*
	// result from convolution along rows
	Matrix row_conv( sub_row_size(), sub_col_size() );

   // determine center of filter kernel
	size_type half_row_size = row_filter.size() / 2;

	// filter along each row
	for( size_type i = row_conv.row_begin(); i != row_conv.row_end(); ++i )
	{
		for( size_type j = row_conv.col_begin(); j != row_conv.col_end(); ++j )
		{
			// left border case
			size_type start_k = 0;
			if( j > half_row_size )
			{
				start_k = j - half_row_size;
			}

			// right border case
			size_type end_k = j + half_row_size + 1;
			if( j >= (sub_col_size() - half_row_size) )
			{
				end_k = sub_col_size();
			}

			// convolve at position (i,j)
			for( size_type k = start_k; k != end_k; ++k )
			{
				row_conv(i, j) += at(i, k) * row_filter( k - start_k );
			}
		}
	}

	size_type half_col_size = col_filter.size() / 2;

	// filter along each row
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			// left border case
			size_type start_k = 0;
			if( i > half_col_size )
			{
				start_k = i - half_col_size;
			}

			// right border case
			size_type end_k = i + half_col_size + 1;
			if( i >= (sub_col_size() - half_col_size) )
			{
				end_k = sub_col_size();
			}

			// convolve at position (i,j)
			at(i, j) = 0;
			for( size_type k = start_k; k != end_k; ++k )
			{
				// store result in *this
				at(i, j) += row_conv(k, j) * col_filter( k - start_k );
			}
		}
	}
	 */
	return( *this );
}

/**
	Convolve matrix with given filter.

	TODO Carefully check function as this is mostly copied from older code I
	wrote

	Implementation note:
	We use temporary storage since a convolution uses values
	around a given point, e.g., suppose pixel (3,3) is modified by a
	convolution. When the filter is moved and centered over pixel (3,4), it
	will need the original value of (3,3), not the new value.
	A smarter implementation might save only (filter.row_size() / 2) rows and
	(filter.col_size() / 2) columns from the original matrix. This, of course,
	would be more complicated and possibly incorrect.

	@param[in] filter Filter to convolve matrix with
 */
Matrix&
Matrix::convolve( const Matrix& filter )
{
	// convolution result
	Matrix conv( sub_row_size(), sub_col_size() );

	// determine center of filter matrix
	size_type filter_row_center = filter.real_row_size() / 2;
	size_type filter_col_center = filter.real_col_size() / 2;

	// convolve matrix with filter--do not process border
	for( size_type i = filter_row_center;
			i != (sub_row_size() - filter_row_center);
			++i )
	{
		for( size_type j = filter_col_center;
				j != (sub_col_size() - filter_col_center);
				++j )
		{
			// calculate difference to determine where in matrix we are
			size_type row_diff = i - filter_row_center;
			size_type col_diff = j - filter_col_center;

			// calculate convolution for position (i,j) by iterating through filter
			for( size_type r = 0; r != filter.real_row_size(); ++r )
			{
				for( size_type s = 0; s != filter.real_col_size(); ++s )
				{
					// update sum as pass through filter
					conv(i, j) += at( row_diff + r, col_diff + s ) * filter(r, s);
				}
			}
		}
	}

	// copy convolution result into *this
	for( size_type i = row_begin(), c_i = 0; i != row_end(); ++i, ++c_i )
	{
		for( size_type j = col_begin(), c_j = 0; j != col_end(); ++j, ++c_j )
		{
			at(i, j) = conv(c_i, c_j);
		}
	}

	return( *this );
}

#include "Linear_Algebra.hpp"

/**
	Obtain principal-component vectors and related quantities.
	P = princomp(X, Q) computes the principal-component vectors of the vector
	population contained in the rows of X, a matrix of size K-by-n where K is
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
Matrix::prin_comp( unsigned Q )
{
	if( Q > sub_col_size() )
	{
		err_quit( "Matrix::prin_comp: %u >= %u\n", Q, sub_col_size() );
	}

	// [K, n] = size(X);
	// X = double(X);

	// obtain mean vector and covariance matrix of the vectors in X
	// [P.Cx, P.mx] = covmatrix(X);
	// P.mx = P.mx';
	Matrix Cx;
	Vector mx;
	covariance( Cx, mx );

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
	Matrix Y = A * this->transpose();

	// obtain reconstructed matrix
	// P.X = (P.A' * P.Y)' + Mx;
	Matrix Mx( sub_row_size(), mx.sub_size() );
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
	Computes the covariance matrix 'cov' and the mean vector 'mean' of a vector
	population organized as the rows of this matrix. 'cov' is of size N-by-N and
	'mean' is of size N-by-1, where N is the dimension of the vectors (i.e., the
	number of columns of X).  

	@param[out] cov Covariance matrix
	@param[out] mean Mean/average of each column of this matrix

	@post This matrix is altered upon returning such that the mean is subtracted
	from each column. Hence, a copy should be made prior to calling covariance
	if the original matrix is to be preserved.  
 */
void
Matrix::covariance( Matrix& cov, Vector& mean )
{
	// allocate mean vector if necessary
	if( mean.sub_size() != sub_col_size() )
	{
		mean = Vector( sub_col_size() );
	}

	// compute an unbiased estimate of mean for each column
	for( size_type j = col_begin(), k = mean.vec_begin();
			j != col_end(); ++j, ++k )
	{
		unsigned   item_num = 0;   // current item number
		value_type unused_var = 0; // unused, but we must pass a variance pointer

		for( size_type i = row_begin(); i != row_end(); ++i )
		{
			ws_tools::stat( at(i,j), item_num++, &mean(k), &unused_var );
		}
	}

	// for each column, subtract the column's mean from each component of column
	for( size_type j = col_begin(), k = mean.vec_begin();
			j != col_end(); ++j, ++k )
	{
		for( size_type i = row_begin(); i != row_end(); ++i )
		{
			at(i,j) -= mean(k);
		}
	}

	// compute an unbiased estimate of 'cov': E{X' * X} where X = *this
	// Note that the product is X'*X since the population vectors are the rows of
	// X, and this equivalent to the following inefficient code:
	// cov = (transpose() * (*this)) / (sub_row_size() - 1);
	// which does not make use of the symmetry of X' * X and creates unnecessary
	// temporary objects.

	// allocate covariance matrix if necessary
	if( cov.sub_row_size() != sub_col_size()
			|| cov.sub_col_size() != sub_col_size() )
	{
		cov.resize( sub_col_size(), sub_col_size() );
	}

	// same as above but avoids computations by computing multiplication
	// directly and makes use of the symmetry of the covariance matrix

	// for each row and column
	for( size_type j = col_begin(); j != col_end(); ++j )
	{
		for( size_type i = j; i != col_end(); ++i )
		{
			// compute dot product of a row in X' with a column in X, but using
			// only X (i.e., without explicitly computing X') 
			value_type sum = 0;
			for( size_type k = row_begin(); k != row_end(); ++k )
			{
				sum += at(k,j) * at(k,i);
			}

			// store value (i,j) (note that the covariance matrix is symmetric)
			cov(i,j) = cov(j, i) = sum;
		}
	}

	// finish computation of the covariance matrix (the -1 in the division
	// provides an unbiased estimation) 
	cov /= (sub_row_size() - 1);
}

/**
	Computes the covariance matrix 'cov' and the mean vector 'mean' of a vector
	population organized as the rows of this matrix. 'cov' is of size N-by-N and
	'mean' is of size N-by-1, where N is the dimension of the vectors (i.e., the
	number of columns of X).  

	@param[out] cov Covariance matrix
	@param[out] mean Mean/average of each column of this matrix

	@post This matrix is altered upon returning such that the mean is subtracted
	from each column. Hence, a copy should be made prior to calling covariance
	if the original matrix is to be preserved.  
 */
Matrix
Matrix::covariance( ) const
{
	Vector mean( sub_col_size() );
	Matrix cov( sub_col_size(), sub_col_size() );

	// compute an unbiased estimate of mean for each column
	for( size_type j = col_begin(), k = mean.vec_begin();
			j != col_end(); ++j, ++k )
	{
		unsigned   item_num = 0;   // current item number
		value_type unused_var = 0; // unused, but we must pass a variance pointer

		for( size_type i = row_begin(); i != row_end(); ++i )
		{
			ws_tools::stat( at(i,j), item_num++, &mean(k), &unused_var );
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

	// for each row and column
	size_type mj = mean.vec_begin();
	for( size_type j = col_begin(); j != col_end(); ++j )
	{
		const value_type mean_j = mean(mj);

		size_type mi = mean.vec_begin();
		for( size_type i = j; i != col_end(); ++i )
		{
			const value_type mean_i = mean(mi);

			// compute dot product of a row in X' with a column in X, but using
			// only X (i.e., without explicitly computing X') 
			value_type sum = 0;
			for( size_type k = row_begin(); k != row_end(); ++k )
			{
				sum += (at(k,j) - mean_j) * (at(k,i) - mean_i);
			}

			// store value (i,j)
			// (note that the covariance matrix is symmetric and the -1 in the
			// division yields an unbiased estimation)
			cov(i,j) = cov(j, i) = sum /= (sub_row_size() - 1);
			++mi;
		}
		++mj;
	}

	return( cov );
}

/**
	Resize the matrix to have the given number of rows and columns.

	@param[in] new_row New number of rows
	@param[in] new_col New number of columns
	@param[in] interpolate Whether to perform interpolation of matrix values
		(slower) or just to allocate a larger matrix (faster but data is lost)
	@retval this This matrix
 */
Matrix&
Matrix::resize( size_type new_row, size_type new_col, bool interpolate )
{
	// do nothing if we already have the requested size
	if( new_row == real_row_size() && new_col == real_col_size() )
	{
		return( *this );
	}
	else if( !interpolate )
	{
		allocate_matrix( new_row, new_col );
		return( *this );
	}

	// compute scale factors and return the rescaled matrix
	// double scale_row = new_row / static_cast<double>( real_row_size() );
	// double scale_col = new_col / static_cast<double>( real_col_size() );
	// return( scale_nointerp( scale_row, scale_col ) );

	// rescale the matrix with interpolation
	return( bicubic_scale( new_row, new_col ) );
}

/**
	Scale matrix by given scale factor. Both the width and the height are scaled
	by the same amount.

	@param[in] s Scale factor for both x and y coordinates (columns and rows)
	@retval this This matrix
 */
Matrix&
Matrix::scale( double s )
{
	return( scale( s, s ) );
}

/**
	Scale matrix by given scale factors using bicubic interpolation.

	Code is based off scaling code in the Gimp program in the file
		'gimp-2.2.6/app/paint-funcs/paint-funcs.c'
	or
		'gimp-2.3.6/app/paint-funcs/scale-funcs.c'
	depending on the version of Gimp.

	Data is scaled a row at a time.

	@param[in] s_y Scale factor for y coordinates (number of rows)
	@param[in] s_x Scale factor for x coordinates (number of columns)
	@retval this This matrix scaled by the given amount
 */
Matrix&
Matrix::scale( double s_y, double s_x )
{
	// calculate new size
	size_type new_height = static_cast<size_type>(
			round( s_y * real_row_size() ) );
	size_type new_width  = static_cast<size_type>(
			round( s_x * real_col_size() ) );

	return( bicubic_scale( new_height, new_width ) );
}

/**
	Scale matrix to the given size using bicubic interpolation.

	Code is based off scaling code in the GIMP program in the file
		'gimp-2.2.6/app/paint-funcs/paint-funcs.c'
	or
		'gimp-2.3.6/app/paint-funcs/scale-funcs.c'
	depending on the version of GIMP.

	Bicubic interpolation fits a cubic polynomial first to each set of 4 data
	points along a row and then to 4 points down a column.  

	@param[in] s_y Scale factor for y coordinates (number of rows)
	@param[in] s_x Scale factor for x coordinates (number of columns)
	@retval this This matrix scaled by the given amount
 */
Matrix&
Matrix::bicubic_scale( size_type new_height, size_type new_width )
{
	// output matrix after rescaling
	Matrix scaled_mat( new_height, new_width );

	int old_y = -4;

	// find the (actual) ratio of the old height to the new height
	double y_ratio = static_cast<double>( real_row_size() )
			/ static_cast<double>( scaled_mat.real_row_size() );

	// a single row from the original, unscaled matrix
	// Two points:
	// - The actual matrix data will be filled in by the function
	// get_scaled_row() before being used.
	// - Two additional pixels are added at both the front and the back of this
	// array to handle pixels along the image border; hence, the range of valid
	// indices is [-2, col_size() + 1]. We have to do this since cubic
	// interpolation needs 4 input values per single output value, so
	// calculation of a border value (i.e., located at column position 0) will
	// check positions -2, -1, 0, and 1 to compute the single output at position
	// 0.
	// Another note: We should not use the Vector class here since we use
	// negative indices below. If we used a Vector, then the at() function would
	// be called with say, -1, which will become a very large positive number
	// since at() takes in unsigned integers (hence, an implicit cast occurs).
	// We are sure to get a segfault or strange behavior by doing this.
	double* unscaled_row = new double[ real_col_size() + 4 ];
	unscaled_row += 2;

	// allocate 4 temporary data arrays for storing the scaled output--we
	// need 4 since we are doing bicubic interpolation, which means we need to
	// keep track of multiple rows so we can do interpolation across them
	double* scaled_rows[4];
	for( size_type i = 0; i != 4; ++i )
	{
		scaled_rows[i] = new double[ scaled_mat.real_col_size() ];

		// initialize to 0: it's possible that this doesn't need to be done since
		// I fixed some other bugs, but it doesn't hurt either
		for( size_type x = 0; x != scaled_mat.real_col_size(); ++x )
		{
			scaled_rows[i][x] = 0;
		}
	}

	// scale the selected region
	for( size_type y = 0; y < scaled_mat.real_row_size(); ++y )
	{
		// scaled matrix is smaller
		if( scaled_mat.real_row_size() < real_row_size() )
		{
			const double inv_ratio = 1.0 / y_ratio;

			// if this is the first time through, load the first row 
			if( y == 0 )
			{
				get_scaled_row( unscaled_row, &scaled_rows[0],
						0, scaled_mat.real_col_size() );
			}
			int new_y = static_cast<int>(y * y_ratio);
			double frac = 1.0 - (y * y_ratio - new_y);

			for( size_type x = 0; x < scaled_mat.real_col_size(); ++x )
			{
				const value_type new_value =
						static_cast<value_type>( frac * scaled_rows[3][x] );
				scaled_mat(y,x) = new_value;
			}

			int max = (int) ((y + 1) * y_ratio) - new_y - 1;

			get_scaled_row( unscaled_row, &scaled_rows[0],
					++new_y, scaled_mat.real_col_size() );

			while( max > 0 )
			{
				for( size_type x = 0; x < scaled_mat.real_col_size(); ++x )
				{
					const value_type new_value =
							static_cast<value_type>( scaled_rows[3][x] );
					scaled_mat(y,x) += new_value;
				}
				get_scaled_row( unscaled_row, &scaled_rows[0],
						++new_y, scaled_mat.real_col_size() );
				--max;
			}
			frac = (y + 1) * y_ratio - ((int) ((y + 1) * y_ratio));

			for( size_type x = 0; x < scaled_mat.real_col_size(); ++x )
			{
				const value_type new_value =
						static_cast<value_type>( frac * scaled_rows[3][x] );
				scaled_mat(y,x) += new_value;
				scaled_mat(y,x) *= inv_ratio;
			}
		}

		// scaled matrix is larger
		else if( scaled_mat.real_row_size() > real_row_size() )
		{
			const int new_y = static_cast<int>( floor( y * y_ratio - 0.5 ) );

			while( old_y <= new_y )
			{
				// get the necesary lines from the source matrix, scale them,
				// and put them into scaled_rows[]
				get_scaled_row( unscaled_row, &scaled_rows[0],
						old_y + 2, scaled_mat.real_col_size() );
				++old_y;
			}

			// apply cubic interpolation down each column
			const double dy = (y * y_ratio - 0.5) - new_y;

			const double p0 = cubic_0( dy );  // = cubic( dy, 1, 0, 0, 0 );
			const double p1 = cubic_1( dy );  // = cubic( dy, 0, 1, 0, 0 );
			const double p2 = cubic_2( dy );  // = cubic( dy, 0, 0, 1, 0 );
			const double p3 = cubic_3( dy );  // = cubic( dy, 0, 0, 0, 1 );
			for( size_type x = 0; x != scaled_mat.real_col_size(); ++x )
			{
				const value_type new_value =
						static_cast<value_type>( 
								(p0 * scaled_rows[0][x] + p1 * scaled_rows[1][x]
								 + p2 * scaled_rows[2][x] + p3 * scaled_rows[3][x])
						);
				scaled_mat(y,x) = new_value;
			}
		}
		else // scaled_mat.real_row_size() == real_row_size()
		{
			get_scaled_row( unscaled_row, &scaled_rows[0],
					y, scaled_mat.real_col_size() );

			// copy row into output matrix
			for( size_type x = 0; x != scaled_mat.real_col_size(); ++x )
			{
				const value_type new_value =
						static_cast<value_type>( scaled_rows[3][x] );
				scaled_mat(y,x) = new_value;
			}
		}
	}

	// free temporary arrays
	for( size_type i = 0; i != 4; ++i )
	{
		delete [] scaled_rows[i];
	}

	unscaled_row -= 2;  // correct for previous offset of 2
	delete [] unscaled_row;

	// copy the scaled matrix into this matrix, which is then returned
	return( *this = scaled_mat );
}

/**
	Get the necesary lines from the source matrix, scale them, and put them into
	scaled_rows[].

	@param[in,out] unscaled_row A single row from the original, unscaled matrix
	@param[in,out] scaled_rows A set of rows in the new, scaled matrix
	@param[in,out] scaled_row Pointer to output row of matrix
	@param[in] y Row position
	@param[in] scaled_width Width of scaled matrix
 */
void
Matrix::get_scaled_row( double* unscaled_row, double** scaled_rows,
		int y, size_type scaled_width ) const
{
	// rotate_pointers( scaled_rows, 4 );
	rotate_pointers( scaled_rows );

	if( y < 0 )
	{
		y = 0;
	}
	if( y < (int) real_row_size() )
	{
		// copy row y from the original matrix;
		// the border pixels are set to their nearest neighbor
		// (note: 'unscaled_row' has 2 extra elements at both the front and
		// back; hence, the indices are valid)
		unscaled_row[-2] = at(y,0);
		unscaled_row[-1] = at(y,0);
		for( size_type x = 0; x != real_col_size(); ++x )
		{
			unscaled_row[x] = at(y,x);
		}
		unscaled_row[ real_col_size() ]     = unscaled_row[real_col_size() - 1];
		unscaled_row[ real_col_size() + 1 ] = unscaled_row[real_col_size() - 1];

		// if scaling up, make the row longer
		if( scaled_width > real_col_size() )
		{
			upsample_row( unscaled_row, scaled_rows[3], scaled_width );
		}

		// if scaling down, make the row shorter
		else if( scaled_width < real_col_size() )
		{
			downsample_row( unscaled_row, scaled_rows[3], scaled_width );
		}

		else // scaled_width == real_col_size(), so no scaling is needed
		{
			memcpy( scaled_rows[3], unscaled_row, sizeof(double) * scaled_width );
		}
	}
	else // y >= real_row_size()
	{
		// directly copy the previous row
		memcpy( scaled_rows[3], scaled_rows[2], sizeof(double) * scaled_width );
	}
}

/**
	Upsample the row.

	@param[in] src_row A single row from the original, unscaled matrix
	@param[out] dest_row A single row from the new, scaled matrix
	@param[in] scaled_width Width of scaled matrix
 */
void
Matrix::upsample_row( const double* src, double* dest,
		size_type scaled_width ) const
{
	const double ratio = real_col_size() / static_cast<double>( scaled_width );

	// this could be optimized much more by precalculating the coefficients for
	// each x
	for( size_type x = 0; x != scaled_width; ++x )
	{
		// +2, -2 is there because (int) rounds towards 0 and we need
		// to round down (src_col is allowed to start at -1 since we actually
		// have two elements before the start of 'src' and note that we take
		// s[-1], which is src_col[-2])
		int src_col = ((int) (x * ratio + 2.0 - 0.5)) - 2;
		assert( src_col >= -1 );
		const double* s = &src[src_col];

		double frac = (x * ratio - 0.5) - src_col;

		// apply cubic interpolation
		*dest = cubic( frac, s[-1], s[0], s[1], s[2] );
		++dest;
	}
}

/**
	Downsample the row.

	Calculate the weighted average of pixel data that each output pixel must
	receive, taking into account that it always scales down, i.e., there's
	always more than one input pixel per output pixel.

	@param[in] src_row A single row from the original, unscaled matrix
	@param[out] dest_row A single row from the new, scaled matrix
	@param[in] scaled_width Width of scaled matrix
 */
void
Matrix::downsample_row( const double* src_row, double* dest_row,
		size_type scaled_width ) const
{
	const double avg_ratio = scaled_width
			/ static_cast<double>( real_col_size() );
	const double inv_width = 1.0 / scaled_width;

	// pointers to the start of the given data
	const double* src_ptr  = src_row;
	double*       dest_ptr = dest_row;

	size_type slice_pos = 0;    // slice position relative to scaled_width
	double    accum     = 0.0;  // output accumulator

	// for the current row
	for( size_type x = 0; x < scaled_width; ++x )
	{
		do // accumulate whole pixels
		{
			accum += *src_ptr++;

			slice_pos += scaled_width;
		}
		while( slice_pos < real_col_size() );
		slice_pos -= real_col_size();

		assert( slice_pos < scaled_width );

		if( slice_pos == 0 )
		{
			// Simplest case: we have reached a whole pixel boundary.  Store
			// the average value per channel and reset the accumulators for
			// the next round.

			// The main reason to treat this case separately is to avoid an
			// access to out-of-bounds memory for the first pixel.
			*dest_ptr++ = accum * avg_ratio;
			accum = 0.0;
		}
		else
		{
			// We have accumulated a whole pixel per channel where just a
			// slice of it was needed.  Subtract now the previous pixel's
			// extra slice.
			// Note: we need an int here, so it can be negative, which is valid
			// since src_ptr corresponds to row in scale() with an extra element
			// allocated at its front.
			double slice = slice_pos * inv_width * src_ptr[-1];
			*dest_ptr++ = (accum - slice) * avg_ratio;

			// That slice is the initial value for the next round.
			accum = slice;
		}
	}

	// verify that src_ptr points to the next-to-last position and that slice_pos
	// should be zero.
	assert( static_cast<size_type>( src_ptr - src_row ) == real_col_size()
				&& slice_pos == 0 );
}

/**
	Scale image by given scale factors without interpolation between pixels.

	No interpolation is performed unlike with Matrix::scale(), so pixels are
	simply removed when downsampling or form blocks when upsampling.

	@param[in] s_y Scale factor for y coordinates (rows)
	@param[in] s_x Scale factor for x coordinates (columns)
*/
Matrix&
Matrix::nointerp_scale( double s_y, double s_x )
{
	// allocate new, scaled image with the scaled number of rows and columns
	Matrix scaled_mat(
			static_cast<size_type>( s_y * real_row_size() ),
			static_cast<size_type>( s_x * real_col_size() )
		);

	// for each pixel of the scaled image, calculate the pixel location that
	// maps to it from the original image
	for( size_type i = scaled_mat.row_begin(); i != scaled_mat.row_end(); ++i )
	{
		size_type y = static_cast<size_type>( i / s_y );

		for( size_type j = scaled_mat.col_begin();
				j != scaled_mat.col_end();
				++j )
		{
			size_type x = static_cast<size_type>( j / s_x );

			// if the pixel position is inside original image
			if( y >= row_begin() && y < row_end()
					&& x >= col_begin() && x < col_end() )
			{
				scaled_mat(i,j) = at(y,x);
			}
		}
	}

	return( (*this = scaled_mat) );
}

/**
	Reshape the matrix to have new dimensions but the same elements.

	For example, if A is a 2x2 matrix, then A.reshape(4,1) will return a 4x1
	matrix (column vector) with the same elements. A.reshape(5,1), for instance,
	is invalid since the total number of elements in the output matrix is not
	the same.  

	@param[in] new_rows New number of rows
	@param[in] new_cols New number of columns
 */
Matrix
Matrix::reshape( size_type new_rows, size_type new_cols ) const
{
	if( (new_rows * new_cols) != (sub_row_size() * sub_col_size()) )
	{
		err_quit( "Matrix::reshape: Number of elements (%u) in new matrix "
				"does not equal number of elements in this matrix (%u)\n",
				(new_rows * new_cols), (sub_row_size() * sub_col_size()) );
	}

	Matrix A( new_rows, new_cols );

	// iterators to this matrix
	const_iterator iter     = const_begin();
	const_iterator end_iter = const_end();

	// iterator to output matrix
	iterator A_iter = A.begin();

	// copy matrix
	for( ; iter != end_iter; )
	{
		*A_iter++ = *iter++;
	}

	return( A );
}

/**
	Apply fourier-transform to the this matrix.

	The output is not normalized.

	@param[out] real_data The real components of the transform
	@param[out] imag_data The imaginary components of the transform
 */
void
Matrix::fft( Matrix& real_data, Matrix& imag_data ) const
{
#ifdef HAVE_FFTW

	const unsigned   N = sub_row_size() * sub_col_size();
	complex<double>* in;
	complex<double>* out;
	fftw_plan p;

	in  = new complex<double>[ N ];
	out = new complex<double>[ N ];

	// copy input (there is only a real part)
	unsigned k = 0;
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			in[k] = complex<double>( at(i,j), 0 );
			++k;
		}
	}

	// make a plan and apply the forward FFT
	p = fftw_plan_dft_2d( sub_row_size(), sub_col_size(),
			reinterpret_cast<fftw_complex*>(in),
			reinterpret_cast<fftw_complex*>(out),
			FFTW_FORWARD, FFTW_ESTIMATE
		);
	fftw_execute(p);

	real_data.resize( sub_row_size(), sub_col_size() );
	imag_data.resize( sub_row_size(), sub_col_size() );

	// copy real and imaginary parts
	k = 0;
	for( size_type i = real_data.row_begin(); i != real_data.row_end(); ++i )
	{
		for( size_type j = real_data.col_begin(); j != real_data.col_end(); ++j )
		{
			real_data(i,j) = out[k].real();
			imag_data(i,j) = out[k].imag();
			++k;
		}
	}

	fftw_destroy_plan(p);

	delete [] in;
	delete [] out;

#else
	err_warn( "Matrix::fft: Unable to apply forward Fourier transform: "
			" FFTW (http://www.fftw.org/) is not installed.\n" );
	return;
#endif // HAVE_FFTW
}

/**
	Apply inverse fourier-transform to the given data.

	@param[in] real_data The real components of the transform
	@param[in] imag_data The imaginary components of the transform

	@post This matrix holds the inverse of the normalized FFT.
 */
void
Matrix::ifft( Matrix& real_data, Matrix& imag_data )
{
#ifdef HAVE_FFTW

	// we should have the same number of components
	if( real_data.sub_row_size() != imag_data.sub_row_size()
			|| real_data.sub_col_size() != imag_data.sub_col_size()
		)
	{
		err_quit( "Matrix::ifft: Real component matrix's size (%u x %u)"
				" does not match imaginary component matrix's size (%u x %u)\n",
				real_data.sub_row_size(),
				real_data.sub_col_size(),
				imag_data.sub_row_size(),
				imag_data.sub_col_size() );
	}

	// get the size of the data
	const size_type data_rows = real_data.sub_row_size();
	const size_type data_cols = real_data.sub_col_size();
	const unsigned  N         = data_rows * data_cols;

	// get input and output data arrays
	complex<double>* in  = new complex<double>[ N ];
	complex<double>* out = new complex<double>[ N ];

	// copy input
	unsigned k = 0;
	for( size_type i = real_data.row_begin(); i != real_data.row_end(); ++i )
	{
		for( size_type j = real_data.col_begin(); j != real_data.col_end(); ++j )
		{
			in[k] = complex<double>( real_data(i,j), imag_data(i,j) );
			++k;
		}
	}

	// make a plan and apply the inverse FFT
	fftw_plan p = fftw_plan_dft_2d(
			data_rows, data_cols,
			reinterpret_cast<fftw_complex*>(in),
			reinterpret_cast<fftw_complex*>(out),
			FFTW_BACKWARD, FFTW_ESTIMATE
		);
	fftw_execute(p);

	// ensure this matrix is of the correct size
	resize( data_rows, data_cols );

	// copy only the real part after rescaling
	k = 0;
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			at(i,j) = out[k].real() / N;
			++k;
		}
	}

	fftw_destroy_plan(p);

	delete [] in;
	delete [] out;

#else
	err_warn( "Matrix::ifft: Unable to apply inverse Fourier transform: "
			" FFTW (http://www.fftw.org/) is not installed.\n" );
	return;
#endif // HAVE_FFTW
}

#if 0
/**
	Apply inverse fourier-transform to the given data.

	@param[in] real_data The real components of the transform
	@param[in] imag_data The imaginary components of the transform

	@post This matrix holds the inverse of the normalized FFT.
 */
void
Matrix::ifft( Matrix& real_data, Matrix& imag_data )
{
}
#endif
