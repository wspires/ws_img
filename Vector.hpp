/**
	@file   Vector.hpp
	@author Wade Spires
	@date   2005/10/28
	@brief  Definition/implementation of vector type.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef VECTOR_H
#define VECTOR_H

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include  <algorithm>
#include  <functional>
#include  <stdexcept>
#include  <iomanip>
#include  <iostream>
#include  <iterator>
#include  <sstream>
#include  <string>
#include  <vector>

// c headers
#include  <cfloat>
#include  <climits>
#include  <cmath>

#include "Matrix.hpp"
#include "ws_tools.hpp"

namespace ws_img
{

/**
	@brief Class representing a mathematical vector of data for
	performing numerical computations
 */
class Vector
{

public:

	enum Orient_Type { Row, Col };
	typedef Matrix_Region Vector_Region;

private:

	Matrix  _vector;       // vector stored internally as a matrix

public:

	typedef Matrix::value_type       value_type;
	typedef Matrix::size_type        size_type;
	typedef Matrix::reference        reference;
	typedef Matrix::const_reference  const_reference;
	typedef Matrix::pointer          pointer;
	typedef Matrix::const_pointer    const_pointer;
	typedef Matrix::difference_type  difference_type;

	// vector stored as a single row of a matrix
	typedef Matrix::row_iterator        iterator;
	typedef Matrix::const_row_iterator  const_iterator;

	typedef Matrix::row_col_iterator        row_col_iterator;
	typedef Matrix::const_row_col_iterator  const_row_col_iterator;

	/**
		Allocate single, uninitialized, element vector.
		Note:
			Vector( value_type default_value = 0 ) : vec( 1, default_value )
			is ambiguous with Vector( size_type, value_type = 0 ), 
			e.g., Vector( 5 ) can reference either one.
	 */
	Vector( );

	/**
		Allocate M-dimensional vector. 
		@param[in] M Number of elements in vector
		@param[in] default_value Default value to initialize each element to
	 */
	Vector( size_type, value_type = 0 );

	/**
		Allocate M-dimensional vector. 
		@param[in] M Number of elements in vector
		@param[in] default_value Default value to initialize each element to
	 */
	Vector( const_iterator, const_iterator );

	/**
		Construct vector by reading given matrix file.
		@param[in] matrix_file_name Name of matrix file to read matrix from
	 */
	Vector( const std::string& );

	/**
		Construct vector by using the given range.

		For instance, passing 0, 2, and 10 yields the vector [ 0 2 4 6 8 10 ].

		Warning: passing 1, -.05, .55 only yields [ 1 .95 ... .60 ], i.e., the
		.55 is omitted. This is due to the inherent imprecision of floating-point
		math.
		TODO we may be able to rescale the values s.t. they're >= 1

		@param[in] begin First value in the range
		@param[in] inc Amount to increment each successive value
		@param[in] end Last value in the range
	 */
	Vector( const double, const double, const double );

	/**
		Construct vector by copying elements from matrix.
		@param[in] M Matrix
	 */
	Vector( const Matrix& );

	/**
		Assign vector contents of matrix.
		@param[in] M Row or column matrix
	 */
	Vector& operator=( const Matrix& );

	Vector& operator=( const value_type );

	/**
		Default destructor does nothing since no member variables are
		dynamically allocated.
	 */
	virtual ~Vector( );

	/**
		Access row i in vector.
		@param[in] i Vector index
		@return ref Reference to position i in vector
	 */
	reference at( size_type );

	/**
		Access row i in vector.
		@param[in] i Vector index
		@return ref Reference to position i in vector
	 */
	const_reference at( size_type ) const;

	/**
		Access row i in vector.
		@param[in] i Vector index
		@return ref Reference to position i in vector
	 */
	reference operator()( size_type );

	/**
		Access row i in vector.
		@param[in] i Vector index
		@return ref Reference to position i in vector
	 */
	const_reference operator()( size_type ) const;

	/**
		Access row i in vector.
		@param[in] i Vector index
		@return ref Reference to position i in vector
	 */
	reference operator[]( size_type );

	/**
		Access row i in vector.
		@param[in] i Vector index
		@return ref Reference to position i in vector
	 */
	const_reference operator[]( size_type ) const;

	/**
		Get real size of vector.
		@return size Vector size
	 */
	size_type real_size( ) const;

	/**
		Get size of vector.
		@return size Vector size
	 */
	size_type size( ) const;

	/**
		Resize vector.
		@param[in] new_size New size of vector
		@param[in] rescale Whether to perform rescaling of matrix values,
			i.e., interpolation, (slower) or just to allocate a larger matrix
			(faster but data is lost)

		@retval this This vector
	 */
	Vector& resize( size_type, bool = true );

	/**
		Read vector as a MATLAB matrix from file matrix_file_name.
		@param[in] file_name Name of vector file to read vector from
	 */
	void read( const std::string& );

	/**
		Write contents of matrix to file file_name using specified precision for
		each value.
		@param[in] file_name Name of file to write matrix to
		@param[in] precision Precision to use for floating-point values
	 */
	void write( const std::string& = "", unsigned = 5 ) const;

	/**
		Get iterator to start of vector.
		@retval iter Iterator to row
	 */
	iterator begin( );

	/**
		Get iterator to start of vector.
		@retval iter Iterator to row
	 */
	const_iterator const_begin( ) const;

	/**
		Get iterator to end of vector.
		@retval iter Iterator to row
	 */
	iterator end( );

	/**
		Get iterator to end of vector.
		@retval iter Iterator to row
	 */
	const_iterator const_end( ) const;

	/**
		Clear contents of vector.
	 */
	void clear( );

public:

	Vector& operator+=( const Vector& rhs );
	Vector  operator+( const Vector& ) const;

	Vector& operator-=( const Vector& );
	Vector  operator-( const Vector& ) const;

	Matrix  operator*( const Vector& ) const;
	Matrix  operator*( const Matrix& ) const;

	Vector& operator+=( const value_type );
	Vector  operator+( const value_type ) const;
	Vector& operator-=( const value_type );
	Vector  operator-( const value_type ) const;
	Vector& operator*=( const value_type );
	Vector  operator*( const value_type ) const;
	Vector& operator/=( const value_type );
	Vector  operator/( const value_type ) const;

	/**
		Add vector by value.
		Some older versions of gcc complain if friend functions are not defined
		in a header file.
		@param[in] rhs Value on right-hand side of +
		@retval vector Vector result 
	 */
	friend Vector operator+( const value_type value, const Vector& rhs )
	{
		return( rhs + value );
	}

	/**
		Subtract vector from value.
		Some older versions of gcc complain if friend functions are not defined
		in a header file.
		@param[in] rhs Value on right-hand side of -
		@retval vector Vector result 
	 */
	friend Vector operator-( const value_type value, const Vector& rhs )
	{
		return( rhs - value );
	}

	/**
		Multiply matrix and vector resulting in a vector.
		Some older versions of gcc complain if friend functions are not defined
		in a header file.
		@param[in] lhs Matrix on left-hand side of *
		@param[in] rhs Vector on right-hand side of *
		@retval vector Vector result 
	 */
	friend Vector operator*( const Matrix& lhs, const Vector& rhs )
	{
		// Matrix(rhs) converts rhs into a matrix containing a single column
		// TODO speed this up by implementing the multiplication directly
		return( Vector( lhs * Matrix(rhs) ) );
	}

	/**
		Multiply vector by value.
		Some older versions of gcc complain if friend functions are not defined
		in a header file.
		@param[in] rhs Value on right-hand side of *
		@retval vector Vector result 
	 */
	friend Vector operator*( const value_type value, const Vector& rhs )
	{
		return( rhs * value );
	}

	/**
		Divide vector by value.
		Some older versions of gcc complain if friend functions are not defined
		in a header file.
		@param[in] rhs Value on right-hand side of /
		@retval vector Vector result 
	 */
	friend Vector operator/( const value_type value, const Vector& rhs )
	{
		return( rhs / value );
	}

	Vector operator-( ) const;

	value_type sum( ) const;
	value_type dot_product( const Vector& ) const;
	Vector cross_product( const Vector& ) const;

	Vector transpose( );

	Vector& sort( );
	Vector& reverse_sort( );

	value_type two_norm( ) const;
	value_type two_norm_2( ) const;
	value_type inf_norm( );
	void normalize( );

	Vector convolve( const Vector& ) const;

	static Vector random( size_type, value_type = INT_MAX );

	std::vector<value_type> to_vector( ) const;
	double* to_array( ) const;

	static Vector get_row_vector( const Matrix&, size_type );
	static Vector get_col_vector( const Matrix&, size_type );
	static void set_row_vector( const Vector&, Matrix&, size_type );
	static void set_col_vector( const Vector&, Matrix&, size_type );

	void stat( double*, double*, double* ) const;

	value_type get_min( ) const;
	value_type get_max( ) const;
	void get_min_max( value_type*, value_type* ) const;

	Vector sin( ) const;
	Vector cos( ) const;
	Vector tan( ) const;
	Vector abs( ) const;

	// Fourier transform functions
	void fft( Vector&, Vector& ) const;
	void ifft( Vector&, Vector& );

public:  // region functions

	/**
		Get starting position of the region.
		
		The position actually refers to the current submatrix region, which is the
		first column of the vector (0) unless explicitly changed.

		@retval begin
	 */
	size_type vec_begin( ) const;

	/**
		Return last position of the region.
		@retval end 
	 */
	size_type vec_end( ) const;

	/**
		Return number of elements in subvector marked by region.

		The size actually refers to the current subvector region, which is the
		entire vector unless explicitly changed.

		@retval num_elem Number of elements
	 */
	size_type sub_size( ) const;

	/**
		Get current vector region.
		@retval region Current vector region
	 */
	Vector_Region get_region( ) const;

	/**
		Get current vector region.
		@retval region Current vector region
	 */
	Vector_Region region( ) const;

	/**
		Set vector region to the given region.

		@param[in] region_begin Starting position of the region
		@param[in] region_size Number of elements in the region
	 */
	void set_region( size_type, size_type );

	/**
		Set vector region to the given region.
		@param region New vector region
	 */
	void set_region( const Vector_Region& );

	/**
		Save current region.
		This is useful if only the current region needs to be saved before
		changing to multiple regions that do not need saving.
	 */
	void save_region( );

	/**
		Save current region and set new region.

		@param[in] region_begin Starting position of the region
		@param[in] region_size Number of elements in the region
	 */
	void push_region( size_type, size_type );

	/**
		Save current region and set new region.
		@param[in] region New vector region
	 */
	void push_region( const Vector_Region& );

	/**
		Restore last saved region.
	 */
	void pop_region( );

	/**
		Restore last saved region.
	 */
	void restore_region( );

	/**
		Reset vector region to the entire vector.
	 */
	void reset_region( );

};

inline Vector&
Vector::operator=( const value_type value )
{
	_vector = value;
	return( *this );
}

/**
	Access row i in vector.
	@param[in] i Vector index
	@return ref Reference to position i in vector
 */
inline Vector::reference
Vector::at( size_type i )
{
	return( _vector.at( 0, i ) );
}

/**
	Access row i in vector.
	@param[in] i Vector index
	@return ref Reference to position i in vector
 */
inline Vector::const_reference
Vector::at( size_type i ) const
{
	return( _vector.at( 0, i ) );
}

/**
	Access row i in vector.
	@param[in] i Vector index
	@return ref Reference to position i in vector
 */
inline Vector::reference
Vector::operator()( size_type i )
{
	return( at(i) );
}

/**
	Access row i in vector.
	@param[in] i Vector index
	@return ref Reference to position i in vector
 */
inline Vector::const_reference
Vector::operator()( size_type i ) const
{
	return( at(i) );
}

/**
	Access row i in vector.
	@param[in] i Vector index
	@return ref Reference to position i in vector
 */
inline Vector::reference
Vector::operator[]( size_type i )
{
	return( at(i) );
}

/**
	Access row i in vector.
	@param[in] i Vector index
	@return ref Reference to position i in vector
 */
inline Vector::const_reference
Vector::operator[]( size_type i ) const
{
	return( at(i) );
}

/**
	Get real size of vector.
	@return size Vector size
 */
inline Vector::size_type
Vector::real_size( ) const
{
	return( _vector.real_col_size() );
}

/**
	Get size of vector.
	@return size Vector size
 */
inline Vector::size_type
Vector::size( ) const
{
	return( _vector.sub_col_size() );
}

/**
	Resize vector.
	@param[in] new_size New size of vector
	@param[in] rescale Whether to perform rescaling of matrix values,
		i.e., interpolation, (slower) or just to allocate a larger matrix
		(faster but data is lost)

	@retval this This vector
 */
inline Vector&
Vector::resize( size_type new_size, bool rescale )
{
	_vector.resize( 1, new_size, rescale );
	return( *this );
}

/**
	Get iterator to start of vector.
	@retval iter Iterator to row
 */
inline Vector::iterator
Vector::begin( )
{
	// get matrix row iterator
	return( _vector.row_begin( 0 ) );
}

/**
	Get iterator to start of vector.
	@retval iter Iterator to row
 */
inline Vector::const_iterator
Vector::const_begin( ) const
{
	return( _vector.const_row_begin( 0 ) );
}

/**
	Get iterator to end of vector.
	@retval iter Iterator to row
 */
inline Vector::iterator
Vector::end( )
{
	return( _vector.row_end( 0 ) );
}

/**
	Get iterator to end of vector.
	@retval iter Iterator to row
 */
inline Vector::const_iterator
Vector::const_end( ) const
{
	return( _vector.const_row_end( 0 ) );
}

/**
	Clear contents of vector.
 */
inline void
Vector::clear( )
{
	_vector.clear();
}

/**
	Get starting position of the region.
	
	The position actually refers to the current submatrix region, which is the
	first column of the vector (0) unless explicitly changed.

	@retval begin
 */
inline Vector::size_type
Vector::vec_begin( ) const
{
	return( _vector.col_begin() );
}

/**
	Return last position of the region.
	@retval end 
 */
inline Vector::size_type
Vector::vec_end( ) const
{ 
	return( _vector.col_end() );
}

/**
	Return number of elements in subvector marked by region.

	The size actually refers to the current subvector region, which is the
	entire vector unless explicitly changed.

	@retval num_elem Number of elements
 */
inline Vector::size_type
Vector::sub_size( ) const
{ 
	return( _vector.sub_col_size() );	
}

/**
	Get current vector region.
	@retval region Current vector region
 */
inline Vector::Vector_Region
Vector::get_region( ) const
{
	return( _vector.get_region() );
}

/**
	Get current vector region.
	@retval region Current vector region
 */
inline Vector::Vector_Region
Vector::region( ) const
{
	return( _vector.region() );
}

/**
	Set vector region to the given region.

	@param[in] region_begin Starting position of the region
	@param[in] region_size Number of elements in the region
 */
inline void
Vector::set_region( size_type region_begin, size_type region_size )
{
	// vector physically stored as a row vector store
	_vector.set_region( 0, region_begin, 1, region_size );
}

/**
	Set vector region to the given region.
	@param region New vector region
 */
inline void
Vector::set_region( const Vector_Region& region )
{
	_vector.set_region( region );
}

/**
	Save current region.
	This is useful if only the current region needs to be saved before
	changing to multiple regions that do not need saving.
 */
inline void
Vector::save_region( )
{
	_vector.save_region();
}

/**
	Save current region and set new region.

	@param[in] region_begin Starting position of the region
	@param[in] region_size Number of elements in the region
 */
inline void
Vector::push_region( size_type region_begin, size_type region_size )
{
	// vector physically stored as a row vector store
	_vector.push_region( 0, region_begin, 1, region_size );
}

/**
	Save current region and set new region.
	@param[in] region New vector region
 */
inline void
Vector::push_region( const Vector_Region& region )
{
	_vector.push_region( region );
}

/**
	Restore last saved region.
 */
inline void
Vector::pop_region( )
{
	_vector.pop_region();
}

/**
	Restore last saved region.
 */
inline void
Vector::restore_region( )
{
	_vector.restore_region();
}

/**
	Reset vector region to the entire vector.
 */
inline void
Vector::reset_region( )
{
	_vector.reset_region();
}

} // namespace ws_img

#endif // VECTOR_H
