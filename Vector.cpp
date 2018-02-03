/**
	@file   Vector.cpp
	@author Wade Spires
	@date   2005/10/28
	@brief  Definition/implementation of vector type.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Vector.hpp"

// c++ headers
#include  <string>
#include  <vector>

#include "ws_tools.hpp"

using std::string;
using std::vector;

using namespace ws_img;
using namespace ws_tools;

typedef Vector::value_type value_type;

/**
	Allocate single, uninitialized, element vector.
	Note:
		Vector( value_type default_value = 0 ) : vec( 1, default_value )
		is ambiguous with Vector( size_type, value_type = 0 ), 
		e.g., Vector( 5 ) can reference either one.
 */
Vector::Vector( )
{ }

/**
	Allocate M-dimensional vector. 
	@param[in] M Number of elements in vector
	@param[in] default_value Default value to initialize each element to
 */
Vector::Vector( size_type M, value_type default_value )
: _vector( 1, M, default_value )
{ }

/**
	Allocate M-dimensional vector. 
	@param[in] M Number of elements in vector
	@param[in] default_value Default value to initialize each element to
 */
Vector::Vector( const_iterator vec_begin, const_iterator vec_end )
: _vector( 1, vec_end - vec_begin )
{
	iterator vec_iter = begin();
	for( ; vec_begin != vec_end; ++vec_begin, ++vec_iter )
	{
		*vec_iter = *vec_begin;
	}
}

/**
	Construct vector by reading given matrix file.
	@param[in] matrix_file_name Name of matrix file to read matrix from
 */
Vector::Vector( const std::string& matrix_file_name )
{
	read( matrix_file_name );
}

/**
	Construct vector by using the given range.

	For instance, passing 0, 2, and 10 yields the vector [ 0 2 4 6 8 10 ].

	Warning: passing 1, -.05, .55 only yields [ 1 .95 ... .60 ], i.e., the .55 is
	omitted. This is due to the inherent imprecision of floating-point math.
	TODO we may be able to rescale the values s.t. they're >= 1

	@param[in] begin First value in the range
	@param[in] inc Amount to increment each successive value
	@param[in] end Last value in the range
 */
Vector::Vector( const double begin, const double inc, const double end )
{
	// verify ranges are valid
	if( end < begin && inc >= 0 )
	{
		err_quit( "Vector constructor: "
				"end (%lf) < begin (%lf) && inc (%lf) >= 0\n", end, begin, inc );
	}
	else if( end > begin && inc <= 0 )
	{
		err_quit( "Vector constructor: "
				"end (%lf) > begin (%lf) && inc (%lf) <= 0\n", end, begin, inc );
	}
	else if( end == begin && inc != 0 )
	{
		err_quit( "Vector constructor: "
				"end (%lf) == begin (%lf) && inc (%lf) != 0\n", end, begin, inc );
	}

	/*
	// debugging precision problems by rescaling
	double min;
	min = std::min<double>( fabs(begin), fabs(inc) );
	min = std::min<double>( min, fabs(end) );
	double scale_factor = 1;
	if( min != 0 )
	{
		while( min < 1 )
		{
			scale_factor *= 10.0;
			min *= 10.0;
		}
	}
	 */

	// determine the number of elements in the vector--add 1 since it is
	// a few notes:
	// - we add +1 since we include both the values begin and end
	// - all signs will cancel since we check for a valid inc above
	// - for the case that inc does not evenly divide the
	// difference, we use floor() to decrease the count by 1
	size_type num_elements = 1;
	size_type i = vec_begin();
	if( inc > 0 )
	{
		num_elements += static_cast<size_type>(
				// round( (end - begin) / inc )
				round( floor( (end - begin) / inc ) )
		);

		// allocate space for this vector (resize without interpolation)
		resize( num_elements, false );

		// for( double val = begin; val <= end; ++i, val += inc )
		for( double val = begin; i != vec_end(); ++i, val += inc )
		{
			at(i) = static_cast<value_type>( val );
		}
	}
	else if( inc < 0 )
	{
		/*
		// debugging precision problems by rescaling
		const double rescaled_begin = scale_factor * begin;
		const double rescaled_inc   = scale_factor * inc;
		const double rescaled_end   = scale_factor * end;
		fprintf( stderr, "%lf %lf, %lf %lf\n",
				round( floor( (rescaled_end - rescaled_begin) / rescaled_inc ) ),
				round( (rescaled_end - rescaled_begin) / rescaled_inc ),
				round( floor( (end - begin) / inc ) ),
				round( (begin - end) / -inc )
		);
		num_elements += static_cast<size_type>(
				// round( floor( (rescaled_end - rescaled_begin) / rescaled_inc ) )
				round( floor( (rescaled_begin - rescaled_end) / -rescaled_inc ) )
				// round( (rescaled_end - rescaled_begin) / rescaled_inc )
		);
		 */

		num_elements += static_cast<size_type>(
				round( floor( (end - begin) / inc ) )
				// round( (begin - end) / -inc )
				// round( floor( (begin - end) / -inc ) )
		);

		// allocate space for this vector (resize without interpolation)
		resize( num_elements, false );

		// for( double val = begin; val >= end; ++i, val += inc )
		for( double val = begin; i != vec_end(); ++i, val += inc )
		{
			at(i) = static_cast<value_type>( val );
		}
	}
	else // inc == 0
	{
		// allocate space for this vector (resize without interpolation)
		resize( num_elements, false );

		at(i) = static_cast<value_type>( begin );
		++i;
	}

	if( i != vec_end() )
	{
		set_region( vec_begin(), sub_size() - 1 );
	}
}

/**
	Construct vector by copying elements from matrix.
	@param[in] M Matrix
 */
Vector::Vector( const Matrix& M )
: _vector( 1, M.sub_row_size() * M.sub_col_size() )
{
	// concatenate all rows of the matrix into a single vector
	size_type k = vec_begin(); 
	for( size_type i = M.row_begin(); i != M.row_end(); ++i )
	{
		for( size_type j = M.col_begin(); j != M.col_end(); ++j, ++k )
		{
			at(k) = M(i,j);
		}
	}

	/* // Version 2: matrix must be a row or column vector already
	typedef Matrix::const_col_iterator Col_Iter;

	// if matrix is a column vector
	if( M.sub_col_size() == 1 )
	{
		// _vector = Matrix( 1, M.sub_row_size() );
		_vector.resize( 1, M.sub_row_size(), false );

		// get iterators to vector's and matrix's column
		iterator iter       = begin();
		Col_Iter M_iter     = M.const_col_begin( M.col_begin() );
		Col_Iter M_end_iter = M.const_col_end( M.col_begin() );

		// copy each element of column
		for( ; M_iter != M_end_iter; ++M_iter, ++iter )
		{
			*iter = *M_iter;	
		}
	}

	// if matrix is a row vector
	else if( M.sub_row_size() == 1 )
	{
		// _vector = Matrix( 1, M.sub_col_size() );
		_vector = Matrix( 1, M.sub_col_size() );

		// get iterators to vector's and matrix's row
		iterator iter       = begin();
		Col_Iter M_iter     = M.const_row_begin( M.row_begin() );
		Col_Iter M_end_iter = M.const_row_end( M.row_begin() );

		// copy each element of row
		for( ; M_iter != M_end_iter; ++M_iter, ++iter )
		{
			*iter = *M_iter;	
		}
	}
	else // matrix is neither a row nor a column vector
	{
		err_quit( "Vector(const Matrix&): Matrix is not a vector\n" );
	}
	 */
}

/**
	Assign vector contents of matrix.
	@param[in] M Row or column matrix
 */
Vector&
Vector::operator=( const Matrix& M )
{
	// concatenate all rows of the matrix into a single vector
	_vector.resize( 1, M.sub_row_size() * M.sub_col_size(), true );
	size_type k = vec_begin(); 
	for( size_type i = M.row_begin(); i != M.row_end(); ++i )
	{
		for( size_type j = M.col_begin(); j != M.col_end(); ++j, ++k )
		{
			at(k) = M(i,j);
		}
	}

	/* // Version 2: matrix must be a row or column vector already
	typedef Matrix::const_col_iterator Col_Iter;

	// if matrix is a column vector
	if( M.sub_col_size() == 1 )
	{
		_vector = Matrix( 1, M.sub_row_size() );

		// get iterators to vector's and matrix's column
		iterator iter       = begin();
		Col_Iter M_iter     = M.const_col_begin( M.col_begin() );
		Col_Iter M_end_iter = M.const_col_end( M.col_begin() );

		// copy each element of column
		for( ; M_iter != M_end_iter; ++M_iter, ++iter )
		{
			*iter = *M_iter;	
		}
	}

	// if matrix is a row vector
	else if( M.sub_row_size() == 1 )
	{
		_vector = Matrix( 1, M.sub_col_size() );

		// get iterators to vector's and matrix's row
		iterator iter       = begin();
		Col_Iter M_iter     = M.const_row_begin( M.row_begin() );
		Col_Iter M_end_iter = M.const_row_end( M.row_begin() );

		// copy each element of row
		for( ; M_iter != M_end_iter; ++M_iter, ++iter )
		{
			*iter = *M_iter;	
		}
	}
	else // matrix is neither a row nor a column vector
	{
		err_quit( "Matrix is not a vector\n" );
	}
	*/

	return( *this );
}

/**
	Read vector as a MATLAB matrix from file matrix_file_name.
	@param[in] file_name Name of vector file to read vector from
 */
void
Vector::read( const std::string& file_name )
{
	_vector.read( file_name );

	// store vector as a matrix row
	if( _vector.real_col_size() == 1 )
	{
		_vector = _vector.transpose();
	}

	// verify that vector is a matrix
	if( _vector.real_row_size() != 1 )
	{
		err_quit( "Vector in file '%s' has invalid dimensions (%u x %u)\n",
			file_name.c_str(), _vector.real_row_size(),
			_vector.real_col_size() );
	}
}

/**
	Write contents of matrix to file file_name using specified precision for
	each value.
	@param[in] file_name Name of file to write matrix to
	@param[in] precision Precision to use for floating-point values
 */
void
Vector::write( const std::string& file_name, unsigned precision ) const
{
	_vector.write( file_name, precision );
}

/**
	Default destructor does nothing since no member variables are
	dynamically allocated.
 */
Vector::~Vector( )
{ }

Vector&
Vector::operator+=( const Vector& rhs )
{
	_vector += rhs._vector;
	return( *this );
}

Vector
Vector::operator+( const Vector& rhs ) const
{
	return( Vector( *this ) += rhs );
}

Vector&
Vector::operator-=( const Vector& rhs )
{
	_vector -= rhs._vector;
	return( *this );
}

Vector
Vector::operator-( const Vector& rhs ) const
{
	return( Vector( *this ) -= rhs );
}

/**
  	Multiply two vectors: this vector is a column vector and rhs is a row vector,
	which yields a matrix, e.g., w'v where w and v are column vectors.
	To perform element by element multiplication, see dot_product().
 */
Matrix
Vector::operator*( const Vector& rhs ) const
{
	// convert this vector into a column vector to multiply with the given
	// row vector
	return( Matrix(_vector) * rhs._vector );
}

/**
	Multiply this vector times a matrix, yielding a matrix,

	If we denote this vector by v and the given matrix by A, we compute vA only
	if A has a single row; otherwise, we compute v'A.

	@param[in] rhs Matrix on the right-hand side of this vector
	@retval result Result of operation
 */
Matrix
Vector::operator*( const Matrix& rhs ) const
{
	// if A is a row matrix, then do vA for column vector v
	if( rhs.sub_row_size() == 1 )
	{
		return( _vector.transpose() * rhs );
	}

	// perform v'A (we leave the vector as row vector/matrix for multiplication)
	return( _vector * rhs );
}

Vector&
Vector::operator+=( const value_type rhs )
{
	_vector += rhs;
	return( *this );
}

Vector&
Vector::operator-=( const value_type rhs )
{
	_vector -= rhs;
	return( *this );
}

Vector&
Vector::operator*=( const value_type rhs )
{
	_vector *= rhs;
	return( *this );
}

Vector
Vector::operator*( const value_type rhs ) const
{
	return( Vector( *this ) *= rhs );
}

Vector&
Vector::operator/=( const value_type rhs )
{
	_vector /= rhs;
	return( *this );
}

Vector
Vector::operator/( const value_type rhs ) const
{
	return( Vector( *this ) /= rhs );
}

Vector
Vector::operator-( ) const
{
	Vector neg;
	neg._vector      = -_vector;

	return( neg );
}

/**
	Compute sum of all elements in the vector.
	@retval vec_sum Vector sum
 */
value_type
Vector::sum( ) const
{
	// we must compute the sum directly instead of relying on the Matrix class's
	// sum() since Matrix computes a Vector of column sums and this vector is
	// actually stored as a row (i.e., multiple columns)
	value_type vec_sum = 0;
	for( size_type i = vec_begin(); i != vec_end(); ++i )
	{
		vec_sum += at(i);
	}
	return( vec_sum );
}

value_type
Vector::dot_product( const Vector& rhs ) const
{
	return( _vector.dot_product( rhs._vector ) );
}

/**
	Perform cross product between the two vectors.
 */
Vector
Vector::cross_product( const Vector& rhs ) const
{
	if( size() != 3 && size() != rhs.size() )
	{
		err_quit( "Vector::cross_product: size mismatch %u versus %u\n",
				size(), rhs.size() );
	}

	Vector cross_vec(3);
	cross_vec(0) = (at(vec_begin() + 1) * rhs(rhs.vec_begin() + 2))
			- (at(vec_begin() + 2) * rhs(rhs.vec_begin() + 1));
	cross_vec(1) = (at(vec_begin() + 2) * rhs(rhs.vec_begin() + 0))
			- (at(vec_begin() + 0) * rhs(rhs.vec_begin() + 2));
	cross_vec(2) = (at(vec_begin() + 0) * rhs(rhs.vec_begin() + 1))
			- (at(vec_begin() + 1) * rhs(rhs.vec_begin() + 0));

	// In simpler terms, the above is equivalent to this. This, however,
	// doesn't respect vector regions.
	// cross_vec(0) = (at(1) * rhs(2)) - (at(2) * rhs(1));
	// cross_vec(1) = (at(2) * rhs(0)) - (at(0) * rhs(2));
	// cross_vec(2) = (at(0) * rhs(1)) - (at(1) * rhs(0));

	return( cross_vec );
}

Vector&
Vector::sort( )
{
	_vector.sort();
	return( *this );
}

Vector&
Vector::reverse_sort( )
{
	_vector.reverse_sort();
	return( *this );
}

value_type
Vector::two_norm( ) const
{
	return( _vector.two_norm() );
}

value_type
Vector::two_norm_2( ) const
{
	return( _vector.two_norm_2() );
}

value_type
Vector::inf_norm( )
{
	return( _vector.inf_norm() );
}

void
Vector::normalize( )
{
	_vector.normalize();
}

/**
	Convolve this vector with right-hand side vector.

	Subregions are correctly handled.
	@param[in] rhs Right-hand side of convolution operation
	@retval y Output of convolution
 */
Vector
Vector::convolve( const Vector& rhs ) const
{
	if( sub_size() == 0 || rhs.sub_size() == 0 )
	{
		Vector empty;
		return( empty );
	}

	// perform convolution x*h where the length of x < length of h
	// order doesn't matter since convolution is commutative
	const Vector* x = 0;
	const Vector* h = 0;
	if( sub_size() > rhs.sub_size() )
	{
		x = &rhs;
		h = this;
	}
	else // sub_size() <= rhs.sub_size()
	{
		x = this;
		h = &rhs;
	}

	// y = x * h where * is the convolution operator
	Vector y( x->sub_size() + h->sub_size() - 1, 0 );

	// first section of convolution: not all of h is used
	for( size_type n = 0; n != x->sub_size() - 1; ++n )
	{
		for( size_type k = 0; k != (n + 1); ++k )
		{
			y(n) += x->at(x->vec_begin() + k) * h->at(h->vec_begin() + n - k);

			// fprintf( stderr, "Start %u (k): %lf <- %lf * %lf\n",
					// k, y(n), x->at(x->vec_begin() + k),
					// h->at(h->vec_begin() + n - k) );
		}
		// fprintf( stderr, "Start %u: %lf\n", n, y(n) );
	}

	// middle: all of x and h are used
	for( size_type n = x->sub_size() - 1; n != h->sub_size(); ++n )
	{
		for( size_type k = 0; k != x->sub_size(); ++k )
		{
			y(n) += x->at(x->vec_begin() + k) * h->at( h->vec_begin() + n - k);

			// fprintf( stderr, "Mid %u (k): %lf <- %lf * %lf\n",
					// k, y(n), x->at(x->vec_begin() + k),
					// h->at(h->vec_begin() + n - k) );
		}
		// fprintf( stderr, "Mid %u: %lf\n", n, y(n) );
	}

	// last section of convolution: not all of x is used
	size_type i = 0;
	for( size_type n = h->sub_size(); n != y.sub_size(); ++n )
	{
		++i;  // shift start of x one position to the right

		for( size_type k = i; k != x->sub_size(); ++k )
		{
			y(n) += x->at(x->vec_begin() + k) * h->at( h->vec_begin() + n - k);

			// fprintf( stderr, "End %u (k): %lf <- %lf * %lf\n",
					// k, y(n), x->at(x->vec_begin() + k),
					// h->at( h->vec_begin() + n - k) );
		}
		// fprintf( stderr, "End %u: %lf\n", n, y(n) );
	}

	// TODO Can't push_region since copying will only copy the current region
	// give the region of the output signal that had all values from both
	// signals x and h involved (i.e., no border values); we actually
	// just push a new region, so a user can pop it if they need the complete
	// output
	// const size_type filter_center = x->sub_size() / 2 + 1;
	// const size_type filter_size = h->sub_size() - x->sub_size() + 1;
	// y.push_region( filter_center, filter_size );
	
	return( y );
}

/**
	Convolve this vector with right-hand side vector.
	This version correctly does not handle vector subregions but is a little
	easier to understand as it follows the definition of convolution strictly.  
	@retval y Output of convolution
Vector
Vector::convolve( const Vector& rhs ) const
{
	if( size() == 0 || rhs.size() == 0 )
	{
		Vector empty;
		return( empty );
	}

	// perform convolution x*h where the length of x < length of h
	// order doesn't matter since convolution is commutative
	const Vector* x = 0;
	const Vector* h = 0;
	if( size() > rhs.size() )
	{
		x = &rhs;
		h = this;
	}
	else
	{
		x = this;
		h = &rhs;
	}

	// y = x * h where * is the convolution operator
	Vector y( x->size() + h->size() - 1 );

	// first section of convolution: not all of h is used
	for( size_type n = 0; n != x->size() - 1; ++n )
	{
		for( size_type k = 0; k != (n + 1); ++k )
		{
			y(n) += x->at(k) * h->at(n - k);
		}
		// fprintf( stderr, "Start %u: %lf\n", n, y(n) );
	}

	// middle: all of x and h are used
	for( size_type n = x->size() - 1; n != h->size(); ++n )
	{
		for( size_type k = 0; k != x->size(); ++k )
		{
			y(n) += x->at(k) * h->at(n - k);
		}
		// fprintf( stderr, "Mid %u: %lf\n", n, y(n) );
	}

	// last section of convolution: not all of x is used
	size_type i = 0;
	for( size_type n = h->size(); n != y.size(); ++n )
	{
		++i;  // shift start of x one position to the right

		for( size_type k = i; k != x->size(); ++k )
		{
			y(n) += x->at(k) * h->at(n - k);
		}
		// fprintf( stderr, "End %u: %lf\n", n, y(n) );
	}
	
	return( y );
}
 */

Vector
Vector::random( size_type size, value_type max_value )
{
	return( Vector( Matrix::random( 1, size, max_value ) ) );
}

std::vector<value_type>
Vector::to_vector( ) const
{
	std::vector<value_type> vec;

	const_iterator iter     = const_begin();
	const_iterator end_iter = const_end();
	for( ; iter != end_iter; ++iter )
	{
		vec.push_back( *iter );
	}

	return( vec );
}

double*
Vector::to_array( ) const
{
	return( _vector.to_array() );
}

Vector
Vector::get_row_vector( const Matrix& matrix, size_type row )
{
	Vector v( matrix.sub_col_size() );

	size_type i = v.vec_begin();
	for( size_type j = matrix.col_begin(); j != matrix.col_end(); ++i, ++j )
	{
		v(i) = matrix(row, j);
	}

	return( v );
}

Vector
Vector::get_col_vector( const Matrix& matrix, size_type col )
{
	Vector v( matrix.sub_col_size() );

	size_type j = v.vec_begin();
	for( size_type i = matrix.row_begin(); i != matrix.row_end(); ++i, ++j )
	{
		v(j) = matrix(i, col);
	}

	return( v );
}

void
Vector::set_row_vector( const Vector& vector, Matrix& matrix,
		size_type row )
{
	if( vector.size() != matrix.sub_col_size() )
	{
		return;
	}

	size_type i = vector.vec_begin();
	for( size_type j = matrix.col_begin(); j != matrix.col_end(); ++i, ++j )
	{
		matrix(row, j) = vector(i);
	}
}

void
Vector::set_col_vector( const Vector& vector, Matrix& matrix,
		size_type col )
{
	if( vector.size() != matrix.sub_row_size() )
	{
		return;
	}

	size_type j = vector.vec_begin();
	for( size_type i = matrix.row_begin(); i != matrix.row_end(); ++i, ++j )
	{
		matrix(i, col) = vector(j);
	}
}

/**
	Compute vector statistics--average, variance, and standard deviation.

	@param[out] avg Average value of all pixels in vector
	@param[out] var Variance of all pixels in vector
	@param[out] std Standard deviation of all pixels in vector
 */
void
Vector::stat( double* avg, double* var, double* std ) const
{
	_vector.stat( avg, var, std );
}

/**
	Find minimum value in vector.
	@retval max Minimum value in vector.
 */
value_type
Vector::get_min( ) const
{
	return( _vector.get_min() );
}

/**
	Find maximum value in vector.
	@retval max Maximum value in vector.
 */
value_type
Vector::get_max( ) const
{
	return( _vector.get_max() );
}

/**
	Find minimum and maximum value in vector.
	@retval max Maximum value in vector.
 */
void
Vector::get_min_max( value_type* min, value_type* max ) const
{
	_vector.get_min_max( min, max );
}

/**
	Compute the sine of each element in this vector.

	@retval v Result vector
 */
Vector
Vector::sin( ) const
{
	return( Vector( _vector.sin() ) );
}

/**
	Compute the cosine of each element in this vector.

	@retval v Result vector
 */
Vector
Vector::cos( ) const
{
	return( Vector( _vector.cos() ) );
}

/**
	Compute the tangent of each element in this vector.

	@retval v Result vector
 */
Vector
Vector::tan( ) const
{
	return( Vector( _vector.tan() ) );
}

/**
	Compute the absolute value of each element in this vector.

	@retval v Result vector
 */
Vector
Vector::abs( ) const
{
	return( Vector( _vector.abs() ) );
}

/**
	Apply fourier-transform to the this vector.

	The output is not normalized.

	@param[out] real_data The real components of the transform
	@param[out] imag_data The imaginary components of the transform
 */
void
Vector::fft( Vector& real_data, Vector& imag_data ) const
{
	_vector.fft( real_data._vector, imag_data._vector );
}

/**
	Apply inverse fourier-transform to the this vector.

	@param[in] real_data The real components of the transform
	@param[in] imag_data The imaginary components of the transform

	@post This vector holds the inverse of the normalized FFT.
 */
void
Vector::ifft( Vector& real_data, Vector& imag_data )
{
	_vector.ifft( real_data._vector, imag_data._vector );
}
