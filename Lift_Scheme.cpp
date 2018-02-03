/**
	@file   Lift_Scheme.cpp
	@author Wade Spires
	@date   2005/10/16
	@brief  Class Lift_Scheme.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Lift_Scheme.hpp"

using namespace ws_img;

using std::max;
using std::vector;

// define this to use matrix iterator
// #define USE_ITER

// typedef Matrix::Lift_Filter  Lift_Filter;
// typedef Matrix::Moments      Moments;
// typedef Matrix::Lifts        Lifts;

typedef Matrix::value_type           value_type;
typedef Lift_Scheme::size_type       size_type;
typedef Lift_Scheme::iterator        iterator;
typedef Lift_Scheme::const_iterator  const_iterator;

/**
	Apply lifting scheme to 1D signal.

	@param[in,out] signal
	@param[in] num_levels
 */
void
Lift_Scheme::lift( Vector& signal, unsigned num_levels )
{
	// if N or N_tilde is 1, then no corresponding biorthogonal wavelet exists,
	// so use Haar wavelet instead
	if( ws_tools::is_odd( get_N() ) == 1
			|| ws_tools::is_odd( get_N_tilde() ) == 1 )
	{
		lift_Haar( signal, num_levels );
		return;
	}

	// create initial lift and moment matrix sets
	// Lifts   lifts( signal.sub_size(), 1, get_N(), get_N_tilde(), num_levels );
	// Moments moments( signal.sub_size(), 1, get_N_tilde() );

	// calculate values for lifts and moments
	// calc_moments_lifts( signal.sub_size(), 1, num_levels, moments, lifts );

	// only recalculate lifts and moments if the matrix dimensions change
	// otherwise, reuse previous values
	if( _current_row_size != signal.sub_size() )
	{
		_current_row_size = signal.sub_size();
		_current_col_size = 1;

		// create initial lift and moment matrix sets
		// fprintf( stderr, "%u\n", get_N_tilde() );
		_lifts.set( signal.sub_size(), 1, get_N(), get_N_tilde(), num_levels );
		_moments.set( signal.sub_size(), 1, get_N_tilde() );

		_moments.X.write();

		// calculate values for lifts and moments
		calc_moments_lifts( signal.sub_size(), 1,
				num_levels, _moments, _lifts );
	}

	print_filter_moments_lifts( _filter, _moments, _lifts );

	// lift( signal, lifts, moments, num_levels );
	lift( signal, _lifts, _moments, num_levels );
}

/**
 */
void
Lift_Scheme::lift( Vector& signal, Lifts& lifts, Moments& moments,
	unsigned num_levels )
{
	// number of iterations
	num_levels = get_num_levels( num_levels, signal.sub_size(), get_N(),
			get_N_tilde() );

	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	// forward transform
	unsigned row = num_levels - 1;
	for( unsigned step = 1; step < max_step; step *= 2 )
	{
		Predict( signal.begin(), step, signal.sub_size() );

		Update( signal.begin(), step, lifts.X.const_row_begin(row),
			signal.sub_size() );

		--row;
	}
}

// bool DBG = false;

/**
	Apply lifting scheme to 2D signal.

	@param[in,out] signal
	@param[in] num_levels
 */
void
Lift_Scheme::lift( Matrix& signal, unsigned num_levels )
{
	// if N or N_tilde is 1, then there is no corresponding biorthogonal
	// wavelet, so use Haar wavelet instead
	// if( ws_tools::is_odd( get_N() ) == 1 || ws_tools::is_odd( get_N_tilde() ) == 1 )
	if( get_N() == 1 || get_N_tilde() == 1 )
	{
		lift_Haar( signal, num_levels );
		return;
	}

	// only recalculate lifts and moments if the matrix dimensions change
	// otherwise, reuse previous values
	if( _current_row_size != signal.sub_row_size()
			|| _current_col_size != signal.sub_col_size()
			|| _current_row_size == 0 || _current_col_size == 0 )
	{
		_current_row_size = signal.sub_row_size();
		_current_col_size = signal.sub_col_size();

		// create initial lift and moment matrix sets
		_lifts.set( signal.sub_col_size(), signal.sub_row_size(),
				get_N(), get_N_tilde(), num_levels );

		_moments.set( signal.sub_col_size(), signal.sub_row_size(),
				get_N_tilde() );

		// calculate values for lifts and moments
		calc_moments_lifts( signal.sub_col_size(), signal.sub_row_size(),
				num_levels, _moments, _lifts );

	// print_filter_moments_lifts( _filter, _moments, _lifts );
	}

	lift( signal, _lifts, _moments, num_levels );
}

/**
 */
void
Lift_Scheme::print_filter_moments_lifts( const Lift_Filter& filter,
	const Moments& moments, const Lifts& lifts )
{
	// print filter
	filter.write();
	fprintf( stdout, "\n" );

	// print updated moments
	fprintf( stdout, "X moment:\n" );
	moments.X.write();
	fprintf( stdout, "Y moment:\n" );
	moments.Y.write();

	fprintf( stdout, "\n" );

	// print calculated lifts
	fprintf( stdout, "X lifts:\n" );
	lifts.X.write( "", 5 );

	fprintf( stdout, "Y lifts:\n" );
	lifts.Y.write( "", 5 );

	fprintf( stdout, "\n" );
}

/**
	Apply lifting scheme to 2D signal.

	@param[in,out] signal
	@param[in] num_levels
void
Lift_Scheme::lift( Matrix& signal, unsigned num_levels )
{
	// if N or N_tilde is 1, then there is no corresponding biorthogonal
	// wavelet, so use Haar wavelet instead
	if( ws_tools::is_odd( get_N() ) == 1 || ws_tools::is_odd( get_N_tilde() ) == 1 )
	{
		lift_Haar( signal, num_levels );
		return;
	}

	// create initial lift and moment matrix sets
	Lifts    lifts( signal.sub_col_size(), signal.sub_row_size(), get_N(), get_N_tilde(),
			num_levels );
	Moments  moments( signal.sub_col_size(), signal.sub_row_size(), get_N_tilde() );

	// calculate values for lifts and moments
	calc_moments_lifts( signal.sub_col_size(), signal.sub_row_size(), num_levels,
			moments, lifts );

	lift( signal, lifts, moments, num_levels );
}
 */

/**
	Apply lifting scheme to 2D signal using given filter, lifts, and moments.

	@param[in,out] signal
	@param[in] num_levels
 */
void
Lift_Scheme::lift( Matrix& signal, Lifts& lifts, Moments& moments,
	unsigned num_levels )
{
	// get number of iterations in each direction (and total iterations)
	// note that to compute row_levels we use the # of columns (not rows) since
	// we filter _along_ each row for row_levels iterations
	unsigned row_levels = get_num_levels( num_levels, signal.sub_col_size(),
			get_N(), get_N_tilde() );

	unsigned col_levels = get_num_levels( num_levels, signal.sub_row_size(),
			get_N(), get_N_tilde() );

	// calculate the total number of iterations;
	// we set max_step = 2^num_levels since we use step sizes that are
	// progressively larger powers of 2
	num_levels = max<unsigned>( row_levels, col_levels );
	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	// for each step size separating s components at this level
	for( unsigned step = 1; step < max_step; step *= 2 )
	{
		// transform rows
		if( row_levels != 0 )
		{

// DBG = true;

			// transform each row separately
			for( unsigned i = signal.row_begin(); i < signal.row_end(); i += step )
			{
#ifdef USE_ITER
				Predict( signal.row_begin(i), step, signal.sub_col_size() );

				Update( signal.row_begin(i), step,
						lifts.X.const_row_begin(row_levels - 1), signal.sub_col_size()
				);
#else
				// get number of cases for each type
				const vector<unsigned> predict_cases
						= get_coef_cases( signal.sub_col_size(), step, get_N() );

				const vector<unsigned> update_cases
						= get_coef_cases( signal.sub_col_size(), step,
								get_N_tilde() );

				/*
if( DBG )
{
fprintf( stderr, "Lift row %u, step = %u, size = %u\n",
	i, step, signal.sub_col_size() );
}
				 */
				Predict_row( signal, i, step, predict_cases );

				Update_row( signal, i, step, lifts.X, row_levels - 1,
						update_cases );
// DBG = false;

#endif // USE_ITER
			}
			--row_levels;
		}
// DBG = false;

		// transform columns
		if( col_levels != 0 )
		{

			// transform each column separately
			for( unsigned j = signal.col_begin(); j < signal.col_end(); j += step )
			{
#ifdef USE_ITER

				Predict( signal.col_begin(j), step, signal.sub_row_size() );

				Update( signal.col_begin(j), step,
					lifts.Y.const_row_begin(col_levels - 1), signal.sub_row_size()
				);
#else
				// get number of cases for each type
				const vector<unsigned> predict_cases
						= get_coef_cases( signal.sub_row_size(), step, get_N() );

				const vector<unsigned> update_cases
						= get_coef_cases( signal.sub_row_size(), step,
								get_N_tilde() );

				Predict_col( signal, j, step, predict_cases );

				Update_col( signal, j, step, lifts.Y, col_levels - 1,
						update_cases );
#endif // USE_ITER
			}
			--col_levels;
		}
	}
}

/**
 */
void
Lift_Scheme::lift_Haar( Vector& signal, unsigned num_levels )
{
	// number of iterations
	num_levels = get_num_levels_Haar( num_levels, signal.sub_size() );
	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	for( unsigned step = 1; step < max_step; step *= 2 )
	{
		lift_Haar_1D( signal.begin(), step, signal.sub_size() );
	}
}

/**
	Lifting scheme.
	TODO: fix column transform to not create and copy vector

	@param[in,out] vec Image to lift
	@param[in] num_levels Number of lifting steps
 */
void
Lift_Scheme::lift_Haar( Matrix& signal, unsigned num_levels )
{
	// number of iterations: take image's size into account
	unsigned row_levels = get_num_levels_Haar( num_levels,
			signal.sub_col_size() );

	unsigned col_levels = get_num_levels_Haar( num_levels,
			signal.sub_row_size() );

	num_levels = max<unsigned>( row_levels, col_levels );

	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	// apply forward transformation
	for( unsigned step = 1; step < max_step; step *= 2 )
	{
		// transform rows
		if( row_levels != 0 )
		{
			// transform each row separately
			for( unsigned i = signal.row_begin(); i < signal.row_end(); i += step )
			{
				// use 1D Haar on row
				lift_Haar_1D( signal.row_begin(i), step, signal.sub_col_size() );
			}
			--row_levels;
		}

		// transform columns
		if( col_levels != 0 )
		{
			// transform each column separately
			for( unsigned j = signal.col_begin(); j < signal.col_end(); j += step )
			{
				// use 1D Haar on column
				lift_Haar_1D( signal.col_begin(j), step, signal.sub_row_size() );
			}
			--col_levels;
		}
	}
}

/**
	Find Haar and Lambda coefficients for the given vector. Boundary cases are
	not handled, so all signals must have even length.

	TODO: Does not seem to place integers on the border consistently (i.e., case
	of 18x1)
	Need to fix this.

	s_* => average signal
	d_* => difference signal

	@param[in,out] signal
	@param[in] step
 */
void
Lift_Scheme::lift_Haar_1D( iterator signal, unsigned step,
	unsigned size )
{
	// number of coefficients at this level
	unsigned num_coef = get_num_coef( size, step );

	// number of values between coefficients
	unsigned coef_offset = step * 2;

	// number of coefficients for averages and differences
	unsigned d_num = num_coef / 2;
	unsigned s_num = num_coef - d_num;

	// must calculate signal average to handle boundary case (since an odd number
	// of coefficients => last s has no d to its right)
	double s_average = 0;
	if( ws_tools::is_odd( num_coef ) )
	{
		unsigned offset = coef_offset / 2;

		// calculate average of all current coefficients at this level
		iterator signal_iter = signal;
		for( unsigned i = 0; i != num_coef; ++i, signal_iter += offset )
		{
			//fprintf( stderr, "i: %d, pos: %d\n", i, pos );
			assert( pos < size );
			s_average += *signal_iter;
		}
		s_average /= num_coef;
	}

	double s_sum = 0;  // average of s signal

	// starting positions in vector for difference and signal
	iterator avg_signal  = signal;
	iterator diff_signal = signal + step;

	// calculate each coefficient
	for( ; d_num != 0; --d_num )
	{
		// predict and update new coefficients
		*diff_signal -= *avg_signal;
		*avg_signal  += .5 * (*diff_signal);

		// update sum of s coefficients in case of boundary value
		s_sum += *avg_signal;

		// next coefficient positions
		diff_signal += coef_offset;
		avg_signal  += coef_offset;
	}

	// handle boundary case with different filter
	if( ws_tools::is_odd( num_coef ) )
	{
		*avg_signal = (s_average * double(s_num)) - s_sum;
		s_sum += *avg_signal;
	}
}

/**
	Apply inverse lifting scheme to signal.

	@param[in,out] signal
	@param[in] num_levels
 */
void
Lift_Scheme::inverse_lift( Vector& signal, unsigned num_levels )
{
	if( get_N() == 1 || get_N_tilde() == 1 )
	{
		inverse_lift_Haar( signal, num_levels );
		return;
	}

	// create initial lift and moment matrix sets
	// Lifts        lifts( signal.sub_size(), 1, get_N(), get_N_tilde(), num_levels );
	// Moments      moments( signal.sub_size(), 1, get_N_tilde() );

	// calculate values for lifts and moments
	// calc_moments_lifts( signal.sub_size(), 1, num_levels, moments, lifts );

	// only recalculate lifts and moments if the matrix dimensions change
	// otherwise, reuse previous values
	if( _current_row_size != signal.sub_size() )
	{
		_current_row_size = signal.sub_size();
		_current_col_size = 1;

		// create initial lift and moment matrix sets
		_lifts   = Lifts( signal.sub_size(), 1, get_N(), get_N_tilde(),
				num_levels );
		_moments = Moments( signal.sub_size(), 1, get_N_tilde() );

		// calculate values for lifts and moments
		calc_moments_lifts( signal.sub_size(), 1,
				num_levels, _moments, _lifts );
	}

	inverse_lift( signal, _lifts, _moments, num_levels );
}

/**
	Inverse lifting.

	@param[in,out] signal
	@param[in] lifts
	@param[in] moments
	@param[in] param
 */
void
Lift_Scheme::inverse_lift( Vector& signal, Lifts& lifts, Moments& moments,
	unsigned num_levels )
{
	// number of iterations
	num_levels = get_num_levels( num_levels, signal.sub_size(), get_N(),
			get_N_tilde() );
	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	// change signs of filter and lifts
	invert_filter_lift( lifts );

	// inverse transform
	unsigned row = num_levels - 1;
	for( unsigned step = max_step / 2; step >= 1; step /= 2 )
	{
		Update( signal.begin(), step, lifts.X.const_row_begin(row),
			signal.sub_size() );

		Predict( signal.begin(), step, signal.sub_size() );

		--row;
	}

	// revert signs of filter/lifts back
	invert_filter_lift( lifts );
}

/**
	Apply inverse lifting scheme to signal.

	@param[in,out] signal
	@param[in] num_levels
 */
void
Lift_Scheme::inverse_lift( Matrix& signal, unsigned num_levels )
{
	if( get_N() == 1 || get_N_tilde() == 1 )
	// if( ws_tools::is_odd( get_N() ) == 1 || ws_tools::is_odd( get_N_tilde() ) == 1 )
	{
		inverse_lift_Haar( signal, num_levels );
		return;
	}

	// only recalculate lifts and moments if the matrix dimensions change
	// otherwise, reuse previous values
	if( _current_row_size != signal.sub_row_size()
			|| _current_col_size != signal.sub_col_size()
			|| _current_row_size == 0 || _current_col_size == 0 )
	{
		// save size of signal, so we can avoid recalculation possibly
		_current_row_size = signal.sub_row_size();
		_current_col_size = signal.sub_col_size();

		// create initial lift and moment matrix sets
		_lifts = Lifts( signal.sub_col_size(), signal.sub_row_size(),
				get_N(), get_N_tilde(), num_levels);

		_moments = Moments( signal.sub_col_size(), signal.sub_row_size(),
				get_N_tilde() );

		// calculate values for lifts and moments
		calc_moments_lifts( signal.sub_col_size(), signal.sub_row_size(),
				num_levels, _moments, _lifts );
	}

	// print_filter_moments_lifts( _filter, _moments, _lifts );

	inverse_lift( signal, _lifts, _moments, num_levels );
}

/**
	Apply inverse lifting scheme to signal.

	@param[in,out] signal
	@param[in] num_levels
void
Lift_Scheme::inverse_lift( Matrix& signal, unsigned num_levels )
{
	if( get_N() == 1 || get_N_tilde() == 1 )
	// if( ws_tools::is_odd( get_N() ) == 1 || ws_tools::is_odd( get_N_tilde() ) == 1 )
	{
		inverse_lift_Haar( signal, num_levels );
		return;
	}

	// create initial lift and moment matrix sets
	Lifts       lifts( signal.sub_col_size(), signal.sub_row_size(), get_N(), get_N_tilde(),
		num_levels);
	Moments     moments( signal.sub_col_size(), signal.sub_row_size(), get_N_tilde() );

	// calculate values for lifts and moments
	calc_moments_lifts( signal.sub_col_size(), signal.sub_row_size(), num_levels,
			moments, lifts );

	inverse_lift( signal, lifts, moments, num_levels );
}
 */

/**
	Inverse lifting.

	@param[in,out] signal
	@param[in] lifts
	@param[in] moments
	@param[in] num_levels Number of lifting steps
 */
void
Lift_Scheme::inverse_lift( Matrix& signal, Lifts& lifts, Moments& moments,
	unsigned num_levels )
{
	// get number of iterations in each direction (and total iterations)
	// note that to compute row_levels we use the # of columns (not rows) since
	// we filter _along_ each row for row_levels iterations
	unsigned row_levels = get_num_levels( num_levels, signal.sub_col_size(),
			get_N(), get_N_tilde());

	unsigned col_levels = get_num_levels( num_levels, signal.sub_row_size(),
			get_N(), get_N_tilde() );

	num_levels = max<unsigned>( row_levels, col_levels );
	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	// change signs of filter and lifts
	invert_filter_lift( lifts );

	// for each step size separating s components at this level
	for( unsigned step = max_step / 2; step >= 1; step /= 2 )
	{
		// transform columns
		if( num_levels <= col_levels )
		{
			// get number of cases for each type
			const vector<unsigned> update_cases
					= get_coef_cases( signal.sub_row_size(), step, get_N_tilde() );

			const vector<unsigned> predict_cases
					= get_coef_cases( signal.sub_row_size(), step, get_N() );

			// transform each column separately
			for( unsigned j = signal.col_begin(); j < signal.col_end();
					j += step )
			{
#ifdef USE_ITER
				Update( signal.col_begin(j), step,
					lifts.Y.const_row_begin(col_levels - num_levels),
					signal.sub_row_size() );

				Predict( signal.col_begin(j), step, signal.sub_row_size() );
#else
				Update_col( signal, j, step, lifts.Y, col_levels - num_levels,
						update_cases );
				// Update( signal.col_begin(j), step,
					// lifts.Y.const_row_begin(col_levels - num_levels),
					// signal.sub_row_size() );

				Predict_col( signal, j, step, predict_cases );
#endif // USE_ITER
			}
		}

		// transform rows
		if( num_levels <= row_levels )
		{
			// get number of cases for each type
			const vector<unsigned> update_cases
					= get_coef_cases( signal.sub_col_size(), step, get_N_tilde() );

			const vector<unsigned> predict_cases
					= get_coef_cases( signal.sub_col_size(), step, get_N() );

			// transform each row separately
			for( unsigned i = signal.row_begin(); i < signal.row_end(); i += step )
			{
#ifdef USE_ITER
				Update( signal.row_begin(i), step,
					lifts.X.const_row_begin(row_levels - num_levels),
					signal.sub_col_size() );

				Predict( signal.row_begin(i), step, signal.sub_col_size() );
#else
				Update_row( signal, i, step, lifts.X, row_levels - num_levels,
						update_cases );

				Predict_row( signal, i, step, predict_cases );
#endif // USE_ITER
			}
		}

		// move up one level
		--num_levels;
	}

	// revert signs of filter/lifts back
	invert_filter_lift( lifts );
}

/**
	Inverse lifting scheme.
	@param[in,out] signal Vector to apply inverse lifting
	@param[in] num_levels Number of lifting steps
 */
void
Lift_Scheme::inverse_lift_Haar( Vector& signal, unsigned num_levels )
{
	// number of iterations
	num_levels = get_num_levels_Haar( num_levels, signal.sub_size() );

	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	for( unsigned step = max_step / 2; step >= 1; step /= 2 )
	{
		inverse_lift_Haar_1D( signal.begin(), step, signal.sub_size() );
	}
}

/**
	Inverse lifting scheme.
	@param[in,out] signal Signal to apply inverse lifting
	@param[in] num_levels Number of lifting steps
 */
void
Lift_Scheme::inverse_lift_Haar( Matrix& signal, unsigned num_levels )
{
	// number of iterations: take image's size into account
	unsigned row_levels = get_num_levels_Haar( num_levels,
			signal.sub_col_size() );

	unsigned col_levels = get_num_levels_Haar( num_levels,
			signal.sub_row_size() );

	num_levels = max<unsigned>( row_levels, col_levels );

	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	// apply forward transformation
	//for( unsigned step = max_step / 2; step >= 1; step /= 2 )
	for( unsigned step = max_step / 2; step >= 1; step /= 2 )
	{
		// transform columns
		if( num_levels <= col_levels )
		{
			// transform each column separately
			for( unsigned j = signal.col_begin(); j < signal.col_end(); j += step )
			{
				// use 1D Haar on column
				inverse_lift_Haar_1D( signal.col_begin(j), step,
					signal.sub_row_size() );
			}
		}

		// transform rows
		if( num_levels <= row_levels )
		{
			// transform each row separately
			for( unsigned i = signal.row_begin(); i < signal.row_end(); i += step )
			{
				// use 1D Haar on row
				inverse_lift_Haar_1D( signal.row_begin(i), step,
					signal.sub_col_size() );
			}
		}

		// move up one level
		--num_levels;
	}
}

/**
	Find Haar and Lambda coefficients for the given vector. Boundary cases are
	not handled, so all signals must have even length.

	TODO: Only works for signals with length that are power of 2: 4, 8, 16, ...
	Need to fix this.
	--Appears to work now, but should probably test a lot

	Note: Can replace "* 2" and "/ 2" with "<< 1" and ">> 1" for speed.

	s_* => average signal
	d_* => difference signal

	@param[in,out] signal
	@param[in] step
 */
void
Lift_Scheme::inverse_lift_Haar_1D( iterator signal, unsigned step,
	unsigned size )
{
	// number of coefficients at this level
	unsigned num_coef = get_num_coef( size, step );

	// number of values between coefficients
	unsigned coef_offset = step * 2;

	// number of coefficients for d and s
	unsigned d_num = num_coef / 2;
	unsigned s_num = num_coef - d_num;

	// must calculate signal average to handle boundary case (since an odd number
	// of coefficients => last s has no d to its right)
	double s_average = 0;
	if( ws_tools::is_odd( num_coef ) )
	{
		unsigned offset = coef_offset;

		// calculate average of only s coefficients
		iterator signal_iter = signal;
		for( unsigned i = 0; i != s_num; ++i, signal_iter += offset )
		{
			s_average += *signal_iter;
		}
		s_average /= s_num;
	}

	double ds_sum = 0;  // average of s signal

	// starting positions in vector for difference and signal
	iterator avg_signal  = signal;
	iterator diff_signal = signal + (coef_offset / 2);

	// calculate each coefficient
	for( ; d_num != 0; --d_num )
	{
		// update and predict new coefficients
		*avg_signal  -= .5 * (*diff_signal);
		*diff_signal += *avg_signal;

		// update sum of s and d signal in case of boundary value
		ds_sum += *avg_signal;
		ds_sum += *diff_signal;

		// next coefficient positions
		diff_signal += coef_offset;
		avg_signal  += coef_offset;
	}

	// handle boundary case with different filter
	if( ws_tools::is_odd( num_coef ) )
	{
		*avg_signal = (s_average * double(num_coef)) - ds_sum;
		ds_sum += *avg_signal;
	}
}

/**
	Calculate the moments and lifting coefficients to preserve the moments of the
	wavelet.

	@param[in] width
	@param[in] height
	@param[in] num_levels
	@param[in,out] moments
	@param[in,out] lifts
 */
void
Lift_Scheme::calc_moments_lifts( unsigned width, unsigned height,
	unsigned num_levels, Moments& moments, Lifts& lifts )
{
	// maximum vanishing moments
	unsigned max_N = max<unsigned>(get_N(), get_N_tilde()) - 1;
	unsigned row_levels = get_num_levels( num_levels, (width - 1) / max_N );
	unsigned col_levels = get_num_levels( num_levels, (height - 1) / max_N );
	num_levels = max<unsigned>( row_levels, col_levels );

	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	// find lifting coefficients for each level
	for( unsigned step = 1; step < max_step; step *= 2 )
	{
		// calculate moments and lifts in x direction
		if( row_levels != 0 )
		{
			calc_moment( moments.X, width, step );

			calc_lift( lifts.X.row_begin(row_levels - 1),
					moments.X, width, step );

			--row_levels;
		}

		// calculate moments and lifts in y direction
		if( col_levels != 0 )
		{
			calc_moment( moments.Y, height, step );

			calc_lift( lifts.Y.row_begin(col_levels - 1),
				moments.Y, height, step );

			--col_levels;
		}
	}
}

/**
	Calculate the integral moment tuple for the current level of calculations.

	@param[in,out] moment
	@param[in] size
	@param[in] step
 */
void
Lift_Scheme::calc_moment( Matrix& moment, unsigned size, unsigned step )
{
	assert( get_N() == _filter.sub_col_size() );

	vector<unsigned> cases = get_coef_cases( size, step, get_N() );
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	// distance between next successive s's (or next d's)
	unsigned coef_offset = 2 * step;

	// case 1: # left s < # right s
	unsigned filter_row    = 1;     // skip first row's coefficients
	unsigned diff_signal   = step;  // position of first d
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		unsigned avg_signal = 0;  // first s is the first coefficient

		for( unsigned j = 0; j != get_N(); ++j )
		{
			// update (int, mom_1, mom_2, ...)
			for( unsigned k = 0; k != get_N_tilde(); ++k )
			{
				assert( avg_signal >= 0 && avg_signal < moment.sub_row_size() );
				/*
				fprintf( stderr, "%d\n", avg_signal );
				if( avg_signal < 0 || avg_signal >= moment.sub_row_size() )
				{
					fprintf( stderr, "1 avg_signal: %d i(teration): %d row_size: %d ",
						avg_signal, i, moment.sub_row_size() );
					fprintf( stderr, "step: %d\n", coef_offset );
				}
				 */

				moment(avg_signal, k) +=
					_filter(filter_row, j) * moment(diff_signal, k);
			}
			avg_signal += coef_offset;  // next s coefficient position
		}

		++filter_row;                // next filter row
		diff_signal += coef_offset;  // next d coefficient position
	}

	// case 2: # left s = # right s
	// last filter row (and stay in this row for iteration)
	filter_row = _filter.sub_row_size() - 1;
	unsigned start_pos = 0;
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		//unsigned avg_signal = i * coef_offset;
		unsigned avg_signal = start_pos;

		for( unsigned j = 0; j != get_N(); ++j )
		{
			// update (int, mom_1, mom_2, ...)
			for( unsigned k = 0; k != get_N_tilde(); ++k )
			{
				assert( avg_signal >= 0 && avg_signal < moment.sub_row_size() );

				moment(avg_signal, k) +=
					_filter(filter_row, j) * moment(diff_signal, k);
			}
			avg_signal += coef_offset;  // next s coefficient position
		}

		start_pos   += coef_offset;
		diff_signal += coef_offset;  // next d coefficient position (in same row)
	}

	// case 3: # left s > # right s
	filter_row = _filter.sub_row_size() - 2;  // start at next to last

	// skip over equals to get to greater thans
	unsigned saved_pos = start_pos - coef_offset;
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		unsigned avg_signal = saved_pos;

		for( int j = (get_N() - 1); j >= 0; --j )
		{
			for( unsigned k = 0; k != get_N_tilde(); ++k )
			{
				assert( avg_signal >= 0 && avg_signal < moment.sub_row_size() );

				moment(avg_signal, k) +=
					_filter(filter_row, j) * moment(diff_signal, k);
			}
			avg_signal += coef_offset;  // next s coefficient position
		}

		diff_signal += coef_offset;  // next d coefficient position
		--filter_row;                // previous filter row
	}
}

/**
	@param[in,out] lift_vec
	@param[in] moment
	@param[in] size
	@param[in] step
 */
void
Lift_Scheme::calc_lift( iterator lift_vec, Matrix& moment, unsigned size,
	unsigned step )
{
	assert( get_N_tilde() == moment.sub_col_size() );

	vector<unsigned> cases = get_coef_cases( size, step, get_N_tilde() );
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	// distance between next successive s's (or next d's)
	unsigned coef_offset = 2 * step;

	Matrix tmp_lift( get_N_tilde(), get_N_tilde() );
	Vector b( get_N_tilde() );

	// case 1: # left s < # right s
	unsigned diff_signal = step;  // position of first d
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		unsigned avg_signal = 0;  // first s is the first coefficient

		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			// update (int, mom_1, mom_2, ...)
			for( unsigned k = 0; k != get_N_tilde(); ++k )
			{
				assert( avg_signal >= 0 && avg_signal < moment.sub_row_size() );
			/*

				if( 0 || avg_signal >= moment.sub_row_size() )
				{
					fprintf( stderr, "1 avg_signal: %d, mom.size: %d\n",
						avg_signal, moment.sub_row_size() );
				}

				if( k > tmp_lift.sub_row_size() )
				{
					fprintf( stderr, "1 k: %d > row_size: %d\n",
						k, tmp_lift.sub_row_size() );
				}

				if( j > tmp_lift.sub_col_size() )
				{
					fprintf( stderr, "1 j: %d > col_size: %d\n",
						j, tmp_lift.sub_col_size() );
				}
			 */

				tmp_lift(k, j) = moment(avg_signal, k);
			}
			avg_signal += coef_offset;  // next s coefficient position
		}

		// solve system using row at given d position in moment
		Vector x( tmp_lift.sub_col_size() );
		// Lin_Alg::solve_system( tmp_lift, moment(diff_signal), x );

		Vector v( moment.const_row_begin(diff_signal),
			moment.const_row_end(diff_signal) );

		Lin_Alg::solve_system( tmp_lift, v, x );
		//Vector x = moment(diff_signal);
		//Lin_Alg::solve_system( tmp_lift, x );

		// copy solution into lift_vec
		/*
		for( unsigned j = x.sub_size(); j != x.sub_size(); ++j )
		{
			*lift_vec = x(j);
		}
		 */
		for( unsigned j = x.vec_begin(); j != x.vec_end(); ++j )
		{
			*lift_vec = x(j);
			++lift_vec;
		}

		diff_signal += coef_offset;  // next d coefficient position
	}

	// case 2: # left s = # right s
	unsigned start_pos = 0;
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		//unsigned avg_signal = i * coef_offset;
		unsigned avg_signal = start_pos;

		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			// update (int, mom_1, mom_2, ...)
			for( unsigned k = 0; k != get_N_tilde(); ++k )
			{
				assert( avg_signal >= 0 && avg_signal < moment.sub_row_size() );

				// should we maybe offset k by offset between averages? 
				tmp_lift(k, j) = moment(avg_signal, k);
			}
			avg_signal += coef_offset;  // next s coefficient position
		}

		// solve system using column at given d position in moment
		Vector x( tmp_lift.sub_row_size() );
		// Lin_Alg::solve_system( tmp_lift, moment(diff_signal), x );
		Vector v( moment.const_row_begin(diff_signal),
			moment.const_row_end(diff_signal) );

		Lin_Alg::solve_system( tmp_lift, v, x );
		//Vector x = moment(diff_signal);
		//Lin_Alg::solve_system( tmp_lift, x );

		// copy solution into lift_vec
		for( unsigned j = x.vec_begin(); j != x.vec_end(); ++j )
		{
			*lift_vec = x(j);
			++lift_vec;
		}

		start_pos   += coef_offset;
		diff_signal += coef_offset;  // next d coefficient position
	}

	// case 3: # left s > # right s
	unsigned saved_pos = start_pos - coef_offset;
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		//unsigned avg_signal = (num_s_equal - 1) * coef_offset;
		unsigned avg_signal = saved_pos;  // return s to the saved starting position

		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			// update (int, mom_1, mom_2, ...)
			for( unsigned k = 0; k != get_N_tilde(); ++k )
			{
				assert( avg_signal >= 0 && avg_signal < moment.sub_row_size() );

				tmp_lift(k, j) = moment(avg_signal, k);
			}
			avg_signal += coef_offset;  // next s coefficient position
		}

		// solve system using row at given d position in moment
		Vector x( tmp_lift.sub_col_size() );
		// Lin_Alg::solve_system( tmp_lift, moment(diff_signal), x );
		Vector v( moment.const_row_begin(diff_signal),
				moment.const_row_end(diff_signal) );

		Lin_Alg::solve_system( tmp_lift, v, x );
		//Vector x = moment(diff_signal);
		//Lin_Alg::solve_system( tmp_lift, x );

		// copy solution into lift_vec
		/*
		for( unsigned j = x.sub_size(); j != x.sub_size(); ++j )
		{
			*lift_vec = x(j);
		}
		 */
		// copy solution into lift_vec
		for( unsigned j = x.vec_begin(); j != x.vec_end(); ++j )
		{
			*lift_vec = x(j);
			++lift_vec;
		}

		diff_signal += coef_offset;  // next d coefficient position
	}
}

/**
	Predict difference coefficients d. Overwrite signal.

	@param[in,out] signal
	@param[in] step
 */
void
Lift_Scheme::Predict( iterator signal, unsigned step, size_type size )
{
	// vector<unsigned> cases = get_coef_cases( signal.sub_size(), step, get_N() );
	vector<unsigned> cases = get_coef_cases( size, step, get_N() );
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	// distance between next successive s's (or next d's)
	unsigned coef_offset = 2 * step;

	unsigned filter_row = 1; // start at row 1 (skip row 0)

	iterator diff_signal = signal + step;  // position of d's
	iterator avg_signal  = signal;         // position of s's

	// case 1: # left s < # right s
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		// avg_signal = 0;
		avg_signal = signal;

		// d = d - prediction
		for( unsigned j = 0; j != get_N(); ++j )
		{
			(*diff_signal) -= (*avg_signal) * _filter(filter_row, j );
			avg_signal += coef_offset;
		}

		diff_signal += coef_offset;
		++filter_row;
	}

	// TODO make separate function where we unroll the inner 'for' loop for
	// various values of N to optimize since this is the inner-most computation
	// TODO We could also use template meta-programming by using N and N~ as
	// template parameters so that this happens at compile time.

	// case 2: # left s = # right s
	filter_row = _filter.sub_row_size() - 1;  // use row for equal #s
	iterator start_pos = signal;
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		avg_signal = start_pos;  // return to

		// d = d - prediction
		for( unsigned j = 0; j != get_N(); ++j )
		{
			(*diff_signal) -= (*avg_signal) * _filter(filter_row, j);
			avg_signal += coef_offset;
		}

		diff_signal += coef_offset;
		start_pos   += coef_offset;  // change starting position of s
	}

	// case 3: # left s > # right s
	// first row with unequal # of s's (relative to end)
	filter_row = _filter.sub_row_size() - 2;
	iterator saved_pos = (start_pos - coef_offset);
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		avg_signal = saved_pos;  // position to return for s

		// traverse opposite order
		// d = d - prediction
		for( int j = (get_N() - 1); j >= 0; --j )
		{
			(*diff_signal) -= (*avg_signal) * _filter(filter_row, j);
			avg_signal += coef_offset;
		}
		--filter_row; // reverse the order of the row traversal

		diff_signal += coef_offset;
	}
}

/**
	Predict difference coefficients d. Overwrite signal.

	TODO we can speed this up even more by replacing loops that involve get_N()
	to precomputed versions for specific Ns, or using template meta-programming.

	@param[in,out] signal
	@param[in] step
	@param[in] cases The number of cases to examine, which is the result to a
		call to get_coef_cases().
 */
void
Lift_Scheme::Predict_row( Matrix& signal, const size_type row_pos,
		const unsigned step, const vector<unsigned>& cases )
{
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	// distance between next successive s's (or next d's)
	unsigned coef_offset = 2 * step;

	unsigned filter_row = 1; // start at row 1 (skip row 0)

	size_type diff_pos = signal.col_begin() + step;  // position of d's
	size_type avg_pos  = signal.col_begin();         // position of s's

	// case 1: # left s < # right s
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		avg_pos = signal.col_begin();

		// d = d - prediction
		for( unsigned j = 0; j != get_N(); ++j )
		{
			signal(row_pos, diff_pos) -=
					signal(row_pos, avg_pos) * _filter(filter_row, j);

			avg_pos += coef_offset;
		}

		diff_pos += coef_offset;
		++filter_row;
	}

	// case 2: # left s = # right s
	filter_row = _filter.sub_row_size() - 1;  // use row for equal #s
	size_type start_pos = signal.col_begin();
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		avg_pos = start_pos;  // return to starting position

		// if( get_N() == 4 )
		{
			/*
			value_type* filter = new value_type[ _filter.col_size() ];
			size_type jj = 0; 
			for( size_type j = _filter.col_begin();
					j != _filter.col_end();
					++j, ++jj )
			{
				filter[jj] = _filter(filter_row, j);
			}

			signal(row_pos, diff_pos) -=
					signal(row_pos, avg_pos) * filter[0]
					+ signal(row_pos, avg_pos + coef_offset) * filter[1]
					+ signal(row_pos, avg_pos + 2 * coef_offset) * filter[2]
					+ signal(row_pos, avg_pos + 3 * coef_offset) * filter[3];

			delete [] filter;
			 */
		}
		// else
		{
			// d = d - prediction
			for( unsigned j = 0; j != get_N(); ++j )
			{
				signal(row_pos, diff_pos) -=
						signal(row_pos, avg_pos) * _filter(filter_row, j);

				avg_pos += coef_offset;
			}
		}

		diff_pos  += coef_offset;
		start_pos += coef_offset;  // change starting position of s
	}

	// case 3: # left s > # right s
	// first row with unequal # of s's (relative to end)
	filter_row = _filter.sub_row_size() - 2;
	size_type saved_pos = (start_pos - coef_offset);
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		avg_pos = saved_pos;  // position to return for s

		// traverse opposite order
		// d = d - prediction
		for( int j = (get_N() - 1); j >= 0; --j )
		{
			signal(row_pos, diff_pos) -=
					signal(row_pos, avg_pos) * _filter(filter_row, j);

			avg_pos += coef_offset;
		}
		--filter_row; // reverse the order of the row traversal

		diff_pos += coef_offset;
	}
}

/**
	Predict difference coefficients d. Overwrite signal.

	@param[in,out] signal
	@param[in] step
 */
void
Lift_Scheme::Predict_col( Matrix& signal, const size_type col_pos,
		const unsigned step, const vector<unsigned>& cases )
{
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	// distance between next successive s's (or next d's)
	unsigned coef_offset = 2 * step;

	unsigned filter_row = 1; // start at row 1 (skip row 0)

	size_type diff_pos = signal.row_begin() + step;  // position of d's
	size_type avg_pos  = signal.row_begin();  // position of s's

	// case 1: # left s < # right s
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		avg_pos = signal.row_begin();

		// d = d - prediction
		for( unsigned j = 0; j != get_N(); ++j )
		{
			signal(diff_pos, col_pos) -=
					signal(avg_pos, col_pos) * _filter(filter_row, j);

			avg_pos += coef_offset;
		}

		diff_pos += coef_offset;
		++filter_row;
	}

	// case 2: # left s = # right s
	filter_row = _filter.sub_row_size() - 1;  // use row for equal #s
	size_type start_pos = signal.row_begin();
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		avg_pos = start_pos;  // return to starting position

		// d = d - prediction
		for( unsigned j = 0; j != get_N(); ++j )
		{
			signal(diff_pos, col_pos) -=
					signal(avg_pos, col_pos) * _filter(filter_row, j);

			avg_pos += coef_offset;
		}

		diff_pos  += coef_offset;
		start_pos += coef_offset;  // change starting position of s
	}

	// case 3: # left s > # right s
	// first row with unequal # of s's (relative to end)
	filter_row = _filter.sub_row_size() - 2;
	size_type saved_pos = (start_pos - coef_offset);
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		avg_pos = saved_pos;  // position to return for s

		// traverse opposite order
		// d = d - prediction
		for( int j = (get_N() - 1); j >= 0; --j )
		{
			signal(diff_pos, col_pos) -=
					signal(avg_pos, col_pos) * _filter(filter_row, j);

			avg_pos += coef_offset;
		}
		--filter_row; // reverse the order of the row traversal

		diff_pos += coef_offset;
	}
}

/**
	Update stage.

	@param[in,out] signal
	@param[in] step
	@param[in] lift_vec
 */
void
Lift_Scheme::Update( iterator signal, unsigned step,
	const_iterator lift_vec, size_type size )
{
	vector<unsigned> cases = get_coef_cases( size, step, get_N_tilde() );
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	// distance between next successive s's (or next d's)
	unsigned coef_offset = 2 * step;

	iterator diff_signal = signal + step;  // position of d's
	iterator avg_signal  = signal;         // position of s's

	// case 1: # left s < # right s
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		// avg_signal = 0;  // first s is always the first element of signal
		avg_signal = signal;  // first s is always the first element of signal

		// s = d + update
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			*avg_signal += (*diff_signal) * (*lift_vec);

			avg_signal += coef_offset;
			++lift_vec;
		}
		diff_signal += coef_offset;
	}

	// case 2: # left s = # right s
	iterator start_pos = signal;
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		avg_signal = start_pos;  // return to

		// s = d + update
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			*avg_signal += (*diff_signal) * (*lift_vec);

			avg_signal += coef_offset;
			++lift_vec;
		}

		diff_signal += coef_offset;
		start_pos   += coef_offset;
	}

	// case 3: # left s > # right s
	iterator saved_pos = start_pos - coef_offset;
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		avg_signal = saved_pos;  // return s to the saved starting position

		// s = d + update
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			*(avg_signal) += (*diff_signal) * (*lift_vec);

			avg_signal += coef_offset;
			++lift_vec;
		}

		diff_signal += coef_offset;
	}
}

/**
	Update stage for rows.

	The row row_pos of the signal is updated. The lifting coefficients in
	lift_coef at row lift_row are used. The given step size separates values.

	@param[in,out] signal Signal
	@param[in] row_pos Row position
	@param[in] step Step-size between values--the number of values separating
		the average and signal coefficients. Successive values for the average
		part of the signal are spaced 2 * step apart.  
	@param[in] lift_coef Lifting coefficients
	@param[in] lift_row Row position of lifting coefficients
 */
void
Lift_Scheme::Update_row( Matrix& signal, const size_type row_pos,
		const unsigned step,
		const Matrix& lift_coef, const size_type lift_row,
		const vector<unsigned>& cases )
{
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	/*
if( DBG )
{
fprintf( stderr, "\n" );
fprintf( stderr, "  Update\n" );
}
	 */

	// distance between next successive s's (or next d's)
	const unsigned coef_offset = 2 * step;

	size_type diff_pos = signal.col_begin() + step;  // position of d's
	size_type avg_pos  = signal.col_begin();         // position of s's
	size_type lift_col = lift_coef.col_begin();      // position of lifting coef.

	/*
if( DBG )
{
	fprintf( stderr, "    Num Less: %u\n", cases[0] );
}
	 */

	// case 1: # left s < # right s
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		// first s is always the first element of signal
		avg_pos = signal.col_begin();

		// s = d + update
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{

			/*
if( DBG )
{
	const double new_val = signal(row_pos, avg_pos)
			+ signal(row_pos, diff_pos) * lift_coef(lift_row, lift_col);

	fprintf( stderr, "    signal(%u,%u) += "
			"signal(%u,%u) * lift_coef(%u,%u)\n",
			row_pos, avg_pos, row_pos, diff_pos, lift_row, lift_col );
	fprintf( stderr, "      => %lf += %lf * %lf => %lf\n",
			signal(row_pos, avg_pos),
			signal(row_pos, diff_pos), lift_coef(lift_row, lift_col),
			new_val );
}
			 */

			signal(row_pos, avg_pos) +=
					signal(row_pos, diff_pos) * lift_coef(lift_row, lift_col);

			avg_pos += coef_offset;
			++lift_col;
		}
/*
if( DBG )
{
	fprintf( stderr, "\n" );
}
 */
		diff_pos += coef_offset;
	}

/*
if( DBG )
{
	fprintf( stderr, "    Num Equal: %u\n", cases[1] );
}
 */

	// case 2: # left s = # right s
	size_type start_pos = signal.col_begin();
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		avg_pos = start_pos;  // return to starting position

		// s = d + update
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{

			/*
if( DBG )
{
	const double new_val = signal(row_pos, avg_pos)
			+ signal(row_pos, diff_pos) * lift_coef(lift_row, lift_col);

	fprintf( stderr, "    signal(%u,%u) += "
			"signal(%u,%u) * lift_coef(%u,%u)\n",
			row_pos, avg_pos, row_pos, diff_pos, lift_row, lift_col );
	fprintf( stderr, "      => %lf += %lf * %lf => %lf\n",
			signal(row_pos, avg_pos),
			signal(row_pos, diff_pos), lift_coef(lift_row, lift_col),
			new_val );
}
			 */

			signal(row_pos, avg_pos) +=
					signal(row_pos, diff_pos) * lift_coef(lift_row, lift_col);

			avg_pos += coef_offset;
			++lift_col;
		}
/*
if( DBG )
{
	fprintf( stderr, "\n" );
}
 */

		diff_pos  += coef_offset;
		start_pos += coef_offset;  // move starting position to next coef.
	}

/*
if( DBG )
{
	fprintf( stderr, "    Num Greater: %u\n", cases[2] );
}
 */

	// case 3: # left s > # right s
	size_type saved_pos = start_pos - coef_offset;
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		avg_pos = saved_pos;  // return s to the saved starting position

		// s = d + update
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{

			/*
if( DBG )
{
	const double new_val = signal(row_pos, avg_pos)
			+ signal(row_pos, diff_pos) * lift_coef(lift_row, lift_col);

	fprintf( stderr, "    signal(%u,%u) += "
			"signal(%u,%u) * lift_coef(%u,%u)\n",
			row_pos, avg_pos, row_pos, diff_pos, lift_row, lift_col );
	fprintf( stderr, "      => %lf += %lf * %lf => %lf\n",
			signal(row_pos, avg_pos),
			signal(row_pos, diff_pos), lift_coef(lift_row, lift_col),
			new_val );
}
			 */

			signal(row_pos, avg_pos) +=
					signal(row_pos, diff_pos) * lift_coef(lift_row, lift_col);

			avg_pos += coef_offset;
			++lift_col;
		}

/*
if( DBG )
{
	fprintf( stderr, "\n" );
}
 */

		diff_pos += coef_offset;
	}
}

/**
	Update stage for columns.

	The row col_pos of the signal is updated. The lifting coefficients in
	lift_coef at row lift_row are used. The given step size separates values.

	@param[in,out] signal Signal
	@param[in] col_pos Column position
	@param[in] step Step-size between values--the number of values separating
		the average and signal coefficients. Successive values for the average
		part of the signal are spaced 2 * step apart.  
	@param[in] lift_coef Lifting coefficients
	@param[in] lift_row Row position of lifting coefficients
 */
void
Lift_Scheme::Update_col( Matrix& signal, const size_type col_pos,
		const unsigned step,
		const Matrix& lift_coef, const size_type lift_row,
		const vector<unsigned>& cases )
{
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	// distance between next successive s's (or next d's)
	const unsigned coef_offset = 2 * step;

	size_type diff_pos = signal.row_begin() + step;  // position of d's
	size_type avg_pos  = signal.row_begin();         // position of s's
	size_type lift_col = lift_coef.col_begin();      // position of lifting coef.

	// case 1: # left s < # right s
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		// first s is always the first element of signal
		avg_pos = signal.row_begin();

		// s = d + update
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			signal(avg_pos, col_pos) +=
					signal(diff_pos, col_pos) * lift_coef(lift_row, lift_col);

			avg_pos += coef_offset;
			++lift_col;
		}
		diff_pos += coef_offset;
	}

	// case 2: # left s = # right s
	size_type start_pos = signal.row_begin();
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		avg_pos = start_pos;  // return to starting position

		// s = d + update
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
/*
	// fprintf( stderr, "\t%u -- (%u, %u) %s %lf",
			// j, lift_row, lift_col,
			// lift_coef.get_region().to_string().c_str(),
			// lift_coef(lift_row, lift_col)
	// );

fprintf( stderr, "  %u -- (%u, %u) -> %lf (%u, %u) -> %lf -- %s",
		j,
		avg_pos, col_pos, signal(avg_pos,col_pos),
		diff_pos, col_pos, signal(diff_pos,col_pos),
		signal.get_region().to_string().c_str()
);
 */

			signal(avg_pos, col_pos) +=
					signal(diff_pos, col_pos) * lift_coef(lift_row, lift_col);

// fprintf( stderr, " -- done\n" );

			avg_pos += coef_offset;
			++lift_col;
		}

		diff_pos  += coef_offset;
		start_pos += coef_offset;  // move starting position
	}

	// case 3: # left s > # right s
	size_type saved_pos = start_pos - coef_offset;
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		avg_pos = saved_pos;  // return s to the saved starting position

		// s = d + update
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			signal(avg_pos, col_pos) +=
					signal(diff_pos, col_pos) * lift_coef(lift_row, lift_col);

			avg_pos += coef_offset;
			++lift_col;
		}

		diff_pos += coef_offset;
	}
}

/**
	Get number of iterations for the 3 cases for the coefficients:
		1. # left s < # right s
		2. # left s = # right s
		3. # left s > # right s

	An average coefficents s is calculated using its neighbors. Due to
	image/matrix borders, we have to use different numbers of neighbors
	depending on the position of the coefficient.
 */
vector<unsigned>
Lift_Scheme::get_coef_cases( unsigned size, unsigned step, unsigned n )
{
	assert( n >= 2 );

	unsigned num_coef = get_num_coef( size, step );
	unsigned d_num    = num_coef / 2;
	unsigned num_odd  = int(ws_tools::is_odd( num_coef ));
// fprintf( stderr, "%u %u %u\n", num_coef, d_num, num_odd );

	vector<unsigned> cases( 3 );

	cases[0] = (n / 2) - 1;
	cases[1] = d_num - n + num_odd + 1;
	cases[2] = (n / 2) - num_odd;

	// verify that the sum of all cases is the same as the number of average
	// coefficients
	assert( (cases[0] + cases[1] + cases[2]) == (num_coef / 2) + num_odd );

	return( cases );
}

/**
	Change sign of all filter and lifting coefficients for use in the inverse
	transform.

	@param[in,out] lifts
 */
void
Lift_Scheme::invert_filter_lift( Lifts& lifts )
{
	_filter *= -1;
	lifts.X *= -1;
	lifts.Y *= -1;
}

/**
	Apply lifting scheme to 1D signal.

	@param[in,out] signal
	@param[in] num_levels
 */
void
Lift_Scheme::lift_int( Vector& signal, unsigned num_levels )
{
	// if N or N_tilde is 1, then no corresponding biorthogonal wavelet exists,
	// so use Haar wavelet instead
	if( get_N() == 1 || get_N_tilde() == 1 )
	{
		//lift_Haar_int( signal, num_levels );
		err_quit( "lift_Haar_int() unimplemented\n" );
		return;
	}

	// create initial lift and moment matrix sets
	Lifts   lifts( signal.sub_size(), 1, get_N(), get_N_tilde(), num_levels );
	Moments moments( signal.sub_size(), 1, get_N_tilde() );

	// calculate values for lifts and moments
	calc_moments_lifts( signal.sub_size(), 1, num_levels, moments, lifts );

	lift_int( signal, lifts, moments, num_levels );
}

/**
 */
void
Lift_Scheme::lift_int( Vector& signal, Lifts& lifts, Moments& moments,
	unsigned num_levels )
{
	// number of iterations
	num_levels = get_num_levels( num_levels, signal.sub_size(), get_N(),
			get_N_tilde() );

	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	// forward transform
	unsigned row = num_levels - 1;
	for( unsigned step = 1; step < max_step; step *= 2 )
	{
		Predict_int( signal, step );

		Update_int( signal, step, lifts.X.const_row_begin(row),
				signal.sub_size() );

		--row;
	}
}

/**
	Apply inverse lifting scheme to signal.

	@param[in,out] signal
	@param[in] num_levels
 */
void
Lift_Scheme::inverse_lift_int( Vector& signal, unsigned num_levels )
{
	if( get_N() == 1 || get_N_tilde() == 1 )
	{
		//inverse_lift_Haar_int( signal, num_levels );
		err_quit( "inverse_lift_Haar_int() unimplemented\n" );
		return;
	}

	// create initial lift and moment matrix sets
	Lifts   lifts( signal.sub_size(), 1, get_N(), get_N_tilde(), num_levels );
	Moments moments( signal.sub_size(), 1, get_N_tilde() );

	// calculate values for lifts and moments
	calc_moments_lifts( signal.sub_size(), 1, num_levels, moments, lifts );

	inverse_lift_int( signal, lifts, moments, num_levels );
}

/**
	Inverse lifting.

	@param[in,out] signal
	@param[in] lifts
	@param[in] moments
	@param[in] param
 */
void
Lift_Scheme::inverse_lift_int( Vector& signal, Lifts& lifts, Moments& moments,
	unsigned num_levels )
{
	// number of iterations
	num_levels = get_num_levels( num_levels, signal.sub_size(), get_N(),
			get_N_tilde() );
	unsigned max_step = (num_levels == 0) ? 0 : 1 << num_levels;

	// change signs of filter and lifts
	invert_filter_lift( lifts );

	// inverse transform
	unsigned row = num_levels - 1;
	for( unsigned step = max_step / 2; step >= 1; step /= 2 )
	{
		Update_int( signal, step, lifts.X.const_row_begin(row),
				signal.sub_size() );

		Predict_int( signal, step );

		--row;
	}

	// revert signs of filter/lifts back
	invert_filter_lift( lifts );
}

/**
	Predict difference coefficients d. Overwrite signal. Integer version

	@param[in,out] signal
	@param[in] step
 */
void
Lift_Scheme::Predict_int( Vector& signal, unsigned step )
{
	vector<unsigned> cases = get_coef_cases( signal.sub_size(), step, get_N() );
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	// distance between next successive s's (or next d's)
	unsigned coef_offset = 2 * step;

	unsigned filter_row  = 1; // start at row 1 (skip row 0)
	unsigned diff_signal = step;   // position of d's
	unsigned avg_signal  = 0;      // position of s's

	// accumulates sums
	Vector sums( signal.sub_size() );

	// case 1: # left s < # right s
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		avg_signal = 0;

		// d = d - prediction
		for( unsigned j = 0; j != get_N(); ++j )
		{
			sums(diff_signal) += signal(avg_signal) * _filter(filter_row, j);
			avg_signal += coef_offset;
		}

		diff_signal += coef_offset;
		++filter_row;
	}

	// case 2: # left s = # right s
	filter_row = _filter.sub_row_size() - 1;  // use row for equal #s
	unsigned start_pos = 0;
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		avg_signal = start_pos;  // return to

		// d = d - prediction
		for( unsigned j = 0; j != get_N(); ++j )
		{
			sums(diff_signal) += signal(avg_signal) * _filter(filter_row, j);
			avg_signal        += coef_offset;
		}

		diff_signal += coef_offset;
		start_pos   += coef_offset;  // change starting position of s
	}

	// case 3: # left s > # right s
	// first row with unequal # of s's (relative to end)
	filter_row = _filter.sub_row_size() - 2;
	unsigned saved_pos = (start_pos - coef_offset);
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		avg_signal = saved_pos;  // position to return for s

		// traverse opposite order
		// d = d - prediction
		for( int j = (get_N() - 1); j >= 0; --j )
		{
			sums(diff_signal) += signal(avg_signal) * _filter(filter_row, j);
			avg_signal += coef_offset;
		}
		--filter_row; // reverse the order of the row traversal

		diff_signal += coef_offset;
	}

	// for each d, subtract floor( sum + 1/2 )
	for( unsigned i = step; i != diff_signal; i += coef_offset )
	{
		// signal(i) -= floor( sums(i) + .5 );
	}
}

/**
	Update stage integer version.

	@param[in,out] signal
	@param[in] step
	@param[in] lift_vec
 */
void
Lift_Scheme::Update_int( Vector& signal, unsigned step,
	const_iterator lift_vec, size_type size )
{
	vector<unsigned> cases = get_coef_cases( size, step, get_N_tilde() );
	assert( cases.size() >= 3 );
	assert( int(cases[0]) >= 0 && int(cases[1]) >= 0 && int(cases[2]) >= 0 );

	// distance between next successive s's (or next d's)
	unsigned coef_offset = 2 * step;

	unsigned diff_signal = step;   // position of d's
	unsigned avg_signal  = 0;      // position of s's

	// accumulate sums
	Vector sums( size );

	// case 1: # left s < # right s
	for( unsigned num_less = cases[0]; num_less != 0; --num_less )
	{
		avg_signal = 0;  // first s is always the first element of signal

		// s = d + update
		//for( unsigned j = get_N_tilde(); j != 0; --j )
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			//signal(avg_signal) += signal(diff_signal) * lift_vec( lift_pos );
			sums(avg_signal) += signal(diff_signal) * (*lift_vec);

			avg_signal += coef_offset;
			++lift_vec;
		}
		diff_signal += coef_offset;
	}

	// case 2: # left s = # right s
	unsigned start_pos = 0;
	for( unsigned num_equal = cases[1]; num_equal != 0; --num_equal )
	{
		avg_signal = start_pos;  // return to

		// s = d + update
		//for( unsigned j = get_N_tilde(); j != 0; --j )
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			//signal(avg_signal) += signal(diff_signal) * lift_vec( lift_pos );
			sums(avg_signal) += signal(diff_signal) * (*lift_vec);

			avg_signal += coef_offset;
			++lift_vec;
		}

		diff_signal += coef_offset;
		start_pos   += coef_offset;
	}

	// case 3: # left s > # right s
	unsigned saved_pos = start_pos - coef_offset;
	for( unsigned num_greater = cases[2]; num_greater != 0; --num_greater )
	{
		avg_signal = saved_pos;  // return s to the saved starting position

		// s = d + update
		//for( unsigned j = get_N_tilde(); j != 0; --j )
		for( unsigned j = 0; j != get_N_tilde(); ++j )
		{
			//signal(avg_signal) += signal(diff_signal) * lift_vec( lift_pos );
			sums(avg_signal) += signal(diff_signal) * (*lift_vec);

			avg_signal += coef_offset;
			++lift_vec;
		}

		diff_signal += coef_offset;
	}

	// for each d, subtract floor( sum + 1/2 )
	for( unsigned i = 0; i != avg_signal; i += coef_offset )
	{
		signal(i) += floor( sums(i) + .5 );
	}
}

/**
	Reorder the elements of the given signal in octave form.
	@param[in,out] signal 2-D signal to reorder
	@param[in] num_levels Number of levels in the wavelet decomposition
 */
void
Lift_Scheme::octave_form( Matrix& signal, unsigned num_levels ) const
{
	// get number of iterations in each direction (and total iterations)
	// note that to compute row_levels we use the # of columns (not rows) since
	// we filter _along_ each row for row_levels iterations
	unsigned row_levels = 0;
	unsigned col_levels = 0;
	if( ws_tools::is_odd(get_N()) || ws_tools::is_odd(get_N_tilde()) )
	{
		// use get_num_levels() associated with Haar function
		row_levels = get_num_levels( num_levels, signal.sub_col_size() );
		col_levels = get_num_levels( num_levels, signal.sub_row_size() );
	}
	else // levels for biorthogonal transformation
	{
		// number of rows, for example, is a function of the number of columns
		// (not the number of rows) since the loop iteration below inside the if
		// statement is executed only if the coefficients can be further split
		row_levels = get_num_levels( num_levels, signal.sub_col_size(),
				get_N(), get_N_tilde() );
		col_levels = get_num_levels( num_levels, signal.sub_row_size(),
				get_N(), get_N_tilde() );
	}

	num_levels = max<unsigned>( row_levels, col_levels );

	// current number of rows to reorder and corresponding size of each row
	size_type row_end  = signal.row_end();
	size_type col_size = signal.sub_col_size();

	// current number of columns to reorder and corresponding size of each column
	size_type col_end  = signal.col_end();
	size_type row_size = signal.sub_row_size();

	// for each step size separating s components at this level
	for( unsigned level = 0; level != num_levels; ++level )
	{
		// reorder rows
		if( row_levels != 0 )
		{
			// reorder each row separately
			for( size_type i = signal.row_begin(); i != row_end; ++i )
			{
				octave_form_1D( signal.row_begin(i), col_size );
			}

			--row_levels;

			// the next octave size is half the size of the current octave, but
			// but round up if the size is odd (e.g., 6 -> 3 and 7 -> 4)
			row_end  = (row_end + 1) / 2;
			col_size = (col_size + 1) / 2;
		}

		// reorder columns
		if( col_levels != 0 )
		{
			// reorder each column separately
			for( size_type j = signal.col_begin(); j != col_end; ++j )
			{
				octave_form_1D( signal.col_begin(j), row_size );
			}

			--col_levels;

			// the next octave size is half the size of the current octave, but
			// but round up if the size is odd (e.g., 6 -> 3 and 7 -> 4)
			col_end  = (col_end + 1) / 2;
			row_size = (row_size + 1) / 2;
		}
	}
}

/**
	Reorder the elements of the given signal in octave form.
	@param[in,out] signal Iterator to 1-D signal to reorder
	@param[in] size Size of signal
 */
void
Lift_Scheme::octave_form_1D( iterator signal, size_type size ) const
{
	Vector copy_vec( size );  // temporary place to copy data when reordering

	// iterators to copy from: diff starts one after signal's start
	iterator avg_from  = signal;
	iterator diff_from = signal + 1;

	// iterators to copy to: diff starts at middle of vector
	iterator avg_to  = copy_vec.begin();
	iterator diff_to = copy_vec.begin() + (size - (size / 2));

	// copy each element except last
	for( size_type i = 0; i < (size - 2); i += 2 )
	{
		*avg_to  = *avg_from;
		*diff_to = *diff_from;

		avg_from  += 2;
		diff_from += 2;
		++avg_to;
		++diff_to;
	}

	// move last element (or two)
	if( ws_tools::is_odd( size ) )
	{
		*avg_to  = *avg_from;
	}
	else
	{
		*avg_to  = *avg_from;
		*diff_to = *diff_from;
	}

	// copy vector into signal to change signal's order
	iterator vec_iter = copy_vec.begin();
	iterator vec_end  = copy_vec.end();
	for( ; vec_iter != vec_end; ++vec_iter, ++signal )
	{
		*signal = *vec_iter;
	}
}

/**
	Reorder the elements of the given signal in octave form.
	@param[in,out] signal 2-D signal to reorder
	@param[in] num_levels Number of levels in the wavelet decomposition
 */
void
Lift_Scheme::inverse_octave_form( Matrix& signal, unsigned num_levels ) const
{
	// get number of iterations in each direction (and total iterations)
	// note that to compute row_levels we use the # of columns (not rows) since
	// we filter _along_ each row for row_levels iterations
	unsigned row_levels = 0;
	unsigned col_levels = 0;
	if( ws_tools::is_odd(get_N()) || ws_tools::is_odd(get_N_tilde()) )
	{
		// use get_num_levels() associated with Haar function
		row_levels = get_num_levels( num_levels, signal.sub_col_size() );
		col_levels = get_num_levels( num_levels, signal.sub_row_size() );
	}
	else // levels for biorthogonal transformation
	{
		row_levels = get_num_levels( num_levels, signal.sub_col_size(),
				get_N(), get_N_tilde() );
		col_levels = get_num_levels( num_levels, signal.sub_row_size(),
				get_N(), get_N_tilde() );
	}

	num_levels = max<unsigned>( row_levels, col_levels );

	// for each step size separating s components at this level
	for( int level = (num_levels - 1); level >= 0; --level )
	{
		// current number of rows to reorder and corresponding size of each row
		size_type row_end  = signal.row_end();
		size_type col_size = signal.sub_col_size();

		// current number of columns to reorder and corresponding size of each
		// column
		size_type col_end  = signal.col_end();
		size_type row_size = signal.sub_row_size();

		// reorder columns
		if( col_levels != 0 )
		{
			for( size_type k = 0; k != size_type(level); ++k )
			{
				col_end  = (col_end + 1)  / 2;
				row_size = (row_size + 1) / 2;
			}

			// reorder each column separately
			for( size_type j = signal.col_begin(); j != col_end; ++j )
			{
				inverse_octave_form_1D( signal.col_begin(j), row_size );
			}

			--col_levels;
		}

		// reorder rows
		if( row_levels != 0 )
		{
			for( size_type k = 0; k != size_type(level); ++k )
			{
				row_end  = (row_end + 1)  / 2;
				col_size = (col_size + 1) / 2;
			}

			// reorder each row separately
			for( size_type i = signal.row_begin(); i != row_end; ++i )
			{
				inverse_octave_form_1D( signal.row_begin(i), col_size );
			}

			--row_levels;
		}
	}
}

/**
	Reorder the elements of the given signal in octave form.
	@param[in,out] signal Iterator to 1-D signal to reorder
	@param[in] size Size of signal
 */
void
Lift_Scheme::inverse_octave_form_1D( iterator signal, size_type size ) const
{
	Vector copy_vec( size );  // temporary place to copy data when reordering

	// iterators to copy from: diff starts one after signal's start
	iterator avg_from  = signal;
	iterator diff_from = signal + (size - (size / 2));

	// iterators to copy to: diff starts at middle of vector
	iterator avg_to  = copy_vec.begin();
	iterator diff_to = copy_vec.begin() + 1;

	// copy each element except last
	for( size_type i = 0; i < (size - 2); i += 2 )
	{
		*avg_to  = *avg_from;
		*diff_to = *diff_from;

		++avg_from;
		++diff_from;
		avg_to  += 2;
		diff_to += 2;
	}

	// move last element (or two)
	if( ws_tools::is_odd( size ) )
	{
		*avg_to  = *avg_from;
	}
	else
	{
		*avg_to  = *avg_from;
		*diff_to = *diff_from;
	}

	// copy vector into signal to change signal's order
	iterator vec_iter = copy_vec.begin();
	iterator vec_end  = copy_vec.end();
	for( ; vec_iter != vec_end; ++vec_iter, ++signal )
	{
		*signal = *vec_iter;
	}
}

/**
	Get subbands.
	Code is almost identical to octave_form() except the end of the first for
	loop: the separate subbands HH, LH, and HL are created for the current level.
	Matrix can not be the output of octave_form().

	TODO Write function to combine subbands into a single matrix to verify we get
	the same matrix as calling octave_form().

	Algorithm:
	1. Determine number of transform levels for rows/columns/overall
	2. Apply similar algorithm to octave_form():
			For each level
				Reorder rows if possible
				Reorder columns if possible

				If reordered both rows and columns
				If reordered only rows
				If reordered only columns

	Implementation notes:
		To get and save a subband, we use code similar to the following 3 lines:
			signal.set_region( row_start, col_start, row_size, col_size );
			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::HL );
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

		The first line sets the region to that of the subband we are interested in
		after the signal has been reordered for that level.
		The second line records the level and band this subband occurs at.
		The third line saves the subband into the set of all subbands, which is
		indexed by level and band position.
		Note that 'Wavelet_Subband( signal, ... )' only copies the portion of the
		signal that is inside the region set in the first line by set_region().
		The original region is saved and reset.

	@param[in,out] signal
	@param[in] num_levels
	@post signal is put in octave format as though octave_form() were called.
 */
Wavelet_Subbands
Lift_Scheme::get_subbands( Matrix& signal, unsigned num_levels ) const
{
	// set of all wavelet subbands in signal
	Wavelet_Subbands subbands;

	// level and band from decomposition: used as temporary in multiple places
	Wavelet_Level_Band level_band;

	// We change the active region below to the subband we need, so we save
	// the current region since we need to continue partitioning the data
	// afterwards using the current region.
	Matrix_Region saved_region = signal.get_region();

	// get number of iterations in each direction (and total iterations)
	// note that to compute row_levels we use the # of columns (not rows) since
	// we filter _along_ each row for row_levels iterations
	unsigned row_levels = 0;
	unsigned col_levels = 0;
	if( ws_tools::is_odd(get_N()) || ws_tools::is_odd(get_N_tilde()) )
	{
		// use get_num_levels() associated with Haar function
		row_levels = get_num_levels( num_levels, signal.sub_col_size() );
		col_levels = get_num_levels( num_levels, signal.sub_row_size() );
	}
	else // levels for biorthogonal transformation
	{
		// number of rows, for example, is a function of the number of columns
		// (not the number of rows) since the loop iteration below inside the if
		// statement is executed only if the coefficients can be further split
		row_levels = get_num_levels( num_levels, signal.sub_col_size(),
				get_N(), get_N_tilde() );
		col_levels = get_num_levels( num_levels, signal.sub_row_size(),
				get_N(), get_N_tilde() );
	}

	num_levels = max<unsigned>( row_levels, col_levels );

	// current number of rows to reorder and corresponding size of each row
	size_type row_end  = signal.row_end();
	size_type col_size = signal.sub_col_size();

	// current number of columns to reorder and corresponding size of each column
	size_type col_end  = signal.col_end();
	size_type row_size = signal.sub_row_size();

	// for each step size separating s components at this level
	for( unsigned level = 0; level != num_levels; ++level )
	{
	// fprintf( stderr, "Level: %u\n", level );
	// fprintf( stderr, "Row end: %u, Col size: %u\n", row_end, col_size );
	// fprintf( stderr, "Col end: %u, Row size: %u\n", col_end, row_size );

		// reorder rows--note that row_levels is decremented
		if( row_levels != 0 )
		{
		// fprintf( stderr, "Row level: %u\n", row_levels );

			// reorder each row separately
			for( size_type i = signal.row_begin(); i != row_end; ++i )
			{
				octave_form_1D( signal.row_begin(i), col_size );
			}

			// the next octave size is half the size of the current octave, but
			// but round up if the size is odd (e.g., 6 -> 3 and 7 -> 4)
			row_end  = (row_end + 1) / 2;
			col_size = (col_size + 1) / 2;
		}

		// reorder columns
		if( col_levels != 0 )
		{
		// fprintf( stderr, "Col level: %u\n", col_levels );

			// reorder each column separately
			for( size_type j = signal.col_begin(); j != col_end; ++j )
			{
				octave_form_1D( signal.col_begin(j), row_size );
			}

			// the next octave size is half the size of the current octave, but
			// but round up if the size is odd (e.g., 6 -> 3 and 7 -> 4)
			col_end  = (col_end + 1) / 2;
			row_size = (row_size + 1) / 2;
		}

		// row_end and col_end have been updated to the positions we need
		// with col_size and row_size giving the size of each subband

		// if both the rows and columns were reordered, then we have 3 subbands:
		// LH, HL, and HH
		if( row_levels != 0 && col_levels != 0 )
		{
			--row_levels;
			--col_levels;

		// fprintf( stderr, "Both rows and columns\n" );

			// LH subband: top-right band
			signal.set_region( signal.row_begin(), col_end, row_end, col_size );

			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::LH );

			// Note: copying signal into the subband only copies the current region
			// of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

		// fprintf( stderr, "\tLH: %s\n",
				// signal.get_region().to_string().c_str() );
		// subbands[ level_band ].write();

			// we need signal.col_begin() of the original region,
			// so revert back to it
			signal.set_region( saved_region );

			// HL subband: bottom-left band
			signal.set_region( row_end, signal.col_begin(), row_size, col_size );

			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::HL );

			// Note: copying signal into the subband only copies the current region
			// of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

		// fprintf( stderr, "\tHL: %s\n",
				// signal.get_region().to_string().c_str() );
		// subbands[ level_band ].write();

			// HH subband: bottom-right band
			signal.set_region( row_end, col_end, row_size, col_size );

			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::HH );

			// Note: copying signal into the subband only copies the current region
			// of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

		// fprintf( stderr, "\tHH: %s\n",
				// signal.get_region().to_string().c_str() );
		// subbands[ level_band ].write();
		}

		// if only the rows were reordered, we have one H subband on the right
		else if( row_levels != 0 )
		{
			--row_levels;

		// fprintf( stderr, "Rows only\n" );

			// H subband: right-half band
			// Note that row_size gives both the starting region row and the
			// length of the region since row_size
			// TODO Verify that using row_size is correct for non-powers of 2
			signal.set_region( signal.row_begin(), col_size, row_size, col_size );

			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::H );

			// Note: the assignment only copies the current region of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

	// subbands[ level_band ].write();
	// fprintf( stderr, "\tH: %s\n", signal.get_region().to_string().c_str() );

		}

		// if only the columns were reordered, we have one H subband on the bottom
		else if( col_levels != 0 )
		{
			--col_levels;

		// fprintf( stderr, "Columns only\n" );

			// H subband: bottom-half band
			// Note that row_size gives both the starting region row and the
			// length of the region since row_size
			// TODO Verify that using row_size is correct for non-powers of 2
			signal.set_region( row_size, signal.col_begin(), row_size, col_size );


			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::H );

			// Note: the assignment only copies the current region of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

	// fprintf( stderr, "\tH: %s\n", signal.get_region().to_string().c_str() );
	// subbands[ level_band ].write();
		}

		// revert to original region
		signal.set_region( saved_region );

	// signal.write();
	// fprintf( stderr, "Row end: %u, Col size: %u\n", row_end, col_size );
	// fprintf( stderr, "Col end: %u, Row size: %u\n", col_end, row_size );
	// fprintf( stderr, "\n" );
	}

	// LL subband: top-left band
	signal.set_region( signal.row_begin(), signal.col_begin(),
			row_size, col_size );

	level_band = Wavelet_Level_Band( num_levels, Wavelet_Level_Band::LL );

	// Note: copying signal into the subband only copies the current region
	// of the signal
	subbands[ level_band ] = Wavelet_Subband( signal, level_band );

// fprintf( stderr, "\tLL: %s\n", signal.get_region().to_string().c_str() );
// subbands[ level_band ].write();

	// revert to original region
	signal.set_region( saved_region );

	return( subbands );
}
