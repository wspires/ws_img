/**
	@file   Gauss.hpp
	@author Wade Spires
	@date   2006/1/27
	@brief  Functions based on Gaussian.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef GAUSS_HPP
#define GAUSS_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <algorithm>
#include <string>
#include <vector>

// c headers
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// tools headers
#include "ws_tools.hpp"

// local headers
#include "Image.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

namespace ws_img
{
#ifndef PI
#define PI M_PI
#endif

#ifndef TWO_PI
#define TWO_PI (2 * PI)
#endif

	extern double DEFAULT_STD_SPREAD;  //< Number of std. dev. to use

	// Gaussian functions: 1D and derivatives
	extern double gauss( double, double, double = 0 );
	extern double gauss_dx( double, double );
	extern double gauss_dxdx( double, double );

	// Gaussian functions: 2D and derivative
	extern double gauss_2D( double, double, double );
	extern double gauss_2D_dx( double, double, double );

	// Laplacian of Gaussian function
	extern double laplace_gauss( double, double, double );

	// Gaussian filters: 1D and derivatives
	extern Vector make_gauss_filter( double = 1, double = DEFAULT_STD_SPREAD );
	extern Vector make_gauss_1D_filter( double = 1,
			double = DEFAULT_STD_SPREAD );
	extern Vector make_gauss_dx_filter( double = 1,
			double = DEFAULT_STD_SPREAD );
	extern Vector make_gauss_dxdx_filter( double = 1,
			double = DEFAULT_STD_SPREAD );

	// Gaussian filters: 2D and derivative
	extern Matrix make_gauss_2D_filter( double = 1,
			double = DEFAULT_STD_SPREAD );
	extern Matrix make_gauss_2D_dx_filter( double ,
			double = DEFAULT_STD_SPREAD );

	// Laplacian of Gaussian filter
	extern Matrix make_log_filter( double = 1,
			double = DEFAULT_STD_SPREAD );
	extern void make_log_filter( Vector&, Vector&, double = 1,
			double = DEFAULT_STD_SPREAD );

	extern Matrix& filter_log( Matrix&, double = 1,
			double = DEFAULT_STD_SPREAD );
	extern Image& filter_log( Image&, double = 1,
			double = DEFAULT_STD_SPREAD, bool = true );

} // namespace ws_img

#endif // GAUSS_HPP 
