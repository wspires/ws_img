/**
	@file   Gabor.hpp
	@author Wade Spires
	@date   2006/3/30
	@brief  Class Gabor.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef GABOR_HPP
#define GABOR_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <algorithm>
#include <string>
#include <vector>

// c headers
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// tools headers
#include "ws_tools.hpp"

// Image headers
#include "Gauss.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

// local headers

#ifndef PI
#define M_PI PI
#endif

extern ws_img::Matrix make_gabor_filter( ws_img::Matrix::size_type = 5,
		ws_img::Matrix::size_type = 5, double = 0, double = 0 );
extern double gabor( double, double, double = 0, double = 0, double = .5,
		double = 8, double = 1, double = 0, double = 0 ); 

extern void make_gabor_filter( ws_img::Matrix::size_type,
		ws_img::Matrix::size_type, double, double, double,
		ws_img::Matrix&, ws_img::Matrix& );

#endif // GABOR_HPP 
