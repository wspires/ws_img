/**
	@file   Harris.hpp
	@author Wade Spires
	@date   2006/4/8
	@brief  Class Harris.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef HARRIS_HPP
#define HARRIS_HPP

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
#include "Image.hpp"

// local headers

extern void harris( ws_img::Image&, ws_img::Image::size_type = 5,
		ws_img::Image::value_type = 60, const double = .04 );

#endif // HARRIS_HPP 
