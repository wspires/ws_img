/**
	@file   extrema_detection.hpp
	@author Wade Spires
	@date   2006/1/31
	@brief  Functions for finding extrema in scale-scape.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef _EXTREMA_DETECTION_HPP
#define _EXTREMA_DETECTION_HPP

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
#include "Image.hpp"
#include "Image_Coordinate.hpp"

// local headers
#include "Keypoint.hpp"
#include "Scale_Space.hpp"

extern Keypoint_List extrema_detection( Scale_Space&, unsigned = 3 );

extern void find_scale_extrema( const std::vector<const ws_img::Image*>&,
		const ws_img::Image&, const std::vector<const ws_img::Image*>&,
		Keypoint_List&, ws_img::Image::size_type = 3,
		ws_img::Image::size_type = 3 );

extern bool is_extrema_point( ws_img::Image::value_type,
		ws_img::Image::size_type, ws_img::Image::size_type, const ws_img::Image&,
		const std::vector<const ws_img::Image*>&,
		const std::vector<const ws_img::Image*>&, const ws_img::Image_Region& );

extern bool is_extrema_point_DoG( const ws_img::Image::value_type,
		const std::vector<const ws_img::Image*>&,
		const ws_img::Image_Region&, const bool );

#endif // _EXTREMA_DETECTION_HPP 
