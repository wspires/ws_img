/**
	@file   keypoint_localization.hpp
	@author Wade Spires
	@date   2006/2/12
	@brief  Functions for localizing keypoints.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef _KEYPOINT_LOCALIZATION_HPP
#define _KEYPOINT_LOCALIZATION_HPP

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
#include "Matrix.hpp"
#include "Vector.hpp"

// local headers
#include "Keypoint.hpp"
#include "Scale_Space.hpp"

// .03 and 10 come from Lowe's paper
extern Keypoint_List keypoint_localization( const Keypoint_List&,
		const Scale_Space&, double = .03, double = 10 );

extern bool inverse( const ws_img::Matrix&, ws_img::Matrix& );

#endif // _KEYPOINT_LOCALIZATION_HPP
