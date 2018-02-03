/**
	@file   keypoint_descriptor.hpp
	@author Wade Spires
	@date   2006/2/22
	@brief  Step 3 of SIFT algorithm to perform orientation assignment.
	One or more orientations are assigned to each keypoint location based on
	local image gradient directions. All future operations are performed on
	image data that has been transformed relative to the assigned orientation,
	scale, and location for each feature, thereby providing invariance to these
	transformations.  

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef _KEYPOINT_DESCRIPTOR_HPP
#define _KEYPOINT_DESCRIPTOR_HPP

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
#include "Matrix.hpp"
#include "Vector.hpp"

// local headers
#include "Keypoint.hpp"
#include "Scale_Space.hpp"
#include "keypoint_localization.hpp"
#include "orientation_assignment.hpp"

extern const ws_img::Matrix::size_type              DESC_NUM_BINS;
extern const ws_img::Matrix::value_type             DESC_BIN_WIDTH;
extern const std::vector<ws_img::Vector>::size_type DESC_NUM_HISTS;

// extern ws_img::Matrix Gauss;

Keypoint_List keypoint_descriptor( const Scale_Space&,
		const Keypoint_List& );

std::vector<ws_img::Vector> make_histograms( const ws_img::Matrix&,
		const ws_img::Matrix&, ws_img::Matrix::size_type, unsigned );

#endif // _KEYPOINT_DESCRIPTOR_HPP
