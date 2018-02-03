/**
	@file   Sift.hpp
	@author Wade Spires
	@date   2006/4/11
	@brief  Class Sift.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef _SIFT_HPP
#define _SIFT_HPP

#include "Image.hpp"
#include "Keypoint.hpp"

extern Keypoint_List sift( const ws_img::Image& );

extern void draw_keypoints( ws_img::Image&, const Keypoint_List& );

#endif // _SIFT_HPP 
