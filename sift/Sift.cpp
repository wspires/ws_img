/**
	@file   Sift.cpp
	@author Wade Spires
	@date   2006/4/11
	@brief  Class Sift.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Sift.hpp"

// local headers
#include "extrema_detection.hpp"
#include "keypoint_descriptor.hpp"
#include "keypoint_localization.hpp"
#include "orientation_assignment.hpp"
#include "Scale_Space.hpp"

// my namespace
using namespace ws_img;

typedef Image::size_type   size_type;
typedef Image::value_type  value_type;
typedef Image::pixel_type  pixel_type;
typedef Image_Coordinate   Img_Coord;

/**
	Run SIFT algorithm on the given image.
 */
Keypoint_List sift( const Image& img )
{
	// construct scale-space representation of the image
	Scale_Space scale_space( img );
	// scale_space.write( base_name );

	// detect extrema in scale-space
	Keypoint_List keypoints = extrema_detection( scale_space );

	// localize keypoints
	keypoints = keypoint_localization( keypoints, scale_space );

	// assign orientation to each keypoint
	keypoints = orientation_assignment( scale_space, keypoints );

	// create descriptor for each keypoint
	keypoints = keypoint_descriptor( scale_space, keypoints );

	return( keypoints );
}

/**
	Draw keypoints on image.
 */
void draw_keypoints( Image& img, const Keypoint_List& keypoints )
{
	img.to_gray();
	Image_Region box( 0, 0, 5, 5 );

	for( unsigned i = 0; i != keypoints.size(); ++i )
	{
		Keypoint point = keypoints[i];
		box.set_row_begin( point.img_coord.row );
		box.set_col_begin( point.img_coord.col );
		img.draw( box );
	}
}
