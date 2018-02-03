/**
	@file   extrema_detection.cpp
	@author Wade Spires
	@date   2006/1/31
	@brief  Functions for finding extrema in scale-space.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "extrema_detection.hpp"

using std::string;
using std::vector;

using namespace ws_img;

typedef Image::size_type   size_type;
typedef Image::value_type  value_type;
typedef Image::pixel_type  pixel_type;
typedef Image_Coordinate   Img_Coord;

/**
	1. Scale-space extrema detection: The first stage of computation searches
	over all scales and image locations. It is implemented efficiently by using a
	difference-of-Gaussian function to identify potential interest points that
	are invariant to scale and orientation.  

	@param[in] scale_space
	@param[in] num_scales_to_check
 */
Keypoint_List
extrema_detection( Scale_Space& scale_space,
		unsigned num_scales_to_check )
{
	log_msg( "\nFinding extrema\n" );

	Keypoint_List keypoints;  // keypoints found

	// verify that we have a center image to check for extrema
	if( ws_tools::is_even( num_scales_to_check ) )
	{
		err_quit( "Number of scales to check for extrema must be odd:"
				" %u scales given", num_scales_to_check );
	}

	// pointers to images in scales below and above current scale in an octave
	// to check for extrema
	// note: use pointers to avoid copying full data (can/should we do this with
	// references instead?)
	vector<const Image*> bottom_scales( num_scales_to_check / 2 );
	vector<const Image*> top_scales( num_scales_to_check / 2 );

	// for each octave
	for( unsigned i = 0; i != scale_space.size(); ++i )
	{
		log_msg( "Searching octave %u / %u\n", i, (scale_space.size() - 1) );

		const Octave& octave = scale_space[i];

		// verify that we can check the requested number of scales
		if( octave.size() < num_scales_to_check )
		{
			err_quit( "Requested to check %u scales, but the octave has %u\n",
					octave.size(), num_scales_to_check );
		}

		// number of sets of scale images to check, e.g., if we are checking 3
		// scales at one time with 6 scales per octave, then we have 3 sets to
		// check: [0,2], [1,3], and [2,4] (we exclude [3,5] since 5 has no DoG)
		unsigned num_scale_sets = octave.size() - num_scales_to_check; 

		// for each possible scale starting position
		for( unsigned j = 0; j != num_scale_sets; ++j )
		{
			log_msg( "   Searching scales %u to %u (%u / %u)\n", j,
					j + bottom_scales.size() + top_scales.size(),
					j, (num_scale_sets - 1) );

			// get images at bottom scale
			for( unsigned k = 0; k != bottom_scales.size(); ++k )
			{
				// order the images in bottom_scales s.t. the images closest to the
				// current scale are at the front of the list (this ordering is not
				// used in the traditional SIFT algorithm but may be useful later)
				bottom_scales[ bottom_scales.size() - k - 1 ] = &octave[j + k].D;
			}

			// get image at middle scale (the image to search for extrema in)
			const Image& current_scale = octave[ j + bottom_scales.size() ].D;

			// get images at top scale
			for( unsigned k = 0; k != top_scales.size(); ++k )
			{
				// the index is calculated to skip the bottom and middle scales
				top_scales[ k ] = &octave[j + top_scales.size() + k + 1].D;
			}

			// find extrema in the current scale image
			Keypoint_List scale_points;
			find_scale_extrema( bottom_scales, current_scale, top_scales,
					scale_points );

			// add current points to set of all keypoints
			for( unsigned k = 0; k != scale_points.size(); ++k )
			{
				// save octave/scale position
				Keypoint& point = scale_points[k];
				point.octave_pos = i;
				point.scale_pos  = j + bottom_scales.size();

				keypoints.push_back( point );
			}
		}
	}

	return( keypoints );
}

#if false

/**
	Find extrema in the current set of scale images across position using the
	top and bottom scale image sets.

	@param[in] bottom_scales Set of scale images below the current scale image in
		the pyramid
	@param[in] current_scale Current scale image to find extrema in
	@param[in] top_scales Set of scale images above the current scale image in
		the pyramid
	@param[out] keypoints Set of keypoints detected in the current scale image
	@param[in] row_size Number of rows in neighborhood to check around each
		point (must be odd; default is 3)
	@param[in] col_size Number of columns in neighborhood to check around each
		point (must be odd; default is 3)
 */
void
find_scale_extrema( const vector<const Image*>& bottom_scales,
		const Image& current_scale, const vector<const Image*>& top_scales,
		Keypoint_List& keypoints, size_type row_size, size_type col_size )
{
	// statistics of keypoints found
	static int key_count    = 0;
	static int no_key_count = 0;
	static int total        = 0;

	// verify that we have a center position in the neighborhood we are checking
	if( ws_tools::is_even( row_size ) || ws_tools::is_even( col_size ) )
	{
		err_quit( "Region size must be odd: given (%u x %u)\n",
				row_size, col_size );
	}
	Image::size_type half_row_size = (row_size / 2);
	Image::size_type half_col_size = (col_size / 2);

	// for each position in the current scale image
	for( size_type i = (current_scale.row_begin() + half_row_size);
			i < (current_scale.row_end() - half_row_size);
			++i )
	{
		for( size_type j = (current_scale.col_begin() + half_col_size);
				j < (current_scale.col_end() - half_col_size);
				++j )
		{
			// the pixel to determine whether or not it is an extrema
			Image::value_type value = current_scale(i,j);

			++total;

			// get neighborhood around (i,j)
			Image_Region region(
					i - half_row_size, j - half_col_size, row_size, col_size );

			// add all pixels around a neighborhood to vector--we will compare (i,j)
			// to the min and max of this vector to determine if (i,j) is an extrema
			Vector img_region( region.row_size() * region.col_size() - 1 );
			size_type z = 0;
			for( size_type y = region.row_begin(); y != region.row_end(); y++ )
			{
				for( size_type x = region.col_begin(); x != region.col_end(); x++ )
				{
					// add all neighboring points (not including (i,j)) to vector
					if( y != i && x != j )
					{
						assert( z != img_region.size() );
						img_region(z++) = current_scale(y,x);
					}
				}
			}

			bool is_keypoint = false;  // whether (i,j) is a keypoint

			// determine whether (i,j) is extrema of the current image
			Image::value_type min, max;
			img_region.get_min_max( &min, &max );
			if( value < min  )
			{
				// check if (i,j) is minima in the bottom and top scales (but don't
				// check the top if (i,j) did not pass the bottom scales)
				if( is_extrema_point_DoG( value, bottom_scales, region, true ) )
				{
					is_keypoint =
						is_extrema_point_DoG( value, top_scales, region, true );
				}
			}
			else if( value > max  )
			{
				// check if (i,j) is maxima in the bottom and top scales (but don't
				// check the top if (i,j) did not pass the bottom scales)
				if( is_extrema_point_DoG( value, bottom_scales, region, false ))
				{
					is_keypoint =
						is_extrema_point_DoG( value, top_scales, region, false );
				}
			}

			if( is_keypoint )
			{
				// add new keypoint to list of keypoints found in 'current_scale';
				// the scale and octave position should be set after returning
				++key_count;
				Keypoint keypoint;
				keypoint.img_coord.row = i;
				keypoint.img_coord.col = j;
				keypoints.push_back( keypoint );
			}
			else
			{
				++no_key_count;
			}
		}
	}
	log_msg( "   Keypoints: %u, Non-keypoints: %u, Total: %u\n\n",
			key_count, no_key_count, total );
}

#endif

/**
	Find extrema in the current set of scale images across position using the
	top and bottom scale image sets.

TODO
Causes segfault ws_tools::get_ext_name() due to return() in innermost loop.
I have no idea why this happens since I get the same keypoints either way, but
this method is slightly more efficient.

	@param[in] bottom_scales Set of scale images below the current scale image in
		the pyramid
	@param[in] current_scale Current scale image to find extrema in
	@param[in] top_scales Set of scale images above the current scale image in
		the pyramid
	@param[out] keypoints Set of keypoints detected in the current scale image
	@param[in] row_size Number of rows in neighborhood to check around each
		point (must be odd; default is 3)
	@param[in] col_size Number of columns in neighborhood to check around each
		point (must be odd; default is 3)
 */
void
find_scale_extrema( const vector<const Image*>& bottom_scales,
		const Image& current_scale, const vector<const Image*>& top_scales,
		Keypoint_List& keypoints, size_type row_size, size_type col_size )
{
	// statistics of keypoints found
	static int key_count    = 0;
	static int no_key_count = 0;
	static int total        = 0;

	// verify that we have a center position in the neighborhood we are checking
	if( ws_tools::is_even( row_size ) || ws_tools::is_even( col_size ) )
	{
		err_quit( "Region size must be odd: given (%u x %u)\n",
				row_size, col_size );
	}
	Image::size_type half_row_size = (row_size / 2);
	Image::size_type half_col_size = (col_size / 2);

	// get neighborhood around (i,j)
	Image_Region reg( 0, 0, row_size, col_size );

	// for each position in the current scale image
	for( size_type i = (current_scale.row_begin() + half_row_size);
			i < (current_scale.row_end() - half_row_size);
			++i )
	{
		for( size_type j = (current_scale.col_begin() + half_col_size);
				j < (current_scale.col_end() - half_col_size);
				++j )
		{
			// the pixel to determine whether or not it is an extrema
			Image::value_type value = current_scale(i,j);

			++total;

			// get neighborhood around (i,j)
			reg.set_row_begin( i - half_row_size );
			reg.set_col_begin( j - half_col_size );

			// add point to list of keypoints if it is an extrema
			if( is_extrema_point( value, i, j, current_scale, bottom_scales,
						top_scales, reg ) )
			{
				// add new keypoint to list of keypoints found in 'current_scale';
				// the scale and octave position should be set after returning
				++key_count;
				Keypoint keypoint;
				keypoint.img_coord.row = i;
				keypoint.img_coord.col = j;
				keypoints.push_back( keypoint );
			}
			else
			{
				++no_key_count;
			}
		}
	}
	log_msg( "   Keypoints: %u, Non-keypoints: %u, Total: %u\n\n",
			key_count, no_key_count, total );
}

/**
	Determine if given point is extrema across all images.

	@param[in] value Value to check
	@param[in] i Row position of value
	@param[in] j Column position of value
	@param[in] bottom_scales Set of scale images below the current scale image in
		the pyramid
	@param[in] current_scale Current scale image to find extrema in
	@param[in] top_scales Set of scale images above the current scale image in
		the pyramid
	@param[out] keypoints Set of keypoints detected in the current scale image
	@param[in] region Region inside of each scale image to check
 */
bool
is_extrema_point( value_type value, size_type i, size_type j,
	const Image& current_scale, const vector<const Image*>& bottom_scales,
	const vector<const Image*>& top_scales,
	const Image_Region& reg )
{
	if( value > current_scale(reg.row_begin(), reg.col_begin()) )
	{
		for( size_type y = reg.row_begin(); y != reg.row_end(); ++y )
		{
			for( size_type x = reg.col_begin(); x != reg.col_end(); ++x )
			{
				// if point is not an extrema, stop
				if( y != i && x != j && value <= current_scale(y,x) )
				{
					return( false );
				}
			}
		}

		// check if (i,j) is maxima in the bottom and top scales (but don't
		// check the top if (i,j) did not pass the bottom scales)
		if( is_extrema_point_DoG( value, bottom_scales, reg, false ))
		{
			return( is_extrema_point_DoG( value, top_scales, reg, false ) );
		}
	}
	else if( value < current_scale(reg.row_begin(), reg.col_begin()) )
	{
		for( size_type y = reg.row_begin(); y != reg.row_end(); ++y )
		{
			for( size_type x = reg.col_begin(); x != reg.col_end(); ++x )
			{
				// if point is not an extrema, stop
				if( y != i && x != j && value >= current_scale(y,x) )
				{
					return( false );
				}
			}
		}

		// check if (i,j) is minima in the bottom and top scales (but don't
		// check the top if (i,j) did not pass the bottom scales)
		if( is_extrema_point_DoG( value, bottom_scales, reg, true ) )
		{
			return( is_extrema_point_DoG( value, top_scales, reg, true ) );
		}
	}

	return( false );
}

/**
	Check set of scale images to determine if given value is an extrema of the
	given region in the Difference of Gaussian (DoG).

	@param[in] value Value to determine if it is an extrema
	@param[in] scales Set of scale images to check
	@param[in] region Region inside of each scale image to check
	@param[in] check_for_min Whether we are determining if determining if the
		value is a minima (true) or maxima (false) of the region
	@retval is_extrema Whether value is an extrema
 */
bool
is_extrema_point_DoG( const value_type value, const vector<const Image*>& scales,
		const Image_Region& region, const bool check_for_min )
{
	// if checking to see if value is a minima
	if( check_for_min )
	{
		// check each scale image
		for( unsigned k = 0; k != scales.size(); ++k )
		{
			const Image* img = scales[k];

			// examine neighborhood around each point
			const_cast<Image*>( img )->push_region( region );

			for( size_type i = img->row_begin(); i != img->row_end(); ++i )
			{
				for( size_type j = img->col_begin(); j != img->col_end(); ++j )
				{
					// if value is not a minima, then stop
					if( value >= scales[k]->at(i,j) )
					{
						const_cast<Image*>( scales[k] )->pop_region();
						return( false );
					}
				}
			}

			// reset the region
			const_cast<Image*>( scales[k] )->pop_region();
		}
	}
	else // if checking to see if value is a maxima
	{
		// check each scale image
		for( unsigned k = 0; k != scales.size(); ++k )
		{
			const Image* img = scales[k];

			// examine neighborhood around each point
			const_cast<Image*>( img )->push_region( region );

			for( size_type i = img->row_begin(); i != img->row_end(); ++i )
			{
				for( size_type j = img->col_begin(); j != img->col_end(); ++j )
				{
					// if value is not a maxima, then stop
					if( value <= scales[k]->at(i,j) )
					{
						const_cast<Image*>( scales[k] )->pop_region();
						return( false );
					}
				}
			}

			const_cast<Image*>( scales[k] )->pop_region();
		}
	}

	// value is an extrema for this set of scale images
	return( true );
}
