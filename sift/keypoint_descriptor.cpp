/**
	@file   keypoint_descriptor.cpp
	@author Wade Spires
	@date   2006/2/22
	@brief  Class keypoint_descriptor.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "keypoint_descriptor.hpp"

using std::string;
using std::vector;

// my namespace
using namespace ws_img;

typedef Image::size_type   size_type;
typedef Image::value_type  value_type;
typedef Image::pixel_type  pixel_type;
typedef Image_Coordinate   Img_Coord;

const ws_img::Matrix::size_type      DESC_NUM_BINS  = 8;
const ws_img::Matrix::value_type     DESC_BIN_WIDTH = (2 * PI) / NUM_BINS;
const std::vector<Vector>::size_type DESC_NUM_HISTS = 16;

// filter for all histogram filtering
// ws_img::Matrix Gauss;

/**
	4. Keypoint descriptor: The local image gradients are measured at the
	selected scale in the region around each keypoint. These are transformed
	into a representation that allows for significant levels of local shape
	distortion and change in illumination.  

	@param[in,out] scale_space
	@param[in,out] keypoints
 */
Keypoint_List
keypoint_descriptor( const Scale_Space& scale_space,
	const Keypoint_List& keypoints )
{
	log_msg( "\nKeypoint descriptor\n" );

	Keypoint_List new_keypoints;

	const Matrix::size_type num_bins = DESC_NUM_BINS;

	// dimension of region to construct histogram
	const size_type reg_size = 16;

	Gauss = make_filter( reg_size );

	// examine each keypoint
	for( unsigned k = 0; k != keypoints.size(); ++k )
	{
		const Keypoint&  key    = keypoints[k];
		const Matrix&    mag    = key.grad_mag;
		const Matrix&    dir    = key.grad_dir;

		// create histogram around point (4-by-4 superimposed on 16-by-16 region)
		Keypoint new_key = key;
		vector<Vector> hists = make_histograms( mag, dir, num_bins,
				DESC_NUM_HISTS );

		// set histogram as the keypoint descriptor
		new_key.descriptors = Vector( num_bins * DESC_NUM_HISTS );
		size_type k = 0;
		for( size_type i = 0; i != hists.size(); ++i )
		{
			for( size_type j = hists[i].vec_begin(); j != hists[i].vec_end(); ++j )
			{
				new_key.descriptors(k++) = hists[i](j);
			}
		}
		// new_key.descriptors /= sqrt( new_key.descriptors.sum() );  // normalize
		new_keypoints.push_back( new_key );
	}

	log_msg( "   Keypoints: %u, Old keypoints: %u\n",
			new_keypoints.size(), keypoints.size() );


	return( new_keypoints );
}

/**
	Compute gradient magnitude and direction.

	@param[in] grad_mag Gradient magnitude
	@param[in] grad_dir Gradient direction
	@param[in] num_bins Number of bins in histogram
	@param[in] hist_dim Number of histograms
	@retval hists Histograms
 */
vector<Vector>
make_histograms( const Matrix& grad_mag, const Matrix& grad_dir,
	const Matrix::size_type num_bins, unsigned num_hists )
{
	// number of angles that are mapped to each bin
	double bin_width = (2 * PI) / num_bins;

	// allocate histograms in a grid formation
	unsigned hists_per_row = static_cast<unsigned>( sqrt( num_hists ) );
	vector<Vector> hists( num_hists, Vector(num_bins) );

// gauss.write();
// fprintf( stderr, "%u %u\n", Gauss.row_size(), Gauss.col_size() );
// err_quit( "\n" );

	// current histogram position and region of interest
	unsigned k = 0;
	Image_Region reg( 0, 0, hists_per_row, hists_per_row );

	for( size_type r = 0; r != hists_per_row; ++r )
	{
		// move region for computing histogram
		reg.set_row_begin( r * hists_per_row );

		for( size_type s = 0; s != hists_per_row; ++s )
		{
			reg.set_col_begin( s * hists_per_row );

			// construct histogram of angles in region
			for( size_type i = reg.row_begin(); i != reg.row_end(); ++i )
			{
				for( size_type j = reg.col_begin(); j != reg.col_end(); ++j )
				{
					value_type hist_count = grad_mag(i,j);
					// value_type hist_count = grad_mag(i,j) * Gauss(i,j);
					value_type theta      = grad_dir(i,j);

					// increment histogram bin
					size_type bin =
							static_cast<size_type>( round( theta / bin_width ) );
					hists[k]( bin ) += hist_count;
				}
			}

			++k;
		}
	}

// hist.write();
// err_quit( "\n" );

	return( hists );
}
