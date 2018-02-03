/**
	@file   Scale_Space.cpp
	@author Wade Spires
	@date   2006/1/28
	@brief  Class Scale_Space.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Scale_Space.hpp"

using std::min;
using std::string;
using std::vector;

using namespace ws_img;

// initialize constant, static values
const double    Scale_Space::INITIAL_SIGMA        = .5;
const double    Scale_Space::DEFAULT_SIGMA        = 1.5;
const unsigned  Scale_Space::DEFAULT_NUM_SCALES   = 3;
const unsigned  Scale_Space::DEFAULT_NUM_OCTAVES  = 1;

typedef Image::size_type size_type;

/**
	Construct scale space using given image blurred by sigma.

	@param[in] scales_per_octave (s + 3)
 */
Scale_Space::Scale_Space( const Image& in_img, double sigma,
	unsigned num_octaves, unsigned scales_per_octave )
: _sigma( sigma )
{
	log_msg( "\nConstructing scale-space\n" );

	Image img = in_img;  // make local copy of image to work on

	// convert to gray-scale and scale pixel values to the range [0,1]
	img.to_gray();
	// img.normalize( 1 );
	img.normalize();

	// produce s + 3 images in the stack of blurred images for each octave
	// where scales_per_octave = s + 3
	if( scales_per_octave < 3 )
	{
		err_quit( "scales_per_octave: %u\n", scales_per_octave );
	}

	Vector filter;  // Gaussian filter

	// double the image size and smooth with initial Gaussian
	img.scale( 2 );
	filter = make_gauss_filter( INITIAL_SIGMA );
	// filter = make_gauss_filter( 1 + INITIAL_SIGMA );

// move the boundary of the image inside to avoid spurious edges from
// filtering along the border 
// Image_Region new_reg = img.get_region();
// new_reg.set_row_begin( new_reg.row_begin() + (filter.size() / 2) );
// new_reg.set_col_begin( new_reg.col_begin() + (filter.size() / 2) );
// new_reg.set_row_size( new_reg.row_size() - filter.size() - 1 );
// new_reg.set_col_size( new_reg.col_size() - filter.size() - 1 );
// img.push_region( new_reg );
	img.convolve( filter, filter );
// img.pop_region();

	// for each octave
	// note: we represent both octaves and scale images so as to not break
	// compatibility with previously written code that does not search for extrema
	// across octaves, only across scale images within an octave
	for( unsigned i = 0; i != num_octaves; ++i )
	{
		// add each octave to scale-space
		push_back( Octave() );

		log_msg( "Octave %u / %u:\n", i, size() - 1 );

		// for each pair of images in the octave
		for( unsigned j = 0; j != scales_per_octave; ++j )
		{
			log_msg( "   Scale: %u\n", j );

			// double current_sigma = INITIAL_SIGMA;
			double current_sigma = DEFAULT_SIGMA;

			// first image of octave: use previous image
			Scale_Space_Image first_img;
			first_img.L = img;
			first_img.s = current_sigma;

			// convolve image with Gaussian
			Vector filter = make_gauss_filter( current_sigma );
			if( filter.size() > img.sub_row_size()
					|| filter.size() > img.sub_col_size() )
			{
				break;
			}
			first_img.L.convolve( filter, filter );

			// second image of octave: filtered version of first
			Scale_Space_Image second_img;
			second_img.L = first_img.L;
			second_img.s = 1.44 * current_sigma;  // sqrt(2) * sigma

			// convolve image with Gaussian
			second_img.L.convolve( filter, filter );

			// create difference of Gaussians
			first_img.D = first_img.L - second_img.L;

			// the first image of the next scale set is the second image after
			// downsampling by 2
			img = second_img.L;
			img.scale( .5, .5 );  // use interpolation when scaling
			// scale_img( img, .5, .5 );  // no interpolation

			// add only the first image to the scale-space
			back().push_back( first_img );
		}

		// add empty to coincide with previously written code for extrema detection
		// that expects a final level that has no DoG image
		Scale_Space_Image unused_img;
		back().push_back( unused_img );

		img = back()[ back().size() - 2 ].L; 
	}

	// rescale all images so that they have the same dimension
	for( unsigned i = 0; i != size(); ++i )
	{
		// downsample first image in octave to match the original's size
		double factor = .5;
		Scale_Space_Image& scale_img = at(i)[0];
		scale_img.L.scale(factor, factor);
		scale_img.D.scale(factor, factor);

		// upsample all remaining images to match the original's size
		// (skip the image in octave[1] since it has the same size already)
		Octave& octave = at(i);
		factor = 2;
		for( unsigned j = 2; j != (octave.size() - 1); ++j )
		{
			Scale_Space_Image& scale_img = octave[j];
			scale_img.L.scale(factor, factor);
			scale_img.D.scale(factor, factor);
			factor *= 2;
		}
	}
}

#if false // first method that computes multiple octaves and scale images

/**
	Construct scale space using given image blurred by sigma.

	@param[in] scales_per_octave (s + 3)
 */
Scale_Space::Scale_Space( const Image& in_img, double sigma,
	unsigned num_octaves, unsigned scales_per_octave )
: _sigma( sigma )
{
	Image img = in_img;  // make local copy of image to work on

	// convert to gray-scale and scale pixel values to the range [0,1]
	img.to_gray();
	img.normalize( 1 );

	// produce s + 3 images in the stack of blurred images for each octave
	// where scales_per_octave = s + 3
	if( scales_per_octave < 4 )
	{
		err_quit( "scales_per_octave: %u\n", scales_per_octave );
	}

	// double the image size and smooth with initial Gaussian
	img.scale( 2 );
	Vector filter = make_gauss_filter( INITIAL_SIGMA );
	// Vector filter = make_gauss_filter(
			// sqrt( (_sigma * sigma) - (4 * INITIAL_SIGMA * INITIAL_SIGMA) ) );
	if( filter.size() > img.sub_row_size() || filter.size() > img.sub_col_size() )
	{
		return;
	}
	img.convolve( filter, filter );

	// determine the actual number of octaves to use using input parameters
	num_octaves = get_num_octaves( img, num_octaves, scales_per_octave );

	// allocate scale-space
	for( unsigned i = 0; i != num_octaves; ++i )
	{
		// add each octave to scale-space
		push_back( Octave() );

		// add each scale-space image to the octave
		for( unsigned j = 0; j != scales_per_octave; ++j )
		{
			back().push_back( Scale_Space_Image() );
		}
	}

	// k = 2^(1/s) where scales_per_octave = s + 3
	double k = pow( 2, 1.0 / (scales_per_octave - 3) ); 
	double constant_factor = sqrt( (k * k) - 1.0 );

fprintf( stderr, "Scale factor: %lf\n", k );

	double current_sigma = INITIAL_SIGMA;

	// 1. smooth with larger and larger Gaussians

	/*
	// for each octave
	for( unsigned i = 0; i != size(); ++i )
	{
		// add another octave to scale-space
		Octave& octave = at(i);

fprintf( stderr, "Octave %u / %u:\n   Sigma: %lf\n",
i, size() - 1, current_sigma );

		// add first image of octave
		Scale_Space_Image& first_scale_img = octave.front();

		first_scale_img.s = current_sigma;
		first_scale_img.L = img;

		// initialize previous scale-space image
		Scale_Space_Image* prev_scale_img = &first_scale_img;

		// current_sigma = k;
		current_sigma = _sigma;  // reset sigma

		// for each image in the octave
		for( unsigned j = 1; j != octave.size(); ++j )
		{
fprintf( stderr, "   Sigma: %lf\n", current_sigma );

			// add another image to octave
			Scale_Space_Image& scale_img = octave[j];
			scale_img.s = current_sigma;

			Vector filter = make_gauss_filter( scale_img.s * constant_factor );
			if( filter.size() > img.sub_row_size()
					|| filter.size() > img.sub_col_size() )
			{
				break;
			}

			scale_img.L = img;  // convolve original image with large Gaussian
			// scale_img.L = prev_scale_img->L;  // convolve previous image
			scale_img.L.convolve( filter, filter );

			// create difference of Gaussians
			prev_scale_img->D = scale_img.L - prev_scale_img->L;

			// save previous image and update sigma
			prev_scale_img = &scale_img;
			current_sigma *= k;
		}

		// for the next octave, the initial image is the image with 2*_sigma
		assert( octave.size() >= 2 );
		// current_sigma = octave[ octave.size() - 2 ].s;
		img = octave[ octave.size() - 2 ].L;
		// img = octave[ octave.size() - 3 ].L;

		// downsample image by 2
		scale_img( img, .5, .5 );  // no interpolation
		// img.scale( .5, .5 );  // use interpolation when scaling
	}
	 */

	// 2. smooth with same-sized Gaussian
	// for each octave
	for( unsigned i = 0; i != size(); ++i )
	{
		// add another octave to scale-space
		Octave& octave = at(i);

fprintf( stderr, "Octave %u / %u:\n   Sigma: %lf\n",
i, size() - 1, current_sigma );

		// add first image of octave
		Scale_Space_Image& first_scale_img = octave.front();

		first_scale_img.s = current_sigma;
		first_scale_img.L = img;

		// create single filter
		Vector filter = make_gauss_filter( first_scale_img.s );
		if( filter.size() > first_scale_img.L.sub_row_size()
				|| filter.size() > first_scale_img.L.sub_col_size() )
		{
			break;
		}

		// filter image
		first_scale_img.L.convolve( filter, filter );

		// initialize previous scale-space image
		Scale_Space_Image* prev_scale_img = &first_scale_img;

		// current_sigma = k;
		current_sigma = _sigma;  // reset sigma

		// for each image in the octave
		for( unsigned j = 1; j != octave.size(); ++j )
		{
			// add another image to octave
			Scale_Space_Image& scale_img = octave[j];
			// scale_img.s = pow(2, j) * current_sigma;
			scale_img.s = k * prev_scale_img->s;

fprintf( stderr, "   Sigma: %lf\n", scale_img.s );

			scale_img.L = prev_scale_img->L;  // convolve previous image

			if( filter.size() > scale_img.L.sub_row_size()
					|| filter.size() > scale_img.L.sub_col_size() )
			{
				break;
			}

			scale_img.L.convolve( filter, filter );

			// create difference of Gaussians
			prev_scale_img->D = scale_img.L - prev_scale_img->L;

// move the region inside to remove convolution effects along the borders
size_type half_size = filter.size() / 2;
prev_scale_img->D.set_region( half_size, half_size,
	prev_scale_img->D.row_size() - 2 * half_size,
	prev_scale_img->D.col_size() - 2 * half_size );

			// save previous image and update sigma
			prev_scale_img = &scale_img;
		}

		// for the next octave, the initial image is the image with 2*_sigma
		assert( octave.size() >= 2 );
		// img = octave[ octave.size() - 3 ].L;
		img = octave[ octave.size() - 2 ].L;

		// downsample image by 2
		scale_img( img, .5, .5 );  // no interpolation
		// img.scale( .5, .5 );  // use interpolation when scaling
	}
	/*
	 */
}

#endif

/**
	Determine the number of octaves to use. The number of octaves to actually
	use must allow 'scales_per_octave' scale images to be filled in each
	octave. This must take into account the image dimensions, the number of
	octaves requested, the number of scales in each octave, and the size of the
	Gaussian filters created.

	@param[in] img Input image
	@param[in] num_octaves Number of octaves requested
	@param[in] scales_per_octave Number of scales per octave
	@retval num_octaves Number of octaves to use
 */
unsigned
Scale_Space::get_num_octaves( const Image& img, unsigned num_octaves,
	unsigned scales_per_octave )
{
	// determine the number of octaves by taking into account the image size
	size_type max_img_octaves = static_cast<size_type>(
			log2( min<size_type>( img.sub_row_size(), img.sub_col_size() ) )
		);
	num_octaves = min<size_type>( max_img_octaves, num_octaves );

	size_type row_size = img.sub_row_size();
	size_type col_size = img.sub_col_size();

	double k = pow( 2, 1.0 / (scales_per_octave - 3) ); 

	// simulate creating the octaves: stop when the Gaussian filter mask is
	// larger than the image size since the image is downsampled with each octave
	for( unsigned i = 0; i != num_octaves; ++i )
	{
		double sigma = _sigma;

		// for each image in the octave
		for( unsigned j = 1; j != scales_per_octave; ++j )
		{
			// create filter; compare sizes; quit if filter is larger than image
			Vector filter = make_gauss_filter( sigma );
			if( filter.size() > row_size || filter.size() > col_size )
			{
				return( i );
			}
			sigma *= k;
		}

		// simulate down-sampling by 2
		row_size /= 2;
		col_size /= 2;
	}

	return( num_octaves );
}

/**
	Scale image by given scale factors.

	No interpolation is performed, so pixels are simply removed.
	Image::scale() does perform interpolation.

	@param[in] s_x Scale factor for x coordinates (columns)
	@param[in] s_y Scale factor for y coordinates (rows)
 */
void
Scale_Space::scale_img( Image& img, double s_x, double s_y )
{
	// allocate new, scaled image with the scaled number of rows and columns
	Image scaled_img( static_cast<Image::size_type>( s_y * img.sub_row_size() ),
			static_cast<Image::size_type>( s_x * img.sub_col_size() ) );

	// for each pixel of the scaled image, calculate the pixel location that
	// maps to it from the original image
	for( Image::size_type i = scaled_img.row_begin();
			i != scaled_img.row_end();
			++i )
	{
		Image::size_type y = static_cast<Image::size_type>( i / s_y );

		for( Image::size_type j = scaled_img.col_begin();
				j != scaled_img.col_end();
				++j )
		{
			Image::size_type x = static_cast<Image::size_type>( j / s_x );

			// if the pixel position is inside original image
			if( y >= img.row_begin() && y < img.row_end()
					&& x >= img.col_begin() && x < img.col_end() )
			{
				scaled_img(i,j) = img( y, x );
			}
		}
	}
	img = scaled_img;
}

/**
	Write scale space images using given base name.
	@param[in] base_name Base name to create images from, including directory and
		image name, e.g., 'out_dir/test_img'.
 */
void
Scale_Space::write( const std::string& base_name )
{
	// output scale space
	for( unsigned i = 0; i != size(); ++i )
	{
		Octave& octave = at(i);

		log_msg( "Writing octave %u\n", i );

		for( unsigned j = 0; j != (octave.size() - 1); ++j )
		{
			log_msg( "Writing scale image %u\n", j );

			Scale_Space_Image& scale_img = octave[j];

			string out_name = base_name + ws_tools::format_string( "_%d_%d_", i, j );

			Image L = scale_img.L;
			L.normalize();
			L.write( out_name + "L.pgm" );

			// write normalized image and matrix
			Image D = scale_img.D;
			D.normalize();
			D.write( out_name + "D.pgm" );
			// Matrix m = D;
			// m.write( out_name + "D.m", 5 );
		}

// fprintf( stderr, "Writing scale image %u\n", (octave.size() - 1) );

		// write last image, which has no D
		// string out_name = base_name + ws_tools::format_string( "_%d_%d_",
				// i, (octave.size() - 1) );

		// Image L = octave[ octave.size() - 1 ].L;
		// L.normalize();
		// L.write( out_name + "L.pgm" );
	}
}
