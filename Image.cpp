/**
   @file   Image.cpp
   @author Wade Spires
   @date   2005/08/02
   @brief  Image class for reading, writing, and processing grayscale PGM and
		color PPM images.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Image.hpp"

#include "Pnm_IO.hpp"
#include "Jpeg_IO.hpp"

using namespace ws_img;

using std::min;  using std::max;
using std::vector;
using std::max_element;
using std::min_element;

// initialize class constants
const Image::pixel_type Image::Black_Pixel   = 0;
const Image::pixel_type Image::White_Pixel   = 255;
const Image::pixel_type Image::INVALID_PIXEL = FLT_MAX;
const Image::pixel_type Image::DEFAULT_MAX   = 255;

typedef Image_Coordinate Img_Coord;

typedef Image::value_type  value_type;
typedef Image::uchar       uchar;

typedef Image::Wavelet_Type Wavelet_Type;
typedef Image::Decomp_Type  Decomp_Type;

/**
	Write contents of raster image to file.
	@param[in] image_file_name Name of image file to write image to
	@param[in] mode File mode to write as: binary or ascii (binary is default)
 */
void
Image::write_pnm( const std::string& image_file_name,
	const std::string mode ) const
{
	Pnm_IO::write( *this, image_file_name, mode );
}

/**
	Read image file.
	@param[in] image_file_name Name of image file to read image from
 */
void
Image::read_pnm( const std::string& image_file_name )
{
	Pnm_IO::read( *this, image_file_name );
}

/**
	Write contents of raster image to file.
	@param[in] image_file_name Name of image file to write image to
 */
void
Image::write_jpg( const std::string& image_file_name ) const
{
	int quality = 80;  // quality of JPEG image
	Jpeg_IO::write( *this, image_file_name, quality );
}

/**
	Read contents of JPEG file into raster image.
	@param[in] image_file_name Name of image file to read image from
 */
void
Image::read_jpg( const std::string& image_file_name )
{
	if( !Jpeg_IO::read( *this, image_file_name ) )
	{
		err_quit( "Unable to read JPEG image '%s'\n", image_file_name.c_str() );
	}
}

/**
	Convert image to bilevel (binary).
	Note: 0 -> white, 1 -> black
 */
Image&
Image::to_bilevel( )
{
	if( is_bilevel() )
	{
		return( *this );
	}
	else if( is_color() )
	{
		to_gray();
	}

	reset_region();
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			// 0 -> white, 1 -> black
			if( at(i,j) > 0 )
			{
				at(i,j) = 0;
			}
			else
			{
				at(i,j) = 1;
			}
		}
	}
	_max_value = 1;
	_img_type = Bilevel;

	return( *this );
}

/**
	Convert image to gray-scale.
 */
Image&
Image::to_gray( )
{
	if( is_gray() )
	{
		return( *this );
	}
	else if( is_bilevel() )
	{
		reset_region();
		for( size_type i = row_begin(); i != row_end(); ++i )
		{
			for( size_type j = col_begin(); j != col_end(); ++j )
			{
				// 0 -> white, 1 -> black
				if( at(i,j) == 0 )
				{
					at(i,j) = White_Pixel;
				}
				else // at(i,j) == 1
				{
					at(i,j) = Black_Pixel;
				}
			}
		}
		_max_value = DEFAULT_MAX;
	}
	else if( is_color() )
	{
		// reweight each pixel from each color channel
		reset_region();
		for( size_type i = row_begin(); i != row_end(); ++i )
		{
			for( size_type j = col_begin(); j != col_end(); ++j )
			{
				at(i,j) = .299 * red(i,j) + .587 * green(i,j) + .114 * blue(i,j);
			}
		}
		_green_data.clear();
		_blue_data.clear();
	}
	_img_type = Gray;

	return( *this );
}

/**
	Convert image to color
 */
Image&
Image::to_color( )
{
	if( is_color() )
	{
		return( *this );
	}
	else if( is_bilevel() )
	{
		to_gray();
	}

	// allocate color channels (gray is stored in red, so it needs no allocation)
	_green_data.resize( real_row_size(), real_col_size(), false );
	_blue_data.resize( real_row_size(), real_col_size(), false );

	// copy gray-scale data into other color channels
	reset_region();
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			green(i,j) = red(i,j);
		}
	}
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			blue(i,j) = red(i,j);
		}
	}
	_img_type = Color;

	return( *this );
}

/**
	Change the currently active color channel.
	@param[in] color_type Color to change to
 */
void
Image::set_channel( const Color_Type& color_type )
{
	if( !is_color() )
	{
		return;
	}
	switch( color_type )
	{
		case Red:
			_current_data = &_red_data;
			break;

		case Green:
			_current_data = &_green_data;
			break;

		case Blue:
			_current_data = &_blue_data;
			break;

		default:
			err_quit( "Unknown color channel code: %u\n", color_type );
			break;
	}
}

/**
	Convert image to integral image representation.
	@retval int_img Integral image
 */
Image
Image::to_integral( )
{
	Image int_img( sub_row_size() + 1, sub_col_size() + 1 );
	return( to_integral( int_img ) );
}

/**
	Convert image to integral image representation.
	@param[out] int_img Integral image
 */
Image&
Image::to_integral( Image& int_img )
{
	if( int_img.sub_row_size() != (sub_row_size() + 1)
		|| int_img.sub_col_size() != (sub_col_size() + 1) )
	{
		int_img = Image( sub_row_size() + 1, sub_col_size() + 1 );
	}

	// initialize first row and column to 0
	for( size_type ii = int_img.row_begin(); ii != int_img.row_end(); ++ii )
	{
		int_img(ii,0) = 0;
	}
	for( size_type jj = int_img.col_begin(); jj != int_img.col_end(); ++jj )
	{
		int_img(0,jj) = 0;
	}

	// sum for all positions
	// i and j are input image positions; ii and jj are integral image positions
	for( size_type i = row_begin(), ii = (int_img.row_begin() + 1);
			i != row_end();
			++i, ++ii )
	{
		for( size_type j = col_begin(), jj = (int_img.col_begin() + 1);
				j != col_end();
				++j, ++jj )
		{
			int_img(ii,jj) = int_img(ii, jj - 1) + int_img(ii - 1, jj) + at(i,j)
				- int_img(ii - 1, jj - 1);
		}
	}

	return( int_img );
}

/*
	TODO
	1. Add Gaussian functions for derivatives
	2. Add Gaussian filtering--use symmetry to speed up
	3. Make filtering separable
	4. Normalize Gaussian values (*255)
 */

/*
	Filter image with a 5x5 Gaussian filter using a fixed sigma.
 */
void Image::filter_gauss( )
{
	Matrix gauss( 1, 5 );
	
	gauss( 0, 0 ) = .05; 
	gauss( 0, 1 ) = .25;
	gauss( 0, 2 ) = .4;
	gauss( 0, 3 ) = .25;
	gauss( 0, 4 ) = .05;

	int rows = (int) real_row_size();
	int cols = (int) real_col_size();

	Matrix tmp( rows, cols );

	// filter in x direction
	for( int i = 0; i != rows; i++ )
	{
		for( int j = 0; j != cols; j++ )
		{
			tmp( i, j ) = 0.0;

			int start_k;
			if( j < 2  )  { start_k = 0; }
			else          { start_k = j - 2; }

			int end_k;
			if( j >= cols - 2 )  { end_k = cols - 1; }
			else                 { end_k = j + 2; }
			
			for( int k = start_k; k <= end_k; k++ )
			{
				tmp(i, j) += at(i, k) * gauss( 0, k - start_k );
			}
		}
	}

	// filter in y direction
	for( int i = 0; i != rows; i++ )
	{
		for( int j = 0; j != cols; j++ )
		{
			at(i, j) = 0.0;

			int start_k;
			if( i < 2  )  { start_k = 0; }
			else          { start_k = i - 2; }

			int end_k;
			if( i >= rows - 2 )  { end_k = rows - 1; }
			else                 { end_k = i + 2; }
			
			for( int k = start_k; k <= end_k; k++ )
			{
				at(i, j) += tmp(k, j) * gauss( 0, k - start_k );
			}
		}
	}
}

/**
	Convolve image with given filter.

	@param[in] row_filter Filter to convolve along each row with
 */
Image&
Image::convolve_row( const Vector& row_filter, Image& out )
{
	_red_data.convolve_row( row_filter, out._red_data );
	if( is_color() )
	{
		_green_data.convolve_row( row_filter, out._green_data );
		_blue_data.convolve_row( row_filter, out._blue_data );
	}
	return( out );
}

/**
	Convolve image with given filter.

	@param[in] row_filter Filter to convolve along each row with
 */
Image&
Image::convolve_row( const Vector& row_filter )
{
	_red_data.convolve_row( row_filter );
	if( is_color() )
	{
		_green_data.convolve_row( row_filter );
		_blue_data.convolve_row( row_filter );
	}
	return( *this );
}

/**
	Convolve image with given filter.

	@param[in] col_filter Filter to convolve along each row with
 */
Image&
Image::convolve_col( const Vector& col_filter, Image& out )
{
	_red_data.convolve_col( col_filter, out._red_data );
	if( is_color() )
	{
		_green_data.convolve_col( col_filter, out._green_data );
		_blue_data.convolve_col( col_filter, out._blue_data );
	}
	return( *this );
}

/**
	Convolve image with given filter.

	@param[in] col_filter Filter to convolve along each row with
 */
Image&
Image::convolve_col( const Vector& col_filter )
{
	_red_data.convolve_col( col_filter );
	if( is_color() )
	{
		_green_data.convolve_col( col_filter );
		_blue_data.convolve_col( col_filter );
	}
	return( *this );
}

/**
	Convolve image with given separable filters.

	@param[in] row_filter Filter to convolve along each row with
	@param[in] col_filter Filter to convolve along each column with
 */
Image&
Image::convolve( const Vector& row_filter, const Vector& col_filter )
{
	_red_data.convolve( row_filter, col_filter );
	if( is_color() )
	{
		_green_data.convolve( row_filter, col_filter );
		_blue_data.convolve( row_filter, col_filter );
	}
	return( *this );
}

/**
	Convolve image with given filter.

	@param[in] filter Filter to convolve image with
 */
Image&
Image::convolve( const Matrix& filter )
{
	_red_data.convolve( filter );
	if( is_color() )
	{
		_green_data.convolve( filter );
		_blue_data.convolve( filter );
	}
	return( *this );
}

/**
	Optical flow using Lucas-Kanade method.
	@param[in] frames Sequence of image frames
	@param[out] U_frames Optical flow in X direction for each pair of frames
	@param[out] V_frames Optical flow in Y direction for each pair of frames
void
lucas_kanade( const vector<Image>& frames,
	vector<Matrix>& U_frames, vector<Matrix>& V_frames )
{
	// clear optical flow vectors
	U_frames = vector<Matrix>();
	V_frames = vector<Matrix>();

	if( frames.size() >= 2 )
	{
		for( unsigned i = 0; i != (frames.size() - 1); ++i )
		{
			const Image& frame_1 = frames[i];
			const Image& frame_2 = frames[i + 1];

			Matrix U, V;
			lucas_kanade( frame_1, frame_2, U, V )

			U_frames.push_back( U );
			V_frames.push_back( V );
		}
	}
}
 */

/**
	Optical flow using Lucas-Kanade method.
	@param[in] frame_1 First frame of sequence
	@param[in] frame_2 Second frame of sequence
	@param[out] U Optical flow in X direction
	@param[out] V Optical flow in Y direction
void
lucas_kanade( const Image& frame_1, const Image& frame_2,
	Matrix& U, Matrix& V )
{
	if( frame_1.sub_row_size() != frame_2.sub_row_size() )
	{
		return;
	}
	if( frame_1.sub_col_size() != frame_2.sub_col_size() )
	{
		return;
	}

	U.resize( frame_1.sub_row_size(), frame_1.sub_col_size(), false );
	V.resize( frame_1.sub_row_size(), frame_1.sub_col_size(), false );


	// compute window size

	// convert to gray-scale in case of color image
	// to_gray();

	// calculate Ix, Iy, It
	// Ix, Iy = average of image derivatives after gaussian smoothing
	// It = difference of image pixels without smoothing
	for( size_type i = 0; i != Ix.real_row_size(); ++i )
	{
		for( size_type j = 0; j != Ix.col_size(); ++j )
		{
			Ix(i,j) = (Ix_1(i,j) + Ix_2(i,j)) / 2;
			Iy(i,j) = (Iy_1(i,j) + Iy_2(i,j)) / 2;
			It(i,j) = (frame_1(i,j) - frame_2(i,j));
		}
	}

	size_type patch_row_size = 3;
	size_type patch_col_size = 3;
}
 */

/*
function [u v] = LKOpticalFlow(imfile1, imfile2, wsize )
	im1 = double( imread(imfile1) );
	im2 = double( imread(imfile2) );

	% compute window size
	[height width] = size(imq);
	halfw = floor( wsize / 2 );
	hs = halfw + 1;
	he = height - halfw;
	ws = halfw + 1;
	we = width - halfw;

	% smooth the image and calculate image gradient
	gfilter = fspecial( 'gaussian', [10, 10], 1 );
	imsmooth1 = conv2( im1, gfilter, 'same' );
	imsmooth2 = conv2( im2, gfilter, 'same' );
	[Ix1, Iy1] = gradient( im1 );
	[Ix2, Iy2] = gradient( im2 );

%OR

	% make a kernel for partial derivatives of the gaussian in x and y directions 
	kernel_size = 6 * sigma + 1;
	k = (kernel_size - 1) / 2;
	gauss_kenel_x = zeros( kernel_size, kernel_size );
	for i = 1:kernel_size
		for j = 1:kernel_size
			gauss_kernel_x(i,j) = -( (j - k - 2) / (2 * pi * sigma^3)
				* exp( -( (i - k - 1)^2 + (j - k - 1)^2) / (2 * sigma^2) ));
		end
	end
	gauss_kenel_y = gauss_kenel_x;


	% calculate matrix components
	Ix = (Ix1 + Ix2) / 2;  // average of derivatives
	Ix = (Ix1 + Ix2) / 2;
	It = im2 - im1;

	Ix2  = Ix.*Ix;
	Iy2  = Iy.*Iy;
	IxIy = Ix.*Iy;
	IxIt = Ix.*It;
	IyIt = Iy.*It;

	% calculate optical flow
	u = zeros( height, width );
	v = u;
	for row = hs:he
		for col = ws:we
			a = sum( sum( Ix2( row - halfw:row + halfw, col - halfw:col + halfw ) ) );
			b = sum( sum( Iy2( row - halfw:row + halfw, col - halfw:col + halfw ) ) );
			c = sum( sum( IxIy( row - halfw:row + halfw, col - halfw:col + halfw ) ) );
			A = -sum( sum( IxIt( row - halfw:row + halfw, col - halfw:col + halfw ) ) );
			B = -sum( sum( IyIt( row - halfw:row + halfw, col - halfw:col + halfw ) ) );
			C = [a, c; c, b];
			tmp = pinv(C) * [A; B];
			u(row, col) = tmp(1);
			v(row, col) = tmp(2);
		end
	end

% OR

	nbh
	for i = (1 + floor(nbd_size / 2)):(height - floor(nbd_size / 2) 
		for j = (1 + floor(nbd_size / 2):(width - floor(ndb_size / 2))

			A = zeros(2,2);
			B = zeros(2,1);

			for m = (i - floor(nbd_size / 2)): (i + floor(nbd_size / 2))
				for m = (j - floor(nbd_size / 2)): (i + floor(nbd_size / 2))

					A(1,1) = A(1,1) + Ix(m,n) * Ix(m,n);
					A(1,2) = A(1,2) + Ix(m,n) * Iy(m,n);

					A(2,1) = A(2,1) + Ix(m,n) * Iy(m,n);
					A(2,2) = A(2,2) + Iy(m,n) * Iy(m,n);

					B(1,1) = B(1,1) + Ix(m,n) * It(m,n);
					B(2,1) = B(2,1) + Iy(m,n) * It(m,n);

				end
			end

			A_inv = pinv(A);
			result = A_inv(-B);
			U(i,j) = result(1,1);
			V(i,j) = result(2,1);
		end
	end


	% visualize the optical flow
	figure quiver( u, v, 0 );
*/

/**
	Transpose image--interchange rows and columns.
	@return this This image transposed
 */
Image Image::transpose( )
{
	Image transposed_img = *this;

	transposed_img._red_data = _red_data.transpose();
	if( is_color() )
	{
		transposed_img._green_data = _green_data.transpose();
		transposed_img._blue_data  = _blue_data.transpose();
	}

	return( transposed_img );
}

/**
	Crop image.

	@retval this This image
 */
Image&
Image::crop( )
{
	return( crop( get_region() ) );
}

/**
	Crop image by angle in degrees. The image is rotated around its center.

	@param[in] angle Angle (in degrees) to rotate image by
	@param[in] do_fast_rotate Whether to do a faster rotation (nearest neighbor)
		or do a slower, more accurate rotation (bilinear interpolation)
	@retval this This image
 */
Image&
Image::crop( const Image_Region& region )
{
	// TODO maybe just set_region() and stop, and provide trim() function to
	// actually remove data

	// operator=() only copies the current region
	set_region( region );
	Image I = *this;
	*this = I;

	return( *this );
}

/**
	Resize the matrix to have the given number of rows and columns.

	@param[in] new_row New number of rows in image
	@param[in] new_col New number of columns in image
	@param[in] interpolate Whether to perform interpolation of image values
		(slower) or just to allocate a larger image (faster but data is lost)
	@retval this This image
 */
Image& Image::resize( size_type new_row, size_type new_col, bool interpolate )
{
	_red_data.resize( new_row, new_col, interpolate );
	if( is_color() )
	{
		_green_data.resize( new_row, new_col, interpolate );
		_blue_data.resize( new_row, new_col, interpolate );
	}
	return( *this );
}

/**
	Scale image by given scale factor. Both the width and the height are scaled
	by the same amount.

	@param[in] s Scale factor for both x and y coordinates (columns and rows)
 */
Image& Image::scale( double s )
{
	return( scale( s, s ) );
}

/**
	Scale image by given scale factors using bicubic interpolation.

	Data is scaled a row at a time.

	@param[in] s_y Scale factor for y coordinates (rows)
	@param[in] s_x Scale factor for x coordinates (columns)
	@retval this This image scaled by the given amount
 */
Image& Image::scale( double s_y, double s_x )
{
	_red_data.scale( s_y, s_x );
	if( is_color() )
	{
		_green_data.scale( s_y, s_x );
		_blue_data.scale( s_y, s_x );
	}
	return( *this );
}

/**
	Scale image by given scale factors using no interpolation.

	Data is scaled a row at a time.

	@param[in] s_y Scale factor for y coordinates (rows)
	@param[in] s_x Scale factor for x coordinates (columns)
	@retval this This image scaled by the given amount
 */
Image& Image::nointerp_scale( double s_y, double s_x )
{
	_red_data.nointerp_scale( s_y, s_x );
	if( is_color() )
	{
		_green_data.nointerp_scale( s_y, s_x );
		_blue_data.nointerp_scale( s_y, s_x );
	}
	return( *this );
}

/**
	Rotate image by angle in degrees. The image is rotated around its center.

	@param[in] angle Angle (in degrees) to rotate image by
	@param[in] do_fast_rotate Whether to do a faster rotation (nearest neighbor)
		or do a slower, more accurate rotation (bilinear interpolation)
 */
Image& Image::rotate( const double angle, const bool do_fast_rotate )
{
	// return( rotate( angle,
				// row_begin() + sub_row_size() / 2,
				// col_begin() + sub_col_size() / 2,
				// do_fast_rotate ) );

	return( rotate( angle, real_row_size() / 2, real_col_size() / 2,
				do_fast_rotate ) );
}

/**
	Rotate image by angle in degrees.

	TODO Provide way to switch between keeping original image size (and cropping
	the result) and resizing to show the complete rotated image (with black
	filled in) 

	TODO move to matrix class

	TODO change operator= to for loop to copy (help with color and regions)

	TODO verify how to work for image regions

	TODO support color

	The algorithms come from "Practical Algorithms for Image Analysis"
	by Seul, etal. (vision/prac_alg/ch_2.4/imrotate/imrotate.c)

	@param[in] angle Angle (in degrees) to rotate image by
	@param[in] origin_y Center of rotation (y-coordinate)
	@param[in] origin_x Center of rotation (x-coordinate)
	@param[in] do_fast_rotate Whether to do a faster rotation (nearest neighbor)
		or do a slower, more accurate rotation (bilinear interpolation)
 */
Image& Image::rotate( const double angle, const size_type origin_y,
		const size_type origin_x, const bool do_fast_rotate )
{
	// convert angle to radians and precompute cos/sin of angle
	const double a = angle * (PI / 180);
	const double cos_angle = cos( a );
	const double sin_angle = sin( a );

	const double oy = (double) origin_y;
	const double ox = (double) origin_x;

	double start_y = -ox * sin_angle - oy * cos_angle + oy;
	double start_x = -ox * cos_angle + oy * sin_angle + ox;

	size_type height = sub_row_size() - 1;
	size_type width  = sub_col_size() - 1;

	// allocate rotated image
	Image out_img( sub_row_size(), sub_col_size() );

	// rotate the entire image, not just the current region
	// push_region( 0, 0, real_row_size(), real_col_size() );

	// do a faster rotation (nearest neighbor)
	if( do_fast_rotate )
	{
		// for every pixel (i,j) in the output image, find its corresponding
		// pixel (y,x) in the input image
		double y = start_y;
		for( size_type i = out_img.row_begin(); i != out_img.row_end(); ++i )
		{
			double x = start_x;
			for( size_type j = out_img.col_begin(); j != out_img.col_end(); ++j )
			{
				// convert y and x to image index positions
				const size_type y_pos = (size_type) y;
				const size_type x_pos = (size_type) x;

				// verify the point is inside the image
				if( y_pos >= 1 && y_pos < height && x_pos >= 1 && x_pos < width )
				{
					// get fractional part of y and x
					const double frac_y = y - (double) y_pos;
					const double frac_x = x - (double) x_pos;

					// get new pixel using nearest neighbor
					if( frac_x < .5 )
					{
						out_img(i,j) = (frac_y < 0.5) ?
							at(y_pos, x_pos) : at(y_pos + 1, x_pos);
					}
					else
					{
						out_img(i,j) = (frac_y < 0.5) ?
							at(y_pos, x_pos + 1) : at(y_pos + 1, x_pos + 1);
					}
				}
				else
				{
					out_img(i,j) = Black_Pixel;
				}

				// get next position in input image
				y += sin_angle;
				x += cos_angle;
			}

			start_y += cos_angle;
			start_x -= sin_angle;
			y = start_y;
		}
	}
	else // do a slower, more accurate rotation (bilinear interpolation)
	{
		// for every pixel (i,j) in the output image, find its corresponding
		// pixel (y,x) in the input image
		double y = start_y;
		for( size_type i = out_img.row_begin(); i != out_img.row_end(); ++i )
		{
			double x = start_x;
			for( size_type j = out_img.col_begin(); j != out_img.col_end(); ++j )
			{
				// convert y and x to image index positions
				const size_type y_pos = (size_type) y;
				const size_type x_pos = (size_type) x;

				// verify the point is inside the image
				if( y_pos >= 1 && y_pos < height && x_pos >= 1 && x_pos < width )
				{
					// get fractional part of y and x
					const double frac_y = y - (double) y_pos;
					const double frac_x = x - (double) x_pos;

					const double fx_fy = frac_x * frac_y;

					// apply bilinear interpolation using 4 neighboring pixels
					out_img(i,j) = (
							(1.0 - frac_x - frac_y + fx_fy) * at(y_pos,x_pos)
							+ (frac_x - fx_fy) * at(y_pos, x_pos + 1)
							+ (frac_y - fx_fy) * at(y_pos + 1, x_pos)
							+ fx_fy * at(y_pos + 1, x_pos + 1) + 0.5
						);
				}
				else
				{
					out_img(i,j) = Black_Pixel;
				}

				// get next position in input image
				y += sin_angle;
				x += cos_angle;
			}
			start_y += cos_angle;
			start_x -= sin_angle;
			y = start_y;
		}

	}

	// change region back to its original
	// pop_region();
	
	return( (*this = out_img) );
}

/**
	Rotate image by angle in radians.

	TODO This version doesn't quite work (it's rather old code), but it does
		solve the todo for the version that works above by letting the user use
		either the original size of the image or expanding it to have the entire
		original image inside it. We should probably incorporate it with the
		above code then.  

	@param angle Angle (in radians) to rotate image by
	@param coord Position to rotate around
	@param
 */
Image& Image::rotate( double angle, const Image_Coordinate& rotate_origin,
	bool preserve_size )
{
	double cos_angle = cos( angle );  // precompute cos/sin of angle
	double sin_angle = sin( angle );

	/*
		get size of rotated image:
		rotate each corner of image to find the left-, right-, top-, and
		bottom-most points in rotated image
	 */
	vector<int> row_vec, col_vec;  // place to store
	int row, col;

 	// get each corner point
	Img_Coord top_left( 0, 0 );
	Img_Coord top_right( 0, real_col_size() - 1 );
	Img_Coord bottom_left( real_row_size() - 1, 0 );
	Img_Coord bottom_right( real_row_size() - 1, real_col_size() - 1 );

	// rotate each point
	rotate( top_left, rotate_origin, cos_angle, sin_angle, &row, &col );
	row_vec.push_back( row );  col_vec.push_back( col );

	rotate( top_right, rotate_origin, cos_angle, sin_angle, &row, &col );
	row_vec.push_back( row );  col_vec.push_back( col );

	rotate( bottom_left, rotate_origin, cos_angle, sin_angle, &row, &col );
	row_vec.push_back( row );  col_vec.push_back( col );

	rotate( bottom_right, rotate_origin, cos_angle, sin_angle, &row, &col );
	row_vec.push_back( row );  col_vec.push_back( col );

	// find extent of rotated image
	int row_start = *min_element( row_vec.begin(), row_vec.end() );
	int col_start = *min_element( col_vec.begin(), col_vec.end() );
	int row_end   = *max_element( row_vec.begin(), row_vec.end() );
	int col_end   = *max_element( col_vec.begin(), col_vec.end() );

	fprintf( stderr, "%d %d %d %d\n",
		row_start, col_start, row_end - row_start, col_end - col_start );

	// for each image in input image
	Image new_img( row_end - row_start, col_end - col_start );
	for( unsigned i = 0; i != new_img.real_row_size(); i++ )
	{
		for( unsigned j = 0; j != new_img.real_col_size(); j++ )
		{
			rotate( Img_Coord(i, j), rotate_origin, cos_angle, sin_angle,
				&row, &col );

			int new_i = i - row_start;  // shift point so we have entire image
			int new_j = j - col_start;
			int old_row = row;
			int old_col = col;

			if( old_col >= 0 && old_col < int( real_col_size() )
				&& old_row >= 0 && old_row < int( real_row_size() ) )
 /*
			if( new_i >= 0 && new_i < int( new_img.real_row_size() )
				&& new_j >= 0 && new_j < int( new_img.real_col_size() )
				&& col >= 0 && col < int( real_col_size() )
				&& row >= 0 && row < int( real_row_size() ) )
  */
			{
				if( new_i < 0 || new_i >= int( new_img.real_row_size() )
					|| new_j < 0 || new_j >= int( new_img.real_col_size() ) )
				{
				/*
fprintf( stderr, "(%d %d) (%d %d)\n", i, j,
	new_img.real_row_size(), new_img.real_col_size() );
fprintf( stderr, "(%d %d) (%d %d)\n", new_i, new_j,
	new_img.real_row_size(), new_img.real_col_size() );
fprintf( stderr, "(%d %d) (%d %d)\n\n", row, col,
	real_row_size(), real_col_size() );
				 */
				}
				else
				{
					assert( new_i >= 0 && new_i < int( new_img.real_row_size() )
						&& new_j >= 0 && new_j < int( new_img.real_col_size() )
						&& col >= 0 && col < int( real_col_size() )
						&& row >= 0 && row < int( real_row_size() ) );

					new_img(new_i, new_j) = at(old_row, old_col);
				}
			}
			else
			{
			/*
		fprintf( stderr, "(%d %d) (%d %d)\n", i, j,
			new_img.real_row_size(), new_img.real_col_size() );
		fprintf( stderr, "(%d %d) (%d %d)\n", new_i, new_j,
			new_img.real_row_size(), new_img.real_col_size() );
		fprintf( stderr, "(%d %d) (%d %d)\n\n",
			row, col, real_row_size(), real_col_size() );
			 */
			}
		}
	}
/*
 */


/*
	// for each image in input image
	Image new_img( real_row_size(), real_col_size() );
	for( unsigned i_0 = 0; i_0 != real_row_size(); i_0++ )
	{
		for( unsigned j_0 = 0; j_0 != real_col_size(); j_0++ )
		{
			// shift pixel position in order to rotate around (0,0) instead of
			// rotation point
			int i = i_0 - y_0;
			int j = j_0 - x_0;

			// rotate around (0,0)
			int x_new = int(  (j * cos_angle) + (i * sin_angle) );
			int y_new = int( -(j * sin_angle) + (i * cos_angle) );

			// reshift rotated position so that rotation is around (x_0,y_0)
			int x = int( x_new + x_0 );
			int y = int( y_new + y_0 );

			// only store pixel if its in image's bounds
			if( y >= 0 && unsigned(y) < real_row_size()
				&& x >= 0 && unsigned(x) < real_col_size() )
			{
				new_img(y, x) = at(i_0, j_0);
			}
		}
	}
 */

/*
	if( preserve_size )
	{
		// must take into account all possible rotations though--can't just pick
		// top left corners of original image necessarily
		Image_Region region( left-most point in new_img, right-most point in new_img,
			real_row_size(), real_col_size() );
		*this = new_img.sub_image( region );
		return( );
	}
	else
	{
		*this = new_img;
	}
 */

	*this = new_img;

	return( *this );
}

/**
	Rotate point around given given point.

	@param coord
	@param rotate_origin
	@param cos_angle
	@param sin_angle
	@param rotate_row
	@param rotate_col
 */
void Image::rotate( const Img_Coord& coord, const Img_Coord& rotate_origin,
	double cos_angle, double sin_angle, int* rotate_row, int* rotate_col )
{
	// shift pixel position in order to rotate around (0,0) instead of
	// rotation point
	int i = int(coord.row) - rotate_origin.row;
	int j = int(coord.col) - rotate_origin.col;

	// rotate around (0,0)
	int y_new = int( -(j * sin_angle) + (i * cos_angle) );
	int x_new = int(  (j * cos_angle) + (i * sin_angle) );

	// reshift rotated position so that rotation is around rotation point
	*rotate_row = y_new + int( rotate_origin.row );
	*rotate_col = x_new + int( rotate_origin.col );
}

/**
	Apply transformation matrix to image.
	Rotate image by angle in radians.
	@param angle Angle (in radians) to rotate image by
 */
Image& Image::transform( double angle )
{
	return( *this );
}

/**
	Add corresponding pixels in this image and right-hand side.
	@param rhs Image to add
	@retval this This image
 */
Image& Image::operator+=( const Image& rhs )
{
	_red_data += rhs._red_data;
	if( is_color() )
	{
		_green_data += rhs._green_data;
		_blue_data  += rhs._blue_data;
	}
	return( *this );
}

/**
	Add corresponding pixels in this image and right-hand side.
	@param rhs Image to add
	@retval img Image containing pixel sums
 */
Image Image::operator+( const Image& rhs ) const
{
	return( Image( *this ) += rhs );
}

/**
	Add value to each component of image.
	@param[in] rhs Value on right-hand side of +
	@retval image Image result 
 */
Image&
Image::operator+=( const value_type rhs )
{
	_red_data += rhs;
	if( is_color() )
	{
		_green_data += rhs;
		_blue_data  += rhs;
	}
	return( *this );
}

/**
	Subtract corresponding pixels in this image and right-hand side.
	@param rhs Image to subtract
	@retval this This image
 */
Image& Image::operator-=( const Image& rhs )
{
	_red_data -= rhs._red_data;
	if( is_color() )
	{
		_green_data -= rhs._green_data;
		_blue_data  -= rhs._blue_data;
	}
	return( *this );
}

/**
	Subtract corresponding pixels in this image and right-hand side.
	@param rhs Image to subtract
	@retval img Image containing pixel sums
 */
Image Image::operator-( const Image& rhs ) const
{
	return( Image( *this ) -= rhs );
}

/**
	Subtract value from each component of image.
	@param[in] rhs Value on right-hand side of -
	@retval image Image result 
 */
Image&
Image::operator-=( const value_type rhs )
{
	_red_data -= rhs;
	if( is_color() )
	{
		_green_data -= rhs;
		_blue_data  -= rhs;
	}
	return( *this );
}

/**
	Multiply this image and right-hand side.
	@param rhs Image to multiply
	@retval this This image
 */
Image& Image::operator*=( const Image& rhs )
{
	_red_data *= rhs._red_data;
	if( is_color() && rhs.is_color() )
	{
		_green_data *= rhs._green_data;
		_blue_data  *= rhs._blue_data;
	}
	return( *this );
}

/**
	Multiply this image and right-hand side.
	@param rhs Image to multiply
	@retval this This image
 */
Image Image::operator*( const Image& rhs ) const
{
	return( Image( *this ) *= rhs );
}

/**
	Multiply all pixels in image by given value.
	img * x
	@param x Value to multiply by
	@retval img Image with pixel values scaled by x
 */
Image& Image::operator*=( const value_type rhs )
{
	_red_data *= rhs;
	if( is_color() )
	{
		_green_data *= rhs;
		_blue_data  *= rhs;
	}
	return( *this );
}

/**
	Multiply all pixels in image by given value.
	img * rhs
	@param rhs Value to multiply by
	@retval img Image with pixel values scaled by x
 */
Image Image::operator*( const value_type rhs ) const
{
	return( Image( *this ) *= rhs );
}

/**
	Multiply all pixels in image by given value.
	@param rhs Value to multiply by
	@retval this Image with pixel values scaled by x
 */
Image operator*( const value_type lhs, Image& rhs )
{
	return( rhs * lhs );
}

/**
	Divide all pixels in image by given value.
	img * rhs
	@param rhs Value to divide by
	@retval img Image with pixel values scaled by rhs
 */
Image& Image::operator/=( const value_type rhs )
{
	*this *= (1 / rhs);  // multiply by reciprocal to perform division
	return( *this );
}

/**
	Divide all pixels in image by given value.
	img * rhs
	@param rhs Value to divide by
	@retval img Image with pixel values scaled by rhs
 */
Image Image::operator/( const value_type rhs ) const
{
	return( Image( *this ) /= rhs );
}

/**
	Compute image statistics--average, variance, and standard deviation.

	@param[out] avg Average value of all pixels in image
	@param[out] var Variance of all pixels in image
	@param[out] std Standard deviation of all pixels in image
 */
void Image::stat( double* avg, double* var, double* std ) const
{
	_red_data.stat( avg, var, std );
	if( is_color() )
	{
		// TODO Figure out best way to perform statistics for color image
		// (maybe allow an optional argument to the color channel to use?)
	}
}

/**
	Find minimum value in image.
	@retval max Minimum value in image.
 */
value_type
Image::get_min( ) const
{
	return( _red_data.get_min() );
}

/**
	Find maximum value in image.
	@retval max Maximum value in image.
 */
value_type
Image::get_max( ) const
{
	return( _red_data.get_max() );
}

/**
	Find minimum and maximum value in image.
	@retval max Maximum value in image.
 */
void
Image::get_min_max( value_type* min, value_type* max ) const
{
	_current_data->get_min_max( min, max );
}

/**
	Normalize pixel values to be in the range [0,new_max_value].
	@param[in] new_max_value New maximum pixel value
	@retval this This image
 */
Image&
Image::normalize( const pixel_type new_max_value )
{
	// set I(x,y) = ((I(x,y) - min) / (max - min)) * max_value
	// find minimum and maximum pixel values
	value_type min_pixel, max_pixel;
	get_min_max( &min_pixel, &max_pixel );

	if( min_pixel == max_pixel )
	{
		return( *this );
	}

	// I(x,y) = ((I(x,y) - min) / (max - min)) * max_value
	pixel_type factor = new_max_value / (max_pixel - min_pixel);
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		for( size_type j = col_begin(); j != col_end(); ++j )
		{
			at(i, j) = (at(i, j) - min_pixel) * factor;
		}
	}

	return( *this );
}

/**
	Return subregion from image around region.
	@param[in] region 
	@param[out] region
 */
void Image::sub_image( const Image_Region& region, Image& new_img )
{
	// if region is not inside image, return empty image
	if( region.row_begin() >= real_row_size()
			|| region.col_begin() > real_col_size() )
	{
		return;
	}

	// use end of image if row or column is out of bounds
	unsigned row_end = min<unsigned>( region.row_begin() + region.row_size(),
		real_row_size() );

	unsigned col_end = min<unsigned>( region.col_begin() + region.col_size(),
		real_col_size() );

	unsigned num_rows = row_end - region.row_begin();
	unsigned num_cols = col_end - region.col_begin();

	// allocate image with needed dimension
	Image tmp_region( num_rows, num_cols );

	// copy region from this image into new image
	for( unsigned i = region.row_begin(), y = 0; i != row_end; ++i, ++y )
	{
		for( unsigned j = region.col_begin(), x = 0; j != col_end; ++j, ++x )
		{
			tmp_region(y,x) = at(i,j);
		}
	}

	// overwrite region last in case *this == region
	new_img = tmp_region;
}

/**
	Draw a pixel at coordinate.

	@param[in] intensity Intensity value to make pixels in line
 */
void Image::draw( const Image_Coordinate& coord,
	const Image::pixel_type& intensity )
{
	if( coord.row < real_row_size() && coord.col < real_col_size() )
	{
		at( coord.row, coord.col ) = intensity;
	}
}

/**
	Draw a line between two coordinates. Note that no requirement is made of the
	two coordinates given: either point may lie outside the image's boundary. If
	a point is not inside the image, only the portion of the line that
	intersects the image is drawn.

	@param[in] start_coord Starting coordinate of line
	@param[in] end_coord Ending coordinate of line
	@param[in] intensity Intensity value to make pixels in line
 */
void Image::draw_line( const Image_Coordinate& start_coord,
	const Image_Coordinate& end_coord, const Image::pixel_type& intensity )
{
	// find left- and right-most points
	Image_Coordinate left, right;
	if( start_coord.col < end_coord.col )
	{
		left  = start_coord;
		right = end_coord;
	}
	else
	{
		left  = end_coord;
		right = start_coord;
	}

	// if line is completely outside of image, don't try to draw it
	if( left.col >= real_col_size() )  return;

	// get change in columns and rows: must cast to double since using unsigned
	double x_diff = double(right.col) - double(left.col);
	double y_diff = double(right.row) - double(left.row);

	// if vertical line, then the slope is infinite so handle as special case
	if( x_diff == 0 )
	{
		// find top- and bottom-most rows
		unsigned row_start = min<unsigned>( start_coord.row, end_coord.row );
		unsigned row_end   = max<unsigned>( start_coord.row, end_coord.row );

		// use end of image if bottom-most row goes passed the end
		row_end = min<unsigned>( row_end + 1, real_row_size() );

		// if line is completely outside of image, don't try to draw it
		if( row_start >= real_row_size() )  return;

		// draw each pixel between the two points
		for( unsigned i = row_start; i != row_end; ++i )
		{
			at( i, start_coord.col ) = intensity;
		}
	}
	else  // draw non-vertical line
	{
		// compute slope and y-intercept of line
		double slope = y_diff / x_diff;
		double y_int = left.row - (slope * left.col);

		// use end of image if right-most column goes off the end
		unsigned col_end = min<unsigned>( right.col + 1, real_col_size() );

		// initialize starting column to left-most column, but move it right if
		// the corresponding row is larger than the image's rows since the line
		// may move up into the image (only if slope > 0, though)
		unsigned col_start = left.col;
		for( ; col_start != col_end; ++col_start )
		{
			unsigned y = unsigned( slope * col_start + y_int );  // y = mx + b

			// if line is inside image now, stop
			if( y < real_row_size() )  break;
		}

		// for each x position between left-most point and right-most point,
		// compute the y position using the line equation and draw each pixel
		for( unsigned x = col_start; x != col_end; ++x )
		{
			unsigned y = unsigned( slope * x + y_int );  // y = mx + b
			at( y, x ) = intensity;
		}

		// At this point, we've drawn a line, but the line might look fragmented
		// if the slope is rather steep. We solve this problem by iterating
		// across the y values and solving for x (the reverse of what
		// we did above). This will fill in the line better.

		// if horizontal line, then line is okay and we don't want to divide by 0
		if( slope == 0 )  return;

		// find top- and bottom-most rows
		unsigned row_start = min<unsigned>( start_coord.row, end_coord.row );
		unsigned row_end   = max<unsigned>( start_coord.row, end_coord.row );

		// use end of image if bottom-most row goes passed the end
		row_end = min<unsigned>( row_end + 1, real_row_size() );

		// if line is outside of image, don't try to draw it
		if( row_start >= real_row_size() )  return;

		// for each y position between top-most point and bottom-most point,
		// compute the x position using the line equation and draw each pixel
		for( unsigned y = row_start; y != row_end; ++y )
		{
			unsigned x = unsigned( (y - y_int) / slope );  // x = (y - b) / m
			at( y, x ) = intensity;
		}
	}
}

/**
	Draw rectangle on image.
	@param[in] rect Rectangle to draw inside image
	@param[in] intensity Intensity value to use for drawing rectangle
 */
void
Image::draw( const Image_Region& rect, const Image::pixel_type& intensity )
{
	// Do not draw rectangle that does not lie inside of this image. This is for
	// two reasons:
	// 1. It respects the image's region since we might draw outside it.
	// 2. It addresses a bug in which -1 gets set to a coordinate which gets cast
	// to a really huge unsigned number.
	// TODO Image_Region r;
	// r.row_begin( max<size_type>( rect.row_begin(), row_begin() ) );
	// r.row_begin( min<size_type>( r.row_begin(), row_end() - 1 ? ) );
	if( !rect.inside( region() ) )
	{
		return;
	}

#if 0  // use draw_line() function for drawing rectangle

	// create 4 corner coordinates of rectangle
	Image_Coordinate top_left( rect.row_begin(), rect.col_begin() );

	Image_Coordinate top_right( rect.row_begin(),
		rect.col_begin() + rect.col_size() - 1 );

	Image_Coordinate bottom_left( rect.row_begin() + rect.row_size() - 1,
		rect.col_begin() );

	Image_Coordinate bottom_right( rect.row_begin() + rect.row_size() - 1,
		rect.col_begin() + rect.col_size() - 1 );

	// draw top, bottom, left, and right sides
	draw_line( top_left,    top_right,    intensity );
	draw_line( bottom_left, bottom_right, intensity );
	draw_line( top_left,    bottom_left,  intensity );
	draw_line( top_right,   bottom_right, intensity );

#else  // draw rectangle directly (faster)

	push_region( rect );

	// draw top and bottom sides
	for( size_type i = row_begin(); i != row_end(); ++i )
	{
		at( i, col_begin() )   = intensity;
		at( i, col_end() - 1 ) = intensity;
	}

	// draw left and right sides
	for( size_type j = col_begin(); j != col_end(); ++j )
	{
		at( row_begin(), j )   = intensity;
		at( row_end() - 1, j ) = intensity;
	}

	pop_region();

#endif
}

/**
	Draw outer white, middle black, and inner white boxes on image with
	dimensions and position of given region.
	@param[in] region Region to draw inside image
 */
void Image::draw2( const Image_Region& region )
{
	// outer box: use this region
	Image_Region outer_region = region;

	// middle box: use region inside outer region
	Image_Region middle_region( region.row_begin() + 1, region.col_begin() + 1,
		region.row_size() - 2, region.col_size() - 2 );

	// inner box: use region inside middle region
	Image_Region inner_region( region.row_begin() + 2, region.col_begin() + 2,
		region.row_size() - 4, region.col_size() - 4 );

	// draw outer, middle, and inner boxes as white, black, and white,
	// respectively
	draw( outer_region,  Image::White_Pixel );
	draw( middle_region, Image::Black_Pixel );
	draw( inner_region,  Image::White_Pixel );
}

/**
	Draw circle in image at given position and radius.
	@param[in] center Center coordinate of circle
	@param[in] radius Radius of circle
	@param[in] intensity Intensity value to use for drawing region
 */
void Image::draw_circle( const Image_Coordinate& center, const double radius,
	const Image::pixel_type& intensity )
{
	if( radius <= 0 )  return;

	// get first and last column of circle's position in image
	unsigned col_start = unsigned( max<double>( center.col - radius, 0 ) );
	unsigned col_end   = unsigned( min<double>( center.col + radius + 1,
				real_col_size() ) );

	// for each x position between left-most point and right-most point,
	// compute the y position using the circle equation and draw each pixel
	for( unsigned x = col_start; x != col_end; ++x )
	{
		// (y - y_0)^2 + (x - x_0)^2 = r^2
		// => y = +/- sqrt( r^2 - (x - x_0)^2 ) + y_0
		double x_squared = (x - center.col) * (x - center.col);
		double square_diff = sqrt( (radius * radius) - x_squared );
		double y_1 =  square_diff + center.row;
		double y_2 = -square_diff + center.row;

		if( y_1 >= 0 && y_1 < real_row_size() )
		{
			at( unsigned(y_1), x ) = intensity;
		}
		if( y_2 >= 0 && y_2 < real_row_size() )
		{
			at( unsigned(y_2), x ) = intensity;
		}
	}

	// At this point, we've drawn a circle, but the line might look like a
	// series of dots. We solve this problem by iterating across the y values
	// and solving for x (the reverse of what we did above). This will fill in
	// the line better.  

	// get first and last row of circle's position in image
	unsigned row_start = unsigned( max<double>( center.row - radius, 0 ) );
	unsigned row_end   = unsigned( min<double>( center.row + radius + 1,
				real_row_size() ) );

	// for each y position between top-most point and bottom-most point,
	// compute the x position using the circle equation and draw each pixel
	for( unsigned y = row_start; y != row_end; ++y )
	{
		// (y - y_0)^2 + (x - x_0)^2 = r^2
		// => x = +/- sqrt( r^2 - (y - y_0)^2 ) + x_0
		double y_squared = (y - center.row) * (y - center.row);
		double square_diff = sqrt( (radius * radius) - y_squared );
		double x_1 =  square_diff + center.col;
		double x_2 = -square_diff + center.col;

		if( x_1 >=0 && x_1 < real_col_size() )
		{
			at( y, unsigned(x_1) ) = intensity;
		}
		if( x_2 >= 0 && x_2 < real_col_size() )
		{
			at( y, unsigned(x_2) ) = intensity;
		}
	}
}

/**
	Apply fourier-transform to the this image.

	The output is not normalized. The transform is computed across each color
	channel individually.  

	@param[out] real_data The real components of the transform
	@param[out] imag_data The imaginary components of the transform
 */
void
Image::fft( Image& real_data, Image& imag_data ) const
{
	_red_data.fft( real_data._red_data, imag_data._red_data );
	if( is_color() )
	{
		_green_data.fft( real_data._green_data, imag_data._green_data );
		_blue_data.fft( real_data._blue_data, imag_data._blue_data );
	}
}

/**
	Apply inverse fourier-transform to the this image.

	The transform is computed across each color channel individually.  

	@param[in] real_data The real components of the transform
	@param[in] imag_data The imaginary components of the transform

	@post This image holds the inverse of the normalized FFT.
 */
void
Image::ifft( Image& real_data, Image& imag_data )
{
	_red_data.ifft( real_data._red_data, imag_data._red_data );
	if( is_color() )
	{
		_green_data.ifft( real_data._green_data, imag_data._green_data );
		_blue_data.ifft( real_data._blue_data, imag_data._blue_data );
	}
}

/**
	Apply lifting scheme to image.
	@param[in] N Number of real vanishing moments
	@param[in] N_tilde Number of dual vanishing moments
	@param[in] num_levels Number of transform levels to apply
 */
void Image::lift( unsigned N, unsigned N_tilde, unsigned num_levels )
{
	Lift_Scheme lift_scheme( N, N_tilde );
	lift( lift_scheme, num_levels );
}

/**
	Apply inverse lifting scheme to image.
 */
void Image::inverse_lift( unsigned N, unsigned N_tilde, unsigned num_levels )
{
	Lift_Scheme lift_scheme( N, N_tilde );
	inverse_lift( lift_scheme, num_levels );
}

/**
	Apply lifting scheme to image.
	@param[in] lift_scheme Object to apply lifting scheme with
 */
void Image::lift( Lift_Scheme& lift_scheme, unsigned num_levels )
{
	lift_scheme.lift( _red_data, num_levels );
	if( is_color() )
	{
		lift_scheme.lift( _green_data, num_levels );
		lift_scheme.lift( _blue_data, num_levels );
	}
}

/**
	Apply inverse lifting scheme to image.
	@param[in] lift_scheme Object to apply lifting scheme with
 */
void Image::inverse_lift( Lift_Scheme& lift_scheme, unsigned num_levels )
{
	lift_scheme.inverse_lift( _red_data, num_levels );
	if( is_color() )
	{
		lift_scheme.inverse_lift( _green_data, num_levels );
		lift_scheme.inverse_lift( _blue_data, num_levels );
	}
}

/**
	Reorder the elements of this image in quadrant form.
	@param[in] lift_scheme Object to apply lifting scheme with
 */
void Image::octave_form( const Lift_Scheme& lift_scheme, unsigned num_levels )
{
	lift_scheme.octave_form( _red_data, num_levels );
	if( is_color() )
	{
		lift_scheme.octave_form( _green_data, num_levels );
		lift_scheme.octave_form( _blue_data, num_levels );
	}
}

/**
	Reorder the elements of this image in quadrant form.
	@param[in] lift_scheme Object to apply lifting scheme with
 */
void Image::inverse_octave_form( const Lift_Scheme& lift_scheme,
	unsigned num_levels )
{
	lift_scheme.inverse_octave_form( _red_data, num_levels );
	if( is_color() )
	{
		lift_scheme.inverse_octave_form( _green_data, num_levels );
		lift_scheme.inverse_octave_form( _blue_data, num_levels );
	}
}

/**
	Get wavelet subbands (from lifting).
	@param[in]
 */
Wavelet_Subbands
Image::get_subbands( const Lift_Scheme& lift_scheme, unsigned num_levels )
{
	return( lift_scheme.get_subbands( _red_data, num_levels ) );
	if( is_color() )
	{
		// TODO Figure out how to handle color
		// lift_scheme.get_subbands( _green_data, num_levels );
		// lift_scheme.get_subbands( _blue_data, num_levels );
	}
}

/**
	Get wavelet subbands from traditional wavelet transform.
	@param[in]
 */
Wavelet_Subbands
Image::get_subbands( const Wavelet& wavelet, unsigned num_levels )
{
	return( wavelet.get_subbands( _red_data, num_levels ) );
	if( is_color() )
	{
		// TODO Figure out how to handle color
		// lift_scheme.get_subbands( _green_data, num_levels );
		// lift_scheme.get_subbands( _blue_data, num_levels );
	}
}

/**
	Perform wavelet transformation on image using given wavelet.

	@param[in] wavelet Wavelet to use for the transformation
	@param[in] decomp_type Type of matrix decomposition to perform
		(e.g., standard or pyramid)

 */
void Image::wavelet_transform( const Wavelet& wavelet,
	Decomp_Type decomp_type )
{
	wavelet.forward_transform( _red_data, decomp_type );
	if( is_color() )
	{
		wavelet.forward_transform( _green_data, decomp_type );
		wavelet.forward_transform( _blue_data, decomp_type );
	}
}

/**
	Perform inverse wavelet transformation on image using given wavelet.

	@param[in] wavelet Wavelet to use for the transformation
	@param[in] decomp_type Type of matrix decomposition to perform
		(e.g., standard or pyramid)

 */
void Image::inverse_wavelet_transform( const Wavelet& wavelet,
	Decomp_Type decomp_type )
{
	wavelet.inverse_transform( _red_data, decomp_type );
	if( is_color() )
	{
		wavelet.inverse_transform( _green_data, decomp_type );
		wavelet.inverse_transform( _blue_data, decomp_type );
	}
}

/**
	Perform wavelet transformation on image using given wavelet type.

	@param[in] wavelet_type Type of wavelet to use for the transformation
		(e.g., Haar, Daub_4, etc.)
	@param[in] decomp_type Type of matrix decomposition to perform
		(e.g., standard or pyramid)

 */
void Image::wavelet_transform( Wavelet::Wavelet_Type wavelet_type,
	Wavelet::Decomp_Type decomp_type )
{
	Wavelet wavelet( wavelet_type );
	wavelet_transform( wavelet, decomp_type );
}

/**
	Perform inverse wavelet transformation on image using given wavelet type.

	@param[in] wavelet_type Type of wavelet to use for the transformation
		(e.g., Haar, Daub_4, etc.)
	@param[in] decomp_type Type of matrix decomposition to perform
		(e.g., standard or pyramid)

 */
void Image::inverse_wavelet_transform( Wavelet_Type wavelet_type,
	Decomp_Type decomp_type )
{
	Wavelet wavelet( wavelet_type );
	inverse_wavelet_transform( wavelet, decomp_type );
}
