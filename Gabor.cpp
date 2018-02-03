/**
	@file   Gabor.cpp
	@author Wade Spires
	@date   2006/3/30
	@brief  Class Gabor.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Gabor.hpp"

using std::string;
using std::vector;

using namespace ws_img;

typedef Matrix::size_type size_type;

/**
  	Create Gabor filter using given parameters.
 */
Matrix
make_gabor_filter( size_type row_size, size_type col_size, double phi,
		double theta )
{
	// ensure the filter has odd dimension
	if( ws_tools::is_even( row_size ) )
	{
		++row_size;
	}
	if( ws_tools::is_even( col_size ) )
	{
		++col_size;
	}
	Matrix filter( row_size, col_size );

	// compute each value of the filter
	double row_mid = static_cast<double>( filter.row_size() / 2 );
	double col_mid = static_cast<double>( filter.col_size() / 2 );
	for( size_type i = filter.row_begin(); i != filter.row_end(); ++i )
	{
		for( size_type j = filter.col_begin(); j != filter.col_end(); ++j )
		{
			double row = i - row_mid; 
			double col = j - col_mid; 
			double g = gabor( row, col, phi, theta );

			filter(i,j) = g;
		}
	}

	/*
	// TODO determine if we should do this
	// make derivative mask's coefficients should sum to zero
	Matrix::value_type value = filter.sum().sum()
			/ (filter.row_size() * filter.col_size());
	for( size_type i = filter.row_begin(); i != filter.row_end(); ++i )
	{
		for( size_type j = filter.col_begin(); j != filter.col_end(); ++j )
		{
			filter(i,j) -= value;
		}
	}
	 */

	/*
	// TODO see if this formula gives good results
	[X,Y]=meshgrid(tx,ty);

	g = (alpha * beta) / pi * exp( -alpha^2 * (X*cos(theta)+Y*sin(theta)).^2-...
	beta^2*(-X * sin(theta) + Y * cos(theta)).^2 +...
	j * 2 * pi * f0 * (X * cos(theta) + Y * sin(theta)));
	 */

	return( filter );
}

/**
  	Two-dimensional Gabor function.

  	Descriptions for the input parameters follow and basically come from here: 
	http://matlabserver.cs.rug.nl/edgedetectionweb/web/edgedetection_params.html

	@param[in] x Input coordinate
	@param[in] y Input coordinate

	@param[in] phi Phase offset
		The phase offset φ in the argument of the cosine factor of the Gabor
		function is specified in degrees. Valid values are real numbers between
		-180 and 180. The values 0 and 180 correspond to center-symmetric
		'center-on' and 'center-off' functions, respectively, while -90 and 90
		correspond to anti-symmetric functions. All other cases correspond to
		asymmetric functions.  

		If one single value is specified, one convolution per orientation will be
		computed. If a list of values is given (e.g. 0,90 which is default),
		multiple convolutions per orientation will be computed, one for each value
		in the phase offset list.  

	@param[in] theta Orientation
		This parameter specifies the orientation of the normal of the parallel
		stripes of a Gabor function. Its value is specified in degrees. Valid
		values are real numbers between 0 and 360.  

		For one single convolution, enter one orientation value and set the value
		of the last parameter in the block "number of orientations" to 1.  

		If "number of orientations" is set to an integer value N, N >= 1, then N
		convolutions will be computed. The orientations of the corresponding Gabor
		functions are equidistantly distributed between 0 and 360 degrees in
		increments of 360/N, starting from the value specified in
		"orientation(s)".  

		An alternative way of computing multiple convolutions for different
		orientations is to specify under "orientation(s)" a list of values
		separated by commas (e.g. 0,45,110). In this case, the value of the
		parameter "number of orientations" is ignored.  

	@param[in] gamma Aspect ratio
		This parameter, called more precisely the spatial aspect ratio, specifies
		the ellipticity of the support of the Gabor function. For γ = 1, the
		support is circular. For γ < 1 the support is elongated in orientation of
		the parallel stripes of the function. Default value is γ = 0.5.  

	@param[in] lambda Wavelength
		This is the wavelength of the cosine factor of the Gabor filter kernel and
		herewith the preferred wavelength of this filter. Its value is specified
		in pixels. Valid values are real numbers equal to or greater than 2. The
		value λ=2 should not be used in combination with phase offset φ=-90 or
		φ=90 because in these cases the Gabor function is sampled in its zero
		crossings.  To prevent the occurence of undesired effects at the image
		borders the wavelength value should be smaller than one fifth of the
		input image size.  

	@param[in] b Bandwidth
		The half-response spatial frequency bandwidth b (in octaves) of a Gabor
		filter is related to the ratio σ / λ, where σ and λ are the standard
		deviation of the Gaussian factor of the Gabor function and the preferred
		wavelength.

		The value of σ cannot be specified directly. It can only be changed
		through the bandwidth b.  

		The bandwidth value must specified as a real positive number. Default is
		1, in which case σ and λ are connected as follows: σ = 0.56 λ.  

		The smaller the bandwidth, the larger σ, the support of the Gabor
		function and the number of visible parallel excitatory and inhibitory
		stripe zones.  

	@param[in] xi Center of receptive field
	@param[in] eta Center of receptive field

	@retval g Value of Gabor function at the given position and with the given
		parameters 
 */
double
gabor( double x, double y, double phi, double theta, double gamma,
		double lambda, double b, double xi, double eta ) 
{
	// convert to radians
	theta *= PI / 180;
	phi   *= PI / 180;

	double x_prime = (x - xi) * cos( theta ) - (y - eta) * sin( theta );
	double y_prime = (x - xi) * sin( theta ) + (y - eta) * cos( theta );

	double sigma = lambda * (1 / PI) * sqrt( log(2) / 2 )
			* (pow(2,b) + 1) / (pow(2,b) - 1);

	double exp_num = (x_prime * x_prime) + (y_prime * y_prime) * (gamma * gamma);
	double exp_den = 2 * sigma * sigma;
	double cos_val = 2 * PI * (x_prime / lambda) + phi;

	return( exp( -exp_num / exp_den ) * cos( cos_val ) );
}

/**
	My own derivation of the Gabor function based on the simplified formula in
	the paper.
 */
void
make_gabor_filter( size_type row_size, size_type col_size, double sigma,
		double theta, double phi, Matrix& real, Matrix& imag )
{
	theta *= PI / 180;
	phi   *= PI / 180;

	if( ws_tools::is_even( row_size ) )
	{
		++row_size;
	}
	if( ws_tools::is_even( col_size ) )
	{
		++col_size;
	}
	real.resize( row_size, col_size, false );
	imag.resize( row_size, col_size, false );

	const double PI_2 = -2 * PI;

	// compute each value of the filter using 1D Gaussian (note that the Gaussian
	// is symmetric about its mean)
	double row_mid = static_cast<double>( real.row_size() / 2 );
	double col_mid = static_cast<double>( real.col_size() / 2 );
	for( size_type i = real.row_begin(); i != real.row_end(); ++i )
	{
		for( size_type j = real.col_begin(); j != real.col_end(); ++j )
		{
			double row = i - row_mid; 
			double col = j - col_mid; 

			double cos_theta_x = cos( PI_2 * theta * col );
			double sin_theta_x = sin( PI_2 * theta * col );

			double cos_phi_y = cos( PI_2 * phi * row );
			double sin_phi_y = sin( PI_2 * phi * row );

			double g = gauss_2D( col, row, sigma );

			real(i,j) = g * (cos_theta_x * cos_phi_y) - (sin_theta_x * sin_phi_y);
			imag(i,j) = g * (cos_theta_x * sin_phi_y) - (sin_theta_x * cos_phi_y );
		}
	}
}
