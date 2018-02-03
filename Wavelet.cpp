/**
	@file   Wavelet.cpp
	@author Wade Spires
	@date   2005/11/21
	@brief  Class Wavelet.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Wavelet.hpp"

using std::max;
using std::min;
using std::string;
using std::vector;

using namespace ws_img;

typedef Wavelet::value_type   value_type;
typedef Wavelet::size_type    size_type;

typedef Wavelet::Wavelet_Type Wavelet_Type;
typedef Wavelet::Decomp_Type  Decomp_Type;

/**
	Perform wavelet transformation on given matrix using GSL.
	@param[in,out] matrix Matrix to transform (can be arbitrary dimension)
	@param[in] wavelet_type Type of wavelet to use for the transformation
		(e.g., Haar, Daub_4, etc.)
	@param[in] decomp_type Type of matrix decomposition to perform
		(e.g., standard or pyramid)
 */
void
Wavelet::forward_transform( Matrix& matrix, Decomp_Type decomp_type ) const
{
	if( decomp_type == Pyramid )
	{
		// forward_pyramid( matrix );
		wavelet2d_nstransform( matrix, 3 );
	}
	else if( decomp_type == Standard )
	{
		forward_standard( matrix );
	}
}

/**
	Perform wavelet transformation using standard decomposition.
	@param[in,out] matrix Matrix to transform (can be arbitrary dimension)
 */
void
Wavelet::forward_standard( Matrix& matrix ) const
{
	gsl_matrix* m = matrix_to_gsl( matrix );

	// allocate working space
	gsl_wavelet_workspace* workspace = gsl_wavelet_workspace_alloc( m->size1 );

	// forward standard transformation
	int error_code =
		gsl_wavelet2d_transform_matrix_forward( _wavelet, m, workspace );

	gsl_wavelet_workspace_free( workspace );

	if( error_code != GSL_SUCCESS )
	{
		err_quit( "A problem occurred while performing the wavelet transform\n" );
	}

	// copy GSL matrix into my matrix
	gsl_to_matrix( m, matrix );
}

/**
	Perform wavelet transformation using pyramid decomposition.
	@param[in,out] matrix Matrix to transform (can be arbitrary dimension)
 */
void
Wavelet::forward_pyramid( Matrix& matrix ) const
{
	gsl_matrix* m = matrix_to_gsl( matrix );

	// allocate working space
	gsl_wavelet_workspace* workspace = gsl_wavelet_workspace_alloc( m->size1 );

	// forward non-standard (pyramid) transformation
	int error_code =
		gsl_wavelet2d_nstransform_matrix_forward( _wavelet, m, workspace );

	gsl_wavelet_workspace_free( workspace );

	if( error_code != GSL_SUCCESS )
	{
		err_quit( "A problem occurred while performing the wavelet transform\n" );
	}

	// copy GSL matrix into my matrix
	gsl_to_matrix( m, matrix );
}

/**
	Perform wavelet transformation on given matrix using GSL.
	@param[in,out] matrix Matrix to transform (can be arbitrary dimension)
	@param[in] wavelet_type Type of wavelet to use for the transformation
		(e.g., Haar, Daub_4, etc.)
	@param[in] decomp_type Type of matrix decomposition to perform
		(e.g., standard or pyramid)
 */
void
Wavelet::inverse_transform( Matrix& matrix, Decomp_Type decomp_type ) const
{
	if( decomp_type == Pyramid )
	{
		inverse_pyramid( matrix );
	}
	else if( decomp_type == Standard )
	{
		inverse_standard( matrix );
	}
}

/**
	Perform wavelet transformation using standard decomposition.
	@param[in,out] matrix Matrix to transform (can be arbitrary dimension)
 */
void
Wavelet::inverse_standard( Matrix& matrix ) const
{
	gsl_matrix* m = matrix_to_gsl( matrix );

	// allocate working space
	gsl_wavelet_workspace* workspace = gsl_wavelet_workspace_alloc( m->size1 );

	// inverse standard transformation
	int error_code =
		gsl_wavelet2d_transform_matrix_inverse( _wavelet, m, workspace );

	gsl_wavelet_workspace_free( workspace );

	if( error_code != GSL_SUCCESS )
	{
		err_quit( "A problem occurred while performing the wavelet transform\n" );
	}

	// copy GSL matrix into my matrix
	gsl_to_matrix( m, matrix );
}

/**
	Perform wavelet transformation using pyramid decomposition.
	@param[in,out] matrix Matrix to transform (can be arbitrary dimension)
 */
void
Wavelet::inverse_pyramid( Matrix& matrix ) const
{
	gsl_matrix* m = matrix_to_gsl( matrix );

	// allocate working space
	gsl_wavelet_workspace* workspace = gsl_wavelet_workspace_alloc( m->size1 );

	// inverse non-standard (pyramid) transformation
	int error_code =
		gsl_wavelet2d_nstransform_matrix_inverse( _wavelet, m, workspace );

	gsl_wavelet_workspace_free( workspace );

	if( error_code != GSL_SUCCESS )
	{
		err_quit( "A problem occurred while performing the wavelet transform\n" );
	}

	// copy GSL matrix into my matrix
	gsl_to_matrix( m, matrix );
}

/**
	Convert my matrix type to a GSL matrix. The returned GSL matrix is square 
	with a dimension that is a power of 2.
	@param[out] matrix My matrix
	@retval m GSL matrix
 */
gsl_matrix*
Wavelet::matrix_to_gsl( const Matrix& matrix ) const
{
	// new dimension for matrix is the smallest power of 2 >= current dimension;
	// current dimension is taken to be the largest size between the number of
	// rows and columns since the new matrix must be square
	size_type max_dimension =
			std::max<size_type>( matrix.sub_row_size(), matrix.sub_col_size() );
	size_type new_dimension =
			(size_type) pow( 2, ceil( log2( max_dimension ) ) );

	// print new dimension and padding amount for testing
	// size_type row_padding = new_dimension - matrix.sub_row_size();
	// size_type col_padding = new_dimension - matrix.sub_col_size();
	// fprintf( stderr, "%u %u %u\n", new_dimension, row_padding, col_padding );

	// allocate square GSL matrix having dimension that is a power of 2
	gsl_matrix* m = gsl_matrix_alloc( new_dimension, new_dimension );
	if( m == NULL )
	{
		err_quit( "Unable to allocate GSL matrix\n" );
	}

	// copy my matrix into GSL matrix: row index needed outside of this loop
	size_t m_i = 0; 
	for( size_type i = matrix.row_begin(); i != matrix.row_end(); ++i, ++m_i )
	{
		// column index needed across 3 loops: m_j in [0, new_dimension)
		size_t m_j = 0;

		// copy each element from my matrix into GSL matrix
		for( size_type j = matrix.col_begin(); j != matrix.col_end(); ++j, ++m_j )
		{
			gsl_matrix_set( m, m_i, m_j, double( matrix(i, j) ) );
		}

		// pad columns at end of row: fill in by copying my matrix in reverse
		// order to reduce the magnitude of the high-pass coefficients
		for( size_type j = matrix.col_end();
				j != 0 && m_j != new_dimension;
				--j, ++m_j )
		{
			gsl_matrix_set( m, m_i, m_j, double( matrix(i, j - 1) ) );
		}

		// pad with 0s if we have exhausted the number of columns in my matrix in
		// the reverse iteration loop above
		for( ; m_j != new_dimension; ++m_j )
		{
			gsl_matrix_set( m, m_i, m_j, 0 );
		}
	}

	// pad rows at end of row: fill in by copying my matrix in reverse order
	// to reduce the magnitude of the high-pass coefficients
	for( size_t i = m_i; i != 0 && m_i != new_dimension; --i, ++m_i )
	{
		for( size_t m_j = 0; m_j != new_dimension; ++m_j )
		{
			// copy from GSL matrix since its width is correct
			gsl_matrix_set( m, m_i, m_j, gsl_matrix_get( m, i - 1, m_j ) );
		}
	}

	// pad with 0s if we have exhausted the number of rows in the loop
	// above that iterates in reverse
	for( ; m_i != new_dimension; ++m_i )
	{
		for( size_t m_j = 0; m_j != new_dimension; ++m_j )
		{
			// copy from GSL matrix since its width is correct
			gsl_matrix_set( m, m_i, m_j, 0 );
		}
	}

	/*
	// print resultant GSL matrix for testing
	for( size_t i = 0; i != m->size1; ++i )
	{
		for( size_t j = 0; j != m->size2; ++j )
		{
			fprintf( stderr, "%.1lf ", gsl_matrix_get( m, i, j ) );
		}
		fprintf( stderr, "\n" );
	}
	fprintf( stderr, "\n" );
	 */

	return( m );
}

/**
	Copy GSL matrix into my matrix. Any extra padding that was added to the GSL
	matrix to make the dimension a power of 2 is not copied.
	The GSL matrix is deallocated after this routine.
	This implies that the transform is not perfectly reconstructible as we are
	losing coefficients along the edges.

	@param[in,out] m GSL matrix to copy from
	@param[out] matrix My matrix to copy to
	@post GSL matrix m is deallocated.
 */
void
Wavelet::gsl_to_matrix( gsl_matrix* m, Matrix& matrix ) const
{
	size_t m_i = 0, m_j = 0;  // indices for GSL matrix

	for( size_type i = matrix.row_begin(); i != matrix.row_end(); ++i, ++m_i )
	{
		m_j = 0;
		for( size_type j = matrix.col_begin(); j != matrix.col_end(); ++j, ++m_j )
		{
			matrix(i, j) = value_type( gsl_matrix_get( m, m_i, m_j ) );
		}
	}
	gsl_matrix_free( m );
}

/**
	Get the GSL wavelet type and size.
	@param[in] wavelet_type Type of wavelet transform to perform
	@retval wavelet Pointer to GSL wavelet object
 */
gsl_wavelet*
Wavelet::get_wavelet( const Wavelet_Type wavelet_type ) const
{
	const gsl_wavelet_type* gsl_type;
	size_t gsl_size;

	switch( wavelet_type )
	{
		case Haar:
			gsl_type = gsl_wavelet_haar;
			gsl_size = 2;
			break;

		case Daub_4:
			gsl_type = gsl_wavelet_daubechies;
			gsl_size = 4;
			break;

		case Daub_6:
			gsl_type = gsl_wavelet_daubechies;
			gsl_size = 6;
			break;

		case Daub_8:
			gsl_type = gsl_wavelet_daubechies;
			gsl_size = 8;
			break;

		case Daub_10:
			gsl_type = gsl_wavelet_daubechies;
			gsl_size = 10;
			break;

		case Daub_12:
			gsl_type = gsl_wavelet_daubechies;
			gsl_size = 12;
			break;

		case Daub_14:
			gsl_type = gsl_wavelet_daubechies;
			gsl_size = 14;
			break;

		case Daub_16:
			gsl_type = gsl_wavelet_daubechies;
			gsl_size = 16;
			break;

		case Daub_18:
			gsl_type = gsl_wavelet_daubechies;
			gsl_size = 18;
			break;

		case Daub_20:
			gsl_type = gsl_wavelet_daubechies;
			gsl_size = 20;
			break;

		case Bspline_103:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 103;
			break;

		case Bspline_105:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 105;
			break;

		case Bspline_202:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 202;
			break;

		case Bspline_204:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 204;
			break;

		case Bspline_206:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 206;
			break;

		case Bspline_208:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 208;
			break;

		case Bspline_301:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 301;
			break;

		case Bspline_303:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 303;
			break;

		case Bspline_305:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 305;
			break;

		case Bspline_307:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 307;
			break;

		case Bspline_309:
			gsl_type = gsl_wavelet_bspline;
			gsl_size = 309;
			break;

		default:
			gsl_type = gsl_wavelet_haar;
			gsl_size = 2;
			break;
	}

	gsl_wavelet* wavelet = gsl_wavelet_alloc( gsl_type, gsl_size );
	if( wavelet == NULL )
	{
		err_quit( "Unable to allocate wavelet\n" );
	}

	return( wavelet );
}

/**
	Set the wavelet coefficients and size.
	@param[in] wavelet_type Type of wavelet transform to perform
	@param[in] is_centered Whether wavelet transform to perform
 */
void
Wavelet::set_wavelet( const Wavelet_Type wavelet_type, bool is_centered )
{
	switch( wavelet_type )
	{
		case Haar:
			_nc = 2;
			_h1 = Haar_Coefficients::ch_2;
			_g1 = Haar_Coefficients::cg_2;
			_h2 = Haar_Coefficients::ch_2;
			_g2 = Haar_Coefficients::cg_2;
			break;

		case Daub_4:
			_nc = 4;
			_h1 = Daub_Coefficients::h_4;
			_g1 = Daub_Coefficients::g_4;
			_h2 = Daub_Coefficients::h_4;
			_g2 = Daub_Coefficients::g_4;
			break;

		case Daub_6:
			_nc = 6;
			_h1 = Daub_Coefficients::h_6;
			_g1 = Daub_Coefficients::g_6;
			_h2 = Daub_Coefficients::h_6;
			_g2 = Daub_Coefficients::g_6;
			break;

		case Daub_8:
			_nc = 8;
			_h1 = Daub_Coefficients::h_8;
			_g1 = Daub_Coefficients::g_8;
			_h2 = Daub_Coefficients::h_8;
			_g2 = Daub_Coefficients::g_8;
			break;

		case Daub_10:
			_nc = 10;
			_h1 = Daub_Coefficients::h_10;
			_g1 = Daub_Coefficients::g_10;
			_h2 = Daub_Coefficients::h_10;
			_g2 = Daub_Coefficients::g_10;
			break;

		case Daub_12:
			_nc = 12;
			_h1 = Daub_Coefficients::h_12;
			_g1 = Daub_Coefficients::g_12;
			_h2 = Daub_Coefficients::h_12;
			_g2 = Daub_Coefficients::g_12;
			break;

		case Daub_14:
			_nc = 14;
			_h1 = Daub_Coefficients::h_14;
			_g1 = Daub_Coefficients::g_14;
			_h2 = Daub_Coefficients::h_14;
			_g2 = Daub_Coefficients::g_14;
			break;

		case Daub_16:
			_nc = 16;
			_h1 = Daub_Coefficients::h_16;
			_g1 = Daub_Coefficients::g_16;
			_h2 = Daub_Coefficients::h_16;
			_g2 = Daub_Coefficients::g_16;
			break;

		case Daub_18:
			_nc = 18;
			_h1 = Daub_Coefficients::h_18;
			_g1 = Daub_Coefficients::g_18;
			_h2 = Daub_Coefficients::h_18;
			_g2 = Daub_Coefficients::g_18;
			break;

		case Daub_20:
			_nc = 20;
			_h1 = Daub_Coefficients::h_20;
			_g1 = Daub_Coefficients::g_20;
			_h2 = Daub_Coefficients::h_20;
			_g2 = Daub_Coefficients::g_20;
			break;

		case Bspline_103:
			_nc = 6;
			_h1 = Bspline_Coefficients::h1_103;
			_g1 = &Bspline_Coefficients::g1_1[2];
			_h2 = &Bspline_Coefficients::h2_1[2];
			_g2 = Bspline_Coefficients::g2_103;
			break;

		case Bspline_105:
			_nc = 10;
			_h1 = Bspline_Coefficients::h1_105;
			_g1 = Bspline_Coefficients::g1_1;
			_h2 = Bspline_Coefficients::h2_1;
			_g2 = Bspline_Coefficients::g2_105;
			break;

		case Bspline_202:
			_nc = 6;
			_h1 = Bspline_Coefficients::h1_202;
			_g1 = &Bspline_Coefficients::g1_2[6];
			_h2 = &Bspline_Coefficients::h2_2[6];
			_g2 = Bspline_Coefficients::g2_202;
			break;

		case Bspline_204:
			_nc = 10;
			_h1 = Bspline_Coefficients::h1_204;
			_g1 = &Bspline_Coefficients::g1_2[4];
			_h2 = &Bspline_Coefficients::h2_2[4];
			_g2 = Bspline_Coefficients::g2_204;
			break;

		case Bspline_206:
			_nc = 14;
			_h1 = Bspline_Coefficients::h1_206;
			_g1 = &Bspline_Coefficients::g1_2[2];
			_h2 = &Bspline_Coefficients::h2_2[2];
			_g2 = Bspline_Coefficients::g2_206;
			break;

		case Bspline_208:
			_nc = 18;
			_h1 = Bspline_Coefficients::h1_208;
			_g1 = Bspline_Coefficients::g1_2;
			_h2 = Bspline_Coefficients::h2_2;
			_g2 = Bspline_Coefficients::g2_208;
			break;

		case Bspline_301:
			_nc = 4;
			_h1 = Bspline_Coefficients::h1_301;
			_g1 = &Bspline_Coefficients::g1_3[8];
			_h2 = &Bspline_Coefficients::h2_3[8];
			_g2 = Bspline_Coefficients::g2_301;
			break;

		case Bspline_303:
			_nc = 8;
			_h1 = Bspline_Coefficients::h1_303;
			_g1 = &Bspline_Coefficients::g1_3[6];
			_h2 = &Bspline_Coefficients::h2_3[6];
			_g2 = Bspline_Coefficients::g2_303;
			break;

		case Bspline_305:
			_nc = 12;
			_h1 = Bspline_Coefficients::h1_305;
			_g1 = &Bspline_Coefficients::g1_3[4];
			_h2 = &Bspline_Coefficients::h2_3[4];
			_g2 = Bspline_Coefficients::g2_305;
			break;

		case Bspline_307:
			_nc = 16;
			_h1 = Bspline_Coefficients::h1_307;
			_g1 = &Bspline_Coefficients::g1_3[2];
			_h2 = &Bspline_Coefficients::h2_3[2];
			_g2 = Bspline_Coefficients::g2_307;
			break;

		case Bspline_309:
			_nc = 20;
			_h1 = Bspline_Coefficients::h1_309;
			_g1 = Bspline_Coefficients::g1_3;
			_h2 = Bspline_Coefficients::h2_3;
			_g2 = Bspline_Coefficients::g2_309;
			break;

		default:
			break;
	}

	// recenter filter position
	_offset = is_centered ? (_nc >> 1) : 0 ;
}

/**
	Get wavelet subbands.
	signal is passed as non-const reference since we change its region, but the
	region and data are the same upon exit.

	@param[in,out] signal
	@param[in] num_levels
	@post signal is put in octave format as though octave_form() were called.
 */
Wavelet_Subbands
Wavelet::get_subbands( Matrix& signal, unsigned num_levels ) const
{
	// set of all wavelet subbands in signal
	Wavelet_Subbands subbands;

	// level and band from decomposition: used as temporary in multiple places
	Wavelet_Level_Band level_band;

	// We change the active region below to the subband we need, so we save
	// the current region since we need to continue partitioning the data
	// afterwards using the current region.
	Matrix_Region saved_region = signal.get_region();

	// number of rows, for example, is a function of the number of columns
	// (not the number of rows) since the loop iteration below inside the if
	// statement is executed only if the coefficients can be further split
	unsigned row_levels = get_num_levels( num_levels, signal.sub_col_size() );
	unsigned col_levels = get_num_levels( num_levels, signal.sub_row_size() );
	num_levels = max<unsigned>( row_levels, col_levels );

	// current number of rows to reorder and corresponding size of each row
	size_type row_end  = signal.row_end();
	size_type col_size = signal.sub_col_size();

	// current number of columns to reorder and corresponding size of each column
	size_type col_end  = signal.col_end();
	size_type row_size = signal.sub_row_size();

	for( unsigned level = 0; level != num_levels; ++level )
	{
		// the next octave size is half the size of the current octave, but
		// but round up if the size is odd (e.g., 6 -> 3 and 7 -> 4)
		row_end  = (row_end + 1)  / 2;
		col_size = (col_size + 1) / 2;

		// the next octave size is half the size of the current octave, but
		// but round up if the size is odd (e.g., 6 -> 3 and 7 -> 4)
		col_end  = (col_end + 1)  / 2;
		row_size = (row_size + 1) / 2;

		// row_end and col_end have been updated to the positions we need
		// with col_size and row_size giving the size of each subband

		// if both the rows and columns were reordered, then we have 3 subbands:
		// LH, HL, and HH
		if( row_levels != 0 && col_levels != 0 )
		{
			--row_levels;
			--col_levels;

			// fprintf( stderr, "Both rows and columns\n" );

			// LH subband: top-right band
			signal.set_region( signal.row_begin(), col_end, row_end, col_size );

			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::LH );

			// Note: copying signal into the subband only copies the current region
			// of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

			// fprintf( stderr, "\tLH: %s\n",
			// signal.get_region().to_string().c_str() );
			// subbands[ level_band ].write();

			// we need signal.col_begin() of the original region,
			// so revert back to it
			signal.set_region( saved_region );

			// HL subband: bottom-left band
			signal.set_region( row_end, signal.col_begin(), row_size, col_size );

			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::HL );

			// Note: copying signal into the subband only copies the current region
			// of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

			// fprintf( stderr, "\tHL: %s\n",
			// signal.get_region().to_string().c_str() );
			// subbands[ level_band ].write();

			// HH subband: bottom-right band
			signal.set_region( row_end, col_end, row_size, col_size );

			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::HH );

			// Note: copying signal into the subband only copies the current region
			// of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

			// fprintf( stderr, "\tHH: %s\n",
			// signal.get_region().to_string().c_str() );
			// subbands[ level_band ].write();
		}

		// if only the rows were reordered, we have one H subband on the right
		else if( row_levels != 0 )
		{
			--row_levels;

			// fprintf( stderr, "Rows only\n" );

			// H subband: right-half band
			// Note that row_size gives both the starting region row and the
			// length of the region since row_size
			// TODO Verify that using row_size is correct for non-powers of 2
			signal.set_region( signal.row_begin(), col_size, row_size, col_size );

			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::H );

			// Note: the assignment only copies the current region of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

		// subbands[ level_band ].write();
		// fprintf( stderr, "\tH: %s\n", signal.get_region().to_string().c_str() );
		}

		// if only the columns were reordered, we have one H subband on the bottom
		else if( col_levels != 0 )
		{
			--col_levels;

			// fprintf( stderr, "Columns only\n" );

			// H subband: bottom-half band
			// Note that row_size gives both the starting region row and the
			// length of the region since row_size
			// TODO Verify that using row_size is correct for non-powers of 2
			signal.set_region( row_size, signal.col_begin(), row_size, col_size );


			level_band = Wavelet_Level_Band( level + 1, Wavelet_Level_Band::H );

			// Note: the assignment only copies the current region of the signal
			subbands[ level_band ] = Wavelet_Subband( signal, level_band );

		// fprintf( stderr, "\tH: %s\n", signal.get_region().to_string().c_str() );
		// subbands[ level_band ].write();
		}

		// revert to original region
		signal.set_region( saved_region );

		// signal.write();
		// fprintf( stderr, "Row end: %u, Col size: %u\n", row_end, col_size );
		// fprintf( stderr, "Col end: %u, Row size: %u\n", col_end, row_size );
	}

	// LL subband: top-left band
	signal.set_region( signal.row_begin(), signal.col_begin(),
			row_size, col_size );

	level_band = Wavelet_Level_Band( num_levels, Wavelet_Level_Band::LL );

	// Note: copying signal into the subband only copies the current region
	// of the signal
	subbands[ level_band ] = Wavelet_Subband( signal, level_band );

	// fprintf( stderr, "\tLL: %s\n", signal.get_region().to_string().c_str() );
	// subbands[ level_band ].write();

	// revert to original region
	signal.set_region( saved_region );


	return( subbands );
}

/**
	Get the wavelet type from given string.
	@param[in] wavelet_str Name of wavelet 
	@retval wavelet_type Type of wavelet
 */
Wavelet_Type
Wavelet::get_wavelet_type( const string& wavelet_str )
{
	Wavelet_Type wavelet_type = Haar;

	if( wavelet_str == "Haar" )
	{
		wavelet_type = Haar;
	}
	else if( wavelet_str == "Daub_4" )
	{
		wavelet_type = Daub_4;
	}
	else if( wavelet_str == "Daub_6" )
	{
		wavelet_type = Daub_6;
	}
	else if( wavelet_str == "Daub_8" )
	{
		wavelet_type = Daub_8;
	}
	else if( wavelet_str == "Daub_10" )
	{
		wavelet_type = Daub_10;
	}
	else if( wavelet_str == "Daub_12" )
	{
		wavelet_type = Daub_12;
	}
	else if( wavelet_str == "Daub_14" )
	{
		wavelet_type = Daub_14;
	}
	else if( wavelet_str == "Daub_16" )
	{
		wavelet_type = Daub_16;
	}
	else if( wavelet_str == "Daub_18" )
	{
		wavelet_type = Daub_18;
	}
	else if( wavelet_str == "Daub_20" )
	{
		wavelet_type = Daub_20;
	}
	else if( wavelet_str == "Bspline_103" )
	{
		wavelet_type = Bspline_103;
	}
	else if( wavelet_str == "Bspline_105" )
	{
		wavelet_type = Bspline_105;
	}
	else if( wavelet_str == "Bspline_202" )
	{
		wavelet_type = Bspline_202;
	}
	else if( wavelet_str == "Bspline_204" )
	{
		wavelet_type = Bspline_204;
	}
	else if( wavelet_str == "Bspline_206" )
	{
		wavelet_type = Bspline_206;
	}
	else if( wavelet_str == "Bspline_208" )
	{
		wavelet_type = Bspline_208;
	}
	else if( wavelet_str == "Bspline_301" )
	{
		wavelet_type = Bspline_301;
	}
	else if( wavelet_str == "Bspline_303" )
	{
		wavelet_type = Bspline_303;
	}
	else if( wavelet_str == "Bspline_305" )
	{
		wavelet_type = Bspline_305;
	}
	else if( wavelet_str == "Bspline_307" )
	{
		wavelet_type = Bspline_307;
	}
	else if( wavelet_str == "Bspline_309" )
	{
		wavelet_type = Bspline_309;
	}
	else
	{
		err_quit( "Wavelet type '%s' is unsupported\n", wavelet_str.c_str() );
	}

	return( wavelet_type );
}

/**
	Get the string representation for the given wavelet type.
	@param[in] wavelet_type Type of wavelet
	@retval wavelet_str Name of wavelet 
 */
string
Wavelet::get_wavelet_type( const Wavelet_Type& wavelet_type )
{
	string wavelet_str;

	switch( wavelet_type )
	{
		case Haar:
			wavelet_str = "Haar";
			break;

		case Daub_4:
			wavelet_str = "Daub_4";
			break;

		case Daub_6:
			wavelet_str = "Daub_6";
			break;

		case Daub_8:
			wavelet_str = "Daub_8";
			break;

		case Daub_10:
			wavelet_str = "Daub_10";
			break;

		case Daub_12:
			wavelet_str = "Daub_12";
			break;

		case Daub_14:
			wavelet_str = "Daub_14";
			break;

		case Daub_16:
			wavelet_str = "Daub_16";
			break;

		case Daub_18:
			wavelet_str = "Daub_18";
			break;

		case Daub_20:
			wavelet_str = "Daub_20";
			break;

		case Bspline_103:
			wavelet_str = "Bspline_103";
			break;

		case Bspline_105:
			wavelet_str = "Bspline_105";
			break;

		case Bspline_202:
			wavelet_str = "Bspline_202";
			break;

		case Bspline_204:
			wavelet_str = "Bspline_204";
			break;

		case Bspline_206:
			wavelet_str = "Bspline_206";
			break;

		case Bspline_208:
			wavelet_str = "Bspline_208";
			break;

		case Bspline_301:
			wavelet_str = "Bspline_301";
			break;

		case Bspline_303:
			wavelet_str = "Bspline_303";
			break;

		case Bspline_305:
			wavelet_str = "Bspline_305";
			break;

		case Bspline_307:
			wavelet_str = "Bspline_307";
			break;

		case Bspline_309:
			wavelet_str = "Bspline_309";
			break;

		default:
			err_quit( "Wavelet type is unsupported\n" );
			break;
	}

	return( wavelet_str );
}

/**
	Non-standard wavelet transform.
void
Wavelet::wavelet2d_nstransform( Matrix& matrix, unsigned num_levels ) const
{
	if( matrix.row_size() != matrix.col_size() )
	{
		// err_quit( "2d dwt works only with square matrix\n" );
	}

	if( binary_logn( matrix.row_size() ) == -1 )
	{
		size_type max_dimension =
				std::max<size_type>( matrix.sub_row_size(), matrix.sub_col_size() );
		size_type new_dimension =
				(size_type) pow( 2, ceil( log2( max_dimension ) ) );
		Matrix tmp( max_dimenstion );
		for( size_type i = matrix.row_begin(); i != matrix.row_end(); ++i )
		{
			for( size_type j = matrix.col_begin(); j != matrix.col_end(); ++j )
			{
				tmp(i,j) = matrix(i,j);
			}
		}
		matrix = tmp;
			// (size_type) pow( 2, ceil( log2( max_dimension ) ) );
			
		// err_quit( "n is not a power of 2\n" );
	}

	if( matrix.row_size() < 2 )
	{
		return;
	}

	// matrix.set_region( 0, 0, 256 , 256 );
	// matrix = matrix;

	// int n = (size == 0) ? 0 : int( floor( log2(size) ) );
	// return( (n < 0) ? 0 : std::min<int>( n, num_levels ) );

	Matrix::size_type size = min<Matrix::size_type>( matrix.row_size(),
			matrix.col_size());
	// size = min<Matrix::size_type>( size, num_levels );

print_filters();

	for( unsigned step = size, level = 0; step >= 2 && level != num_levels;
			step >>= 1, ++level )
	// for( unsigned step = size; step >= 2; step >>= 1 )
	{

	// new dimension for matrix is the smallest power of 2 >= current dimension;
	// current dimension is taken to be the largest size between the number of
	// rows and columns since the new matrix must be square
	// size_type max_dimension =
			// std::max<size_type>( matrix.sub_row_size(), matrix.sub_col_size() );
	// size_type new_dimension =
			// (size_type) pow( 2, ceil( log2( max_dimension ) ) );
	size_type row_work_size =
			(size_type) pow( 2, ceil( log2( matrix.sub_col_size() ) ) );
	fprintf( stderr, "level: %u %u %u\n", level, num_levels, step, size );

	// fprintf( stderr, "level: %u\n", step, size );

		// for every row j
		for( Matrix::size_type row = 0; row < step; ++row )
		{
			// Vector row_work( matrix.col_size() );  // workspace for rows

	// fprintf( stderr, "you should change the size to work with 2^k\n" );
			Vector row_work( matrix.col_size() );  // workspace for rows
			// Vector row_work( row_work_size );  // workspace for rows

			unsigned nmod = _nc * step;  // positive constant = zero mod step
			// fprintf( stderr, "%u\n", nmod );
			nmod -= _offset;  // center support

			unsigned nl = step - 1;
			unsigned nh = step >> 1;

			// for( Vector::size_type ii = 0, i = 0;
					// i < (matrix.col_size() - _nc + 1); i += 2, ++ii )
			for( Vector::size_type ii = 0, i = 0; i < step; i += 2, ++ii )
			{
				unsigned ni = i + nmod;  // pointer to be incremented and wrapped

				// for each wavelet coefficient
				for( Vector::size_type k = 0; k < _nc; ++k )
				{
					// wrap pointer around
					// unsigned jf = nl & (ni + k);
					row_work(ii)      += _h1[k] * matrix(row, i);
					row_work(ii + nh) += _g1[k] * matrix(row, i);
				}
			}

			// copy row into matrix
			for( Vector::size_type j = 0; j < step; ++j )
			{
				matrix(row, j) = row_work(j);
			}
		}

		// for every column j
		for( Matrix::size_type col = 0; col < step; ++col )
		{
			Vector col_work( matrix.row_size() );  // workspace for columns

			unsigned nmod = _nc * step;
			nmod -= _offset;  // center support

			unsigned nl = step - 1;
			unsigned nh = step >> 1;

			for( Vector::size_type ii = 0, i = 0;
					i < (matrix.row_size() - _nc + 1); i += 2, ++ii )
			// for( Vector::size_type ii = 0, i = 0; i < step; i += 2, ++ii )
			{
				unsigned ni = i + nmod;

				// for each wavelet coefficient
				for( Vector::size_type k = 0; k < _nc; ++k )
				{
					// unsigned jf = nl & (ni + k);
					unsigned jf = i;
					col_work(ii)      += _h1[k] * matrix(jf, col);
					col_work(ii + nh) += _g1[k] * matrix(jf, col);
				}
			}

			for( Vector::size_type i = 0; i < step; ++i )
			{
				matrix(i, col) = col_work(i);
			}
		}
	}
}
 */

/**
	Non-standard wavelet transform.
	Pads to power of 2 just to get it working.
 */
void
Wavelet::wavelet2d_nstransform( Matrix& matrix, unsigned num_levels ) const
{
	Matrix* data = &matrix;  // data to use for transform

	if( matrix.row_size() != matrix.col_size() )
	{
		// err_quit( "2d dwt works only with square matrix\n" );
	}
	if( matrix.row_size() < 2 )
	{
		return;
	}

	size_type max_dimension =
		std::max<size_type>( matrix.sub_row_size(), matrix.sub_col_size() );
	size_type new_dimension =
		(size_type) pow( 2, ceil( log2( max_dimension ) ) );

	// copy matrix into temporary matrix with dimension that is a power of 2
	Matrix tmp;
	// if( binary_logn( matrix.row_size() ) == -1 )
	if( new_dimension != matrix.row_size()
			|| new_dimension != matrix.col_size() )
	{
		tmp.resize( new_dimension, new_dimension, false );

		for( size_type i = matrix.row_begin(), ii = tmp.row_begin();
				i != matrix.row_end(); ++i, ++ii )
		{
			for( size_type j = matrix.col_begin(), jj = tmp.col_begin();
					j != matrix.col_end(); ++j, ++jj )
			{
				tmp(ii,jj) = matrix(i,j);
			}
		}
		data = &tmp;
	}

	Matrix::size_type size = min<Matrix::size_type>( data->row_size(),
			data->col_size());

	fprintf( stderr, "%s\n", data->get_region().to_string().c_str() );
// num_levels = 65000;
	for( unsigned step = size, level = 0; step >= 2 && level != num_levels;
			step >>= 1, ++level )
	{
		// for every row j
		for( Matrix::size_type row = 0; row < step; ++row )
		{
			Vector row_work( data->col_size() );  // workspace for rows

			unsigned nmod = _nc * step;  // positive constant = zero mod step
			nmod -= _offset;  // center support

			unsigned nl = step - 1;
			unsigned nh = step >> 1;

			// for( Vector::size_type ii = 0, i = 0;
					// i < (matrix.col_size() - _nc + 1); i += 2, ++ii )
			for( Vector::size_type ii = 0, i = 0; i < step; i += 2, ++ii )
			{
				unsigned ni = i + nmod;  // pointer to be incremented and wrapped

				// for each wavelet coefficient
				for( Vector::size_type k = 0; k < _nc; ++k )
				{
					// wrap pointer around
					unsigned jf = nl & (ni + k);
					row_work(ii)      += _h1[k] * data->at(row, jf);
					row_work(ii + nh) += _g1[k] * data->at(row, jf);
				}
			}

		// TODO do not copy up to step; only copy real data, not any padding
			// copy row into matrix
			for( Vector::size_type j = 0; j < step; ++j )
			{
				data->at(row, j) = row_work(j);
			}
		}

		// for every column j
		for( Matrix::size_type col = 0; col < step; ++col )
		{
			Vector col_work( data->row_size() );  // workspace for columns

			unsigned nmod = _nc * step;
			nmod -= _offset;  // center support

			unsigned nl = step - 1;
			unsigned nh = step >> 1;

			// for( Vector::size_type ii = 0, i = 0;
					// i < (matrix.row_size() - _nc + 1); i += 2, ++ii )
			for( Vector::size_type ii = 0, i = 0; i < step; i += 2, ++ii )
			{
				unsigned ni = i + nmod;

				// for each wavelet coefficient
				for( Vector::size_type k = 0; k < _nc; ++k )
				{
					unsigned jf = nl & (ni + k);
					// unsigned jf = i;
					col_work(ii)      += _h1[k] * data->at(jf, col);
					col_work(ii + nh) += _g1[k] * data->at(jf, col);
				}
			}

			for( Vector::size_type i = 0; i < step; ++i )
			{
				data->at(i, col) = col_work(i);
			}
		}
	}

	if( data != &matrix )
	{
		data->set_region( matrix.get_region() );
		matrix = *data;
		// for( size_type i = matrix.row_begin(), ii = data->row_begin();
				// i != matrix.row_end(); ++i, ++ii )
		// {
			// for( size_type j = matrix.col_begin(), jj = data->col_begin();
					// j != matrix.col_end(); ++j, ++jj )
			// {
				// matrix(i,j) = data->at(ii,jj);
			// }
		// }
	}

	fprintf( stderr, "%s\n", matrix.get_region().to_string().c_str() );
}

/**
	Non-standard wavelet transform.
void
Wavelet::wavelet2d_nstransform( Matrix& matrix, unsigned num_levels ) const
{
	if( matrix.row_size() != matrix.col_size() )
	{
		// err_quit( "2d dwt works only with square matrix\n" );
	}

	if( binary_logn( matrix.row_size() ) == -1 )
	{
		// err_quit( "n is not a power of 2\n" );
	}

	if( matrix.row_size() < 2 )
	{
		return;
	}

	// matrix.set_region( 0, 0, 256 , 256 );
	// matrix = matrix;

	// int n = (size == 0) ? 0 : int( floor( log2(size) ) );
	// return( (n < 0) ? 0 : std::min<int>( n, num_levels ) );

	Matrix::size_type size = min<Matrix::size_type>( matrix.row_size(),
			matrix.col_size());
	// size = min<Matrix::size_type>( size, num_levels );

print_filters();

	for( unsigned step = size, level = 0; step >= 2 && level != num_levels;
			step >>= 1, ++level )
	// for( unsigned step = size; step >= 2; step >>= 1 )
	{
	fprintf( stderr, "level: %u %u %u\n", level, num_levels, step, size );
	// fprintf( stderr, "level: %u\n", step, size );

		// for every row j
		for( Matrix::size_type row = 0; row < step; ++row )
		{
			Vector row_work( matrix.col_size() );  // workspace for rows

			unsigned nmod = _nc * step;  // positive constant = zero mod step
			// fprintf( stderr, "%u\n", nmod );
			nmod -= _offset;  // center support

			unsigned nl = step - 1;
			unsigned nh = step >> 1;

			for( Vector::size_type ii = 0, i = 0; i < step; i += 2, ++ii )
			{
				unsigned ni = i + nmod;  // pointer to be incremented and wrapped

				// for each wavelet coefficient
				for( Vector::size_type k = 0; k < _nc; ++k )
				{
					// wrap pointer around
					unsigned jf = nl & (ni + k);
					row_work(ii)      += _h1[k] * matrix(row, jf);
					row_work(ii + nh) += _g1[k] * matrix(row, jf);
				}
			}

			// copy row into matrix
			for( Vector::size_type j = 0; j < step; ++j )
			{
				matrix(row, j) = row_work(j);
			}
		}

		// for every column j
		for( Matrix::size_type col = 0; col < step; ++col )
		{
			Vector col_work( matrix.row_size() );  // workspace for columns

			unsigned nmod = _nc * step;
			nmod -= _offset;  // center support

			unsigned nl = step - 1;
			unsigned nh = step >> 1;

			for( Vector::size_type ii = 0, i = 0; i < step; i += 2, ++ii )
			{
				unsigned ni = i + nmod;

				// for each wavelet coefficient
				for( Vector::size_type k = 0; k < _nc; ++k )
				{
					unsigned jf = nl & (ni + k);
					col_work(ii) += _h1[k] * matrix(jf, col);
					col_work(ii + nh) += _g1[k] * matrix(jf, col);
				}
			}

			for( Vector::size_type i = 0; i < step; ++i )
			{
				matrix(i, col) = col_work(i);
			}
		}
	}
}
 */

/*
#define ELEMENT(a,stride,i) ((a)[(stride)*(i)])

void
Wavelet::wavelet2d_nstransform( const my_wavelet_type* w,
	double* data, size_t tda, size_t size1, size_t size2, Vector& work )
{
	size_t i, j;

	if (size1 != size2)
	{
		// GSL_ERROR ("2d dwt works only with square matrix", GSL_EINVAL);
	}

	if( work.size() < size1 )
	{
		// GSL_ERROR ("not enough workspace provided", GSL_EINVAL);
	}

	if (binary_logn (size1) == -1)
	{
		// GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
	}

	if (size1 < 2)
	{
		// return GSL_SUCCESS;
	}

	for( i = size1; i >= 2; i >>= 1 )
	{
		// for every row j
		for (j = 0; j < i; j++)
		{
			dwt_step( w, &ELEMENT(data, tda, j), 1, i, work );
		}

		// for every column j
		for (j = 0; j < i; j++)
		{
			dwt_step( w, &ELEMENT(data, 1, j), tda, i, work );
		}
	}
}
 */

/*
void
Wavelet::dwt_step( const my_wavelet_type* w, double* a, size_t stride, size_t n,
	Vector& work )
{
	size_t i, ii;
	size_t jf;
	size_t k;
	size_t nl, ni, nh, nmod;

	for( i = 0; i < work.size(); i++ )
	{
		work(i) = 0.0;
	}

	nmod  = w->nc * n;
	nmod -= w->offset;  // center support

	nl = n - 1;
	nh = n >> 1;

	for (ii = 0, i = 0; i < n; i += 2, ii++)
	{
		ni = i + nmod;
		for (k = 0; k < w->nc; k++)
		{
			jf = nl & (ni + k);
			work(ii) += w->h1[k] * ELEMENT (a, stride, jf);
			work(ii + nh) += w->g1[k] * ELEMENT (a, stride, jf);
		}
	}

	for (i = 0; i < n; i++)
	{
		ELEMENT( a, stride, i ) = work(i);
	}
}
 */

/**
	Determine if n is power of 2.
 */
int
Wavelet::binary_logn( const unsigned n ) const
{
	unsigned ntest;
	unsigned logn = 0;
	unsigned k = 1;

	while( k < n )
	{
		// k <<= 1;
		k *= 2;
		++logn;
	}

	ntest = (1 << logn);

	if( n != ntest )
	{
		return( -1 );  // n is not a power of 2
	}

	return( logn );
}

/*
void
Wavelet::inverse_dwt_step( const my_wavelet_type* w, double* a, size_t stride,
		size_t n, Vector& work )
{
	double ai, ai1;
	size_t i, ii;
	size_t jf;
	size_t k;
	size_t nl, ni, nh, nmod;

	for (i = 0; i < work.size(); i++)
	{
		work(i) = 0.0;
	}

	nmod  = w->nc * n;
	nmod -= w->offset;  // center support

	nl = n - 1;
	nh = n >> 1;

	for (ii = 0, i = 0; i < n; i += 2, ii++)
	{
		ai = ELEMENT (a, stride, ii);
		ai1 = ELEMENT (a, stride, ii + nh);
		ni = i + nmod;
		for (k = 0; k < w->nc; k++)
		{
			jf = (nl & (ni + k));
			work(jf) += (w->h2[k] * ai + w->g2[k] * ai1);
		}
	}

	for (i = 0; i < n; i++)
	{
		ELEMENT (a, stride, i) = work(i);
	}
}
 */

/**
	Print filters.
 */
void
Wavelet::print_filters( ) const
{
	if( _nc == 0 )
	{
		return;
	}

	// print low-pass
	fprintf( stderr, "%lf", _h1[0] );
	for( Vector::size_type i = 1; i < _nc; ++i )
	{
		fprintf( stderr, " %lf", _h1[i] );
	}
	fprintf( stderr, "\n" );

	// print high-pass
	fprintf( stderr, "%lf", _g1[0] );
	for( Vector::size_type i = 1; i < _nc; ++i )
	{
		fprintf( stderr, " %lf", _g1[i] );
	}
	fprintf( stderr, "\n" );

	// print inverse low-pass
	fprintf( stderr, "%lf", _h2[0] );
	for( Vector::size_type i = 1; i < _nc; ++i )
	{
		fprintf( stderr, " %lf", _h2[i] );
	}
	fprintf( stderr, "\n" );

	// print inverse high-pass
	fprintf( stderr, "%lf", _g2[0] );
	for( Vector::size_type i = 1; i < _nc; ++i )
	{
		fprintf( stderr, " %lf", _g2[i] );
	}
	fprintf( stderr, "\n" );
}
