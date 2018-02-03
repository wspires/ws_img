/**
	@file   K_Means.cpp
	@author Wade Spires
	@date   2006/10/20
	@brief  Class K_Means.

	Copyright (c) 2006 Wade Spires. All rights reserved.
 */

#include "K_Means.hpp"

using std::string;
using std::vector;

using namespace ws_tools;
using namespace ws_img;
typedef Matrix::size_type  size_type;
typedef Matrix::value_type value_type;

/**
	K-means clustering.

	Iteratively shift k different mean.
	-Start with random set
	-Classify each vector into K sets
	-Compute next set of K means using the means of each class
 */
Matrix
K_Means::k_means( const unsigned K, const Matrix& data )
{
	// allocate K mean vectors along each row
	Matrix means( K, data.sub_col_size() );

	// labels assigned to each data vector
	vector<size_type> labels( data.sub_row_size(), 0 );

	select_initial_means( data, means );

	const unsigned Num_Iterations = 1000;
	// const double min_error = Matrix::MAX_VALUE;
	const double tol   = .00000001;

	double diff_avg = 0;
	double diff_var = 0;

	unsigned count = 0;
	for( unsigned i = 0; i != Num_Iterations; ++i )
	{
		++count;
		classify_samples( data, means, labels );

		Matrix previous_means = means;

		select_new_means( data, labels, means );

		double difference = compute_difference( previous_means, means );

		stat( difference, i, &diff_avg, &diff_var );

// matrix_to_image( means, 24, 24, format_string( "out/%u", i ) );

		if( difference < tol )
		{
			break;
		}
	}
fprintf( stderr, "%u iterations\n", count );

matrix_to_image( means, 24, 24, format_string( "out/final" ) );

	return( means );
}

/**
	Select the initial set of means.
	TODO just return a vector of indices instead of copying into a matrix of
	means.  
 */
void
K_Means::select_initial_means( const Matrix& data, Matrix& means )
{
	const unsigned K = means.sub_row_size();

	// data vectors that are initially assigned as means
	vector<size_type> initial_means;

	// select the initial K means at random from the data set
	Uniform_Number u( 0, data.sub_row_size() );
	size_type mi = means.row_begin();
	for( unsigned k = 0; k != K; ++k, ++mi )
	{
		unsigned i = 0;  // randomly selected data vector

		bool is_duplicate = false;
		do // find an unused data vector
		{
			is_duplicate = false;

			// choose uniformly dist. number in the range [0, data.sub_row_size() )
			i = static_cast<size_type>( floor( u() ) );

			// skip vector i if we have already used it as a mean
			for( unsigned i = 0; i != initial_means.size(); ++i )
			{
				if( i == initial_means[i] )
				{
					is_duplicate = true;
					break;

				// TODO actually compare the vectors themselves, but do not loop
				// forever since we are selecting vectors at random (maybe keep
				// track of indices to know when we have tried them all)
				}
			}
		}
		while( is_duplicate );

		// keep track of which data item we selected
		initial_means.push_back(i);

		// copy data vector into mean vector
		for( size_type j = data.col_begin(), mj = means.col_begin();
				j != data.col_end();
				++j, ++mj )
		{
			means(mi,mj) = data(i,j);
		}
	}
}

void
K_Means::classify_samples( const Matrix& data, const Matrix& means,
		vector<size_type>& labels )
{
	// compare each data item i to each mean vector k
	unsigned li = 0;
	for( size_type i = data.row_begin(); i != data.row_end(); ++i, ++li )
	{
		// the index to the mean yielding the minimum distance
		size_type  min_idx  = 0;
		value_type min_dist = Matrix::MAX_VALUE; 

		for( size_type k = means.row_begin(); k != means.row_end(); ++k )
		{
			value_type dist = 0; 

			// calculate distance
			size_type mj = means.col_begin();
			for( size_type j = data.col_begin(); j != data.col_end(); ++j, ++mj )
			{
				// add the difference squared to the distance
				// const value_type diff = data(i,j) - means(k,mj);
				// dist += diff * diff;

				// Or can we just add the absolute value of the difference since
				// we don't care about the actual value?
				dist += fabs( data(i,j) - means(k,mj) );
			}

			// check for a new minimum distance
			if( dist < min_dist )
			{
				min_dist = dist;
				min_idx  = k;
			}
		}

		// save the index to the mean k that yields the min. distance to data i
		labels[li] = min_idx;
	}

	// TODO maybe check for 0 counts here before recomputing the means and maybe
	// override the empty class by splitting a non-empty one
}

/*
	Whether to compute the mean using the slower, traditional method that simply
	sums each item and divides by the total; or whether a running average should
	be used that computes the mean in one pass and is also numerically stable.
 */
#define USE_SLOW_MEAN

/**
	Select new mean vectors.
 */
void
K_Means::select_new_means( const Matrix& data, const vector<size_type>& labels,
		Matrix& means )
{
	means = 0; // reset means

	// counts for each class
	vector<unsigned> counts( means.sub_row_size(), 0 );

	// for each data item i with class label mi
	unsigned li = 0;
	for( size_type i = data.row_begin(); i != data.row_end(); ++i, ++li )
	{
		const size_type mi = labels[li];

		// current number of times a data point has received the label mi
		++counts[mi];

		// accumulate the new mean vector mi
		for( size_type j = data.col_begin(), mj = means.col_begin();
				j != data.col_end();
				++j, ++mj )
		{
#ifdef USE_SLOW_MEAN
			means( mi, mj ) += data(i,j);
#else
			const value_type delta = data(i,j) - means(mi,mj);
			means(mi,mj) += delta / counts[mi];
#endif // USE_SLOW_MEAN
		}
	}

	// check for mean vectors that had no data items in its class
	for( unsigned i = 0; i != counts.size(); ++i )
	{
		if( counts[i] == 0 )
		{
			// TODO check for 0 samples
			// fprintf( stderr, "%u zero\n", i );
		}
	}

#ifdef USE_SLOW_MEAN
	// finish computing each mean mi
	size_type mi = means.row_begin();
	for( unsigned i = 0; i != counts.size(); ++i, ++mi )
	{
		// divide each mean vector mi by the number of data items (count) that
		// mapped to this mean
		const value_type count = counts[i];
		for( size_type mj = means.col_begin(); mj != means.col_end(); ++mj )
		{
			means( mi, mj ) /= count;
		}
	}
#endif // USE_SLOW_MEAN
}

double
K_Means::compute_difference( const Matrix& previous_means, const Matrix& means )
{
	// compute difference in the current set of means to the previous means
	value_type max_dist = 0;
	for( size_type i = means.row_begin(); i != means.row_end(); ++i )
	{
		// sum of squared differences (dot product)
		value_type dist = 0;
		for( size_type j = means.col_begin(); j != means.col_end(); ++j )
		{
			const value_type diff = previous_means(i,j) - means(i,j);
			dist += diff * diff;
		}

		if( dist > max_dist )
		{
			max_dist = dist;
		}
	}
	const double difference = sqrt( max_dist );

	return( difference );
}

void
K_Means::matrix_to_image( const Matrix& M, const size_type rows,
		const size_type cols, const string& out )
{
	string out_dir = out;
	check_dir( out_dir );

	// convert each row to an image
	for( size_type mi = M.row_begin(); mi != M.row_end(); ++mi )
	{
		Image img( rows, cols );

		size_type mj = M.col_begin();
		for( size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); ++j, ++mj )
			{
				img(i,j) = M(mi,mj);
			}
		}

		// write to file
		string out_name = out_dir + format_string( "matrix_%u.pgm", mi );
		img.write( out_name );
	}
}
