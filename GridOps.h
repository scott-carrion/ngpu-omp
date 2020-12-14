// Functions and structures for digital image processing and terrain
// modelling with Grid objects.

#ifndef YOUNG_GRIDOPS
#define YOUNG_GRIDOPS

#include <cmath>
#include <vector>

#include <omp.h>

#include "Grid.h"
#include "GridCell.h"


// Bilinear interpolation.
// -- Arguments --
// i : real-number row of position to interpolate at.
// j : real-number column of position to interpolate at.
// v : a vector of grid cells to interpolate between, ordered
//       {upper left, upper right, bottom left, bottom right}.
// -- Returns --
// The interpolated value (z) at the given point (i,j).
template <typename T>
double bilinearInterpolate ( double i, double j,
    const std::vector<GridCell<T> > & v )
{
    // FIXME: THIS FUNCTION CAUSES MULTITHREADED CRASH
    //std::cout << "thread " << omp_get_thread_num() << " starting bilinear interpolate" << std::endl;
    double xmin = v[0].j + 0.5;
    double xrng = v[3].j - v[0].j;
    double ymin = v[3].i + 0.5;
    double yrng = v[0].i - v[3].i;
    double x = (j - xmin) / xrng;
    double y = (i - ymin) / yrng;
    
    double a0 = 1.0 - x; double a1 = x;
    
    double B00 = v[0].value; double B01 = v[1].value;
    double B10 = v[2].value; double B11 = v[3].value;
    
    double c0 = 1.0 - y;
    double c1 = y;
    
    //std::cout << "thread " << omp_get_thread_num() << " ending bilinear interpolate" << std::endl;
    
    return B00 * a0 * c0
        + B10 * a1 * c0
        + B01 * a0 * c1
        + B11 * a1 * c1;
}


// Get nearest four grid cells in grid (for bilinear interpolation).
// -- Arguments --
// g  : grid to sample.
// i0 : real-number row of position to sample around.
// j0 : real-number column of position to sample around.
// -- Returns --
// A vector of 4 points with the grid values at those locations.
template <typename T>
std::vector<GridCell<T> > getNearest4 ( const Grid<T> & g, double i0,
    double j0 )
{
    // get upper-left cell
    int i = floor(i0 - 0.5);
    int j = floor(j0 - 0.5);
    if ( i < 0 ) i = 0;
    if ( j < 0 ) j = 0;
    if ( i + 1 >= g.nrows() ) i = g.nrows() - 2;
    if ( j + 1 >= g.ncols() ) j = g.ncols() - 2;
    
    // add points to vector
    std::vector<GridCell<T> > v;
    v.resize(4);
    v[0] = GridCell<T>(i,   j,   g.getValue(i,   j));
    v[1] = GridCell<T>(i,   j+1, g.getValue(i,   j+1));
    v[2] = GridCell<T>(i+1, j,   g.getValue(i+1, j));
    v[3] = GridCell<T>(i+1, j+1, g.getValue(i+1, j+1));
    
    return v;
}


// Computes positive openness, or the skyview factor.
// -- Source --
// Michael Bishop (?/?/2016?). Adapted by Brennan Young (2/16/2018).
// Not to be redistributed.
// -- Arguments --
// dem : digital elevation model grid.
// daz : azimuth angle step size in search for horizon [deg].
// r   : distance to search for maximum angle with topography.
template <typename T>
Grid<T> skyview ( const Grid<T> & dem, int daz, double r )
{
    double timer_start = omp_get_wtime();

    static const double pi = 4.0 * atan(1.0);
    static const double deg2rad = pi / 180.0;
    const double interval = static_cast<double>(daz) / 360.0;
    const int dmax = floor(r / dem.dx());
    
    bool flag;
    int i, j, k, kk, a, d;
    double ii, jj, z, cosAz, sinAz, angle, Amax, sum;
    std::vector<GridCell<T> > cells;
    Grid<T> out(dem, 1);
    
    int count = 0;

    // (SIGSEGV) #pragma omp parallel for private(count,flag,angle,k,kk,d,sum,a,sinAz,cosAz,Amax,ii,jj,z) shared(out) collapse(2)
    // (SIGSEGV) #pragma omp parallel for collapse(2)
    //#pragma omp parallel for shared(out, interval, dmax, deg2rad) collapse(2)
    //#pragma omp target teams distribute parallel for private(flag, i, j, k, kk, a, d, sinAz, cosAz, Amax, sum, cells) shared(count, out, interval, dmax, deg2rad) collapse(2)
    #pragma omp parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, cells, count, angle, z) shared(interval, dmax, deg2rad, out) collapse(2)
    for ( i = 0; i < dem.nrows(); ++i ) {
        for ( j = 0; j < dem.ncols(); ++j ) {
	    //std::cout << "skyview num threads: " << omp_get_num_threads() << std::endl;
	    

            std::cout << "\rskyview " << (static_cast<double>(count++) / dem.size()) << "%        ";
            
            k = dem.getIndex(i, j);
            
            // scan in each direction
            sum = 0.0;
            for ( a = 0; a < 360; a += daz) {
                sinAz = sin(a * deg2rad);
                cosAz = cos(a * deg2rad);
                
                // find maximum angle from horizontal
                Amax = 0.0;
                
		for ( d = 1; d <= dmax; ++d ) {
                    // get coordinates at point
                    ii = (i+0.5) - d * cosAz;
                    jj = (j+0.5) + d * sinAz;
                    

                    if ( ii < 0 || ii >= dem.nrows()
                            || jj < 0 || jj >= dem.ncols() )
                        break;
                    // interpolate value at point
                    flag = false;
                    cells = getNearest4(dem, ii, jj);
                    
		for ( kk = 0; !flag && kk < cells.size(); ++kk ){
                        if ( cells[kk].value == dem.noData() ) {
                            flag = true;
			}
                    }
                    if ( flag ) {
                        angle = 0.0;
                        continue;
                    }

                    z = bilinearInterpolate(ii, jj, cells);

                    // track maximum angle from horizontal
                    angle = atan((z - dem[k]) / (d * dem.dx()));
                    if ( angle > Amax ) Amax = angle;
                }
                
                Amax = cos(Amax);
                sum += Amax * Amax; 
            }
            out[k] = sum * interval;
        }
    }
    std::cout << "\rskyview " << (static_cast<double>(count++) / dem.size()) << "%\n";
 
    double timer_end = omp_get_wtime();
    std::cout << "Skyview took this many seconds: " << timer_end - timer_start << std::endl;
    return out;
}


// Computes prominence, which is the average downword enclosure.
// -- Source --
// Brennan Young (2/16/2018).
// Not to be redistributed.
// -- Arguments --
// dem : digital elevation model grid.
// daz : azimuth angle step size in search for horizon [deg].
// r   : distance to search for maximum angle with topography.
template <typename T>
Grid<T> prominence ( const Grid<T> & dem, int daz, double r )
{
    double timer_start = omp_get_wtime();

    static const double pi = 4.0 * atan(1.0);
    static const double deg2rad = pi / 180.0;
    const double interval = static_cast<double>(daz) / 360.0;
    const int dmax = floor(r / dem.dx());
    
    bool flag;
    int i, j, k, kk, a, d;
    double ii, jj, z, cosAz, sinAz, angle, Amax, sum;
    std::vector<GridCell<T> > cells;
    Grid<T> out(dem, 1);
    
    int count = 0;
   
   //#pragma omp parallel for shared(out, interval, dmax, deg2rad) collapse(2) 
    //#pragma omp target teams distribute parallel for shared(out, interval, dmax, deg2rad) collapse(2)
    //#pragma omp target teams distribute parallel for private(flag, i, j, k, kk, a, d, sinAz, cosAz, Amax, sum, cells) shared(count, out, interval, dmax, deg2rad) collapse(2)
    //#pragma omp parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, cells, count, angle, z, out) shared(interval, dmax, deg2rad) collapse(2)
    #pragma omp parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, cells, count, angle, z) shared(interval, dmax, deg2rad, out) collapse(2)
    for ( i = 0; i < dem.nrows(); ++i ) {
        for ( j = 0; j < dem.ncols(); ++j ) {
	    std::cout << "\rprominence " << (static_cast<double>(count++) / dem.size()) << "%        ";
            
            k = dem.getIndex(i, j);
            
            // scan in each direction
            sum = 0.0;
            
	    for ( a = 0; a < 360; a += daz) {
                sinAz = sin(a * deg2rad);
                cosAz = cos(a * deg2rad);
                
                // find maximum angle from horizontal
                Amax = 0.0;
                
		for ( d = 1; d <= dmax; ++d ) {
                    // get coordinates at point
                    ii = (i+0.5) - d * cosAz;
                    jj = (j+0.5) + d * sinAz;
                    
                    if ( ii < 0 || ii >= dem.nrows()
                            || jj < 0 || jj >= dem.ncols() )
                        break;
                    
		    // interpolate value at point
                    flag = false;
                    cells = getNearest4(dem, ii, jj);
                    
		    for ( kk = 0; !flag && kk < cells.size(); ++kk ){
                       if ( cells[kk].value == dem.noData() )
				flag = true;
		    }
                    if ( flag ) {
                        angle = 0.0;
                        continue;
                    }
                    z = bilinearInterpolate(ii, jj, cells);
                    
                    // track maximum angle from horizontal
                    angle = atan((z - dem[k]) / (d * dem.dx()));
                    if ( angle < Amax ) Amax = angle;
                }
                
                Amax = sin(Amax);
                sum += Amax * Amax; 
            }
            out[k] = sum * interval;
        }
    }
    std::cout << "\rprominence " << (static_cast<double>(count++) / dem.size()) << "%\n";
   
    double timer_end = omp_get_wtime();
    std::cout << "Prominence took this many seconds: " << timer_end - timer_start << std::endl;
 
    return out;
}

#endif // YOUNG_GRID_OPERATIONS_yyyymmdd
