// Functions and structures for digital image processing and terrain
// modelling with Grid objects.

#ifndef YOUNG_GRIDOPS
#define YOUNG_GRIDOPS

#include <cmath>
#include <vector>

#include <omp.h>

#include "Grid_tc.h"
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
    
    return B00 * a0 * c0
        + B10 * a1 * c0
        + B01 * a0 * c1
        + B11 * a1 * c1;
}

template <typename T>
double bilinearInterpolate_tc ( double i, double j,
    GridCell<T>* v, bool debug_flag )
{
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

    if (debug_flag) { printf("v[0].value == %f; v[1].value == %f; v[2].value == %f; v[3].value == %f\nv[0].i == %d; v[1].i == %d; v[2].i == %d; v[3].i == %d\nv[0].j == %d; v[1].j == %d; v[2].j == %d; v[3].j == %d\n", v[0].value, v[1].value, v[2].value, v[3].value, v[0].i, v[1].i, v[2].i, v[3].i, v[0].j, v[1].j, v[2].j, v[3].j); }
    
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
std::vector<GridCell<T> > getNearest4 ( const Grid_tc<T> & g, double i0,
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

// The below function is a version of the above that does NOT make use of C++ STL
// so that it is compatible with NVIDIA hardware that doesn't support linking against the C++ STL
template <typename T>
void getNearest4_tc(const Grid_tc<T> & g, double i0, double j0, GridCell<T>* v, int nrows, int ncols, int nbands, bool debug_flag) {

    // FIXME: IN THIS FUNCTION, g.nrows() AND g.ncols() RETURN 0 EVEN WHEN THEY SHOULDN'T. WHY IS THIS?
    //if (debug_flag) {
    //	printf("g.nrows() is %d and g.ncols() is %d\n", g.nrows(), g.ncols());
    //}


    // get upper-left cell /*
    int i = floor(i0 - 0.5);
    int j = floor(j0 - 0.5);
    if ( i < 0 ) i = 0;
    if ( j < 0 ) j = 0;
    //if ( i + 1 >= g.nrows() ) i = g.nrows() - 2;
    //if ( j + 1 >= g.ncols() ) j = g.ncols() - 2;
    if ( i + 1 >= nrows ) i = nrows - 2;
    if ( j + 1 >= ncols ) j = ncols - 2;
    
    // add points to vector
    //GridCell<T>* v = (GridCell<T>*)(std::malloc(sizeof(GridCell<T>) * 4));  // XXX instead of a std::vector, use std::malloc() and native array!
    // Instead of allocating a new v every time, pass a pointer to GridCells (should be of length 4 at least) as a parameter


    v[0] = GridCell<T>(i,   j,   g.getValue(i,   j, nrows, ncols));
    v[1] = GridCell<T>(i,   j+1, g.getValue(i,   j+1, nrows, ncols));
    v[2] = GridCell<T>(i+1, j,   g.getValue(i+1, j, nrows, ncols));
    v[3] = GridCell<T>(i+1, j+1, g.getValue(i+1, j+1, nrows, ncols));
    
    if (debug_flag) { 
	int tmp = -1;
	printf("Calculating g.getIndex(%d, %d, %d, %d, %d)...\n", i, j, nrows, ncols, nbands);
	tmp = g.getIndex(i, j, nrows, ncols, nbands, true);
	printf("got %d\n", tmp);
	
	printf("Calculating g.getIndex(%d, %d, %d, %d, %d)...\n", i, j+1, nrows, ncols, nbands);
	tmp = g.getIndex(i, j+1, nrows, ncols, nbands, true);
	printf("got %d\n", tmp);
	
	printf("Calculating g.getIndex(%d, %d, %d, %d, %d)...\n", i+1, j, nrows, ncols, nbands);
	tmp = g.getIndex(i+1, j, true);
	printf("got %d\n", tmp);
	
	printf("Calculating g.getIndex(%d, %d, %d, %d, %d)...\n", i+1, j+1, nrows, ncols, nbands);
	tmp = g.getIndex(i+1, j+1, true);
	printf("got %d\n", tmp);

	printf("Assigning v[0] = g.getValue(%d, %d) == g.data[%d] == %f\n", i, j, g.getIndex(i, j, nrows, ncols, nbands, false), g.getValue(i, j, nrows, ncols)); 
	printf("Assigning v[1] = g.getValue(%d, %d) == g.data[%d] == %f\n", i, j+1, g.getIndex(i, j+1, nrows, ncols, nbands, false), g.getValue(i, j+1, nrows, ncols)); 
	printf("Assigning v[2] = g.getValue(%d, %d) == g.data[%d] == %f\n", i+1, j, g.getIndex(i+1, j, nrows, ncols, nbands, false), g.getValue(i+1, j, nrows, ncols));
	printf("Assigning v[3] = g.getValue(%d, %d) == g.data[%d] == %f\n", i+1, j+1, g.getIndex(i+1, j+1, nrows, ncols, nbands, false), g.getValue(i+1, j+1, nrows, ncols)); 
    
        printf("\nPrinting v[0] stuff:\nv[0].value == %f; v[1].value == %f; v[2].value == %f; v[3].value == %f\nv[0].i == %d; v[1].i == %d; v[2].i == %d; v[3].i == %d\nv[0].j == %d; v[1].j == %d; v[2].j == %d; v[3].j == %d\n", v[0].value, v[1].value, v[2].value, v[3].value, v[0].i, v[1].i, v[2].i, v[3].i, v[0].j, v[1].j, v[2].j, v[3].j);
    }
    //return v;

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
Grid_tc<T> skyview ( const Grid_tc<T> & dem, int daz, double r )
{
    double timer_start = omp_get_wtime();

    static const double pi = 4.0 * atan(1.0);
    static const double deg2rad = pi / 180.0;
    const double interval = static_cast<double>(daz) / 360.0;
    const int dmax = floor(r / dem.dx());
    
    bool flag;
    int i, j, k, kk, a, d;
    double ii, jj, z, cosAz, sinAz, angle, Amax, sum;
    //std::vector<GridCell<T> > cells;
    // GridCell<T>* cells = (GridCell<T>*)(std::malloc(sizeof(T) * 4));  // XXX using a native C array of size 4 to be compatible with NVIDIA hardware
    // Grid_tc<T> out(dem, 1);
    Grid_tc<T> out;
    // XXX pasting in the copy constructor to keep this type trivial

    /* BEGIN PASTED SEMI COPY CONSTRUCTOR */
    out.nr = dem.nr;
    out.nc = dem.nc;
    out.nb = dem.nb;
    out.sz = dem.sz;
    out.vol = dem.vol;
    out.ilv = dem.ilv;
    out.nul = dem.nul;
    out.cw = 0;

    int it;
    std::stringstream ss;
    
    // initialize array to noData
    out.data = new T[out.vol];

    for ( it = 0; it < out.vol; it++ ) { out.data[it] = out.nul; }
    out.stats = new typename Grid_tc<T>::GridStats[out.nb];
    
    if (out.coordinateSys == nullptr) { out.coordinateSys = new char[10000]; } // XXX just using a very large size for now. Might want to FIXME by reconsidering the size to allocate here
    
    // set default band names
    //bandNames = new std::string[nb];
    out.bandNames = new char*[out.nb];
    
    for (int j = 0; j < out.nb; j++) { out.bandNames[j] = new char[10000]; }  // XXX just using a very large size for now. Might want to FIXME by reconsidering the size to allocate here

    for ( it = 0; it < out.nb; it++ ){
        ss.str(std::string());      // clear stringstream
        ss << "Band " << i + 1;
        //bandNames[i] = ss.str();
        strcpy(out.bandNames[it], ss.str().c_str());
    }
    
    // set geography
    //mapInfo = MapInfo_tc(g.mapInfo);
    out.mapInfo = MapInfo_tc(dem.mapInfo);  // XXX this might fail...


    //coordinateSys = std::string(g.coordinateSys);
    strcpy(out.coordinateSys, dem.coordinateSys);  // XXX using strcpy() instead of std::string as before FIXME Happening here... no coordinateSys valid location
    /* END PASTED SEMI COPY CONSTRUCTOR */
    
    int count = 0;

    // hard-coded values from dem here
    int nrows = dem.nrows(); int ncols = dem.ncols(); int nbands = dem.nbands();
    double dem_dx = dem.dx(); T dem_no_data = dem.noData();
    int dem_nb = dem.nb; int dem_nc = dem.nc;
    int dem_ilv = dem.ilv;

    bool debug_flag = false;  // XXX using this to make certain functions i'm debugging print only once

    // (old) CPU: #pragma omp parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, cells, count, angle, z) shared(interval, dmax, deg2rad, out) collapse(2)
    // GPU: #pragma omp target teams distribute parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, count, angle, z) shared(dem, interval, dmax, deg2rad, out) map(tofrom: dem.data[0:dem.get_volume()]) collapse(2) 
    // (new) CPU #pragma omp parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, count, angle, z) shared(interval, dmax, deg2rad, out) collapse(2) 
#pragma omp target teams distribute parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, count, angle, z) shared(dem, interval, dmax, deg2rad, out, nrows, ncols) map(tofrom: out.data[0:out.get_volume()], dem.data[0:dem.get_volume()]) collapse(2)
    //#pragma omp parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, count, angle, z) shared(interval, dmax, deg2rad, out) collapse(2)
    for ( i = 0; i < nrows; ++i ) {
        for ( j = 0; j < ncols; ++j ) { 
	    GridCell<T>* cells = (GridCell<T>*)(std::malloc(sizeof(GridCell<T>) * 4));  // XXX declaring cells here, since it's used internally only anyway
            //std::cout << "\rskyview " << (static_cast<double>(count++) / dem.size()) << "%        ";
            // Using printf() instead
	    //printf("\rskyview %f%%        \n", (static_cast<double>(count++) / dem.size()));
        
            // k = dem.getIndex(i, j);
	    // Begin transposed getIndex function
	    if ( dem_ilv == Grid_tc<T>::INTERLEAVE_BIP ) { k = dem_nb * (i*dem_nc + j); }
	    if ( dem_ilv == Grid_tc<T>::INTERLEAVE_BIP ) { k = dem_nb * (i*dem_nc + j); }
	    else if ( dem_ilv == Grid_tc<T>::INTERLEAVE_BIL ) { k = i*dem_nc*dem_nb + j; }
	    else { k = i*dem_nc + j; } // assume bsq
	    // end transposed getIndex function

        
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
                    

                    if ( ii < 0 || ii >= nrows
                            || jj < 0 || jj >= ncols )
                        { break; }

                    // interpolate value at point
                    flag = false;
                    if (k == 0 && !debug_flag) { getNearest4_tc(dem, ii, jj, cells, nrows, ncols, nbands, true); debug_flag = true; }
		    else { getNearest4_tc(dem, ii, jj, cells, nrows, ncols, nbands, false); }

		// In getNearest4_tc, cells is set to an array of size 4. It is always size 4, so use that instead of
		// what would normally be cells.size()
		for ( kk = 0; !flag && kk < 4; ++kk ){
                        if ( cells[kk].value == dem_no_data ) {
                            flag = true;
			}
                    } 
                    if ( flag ) {
                       angle = 0.0;
                       continue;
                    }

		    // FIXME: ERROR IN CORRECTNESS IN GPU VERSION CAUSED BY INCORRECT SET OF z
                    if (k == 0 /*&& !debug_flag*/) { z = bilinearInterpolate_tc(ii, jj, cells, false); /*debug_flag = true;*/ }
		    else { z = bilinearInterpolate_tc(ii, jj, cells, false); }

                    // track maximum angle from horizontal
                    angle = atan((z - dem[k]) / (d * dem_dx));
                    if ( angle > Amax ) { Amax = angle; /* if (k == 0) { printf("New Amax set: Amax == atan((%f - %f) / (%d * %f))\n", z, dem[k], d, dem_dx); }*/ }
                }
		/* XXX my printf stuff for debugging setting of sum
		if (k == 0) { printf("(k == 0) Amax, before setting it to its cosine, is %f\n", Amax); }
                Amax = cos(Amax);
                if (k == 0) { printf("sum before addition: %f; ", sum); sum += Amax * Amax; printf("sum after addition: %f\n", sum); }
		else { sum += Amax * Amax; }
		*/
		Amax = cos(Amax);
		sum += Amax * Amax;
            }
	    if (k == 0) { printf("k == 0, so we are assigning out[%d] = %f * %f == %f\n", k, sum, interval, sum * interval); }
            out[k] = sum * interval;
        }
    }
    // std::cout << "\rskyview " << (static_cast<double>(count++) / dem.size()) << "%\n";
    // printf("\rskyview %f%%\n", (static_cast<double>(count++) / dem.size()));
 
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
Grid_tc<T> prominence ( const Grid_tc<T> & dem, int daz, double r )
{
    double timer_start = omp_get_wtime();

    static const double pi = 4.0 * atan(1.0);
    static const double deg2rad = pi / 180.0;
    const double interval = static_cast<double>(daz) / 360.0;
    const int dmax = floor(r / dem.dx());
    
    bool flag;
    int i, j, k, kk, a, d;
    double ii, jj, z, cosAz, sinAz, angle, Amax, sum;
    //std::vector<GridCell<T> > cells;
    GridCell<T>* cells;
    //Grid_tc<T> out(dem, 1);
    Grid_tc<T> out;
    // XXX pasting in the copy constructor to keep this type trivial

    /* BEGIN PASTED SEMI COPY CONSTRUCTOR */
    out.nr = dem.nr;
    out.nc = dem.nc;
    out.nb = dem.nb;
    out.sz = dem.sz;
    out.vol = dem.vol;
    out.ilv = dem.ilv;
    out.nul = dem.nul;
    out.cw = 0;

    int it;
    std::stringstream ss;
    
    // initialize array to noData
    out.data = new T[out.vol];

    for ( it = 0; it < out.vol; it++ ) { out.data[it] = out.nul; }
    out.stats = new typename Grid_tc<T>::GridStats[out.nb];
    
    if (out.coordinateSys == nullptr) { out.coordinateSys = new char[10000]; } // XXX just using a very large size for now. Might want to FIXME by reconsidering the size to allocate here
    
    // set default band names
    //bandNames = new std::string[nb];
    out.bandNames = new char*[out.nb];
    
    for (int j = 0; j < out.nb; j++) { out.bandNames[j] = new char[10000]; }  // XXX just using a very large size for now. Might want to FIXME by reconsidering the size to allocate here

    for ( it = 0; it < out.nb; it++ ){
        ss.str(std::string());      // clear stringstream
        ss << "Band " << i + 1;
        //bandNames[i] = ss.str();
        strcpy(out.bandNames[it], ss.str().c_str());
    }
    
    // set geography
    //mapInfo = MapInfo_tc(g.mapInfo);
    out.mapInfo = MapInfo_tc(dem.mapInfo);  // XXX this might fail...


    //coordinateSys = std::string(g.coordinateSys);
    strcpy(out.coordinateSys, dem.coordinateSys);  // XXX using strcpy() instead of std::string as before FIXME Happening here... no coordinateSys valid location
    /* END PASTED SEMI COPY CONSTRUCTOR */
    
    int count = 0;
   
    // CPU: #pragma omp parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, cells, count, angle, z) shared(interval, dmax, deg2rad, out) collapse(2)
    // GPU: #pragma omp target teams distribute parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, cells, count, angle, z) shared(interval, dmax, deg2rad, out) collapse(2)
    #pragma omp target teams distribute parallel for private(flag, i, ii, j, jj, k, kk, a, d, sinAz, cosAz, Amax, sum, cells, count, angle, z) shared(interval, dmax, deg2rad, out) collapse(2)
    for ( i = 0; i < dem.nrows(); ++i ) {
        for ( j = 0; j < dem.ncols(); ++j ) {
	    //std::cout << "\rprominence " << (static_cast<double>(count++) / dem.size()) << "%        ";
	    printf("\rprominence %f%%        \n", (static_cast<double>(count++) / dem.size()));
            
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
                    // cells = getNearest4_tc(dem, ii, jj);  XXX UPDATE ME WITH NEW SIGNATURE
		    
		    // In getNearest4_tc, cells is set to an array of size 4. It is always size 4, so use that instead of
		    // what would normally be cells.size()
                    
		    for ( kk = 0; !flag && kk < 4; ++kk ){
                       if ( cells[kk].value == dem.noData() )
				flag = true;
		    }
                    if ( flag ) {
                        angle = 0.0;
                        continue;
                    }
                    z = bilinearInterpolate_tc(ii, jj, cells);
                    
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
    //std::cout << "\rprominence " << (static_cast<double>(count++) / dem.size()) << "%\n";
    printf("\rprominence %f%%\n", (static_cast<double>(count++) / dem.size()));
   
    double timer_end = omp_get_wtime();
    std::cout << "Prominence took this many seconds: " << timer_end - timer_start << std::endl;
 
    return out;
}

#endif // YOUNG_GRID_OPERATIONS_yyyymmdd

