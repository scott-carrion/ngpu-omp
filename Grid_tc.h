#ifndef YOUNG_GRID_tc
#define YOUNG_GRID_tc

#include <cmath>     // sqrt
#include <limits>    // std::numeric_limits
#include <sstream>   // std::stringstream
#include <string>    // std::string

#include <omp.h>     // OpenMP header

#include "MapInfo_tc.h"

template <typename T>
struct Grid_tc {
    //public:
        typedef struct {
            T minimum;
            T maximum;
            double mean;
            double stdev;
        } GridStats;
    //private:
        // constants
        static const unsigned char INTERLEAVE_NUL = 0;
        static const unsigned char INTERLEAVE_BSQ = 1;
        static const unsigned char INTERLEAVE_BIP = 2;
        static const unsigned char INTERLEAVE_BIL = 3;
        
        // dimensions
        int         nr;                     // number rows
        int         nc;                     // number columns
        int         nb;                     // number bands
        int         sz;                     // number cells per band
        int         vol;                    // total number cells
        //std::string * bandNames;            // names of each band
        char**       bandNames;		    // XXX names of each band, using char* so it's trivially copyable
        int         cw;                     // collar width of noData
        
        // data
        T *         data;                   // data values
        T           nul;                    // no data (null) value
        unsigned char ilv;                  // cell interleave format
        GridStats *  stats;
        
        // geographic data
        MapInfo_tc     mapInfo;                // geographic information

        //std::string coordinateSys;  // coordinate system details
        char* coordinateSys;          // coordinate system details
        
        // helpers
        unsigned char checkInterleave(
            const std::string & ) const;
        
    //public:
        // constructors
        //Grid_tc(int r=1, int c=1, int b=1);    // default
        //Grid_tc(const Grid_tc<T> &);              // copy  XXX COPY CONSTRUCTOR NOT ALLOWED, THIS HAS TO BE TRIVIALLY COPYABLE!
        //Grid_tc(const Grid_tc<T> &, int);         // copy all but data  XXX NOT ALLOWED, THIS HAS TO BE TRIVIALLY COPYABLE!
        
        //destructor
        //~Grid_tc();  // XXX DESTRUCTOR NOT ALLOWED, THIS HAS TO BE TRIVIALLY COPYABLE!
        
        // operators   XXX COPY ASSIGNMENT NOT ALLOWED, THIS HAS TO BE TRIVIALLY COPYABLE!
        //Grid_tc<T> &   operator=(              // assignment
        //              const Grid_tc<T> &);
        T &         operator[](int);        // element access
        T           operator[](int) const;  // element access
        
        // getters
        int         size() const;
        int         volume() const;
        int         ncols() const;
        int         nrows() const;
        int         nbands() const;
        double      x0() const;
        double      y0() const;
        double      dx() const;
        double      dy() const;
        T           noData() const;
        int         get_size() const;
        int         get_volume() const;
        int         get_ncols() const;
        int         get_nrows() const;
        int         get_nbands() const;
        double      get_x0() const;
        double      get_y0() const;
        double      get_dx() const;
        double      get_dy() const;
        T           get_noData() const;
        std::string get_interleave() const;
        MapInfo_tc     get_mapInfo() const;
        std::string get_coordinateSys() const;
        std::string get_bandName(int) const;
        int         get_collarWidth() const;
        
        // grid navigation
        int         getIndex(int, int, int) // index for row, col, band
                      const;
        int         getIndex(int, int)
                      const;
        T           getValue(int) const;    // get value at index
        T           getValue(int, int)      // get value at row, col
                      const;
        T           getValue(int, int, int) // value at row, col, band
                      const;
        void        getCoordinates(int*,    // get row, col, band at index
                      int*, int*, int) const;
        
        // setters
        void        clearGeography();       // remove geographic info
        void        setGeography(           // set getInfo and coordinateSys
                      const std::string&,
                      const std::string&);
        void        set_nbands(int);        // changes number of bands
        void        set_x0(double);         // set anchor point
        void        set_y0(double);
        void        set_dx(double);
        void        set_dy(double);
        void        set_bandName(int,
                      const std::string &);
        void        set_noData(T);          // changes noData value
        void        set_interleave(         // change interleave
                      const std::string &); //   (doesn't reorganize data!)
        void        set_collarWidth(int);   // set width of the noData band
        
        void        setValue(int, T);       // set value at index
        void        setValue(int, int, T);  // set value at row, col
        void        setValue(int, int, int, // set value at row, col, band
                      T);
        
        // statistics
        void        calculateStatistics();  // compute statistics
        const GridStats & getStats(int)     // get stats for band
            const;
        const GridStats & getStats() const; // get stats for band 0
};

// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////
/*
// Default constructor.
template <typename T>
Grid_tc<T>::Grid_tc ( int c, int r, int b )
: nr(r), nc(c), nb(b), sz(r*c), vol(r*c*b), ilv(INTERLEAVE_BSQ),
    nul(-1), cw(0)
{

    int i;
    std::stringstream ss;
    
    // set default noData to -9999 if possible
    if ( sizeof(T) > 1 ) nul = -9999;
    
    // initialize array to noData
    data = new T[vol];

    for ( i = 0; i < vol; i++ ) { data[i] = nul; }
    stats = new GridStats[nb];

    if (coordinateSys == nullptr) { coordinateSys = new char[10000]; } // XXX just using a very large size for now. Might want to FIXME by reconsidering the size to allocate here

    // set default band names
    //bandNames = new std::string[nb];
    bandNames = new char*[nb];
    for (int j = 0; j < nb; j++) { bandNames[j] = new char[10000]; }  // XXX just using a very large size for now. Might want to FIXME by reconsidering the size to allocate here

    std::stringstream* ss_ptr = &ss;

    for ( i = 0; i < nb; i++ ) {
        ss.str(std::string());      // clear stringstream
        ss << "Band " << i + 1;
        //bandNames[i] = ss.str();
        strcpy(bandNames[i], ss.str().c_str());  // XXX using strcpy() instead of std::string interface for copy triviality
        //ss_ptr->str(std::string());
	//ss_ptr->operator<<("Band ").operator<<(i + 1);  // XXX FIXME I JUST COMMENTED OUT THIS STUFF FOR NOW, I DON'T SEE THE POINT OF IT. IT CAUSES A BAD free() SOMEWHERE I THINK!
	//bandNames[i] = ss_ptr->str();
	//strcpy(bandNames[i+1], ss_ptr->str().c_str());  // XXX using strcpy() instread of std::string interface for copy triviality
    }

}
*/
/*
// Copy constructor.
template <typename T>
Grid_tc<T>::Grid_tc ( const Grid_tc<T> & g )
: nr(g.nr), nc(g.nc), nb(g.nb), sz(g.sz), vol(g.vol), ilv(g.ilv),
    nul(g.nul), cw(g.cw)
{
    int i;

    // initialize array to noData
    data = new T[vol];
    for ( i = 0; i < vol; i++ ) { data[i] = g.data[i]; }
    stats = new GridStats[nb];
    
    // set band names to input grid's band names
    bandNames = new std::string[nb];

    for ( i = 0; i < nb; i++ ) {
        bandNames[i] = g.bandNames[i].c_str();
        stats[i] = g.stats[i];
    }
    
    // set geography
    mapInfo = MapInfo_tc(g.mapInfo);
    coordinateSys = std::string(g.coordinateSys);
}
*/
// Semi-copy constructors; copies geography and dimensions, but
// not data.
/*
template <typename T>
Grid_tc<T>::Grid_tc ( const Grid_tc<T> & g, int k )
: nr(g.nr), nc(g.nc), nb(g.nb), sz(g.sz), vol(g.vol), ilv(g.ilv),
    nul(g.nul), cw(0)
{
    int i;
    std::stringstream ss;
    
    // initialize array to noData
    data = new T[vol];

    for ( i = 0; i < vol; i++ ) { data[i] = nul; }
    stats = new GridStats[nb];
    
    if (coordinateSys == nullptr) { coordinateSys = new char[10000]; } // XXX just using a very large size for now. Might want to FIXME by reconsidering the size to allocate here
    
    // set default band names
    //bandNames = new std::string[nb];
    bandNames = new char*[nb];
    
    for (int j = 0; j < nb; j++) { bandNames[j] = new char[10000]; }  // XXX just using a very large size for now. Might want to FIXME by reconsidering the size to allocate here

    for ( i = 0; i < nb; i++ ){
        ss.str(std::string());      // clear stringstream
        ss << "Band " << i + 1;
        //bandNames[i] = ss.str();
        strcpy(bandNames[i], ss.str().c_str());
    }
    
    // set geography
    mapInfo = MapInfo_tc(g.mapInfo);
    //coordinateSys = std::string(g.coordinateSys);
    strcpy(coordinateSys, g.coordinateSys);  // XXX using strcpy() instead of std::string as before FIXME Happening here... no coordinateSys valid location
}
*/
/*
// Destructor.
template <typename T>
Grid_tc<T>::~Grid_tc ()
{
    delete [] data;
    delete [] bandNames;
    delete [] stats;
}

// OPERATORS //////////////////////////////////////////////////////////////////

template <typename T>
Grid_tc<T> & Grid_tc<T>::operator= ( const Grid_tc<T> & g )
{
    int i;
    
    nr = g.nr;
    nc = g.nc;
    nb = g.nb;
    sz = g.sz;
    vol = g.vol;
    ilv = g.ilv;
    nul = g.nul;
    cw = g.cw;
    
    // set array to input grid's values
    if ( data ) delete [] data;
    data = new T[vol];

    for ( i = 0; i < vol; i++ ) { data[i] = g.data[i]; }
    stats = new GridStats[nb];
    
    // set band names to input grid's band names
    if ( bandNames ){ delete [] bandNames; }
    bandNames = new std::string[nb];

    for ( i = 0; i < nb; i++ ){
        bandNames[i] = g.bandNames[i].c_str();
        stats[i] = g.stats[i];
    }
    
    // set geography
    // mapInfo = MapInfo(g.mapInfo);
    mapInfo = MapInfo_tc(g.mapInfo);  // XXX this is to be compatible (albeit temporarily) with the MapInfo_tc constructor
    coordinateSys = std::string(g.coordinateSys);
    
    return *this;
}
*/

template <typename T>
T & Grid_tc<T>::operator[] ( int idx )
{ return data[idx]; }

template <typename T>
T Grid_tc<T>::operator[] ( int idx ) const
{ return data[idx]; }

// HELPERS //////////////////////////////////////////////////////////

// Get an interleave number from the interleave string. Print a
// warning if not a valid interleave format.
template <typename T>
unsigned char Grid_tc<T>::checkInterleave ( const std::string & s )
    const
{
    if ( s == "bsq" ) return INTERLEAVE_BSQ;
    if ( s == "bip" ) return INTERLEAVE_BIP;
    if ( s == "bil" ) return INTERLEAVE_BIL;
    printf("WARNING: '%s' not recognized interleave format, assuming bsq\n",
        s.c_str());
    return INTERLEAVE_NUL;
}

// GETTERS //////////////////////////////////////////////////////////

template <typename T> int Grid_tc<T>::size() const
{ return sz; }
template <typename T> int Grid_tc<T>::volume() const
{ return vol; }
template <typename T> int Grid_tc<T>::ncols() const
{ return nc; }
template <typename T> int Grid_tc<T>::nrows() const
{ return nr; }
template <typename T> int Grid_tc<T>::nbands() const
{ return nb; }
template <typename T> double Grid_tc<T>::x0() const
{ return mapInfo.x0; }
template <typename T> double Grid_tc<T>::y0() const
{ return mapInfo.y0; }
template <typename T> double Grid_tc<T>::dx() const
{ return mapInfo.dx; }
template <typename T> double Grid_tc<T>::dy() const
{ return mapInfo.dy; }
template <typename T> T      Grid_tc<T>::noData() const
{ return nul; }

// deprecated
template <typename T> int Grid_tc<T>::get_size() const
{ return sz; }
template <typename T> int Grid_tc<T>::get_volume() const
{ return vol; }
template <typename T> int Grid_tc<T>::get_ncols() const
{ return nc; }
template <typename T> int Grid_tc<T>::get_nrows() const
{ return nr; }
template <typename T> int Grid_tc<T>::get_nbands() const
{ return nb; }
template <typename T> double Grid_tc<T>::get_x0() const
{ return mapInfo.x0; }
template <typename T> double Grid_tc<T>::get_y0() const
{ return mapInfo.y0; }
template <typename T> double Grid_tc<T>::get_dx() const
{ return mapInfo.dx; }
template <typename T> double Grid_tc<T>::get_dy() const
{ return mapInfo.dy; }
template <typename T> T Grid_tc<T>::get_noData() const
{ return nul; }

template <typename T> std::string Grid_tc<T>::get_interleave () const
{
    if ( ilv == INTERLEAVE_BSQ ) return "bsq";
    if ( ilv == INTERLEAVE_BIP ) return "bip";
    if ( ilv == INTERLEAVE_BIL ) return "bil";
    return "nul";
}

template <typename T> MapInfo_tc Grid_tc<T>::get_mapInfo () const
{ return mapInfo; }

template <typename T> std::string Grid_tc<T>::get_coordinateSys () const
{ return std::string(coordinateSys); }

template <typename T> std::string Grid_tc<T>::get_bandName ( int band )
    const
{ return std::string(bandNames[band]); }

template <typename T> int Grid_tc<T>::get_collarWidth () const
{ return cw; }

// fetch index for a given i (row), j (column), and k (band).
template <typename T>
int Grid_tc<T>::getIndex ( int i, int j, int k ) const
{
    if ( ilv == INTERLEAVE_BIP ) return nb * (i*nc + j) + k;
    if ( ilv == INTERLEAVE_BIP ) return nb * (i*nc + j) + k;
    else if ( ilv == INTERLEAVE_BIL ) return i*nc*nb + k*nc + j;
    else return k*sz + i*nc + j; // assume bsq
}
template <typename T> int Grid_tc<T>::getIndex(int i, int j) const
{ return getIndex(i, j, 0); }

// fetch values given index (idx) or coordinates (row i, col j, band k)
template <typename T> T Grid_tc<T>::getValue ( int idx ) const
{ return data[idx]; }
template <typename T> T Grid_tc<T>::getValue ( int i, int j ) const
{ return data[getIndex(i, j)]; }
template <typename T> T Grid_tc<T>::getValue ( int i, int j, int k ) const
{ return data[getIndex(i, j, k)]; }

// fetch coordinate (row, col, band) for a given index.
template <typename T>
void Grid_tc<T>::getCoordinates ( int * row, int * col, int * band,
    int idx ) const
{
    if ( ilv == INTERLEAVE_BIP ) {
        *row = idx / (nc * nb);
        *col = idx % (nc * nb) / nb;
        *band = idx % nb;
    }
    else if ( ilv == INTERLEAVE_BIL ) {
        *row = idx / (nc * nb);
        *col = idx % nc;
        *band = (idx / nc) % nb;
    }
    else {
        // assume bsq
        *row = (idx % sz) / nc;
        *col = (idx % sz) % nc;
        *band = idx / sz;
    }
}

// SETTERS ////////////////////////////////////////////////////////////////////

template <typename T>
void Grid_tc<T>::clearGeography ()
{
    mapInfo = MapInfo_tc();
    //coordinateSys = std::string();
    strcpy(coordinateSys, "");  // XXX As before, using strcpy() instead
}

template <typename T>
void Grid_tc<T>::setGeography ( const std::string & info,
    const std::string & sys)
{
    // Basically doing the parameterized constructor's job here, since it isn't
    // allowed in order to be a trivial type
    mapInfo = MapInfo_tc(); mapInfo.parseENVI(info);




    //coordinateSys = std::string(sys);
    coordinateSys = new char[sys.length()+1];
    strcpy(coordinateSys, sys.c_str());  // XXX As before, using strcpy() instead
}

template <typename T>
void Grid_tc<T>::set_nbands(int n)
{
    //int               i;  // iterator
    std::stringstream ss;
    
    // adjust dimensions
    nb = n;
    vol = sz * nb;
    
    // replace value array with new dimensions
    delete data;
    data = new T[vol];
    for ( int i = 0; i < vol; i++ ){ data[i] = nul; }
    
    // replace band name array with new number of bands
    delete bandNames;
    //bandNames = new std::string[nb];
    bandNames = new char*[nb];  // XXX use char* instead for copy triviality

    for (int i = 0; i < nb; i++){
        ss.str(std::string()); // clear stringstream;
        ss << "Band " << i + 1;
        //bandNames[i] = ss.str();
        strcpy(bandNames[i], ss.str().c_str());
    }
}

template <typename T> void Grid_tc<T>::set_x0(double v){ mapInfo.x0 = v; }
template <typename T> void Grid_tc<T>::set_y0(double v){ mapInfo.y0 = v; }
template <typename T> void Grid_tc<T>::set_dx(double v){ mapInfo.dx = v; }
template <typename T> void Grid_tc<T>::set_dy(double v){ mapInfo.dy = v; }

template <typename T>
void Grid_tc<T>::set_bandName ( int k, const std::string & name )
{
    if ( k >= 0 && k < nb ) { strcpy(bandNames[k], name.c_str()); }
}

// change the noData value - does not change any of the data in the grid.
template <typename T>
void Grid_tc<T>::set_noData( T nd )
{
    nul = nd;
}

// change the interleave - does not rearrange data in the grid.
template <typename T>
void Grid_tc<T>::set_interleave( const std::string & s )
{
    ilv = checkInterleave(s);
}

// set the width of the no-data band that surrounds the valid data values in
// the grid.
template <typename T>
void Grid_tc<T>::set_collarWidth ( int w )
{
    cw = w;
}

// set values
template<typename T>
void Grid_tc<T>::setValue ( int idx, T val )
{
    data[idx] = val;
}

template<typename T>
void Grid_tc<T>::setValue ( int i, int j, T val )
{
    setValue(getIndex(i, j), val);
}

template<typename T>
void Grid_tc<T>::setValue ( int i, int j, int k, T val )
{
    setValue(getIndex(i, j, k), val);
}

// STATISTICS /////////////////////////////////////////////////////////////////

// Compute statistics for each band.
template<typename T>
void Grid_tc<T>::calculateStatistics ()
{
    int i, j, k, idx;
    int count;
 
    for ( k = 0; k < nb; ++k ) {
        stats[k].minimum = std::numeric_limits<T>::max();
        stats[k].maximum = std::numeric_limits<T>::min();
        stats[k].mean    = 0.0;
        stats[k].stdev   = 0.0;
        count = 0;
        
        // compute min, max, mean, stdev
        for ( i = 0; i < nr; ++i ){
            for ( j = 0; j < nc; ++j ){
                idx = getIndex(i, j, k);
                if ( data[idx] == nul ) continue;
                if ( data[idx] < stats[k].minimum )
                    stats[k].minimum = data[idx];
                if ( stats[k].maximum < data[idx] )
                    stats[k].maximum = data[idx];
                stats[k].mean += data[idx];                 // sum
                stats[k].stdev += data[idx] * data[idx];    // sum square
                ++count;
            }
        }
        stats[k].mean = stats[k].mean / count;
        stats[k].stdev = sqrt(stats[k].stdev / count
            - (stats[k].mean * stats[k].mean));
    }
}

// Get statistics for band.
template<typename T>
const typename Grid_tc<T>::GridStats & Grid_tc<T>::getStats ( int k ) const
{ return stats[k]; }

template<typename T>
const typename Grid_tc<T>::GridStats & Grid_tc<T>::getStats () const
{ return stats[0]; }

#endif // YOUNG_GRID_tc

