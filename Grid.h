#ifndef YOUNG_GRID
#define YOUNG_GRID

#include <cmath>     // sqrt
#include <limits>    // std::numeric_limits
#include <sstream>   // std::stringstream
#include <string>    // std::string

#include "MapInfo.h"

template <typename T>
class Grid {
    public:
        typedef struct {
            T minimum;
            T maximum;
            double mean;
            double stdev;
        } GridStats;
    private:
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
        std::string * bandNames;            // names of each band
        int         cw;                     // collar width of noData
        
        // data
        T *         data;                   // data values
        T           nul;                    // no data (null) value
        unsigned char ilv;                  // cell interleave format
        GridStats *  stats;
        
        // geographic data
        MapInfo     mapInfo;                // geographic information
        std::string coordinateSys;          // coordinate system details
        
        // helpers
        unsigned char checkInterleave(
            const std::string & ) const;
        
    public:
        // constructors
        Grid(int r=1, int c=1, int b=1);    // default
        Grid(const Grid<T> &);              // copy
        Grid(const Grid<T> &, int);         // copy all but data
        
        //destructor
        ~Grid();
        
        // operators
        Grid<T> &   operator=(              // assignment
                      const Grid<T> &);
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
        MapInfo     get_mapInfo() const;
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

// Default constructor.
template <typename T>
Grid<T>::Grid ( int c, int r, int b )
: nr(r), nc(c), nb(b), sz(r*c), vol(r*c*b), ilv(INTERLEAVE_BSQ),
    nul(-1), cw(0)
{
    int i;
    std::stringstream ss;
    
    // set default noData to -9999 if possible
    if ( sizeof(T) > 1 ) nul = -9999;
    
    // initialize array to noData
    data = new T[vol];
    for ( i = 0; i < vol; i++ ) data[i] = nul;
    stats = new GridStats[nb];
    
    // set default band names
    bandNames = new std::string[nb];
    for ( i = 0; i < nb; i++ ) {
        ss.str(std::string());      // clear stringstream
        ss << "Band " << i + 1;
        bandNames[i] = ss.str();
    }
}

// Copy constructor.
template <typename T>
Grid<T>::Grid ( const Grid<T> & g )
: nr(g.nr), nc(g.nc), nb(g.nb), sz(g.sz), vol(g.vol), ilv(g.ilv),
    nul(g.nul), cw(g.cw)
{
    int i;

    // initialize array to noData
    data = new T[vol];
    for ( i = 0; i < vol; i++ ) data[i] = g.data[i];
    stats = new GridStats[nb];
    
    // set band names to input grid's band names
    bandNames = new std::string[nb];
    for ( i = 0; i < nb; i++ ) {
        bandNames[i] = g.bandNames[i].c_str();
        stats[i] = g.stats[i];
    }
    
    // set geography
    mapInfo = MapInfo(g.mapInfo);
    coordinateSys = std::string(g.coordinateSys);
}

// Semi-copy constructors; copies geography and dimensions, but
// not data.
template <typename T>
Grid<T>::Grid ( const Grid<T> & g, int k )
: nr(g.nr), nc(g.nc), nb(g.nb), sz(g.sz), vol(g.vol), ilv(g.ilv),
    nul(g.nul), cw(0)
{
    int i;
    std::stringstream ss;
    
    // initialize array to noData
    data = new T[vol];
    for ( i = 0; i < vol; i++ ) data[i] = nul;
    stats = new GridStats[nb];
    
    // set default band names
    bandNames = new std::string[nb];
    for ( i = 0; i < nb; i++ ){
        ss.str(std::string());      // clear stringstream
        ss << "Band " << i + 1;
        bandNames[i] = ss.str();
    }
    
    // set geography
    mapInfo = MapInfo(g.mapInfo);
    coordinateSys = std::string(g.coordinateSys);
}

// Destructor.
template <typename T>
Grid<T>::~Grid ()
{
    delete [] data;
    delete [] bandNames;
    delete [] stats;
}

// OPERATORS //////////////////////////////////////////////////////////////////

template <typename T>
Grid<T> & Grid<T>::operator= ( const Grid<T> & g )
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
    for ( i = 0; i < vol; i++ ) data[i] = g.data[i];
    stats = new GridStats[nb];
    
    // set band names to input grid's band names
    if ( bandNames ){ delete [] bandNames; }
    bandNames = new std::string[nb];
    for ( i = 0; i < nb; i++ ){
        bandNames[i] = g.bandNames[i].c_str();
        stats[i] = g.stats[i];
    }
    
    // set geography
    mapInfo = MapInfo(g.mapInfo);
    coordinateSys = std::string(g.coordinateSys);
    
    return *this;
}

template <typename T>
T & Grid<T>::operator[] ( int idx )
{ return data[idx]; }

template <typename T>
T Grid<T>::operator[] ( int idx ) const
{ return data[idx]; }

// HELPERS //////////////////////////////////////////////////////////

// Get an interleave number from the interleave string. Print a
// warning if not a valid interleave format.
template <typename T>
unsigned char Grid<T>::checkInterleave ( const std::string & s )
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

template <typename T> int Grid<T>::size() const
{ return sz; }
template <typename T> int Grid<T>::volume() const
{ return vol; }
template <typename T> int Grid<T>::ncols() const
{ return nc; }
template <typename T> int Grid<T>::nrows() const
{ return nr; }
template <typename T> int Grid<T>::nbands() const
{ return nb; }
template <typename T> double Grid<T>::x0() const
{ return mapInfo.x0; }
template <typename T> double Grid<T>::y0() const
{ return mapInfo.y0; }
template <typename T> double Grid<T>::dx() const
{ return mapInfo.dx; }
template <typename T> double Grid<T>::dy() const
{ return mapInfo.dy; }
template <typename T> T      Grid<T>::noData() const
{ return nul; }

// deprecated
template <typename T> int Grid<T>::get_size() const
{ return sz; }
template <typename T> int Grid<T>::get_volume() const
{ return vol; }
template <typename T> int Grid<T>::get_ncols() const
{ return nc; }
template <typename T> int Grid<T>::get_nrows() const
{ return nr; }
template <typename T> int Grid<T>::get_nbands() const
{ return nb; }
template <typename T> double Grid<T>::get_x0() const
{ return mapInfo.x0; }
template <typename T> double Grid<T>::get_y0() const
{ return mapInfo.y0; }
template <typename T> double Grid<T>::get_dx() const
{ return mapInfo.dx; }
template <typename T> double Grid<T>::get_dy() const
{ return mapInfo.dy; }
template <typename T> T Grid<T>::get_noData() const
{ return nul; }

template <typename T> std::string Grid<T>::get_interleave () const
{
    if ( ilv == INTERLEAVE_BSQ ) return "bsq";
    if ( ilv == INTERLEAVE_BIP ) return "bip";
    if ( ilv == INTERLEAVE_BIL ) return "bil";
    return "nul";
}

template <typename T> MapInfo Grid<T>::get_mapInfo () const
{ return mapInfo; }

template <typename T> std::string Grid<T>::get_coordinateSys () const
{ return coordinateSys.c_str(); }

template <typename T> std::string Grid<T>::get_bandName ( int band )
    const
{ return bandNames[band].c_str(); }

template <typename T> int Grid<T>::get_collarWidth () const
{ return cw; }

// fetch index for a given i (row), j (column), and k (band).
template <typename T>
int Grid<T>::getIndex ( int i, int j, int k ) const
{
    if ( ilv == INTERLEAVE_BIP ) return nb * (i*nc + j) + k;
    if ( ilv == INTERLEAVE_BIP ) return nb * (i*nc + j) + k;
    else if ( ilv == INTERLEAVE_BIL ) return i*nc*nb + k*nc + j;
    else return k*sz + i*nc + j; // assume bsq
}
template <typename T> int Grid<T>::getIndex(int i, int j) const
{ return getIndex(i, j, 0); }

// fetch values given index (idx) or coordinates (row i, col j, band k)
template <typename T> T Grid<T>::getValue ( int idx ) const
{ return data[idx]; }
template <typename T> T Grid<T>::getValue ( int i, int j ) const
{ return data[getIndex(i, j)]; }
template <typename T> T Grid<T>::getValue ( int i, int j, int k ) const
{ return data[getIndex(i, j, k)]; }

// fetch coordinate (row, col, band) for a given index.
template <typename T>
void Grid<T>::getCoordinates ( int * row, int * col, int * band,
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
void Grid<T>::clearGeography ()
{
    mapInfo = MapInfo();
    coordinateSys = std::string();
}

template <typename T>
void Grid<T>::setGeography ( const std::string & info,
    const std::string & sys)
{
    mapInfo = MapInfo(info);
    coordinateSys = std::string(sys);
}

template <typename T>
void Grid<T>::set_nbands(int n)
{
    int               i;  // iterator
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
    bandNames = new std::string[nb];
    for (int i = 0; i < nb; i++){
        ss.str(std::string()); // clear stringstream;
        ss << "Band " << i + 1;
        bandNames[i] = ss.str();
    }
}

template <typename T> void Grid<T>::set_x0(double v){ mapInfo.x0 = v; }
template <typename T> void Grid<T>::set_y0(double v){ mapInfo.y0 = v; }
template <typename T> void Grid<T>::set_dx(double v){ mapInfo.dx = v; }
template <typename T> void Grid<T>::set_dy(double v){ mapInfo.dy = v; }

template <typename T>
void Grid<T>::set_bandName ( int k, const std::string & name )
{
    if ( k >= 0 && k < nb ) bandNames[k] = name.c_str();
}

// change the noData value - does not change any of the data in the grid.
template <typename T>
void Grid<T>::set_noData( T nd )
{
    nul = nd;
}

// change the interleave - does not rearrange data in the grid.
template <typename T>
void Grid<T>::set_interleave( const std::string & s )
{
    ilv = checkInterleave(s);
}

// set the width of the no-data band that surrounds the valid data values in
// the grid.
template <typename T>
void Grid<T>::set_collarWidth ( int w )
{
    cw = w;
}

// set values
template<typename T>
void Grid<T>::setValue ( int idx, T val )
{
    data[idx] = val;
}

template<typename T>
void Grid<T>::setValue ( int i, int j, T val )
{
    setValue(getIndex(i, j), val);
}

template<typename T>
void Grid<T>::setValue ( int i, int j, int k, T val )
{
    setValue(getIndex(i, j, k), val);
}

// STATISTICS /////////////////////////////////////////////////////////////////

// Compute statistics for each band.
template<typename T>
void Grid<T>::calculateStatistics ()
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
const typename Grid<T>::GridStats & Grid<T>::getStats ( int k ) const
{ return stats[k]; }

template<typename T>
const typename Grid<T>::GridStats & Grid<T>::getStats () const
{ return stats[0]; }

#endif // YOUNG_GRID_yyyymmdd