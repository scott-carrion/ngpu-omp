// Input/output operations for the Grid object.

#ifndef YOUNG_GRID_IO_20180528
#define YOUNG_GRID_IO_20180528

#include <cmath>    // abs, fabs, sqrt, cos, sin, tan, acos, asin, atan
#include <cstdlib>  // 
#include <fstream>  // std::ifstream, std::ofstream
#include <iostream> // std::cout
#include <string>   // std::string
#include <ctime>    // time_t, time, difftime
#include <set>      // std::set
#include <vector>   // std::vector

#include "Grid.h"

// HELPER FUNCTIONS ///////////////////////////////////////////////////////////

// Returns true if the operating system is little-endian.
bool osLittleEndian ()
{
    int num = 1;
    if ( *(char*) &num == 1) return true; // first byte contains the "1" bit
    return false;
}

// Returns string without leading and trailing whitespace.
std::string gio_trim ( const std::string & s )
{
    int a = s.find_first_not_of(" \t\n");
    int b = s.find_last_not_of(" \t\r\n"); // \r to handle CRLF
    return s.substr(a, b-a+1);
}

// Swaps bytes in data with n bytes from little-endian to big-endian, or vice-
// versa.
void byteswap(void *data, int n)
{
    unsigned char *start = (unsigned char*) data;
    unsigned char *end = start + n - 1;
    unsigned char swap;
    for ( ; start < end; ++start, --end ){
        swap = *start;
        *start = *end;
        *end = swap;
    }
}

// Swaps bytes in an array of length n.
template <class T>
void byteswaparray(T *data, int n){
    #pragma omp target teams distribute parallel for  // TODO: TEST ME
    for ( int i = 0; i < n; ++i ) byteswap(&(data[i]), sizeof(T));
}

// HEADER INFO ////////////////////////////////////////////////////////////////

// Container for information about the format and geography of grid data.
struct HeaderInfo {
    std::string filename;
    int nlines, ncolumns, ngrids, datatype, byteorder;
    double nodata;
    std::string interleave, mapinfo, coordinatesys, bandnames;
    
    HeaderInfo() :
        filename(""), nlines(0), ncolumns(0), ngrids(0), datatype(0),
        byteorder(0), nodata(-9999.0), interleave("bsq"), mapinfo(""),
        coordinatesys(""), bandnames("")
        {}
};

// Read information in an ENVI-format header file into a Grid object.
// -- Arguments --
// fn : the path and name of the ENVI-format header file (*.hdr). If
//      the filename has no extension, ".hdr" is appended. If the
//      extension is not ".hdr", it is replaced with ".hdr".
// -- Returns --
// A HeaderInfo object containing file metadata.
HeaderInfo read_ENVIheader ( std::string fn )
{    
    // check extension
    int d = fn.find_last_of('.');
    if ( fn.size() < 4 || d < 0 || d >= fn.size() )
        fn = fn + ".hdr";
    else if ( fn.substr(d, fn.size()).compare(".hdr") )
        fn = fn.substr(0, d) + ".hdr";
    
    // initialize metadata container
    HeaderInfo info;
    info.filename = fn.substr(0, fn.find_last_of('.')) + ".dat";
    
    // open file
    std::ifstream file(fn.c_str());
    if (!file){
        std::cout << "file open error: " << fn << "\n";
        return info;
    }
    
    // read file
    std::string line;
    std::string fld;
    std::string fldval;
    while ( getline(file, line) ){
        d = line.find(" = ");
        fld = gio_trim(line.substr(0, d));
        if ( line.length() == d + 3 ) continue;
        fldval = gio_trim(line.substr(d + 3));
        
        if ( fld == "samples" )
            info.ncolumns = atoi(fldval.c_str());
        else if ( fld == "lines" )
            info.nlines = atoi(fldval.c_str());
        else if ( fld == "bands" )
            info.ngrids = atoi(fldval.c_str());
        else if ( fld == "data type" )
            info.datatype = atoi(fldval.c_str());
        else if ( fld == "byte order" )
            info.byteorder = atoi(fldval.c_str());
        else if ( fld == "interleave" ) {
            if ( fldval == "BSQ" ) fldval = "bsq";
            else if ( fldval == "BIP" ) fldval = "bip";
            else if ( fldval == "BIL" ) fldval = "bil";
            else if ( fldval != "bsq" && fldval != "bip"
                    && fldval != "bil" ) {
                std::cout << "WARNING: interleave error, assuming bsq\n";
                fldval = "bsq";
            }
            info.interleave = fldval;
        }
        else if ( fld == "data ignore value" )
            info.nodata = atof(fldval.c_str());
        else if ( fld == "map info" )
            info.mapinfo = fldval;
        else if ( fld == "coordinate system string" )
            info.coordinatesys = fldval;
        else if ( fld == "band names" )
            info.bandnames = fldval;
    }
    
    // clean up
    file.close();
    
    // remove '{' and '}' from band names
    if ( info.bandnames.length() > 2 ){
        int a = info.bandnames.find_first_not_of("{");
        int b = info.bandnames.find_last_not_of("}") + 1;
        info.bandnames = info.bandnames.substr(a, b-a);
    }
    
    return info;
}

// GRID INPUT /////////////////////////////////////////////////////////////////

// Helper to read data from a binary grid file.
// -- Arguments --
// grid  : the grid into which data is being read.
// file  : the filestream object from which data is being read.
// data  : pointer of the type of the data being read.
// n     : the number of elements to read. Clamps to values between 0 and grid
//         volume.
// bytes : the number of bytes of each element being read.
// swap  : true if a byte swap is necessary.
template <class T, class U>
void loadGridHelper (Grid<T> & grid, std::ifstream & file, U * data, int n,
    int bytes, bool swap)
{
    if ( n < 0 ) n = 1;
    else if ( grid.get_volume() < n ) n = grid.get_volume();
    try { data = new U [n]; }
    catch (...) { std::cout << "ERR: memory allocation error\n"; exit(1); }
    file.read ((char *) data, n * bytes);
    if ( swap ) byteswaparray(data, n);
    for ( int i = 0; i < n; ++i ) grid[i] = (T) data[i];
    delete [] data;
}

// Read an ENVI format data file (*.dat and *.hdr).
// -- Arguments --
// fn : the path and name of the data file (*.dat).
// -- Returns --
// A floating point Grid object.
Grid<float> loadGrid (const std::string & fn)
{
    // read header information
    HeaderInfo info = read_ENVIheader(fn);
    
    if ( info.datatype == 0 ){
        std::cout << "file '" << fn << "' has no header file (*.hdr)\n";
        exit(1);
    }
    
    // create grid object
    Grid<float> grid = Grid<float>(info.ncolumns, info.nlines, info.ngrids);
    grid.set_noData(info.nodata);
    grid.set_interleave(info.interleave);
    
    // set geographic information
    grid.setGeography(info.mapinfo, info.coordinatesys);
    
    // set band names
    if ( info.bandnames.length() ){
        int i = 0;
        int n = info.bandnames.find_first_of(",", i) - i;
        int j = 0;
        while ( i < info.bandnames.size() - 1 && j < grid.get_nbands() ){
            grid.set_bandName(j++, gio_trim(info.bandnames.substr(i, n)));
            i += n + 1;
            n = info.bandnames.find_first_of(",", i) - i;
        }
    }
    
    // determine if a byte swap is necessary
    bool swap = (osLittleEndian() && (info.byteorder == 1))
        || ( !osLittleEndian() && (info.byteorder == 0) );
    
    // open file
    std::ifstream disk(fn.c_str(), std::ios::in | std::ios::binary);
    if (!disk) {
        std::cout << "file open error\n";
        exit (1);
    }
    
    //read data from file
    if ( info.datatype == 1 ){       // byte (unsigned integer) [8 bit]
        unsigned char * data;
        loadGridHelper(grid, disk, data, grid.get_volume(), 1, swap);
    }
    else if ( info.datatype == 2 ){  // signed integer [16 bit]
        short * data;
        loadGridHelper(grid, disk, data, grid.get_volume(), 2, swap);
    }
    else if ( info.datatype == 3 ){  // signed integer [32 bit]
        int * data;
        loadGridHelper(grid, disk, data, grid.get_volume(), 4, swap);
    }
    else if ( info.datatype == 4 ){  // float [32 bit]
        float * data;
        loadGridHelper(grid, disk, data, grid.get_volume(), 4, swap);
    }
    else if ( info.datatype == 5 ){  // double [64 bit]
        double * data;
        loadGridHelper(grid, disk, data, grid.get_volume(), 8, swap);
    }
    else if ( info.datatype == 6 ){  // float complex, real-imaginary pair
        //std::complex<float> * data;
        //loadGridHelper(grid, disk, data, grid.get_volume(), 8, swap);
    }
    else if ( info.datatype == 9 ){  // double complex, real-imaginary pair
        //std::complex<double> * data;
        //loadGridHelper(grid, disk, data, grid.get_volume(), 16, swap);
    }
    else if ( info.datatype == 12 ){ // unsigned integer [16 bit]
        unsigned short * data;
        loadGridHelper(grid, disk, data, grid.get_volume(), 2, swap);
    }
    else if ( info.datatype == 13 ){ // unsigned integer [32 bit]
        unsigned int * data;
        loadGridHelper(grid, disk, data, grid.get_volume(), 4, swap);
    }
    else if ( info.datatype == 14 ){ // signed integer [64 bit]
        long long * data;
        loadGridHelper(grid, disk, data, grid.get_volume(), 8, swap);
    }
    else if ( info.datatype == 15 ){ // unsigned integer [64 bit]
        unsigned long long * data;
        loadGridHelper(grid, disk, data, grid.get_volume(), 8, swap);
    }
    else {
        std::cout << "ERROR: data type error\n";
        exit(1);
    }
    disk.close ();
    
    grid.calculateStatistics();
    
    return grid;
}

// Helper to read data from a binary grid file.
// -- Arguments --
// grid   : the grid into which data is being read.
// k      : the band being read.
// ngrids : the number of bands in the data file.
// file   : the filestream object from which data is being read.
// data   : pointer of the type of the data being read.
// bytes  : the number of bytes of each element being read.
// swap   : true if a byte swap is necessary.
template <class T, class U>
void loadBandHelper (Grid<T> & grid, int k, int ngrids, std::ifstream & file,
    U * data, int bytes, bool swap)
{
    // load band-sequential
    if ( grid.get_interleave() == "bsq" ){
        int n = grid.get_size();
        file.ignore(n * k * bytes);
        loadGridHelper(grid, file, data, n, bytes, swap);
    }
    
    // load band-interleave by pixel
    else if ( grid.get_interleave() == "bip" ) {
        try { data = new U [bytes]; }
        catch (...) { std::cout << "ERR: memory allocation error\n"; exit(1); }
        for ( int i = 0; i < grid.get_size(); ++i ){
            file.ignore(k * bytes);                     // skip bands before k
            file.read((char *) data, bytes);
            if ( swap ) byteswap(data, bytes);
            grid[i] = (T) data[0];
            file.ignore((ngrids - k - 1) * bytes);      // skip bands after k
        }
        delete [] data;
    }
    
    // load band-interleave by line
    else if ( grid.get_interleave() == "bil" ){
        int nc = grid.get_ncols();
        int j, m = 0;
        try { data = new U [nc]; }
        catch (...) { std::cout << "ERR: memory allocation error\n"; exit(1); }
        for ( int i = 0; i < grid.get_nrows(); ++i ){
            file.ignore(k * nc * bytes);                // skip bands before k
            file.read((char *) data, nc * bytes);
            if ( swap ) byteswaparray(data, nc);
            for ( j = 0; j < nc; ++j ) grid[m++] = (T) data[j];
            file.ignore((ngrids - k - 1) * nc * bytes); // skip bands after k
        }
        delete [] data;
    }
}

// Read a specific band from an ENVI format data file (*.dat and *.hdr).
// -- Arguments --
// fn : the path and name of the data file (*.dat).
// k  : the band to be returned. Modulates k according to the maximum number
//      of bands.
// -- Returns --
// A floating point Grid object.
Grid<float> loadBand (const std::string & fn, int k)
{
    // read header information
    HeaderInfo info = read_ENVIheader(fn);
    
    if ( info.datatype == 0 ){
        std::cout << "file '" << fn << "' has no header file (*.hdr)\n";
        exit(1);
    }
    
    if ( k < 0 ) k = 0;
    else if ( info.ngrids <= k ) k = k % info.ngrids;
    
    // create grid object
    Grid<float> grid = Grid<float>(info.ncolumns, info.nlines, 1);
    grid.set_noData(info.nodata);
    grid.set_interleave(info.interleave);
    
    // set geographic information
    grid.setGeography(info.mapinfo, info.coordinatesys);
    
    // set band names
    if ( info.bandnames.length() ){
        int i = 0;
        int n = info.bandnames.find_first_of(",", i) - i;
        int j = 0;
        while ( i < info.bandnames.size() - 1 && j < info.ngrids ){
            if ( j == k ) grid.set_bandName(0, gio_trim(info.bandnames.substr(i, n)));
            ++j;
            i += n + 1;
            n = info.bandnames.find_first_of(",", i) - i;
        }
    }
    
    // determine if a byte swap is necessary
    bool swap = (osLittleEndian() && (info.byteorder == 1))
        || ( !osLittleEndian() && (info.byteorder == 0) );
    
    // open file
    std::ifstream disk(fn.c_str(), std::ios::in | std::ios::binary);
    if (!disk) {
        std::cout << "file open error\n";
        exit (1);
    }
    
    //read data from file
    if ( info.datatype == 1 ){       // byte (unsigned integer) [8 bit]
        unsigned char * data;
        loadBandHelper(grid, k, info.ngrids, disk, data, 1, swap);
    }
    else if ( info.datatype == 2 ){  // signed integer [16 bit]
        short * data;
        loadBandHelper(grid, k, info.ngrids, disk, data, 2, swap);
    }
    else if ( info.datatype == 3 ){  // signed integer [32 bit]
        int * data;
        loadBandHelper(grid, k, info.ngrids, disk, data, 4, swap);
    }
    else if ( info.datatype == 4 ){  // float [32 bit]
        float * data;
        loadBandHelper(grid, k, info.ngrids, disk, data, 4, swap);
    }
    else if ( info.datatype == 5 ){  // double [64 bit]
        double * data;
        loadBandHelper(grid, k, info.ngrids, disk, data, 8, swap);
    }
    else if ( info.datatype == 6 ){  // float complex, real-imaginary pair
        //std::complex<float> * data;
        //loadBandHelper(grid, k, info.ngrids, disk, data, 8, swap);
    }
    else if ( info.datatype == 9 ){  // double complex, real-imaginary pair
        //std::complex<double> * data;
        //loadBandHelper(grid, k, info.ngrids, disk, data, 16, swap);
    }
    else if ( info.datatype == 12 ){ // unsigned integer [16 bit]
        unsigned short * data;
        loadBandHelper(grid, k, info.ngrids, disk, data, 2, swap);
    }
    else if ( info.datatype == 13 ){ // unsigned integer [32 bit]
        unsigned int * data;
        loadBandHelper(grid, k, info.ngrids, disk, data, 4, swap);
    }
    else if ( info.datatype == 14 ){ // signed integer [64 bit]
        long long * data;
        loadBandHelper(grid, k, info.ngrids, disk, data, 8, swap);
    }
    else if ( info.datatype == 15 ){ // unsigned integer [64 bit]
        unsigned long long * data;
        loadBandHelper(grid, k, info.ngrids, disk, data, 8, swap);
    }
    else {
        std::cout << "ERROR: data type error\n";
        exit(1);
    }
    disk.close ();
    
    grid.calculateStatistics();
    
    return grid;
}

// GRID OUTPUT ////////////////////////////////////////////////////////////////

template <class T> int ENVI_dataType_code ( const Grid<T> & grid )
{
    return 4; // 32-bit float -- DEFAULT
}
template <> int ENVI_dataType_code ( const Grid<unsigned char> & grid )
{
    return 1; // 8-bit byte (unsigned int)
}
template <> int ENVI_dataType_code ( const Grid<short> & grid )
{
    return 2; // 16-bit int
}
template <> int ENVI_dataType_code ( const Grid<int> & grid )
{
    return 3; // 32-bit int
}
template <> int ENVI_dataType_code ( const Grid<double> & grid )
{
    return 5; // 64-bit float
}
template <> int ENVI_dataType_code ( const Grid<unsigned short> & grid )
{
    return 12; // 16-bit unsigned int
}
template <> int ENVI_dataType_code ( const Grid<unsigned int> & grid )
{
    return 13; // 32-bit unsigned int
}
template <> int ENVI_dataType_code ( const Grid<long> & grid )
{
    return 14; // 64-bit int
}
template <> int ENVI_dataType_code ( const Grid<unsigned long> & grid )
{
    return 15; // 64-bit unsigned int
}

// Saves the Grid object to a file and creates a companion header file.
// <!> hard-coded to write header file for floating point data type <!>
template <class T>
void saveGrid ( Grid<T> & grid, std::string filename )
{
    // open header file
    std::string headerfile;
    headerfile = filename.substr(0, filename.find_last_of(".")) + ".hdr";
    
    std::ofstream file(headerfile.c_str(), std::ios::out | std::ios::binary);
    if (!file){
        std::cout << "file open error: " << headerfile << "\n";
        exit(1);
    }
    
    // get byte order (LSF / little-endian, or MSF / big-endian)
    int byteorder = 1;
    if ( osLittleEndian() ) byteorder = 0;
    
    // write to header file
    file << "ENVI\n"
        << "samples = " << grid.get_ncols() << "\n"
        << "lines   = " << grid.get_nrows() << "\n"
        << "bands   = " << grid.get_nbands() << "\n"
        << "header offset = 0\n"
        << "file type = ENVI Standard\n"
        << "data type = " << ENVI_dataType_code(grid) << "\n" // <!> hard-coded to float <!>
        << "interleave = " << grid.get_interleave() << "\n"
        << "byte order = " << byteorder << "\n"
        << "data ignore value = " << grid.get_noData() << "\n"
        << "map info = " << grid.get_mapInfo().asString() << "\n"
        << "coordinate system string = " << grid.get_coordinateSys() << "\n"
        << "band names = {";
    for (int i = 0; i < grid.get_nbands(); i++){
        if ( i ){ file << ", "; }
        file << grid.get_bandName(i);
    }
    file << "}\n";
    
    // close header file
    file.close();
    
    // save data to file
    std::ofstream disk(filename.c_str(), std::ios::out | std::ios::binary);
    disk.write((char*)&grid[0], grid.get_volume() * sizeof(T));
    disk.close();
}

/////////////////////////////////////////////////////////////////////
// GRID CONVERSION //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// Write grid A to grid B, converting from the data type in A to B.
template <class T, class U>
void convert ( const Grid<T> & A, Grid<U> * B )
{
    *B = Grid<U>(A.ncols(), A.nrows(), A.nbands());
    B->set_interleave(A.get_interleave());
    B->set_noData((U) A.noData());
    B->setGeography(
        A.get_mapInfo().asString(), A.get_coordinateSys());
    
    for ( int i = 0; i < A.nbands(); ++i )
        B->set_bandName(i, A.get_bandName(i));
    
    for ( int i = 0; i < A.size(); ++i )
        B->setValue(i, (U) A[i]);
    
    B->calculateStatistics();
}

#endif // YOUNG_GRID_IO_20180528
