// created 07/24/2017 by Brennan Young
// last modified 07/24/2017 by Brennan Young

#ifndef YOUNG_MAPINFO_tc
#define YOUNG_MAPINFO_tc

#include <cstdlib>   // atoi, atof, exit
#include <iostream>  // std::cout
#include <sstream>   // std::stringstream
#include <string>    // std::string
#include <cstring>  // strcpy()  XXX added by SCC

#ifndef YOUNG_DIVIDESTRINGBYDELIMETER_20170724
#define YOUNG_DIVIDESTRINGBYDELIMETER_20170724
#include <vector>    // std::vector
// Trim a string of its leading and trailing whitespace.
std::string trim ( const std::string & s ){
    int a = s.find_first_not_of(" \t\n");
    int b = s.find_last_not_of(" \t\n");
    if ( a == -1 && b == -1 ){ return ""; }
    return s.substr(a, b-a+1);
}
// Split a string into a vector at a delimeter.
std::vector<std::string> split ( std::string s, std::string delim ) {
    int dd   = 0; // section start
    int dend = s.length();
    int d    = s.find(delim);
    std::vector<std::string> result;
    
    while ( dd > -1 ){
        if ( 0 < d ){ result.push_back( s.substr( dd, d-dd ) ); }
        else { result.push_back( s.substr( dd, dend-dd ) ); }
        // next element
        if ( d < 0 ){ dd = d; }
        else {
            dd = d+1;
            d = s.find(delim, dd);
        }
    }
    return result;
}
std::vector<std::string> split ( std::string s, char delim ) {
    return split ( s, std::string( 1, delim ) );
}
#endif // YOUNG_DIVIDESTRINGBYDELIMETER_20170724

class MapInfo_tc {
public:
    // coordinate system
    char* projectionName; // name of the projected coordinate system
    int         projectionZone; // UTM zone
    char* zoneCode;       // UTM zone code (e.g., N or S)
    char* datum;          // underlying datum
    char* units;          // units of measurement
    
    // reference location (in another grid)
    int         tiePoint_x;     // reference x location in file coordinates
    int         tiePoint_y;     // reference y location in file coordinates
    
    // geographic location
    double      x0;             // UTM easting (top left corner)
    double      y0;             // UTM northing (top left corner)
    double      dx;             // cell size in x direction
    double      dy;             // cell size in y direction
    
    // constructors
    //MapInfo_tc();  // XXX not allowed for trivial
    //MapInfo_tc(std::string);       // parse as component of ENVI header file
    // MapInfo_tc(const MapInfo_tc &);  // XXX trivial copy constructor required!
    
    // destructor
    // ~MapInfo_tc();  // XXX Trivial destructor required!
    
    // operators
    // MapInfo_tc & operator=(const MapInfo_tc &);  // XXX trivial copy assignment(?) required(?)!
    
    // operations
    void        clear();        // resets map info to defaults
    void        parseENVI(std::string); // parses info from ENVI header file
                                        // format string
    
    // I/O
    std::string asString() const;
};

// CONSTRUCTORS / DESTRUCTOR //////////////////////////////////////////////////
/*
MapInfo_tc::MapInfo_tc(){ clear(); }
MapInfo_tc::MapInfo_tc( std::string s ){ parseENVI( s ); }
*/
/*
MapInfo_tc::MapInfo_tc( const MapInfo_tc & mi ) {
    projectionName  = std::string(mi.projectionName);
    projectionZone  = mi.projectionZone;
    zoneCode        = std::string(mi.zoneCode);
    datum           = std::string(mi.datum);
    units           = std::string(mi.units);
    tiePoint_x      = mi.tiePoint_x;
    tiePoint_y      = mi.tiePoint_y;
    x0              = mi.x0;
    y0              = mi.y0;
    dx              = mi.dx;
    dy              = mi.dy;
}

MapInfo_tc::~MapInfo_tc(){}

// OPERATORS //////////////////////////////////////////////////////////////////

MapInfo_tc & MapInfo_tc::operator= ( const MapInfo_tc & mi ) {
    projectionName  = std::string(mi.projectionName);
    projectionZone  = mi.projectionZone;
    zoneCode        = std::string(mi.zoneCode);
    datum           = std::string(mi.datum);
    units           = std::string(mi.units);
    tiePoint_x      = mi.tiePoint_x;
    tiePoint_y      = mi.tiePoint_y;
    x0              = mi.x0;
    y0              = mi.y0;
    dx              = mi.dx;
    dy              = mi.dy;
    return *this;
}
*/
// OPERATIONS /////////////////////////////////////////////////////////////////

void MapInfo_tc::clear(){
    //projectionName  = "undefined";
    projectionZone  = 0;
    //zoneCode        = "N";
    //datum           = "undefined";
    //units           = "undefined";
    
    // Using strcpy() instead...
    projectionName = new char[64]; strcpy(projectionName, "undefined");
    zoneCode = new char[64]; strcpy(zoneCode, "N");
    datum = new char[64]; strcpy(datum, "undefined");
    units = new char[64]; strcpy(units, "undefined");

    tiePoint_x      = 0;
    tiePoint_y      = 0;
    x0              = 0.0;
    y0              = 0.0;
    dx              = 1.0;
    dy              = 1.0;
}

// Parses a string for map info, based on the format of the 'map info' property
// for ENVI header files.
void MapInfo_tc::parseENVI ( std::string s ) {
    if ( !s.length() ){
        clear();
        return;
    }
    
    int i = 0;
    
    // split each element
    std::vector<std::string> ss = split(
        s.substr( s.find_first_not_of( '{' ), s.find_last_not_of( '}' ) ),
        ',' );
    
    // trim white space
    for ( ; i < ss.size(); ++i ){ ss[i] = trim( ss[i] ); }
    
    // assign values from input string
    i = 0;

    // XXX Using strcpy() here instead of std::string interface in order to retain copy triviality

    //if ( ss.size() > i ) projectionName = std::string( ss[i++] );
    if ( ss.size() > i) { projectionName = new char[std::string(ss[i++]).length() + 1]; strcpy(projectionName, std::string(ss[i-1]).c_str()); }
    if ( ss.size() > i ) tiePoint_x     = atoi( ss[i++].c_str() );
    if ( ss.size() > i ) tiePoint_y     = atoi( ss[i++].c_str());
    if ( ss.size() > i ) x0             = atof( ss[i++].c_str() );
    if ( ss.size() > i ) y0             = atof( ss[i++].c_str() );
    if ( ss.size() > i ) dx             = atof( ss[i++].c_str() );
    if ( ss.size() > i ) dy             = atof( ss[i++].c_str() );
    if ( ss.size() > i ) projectionZone = atoi( ss[i++].c_str() );

    if ( ss.size() > i) { zoneCode = new char[std::string(ss[i++]).length() + 1]; strcpy(projectionName, std::string(ss[i-1]).c_str()); }
    if ( ss.size() > i) { datum = new char[std::string(ss[i++]).length() + 1]; strcpy(datum, std::string(ss[i-1]).c_str()); }
    if ( ss.size() > i) { units = new char[std::string(ss[i++]).length() + 1]; strcpy(units, std::string(ss[i-1]).c_str()); }



    //if ( ss.size() > i ) zoneCode       = std::string( ss[i++] );
    //if ( ss.size() > i ) datum          = std::string( ss[i++] );
    //if ( ss.size() > i ) units          = std::string( ss[i++] );
}

// Compiles a string representation of the map info object using the format of
// the 'map info' component of ENVI header files.
std::string MapInfo_tc::asString() const {
    std::stringstream ss;
    ss << std::fixed << "{"
       << projectionName << ", "
       << tiePoint_x     << ", "
       << tiePoint_y     << ", "
       << x0             << ", "
       << y0             << ", "
       << dx             << ", "
       << dy             << ", "
       << projectionZone << ", "
       << zoneCode       << ", "
       << datum          << ", "
       << units          << "}";
    return ss.str();
}

#endif // YOUNG_MAPINFO_tc

