// created 07/24/2017 by Brennan Young
// last modified 07/24/2017 by Brennan Young

#ifndef YOUNG_MAPINFO_20170724
#define YOUNG_MAPINFO_20170724

#include <cstdlib>   // atoi, atof, exit
#include <iostream>  // std::cout
#include <sstream>   // std::stringstream
#include <string>    // std::string

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

class MapInfo {
public:
    // coordinate system
    std::string projectionName; // name of the projected coordinate system
    int         projectionZone; // UTM zone
    std::string zoneCode;       // UTM zone code (e.g., N or S)
    std::string datum;          // underlying datum
    std::string units;          // units of measurement
    
    // reference location (in another grid)
    int         tiePoint_x;     // reference x location in file coordinates
    int         tiePoint_y;     // reference y location in file coordinates
    
    // geographic location
    double      x0;             // UTM easting (top left corner)
    double      y0;             // UTM northing (top left corner)
    double      dx;             // cell size in x direction
    double      dy;             // cell size in y direction
    
    // constructors
    MapInfo();
    MapInfo(std::string);       // parse as component of ENVI header file
    MapInfo(const MapInfo &);
    
    // destructor
    ~MapInfo();
    
    // operators
    MapInfo & operator=(const MapInfo &);
    
    // operations
    void        clear();        // resets map info to defaults
    void        parseENVI(std::string); // parses info from ENVI header file
                                        // format string
    
    // I/O
    std::string asString() const;
};

// CONSTRUCTORS / DESTRUCTOR //////////////////////////////////////////////////

MapInfo::MapInfo(){ clear(); }

MapInfo::MapInfo( std::string s ){ parseENVI( s ); }

MapInfo::MapInfo( const MapInfo & mi ) {
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

MapInfo::~MapInfo(){}

// OPERATORS //////////////////////////////////////////////////////////////////

MapInfo & MapInfo::operator= ( const MapInfo & mi ) {
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

// OPERATIONS /////////////////////////////////////////////////////////////////
void MapInfo::clear(){
    projectionName  = "undefined";
    projectionZone  = 0;
    zoneCode        = "N";
    datum           = "undefined";
    units           = "undefined";
    tiePoint_x      = 0;
    tiePoint_y      = 0;
    x0              = 0.0;
    y0              = 0.0;
    dx              = 1.0;
    dy              = 1.0;
}

// Parses a string for map info, based on the format of the 'map info' property
// for ENVI header files.
void MapInfo::parseENVI ( std::string s ) {
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
    if ( ss.size() > i ) projectionName = std::string( ss[i++] );
    if ( ss.size() > i ) tiePoint_x     = atoi( ss[i++].c_str() );
    if ( ss.size() > i ) tiePoint_y     = atoi( ss[i++].c_str());
    if ( ss.size() > i ) x0             = atof( ss[i++].c_str() );
    if ( ss.size() > i ) y0             = atof( ss[i++].c_str() );
    if ( ss.size() > i ) dx             = atof( ss[i++].c_str() );
    if ( ss.size() > i ) dy             = atof( ss[i++].c_str() );
    if ( ss.size() > i ) projectionZone = atoi( ss[i++].c_str() );
    if ( ss.size() > i ) zoneCode       = std::string( ss[i++] );
    if ( ss.size() > i ) datum          = std::string( ss[i++] );
    if ( ss.size() > i ) units          = std::string( ss[i++] );
}

// Compiles a string representation of the map info object using the format of
// the 'map info' component of ENVI header files.
std::string MapInfo::asString() const {
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

#endif // YOUNG_MAPINFO_yyyymmdd