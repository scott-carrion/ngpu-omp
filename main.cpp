// main.cpp
//
// Skyview code parallelized using OpenMP for the O2 Lab
// Written by Scott Carrion
// DO NOT DISTRIBUTE WITHOUT PERMISSION!!

#include "GridIO.h"
#include "GridOps.h"

int main ()
{
    // Make Grid of floats out of input file, NangaSRTMv3.dat
    Grid<float> dem = loadGrid("NangaSRTMv3.dat");
 
    // Set output Grid using skyview function   
    Grid<float> out = skyview(dem, 1, 999999999);
    saveGrid(out, "skyview.dat");  // Save the grid to skyview.dat
    
    // Set output Grid using prominence function
    out = prominence(dem, 1, 999999999);
    saveGrid(out, "prominence.dat");  // Save the grid to prominence.dat
    
    return 0;
}
