// main.cpp
//
// Skyview code parallelized using OpenMP for the O2 Lab
// Written by Scott Carrion
// DO NOT DISTRIBUTE WITHOUT PERMISSION!!

#include "GridIO.h"
#include "GridOps.h"
#include <cstdio>
#include <omp.h>
#include <string>

// XXX temp XXX
#include <type_traits>
#include <iostream>
#include "Grid.h"

pick_k(const std::string& grid, const std::string& target, int k = 100, int start = 0)
{
	Grid<float> dem = loadGrid(grid);
	std::ofstream writer(target, std::ios::out);
	for (int i = start; i < start+k; ++i) {
		writer << dem[i] << std::endl;
	}

	writer.close();
}


// Debug main for solving trivially copyable Grid ADT problem...
int main(int argc, char** argv)
{
	std::cout << std::is_trivially_copyable<dummy_point>::value << " for specified type" << std::endl;
	return 0;
}

/*
int main (int argc, char** argv)
{
    omp_set_num_threads(16);
    
    std::string input = "NangaSRTMv3.dat";
    std::string skyview_output = "skyview.dat";
    std::string skyview_sample_output = "skyview_sample.txt";
    std::string prominence_output = "prominence.dat";
    std::string prominence_sample_output = "prominence_sample.txt";

    if (argc > 1) {
	std::string inputname = std::string(argv[1]);
	input = inputname + ".dat";
	skyview_output = "skyview_" + inputname + ".dat";
	skyview_sample_output = "skyview_sample_" + inputname + ".txt";
	prominence_output = "prominence_" + inputname + ".dat";
	prominence_sample_output = "prominence_sample_" + inputname + ".txt";
    }

    // Make Grid of floats out of input file, NangaSRTMv3.dat
    Grid<float> dem = loadGrid(input);
 
    // Set output Grid using skyview function   
    Grid<float> out = skyview(dem, 1, 999999999);
    saveGrid(out, skyview_output);  // Save the grid to skyview.dat
    
    // Set output Grid using prominence function
    out = prominence(dem, 1, 999999999);
    saveGrid(out, prominence_output);  // Save the grid to prominence.dat
    
   // Grab some samples of the uncompressed skyview and prominence grids for correctness analysis
   pick_k(skyview_output, skyview_sample_output);
   pick_k(prominence_output, prominence_sample_output);
   
   return 0;
}
*/
