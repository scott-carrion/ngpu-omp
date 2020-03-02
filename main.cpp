#include "GridIO.h"
#include "GridOps.h"

int main ()
{
    Grid<float> dem = loadGrid("NangaSRTMv3.dat");
    
    Grid<float> out = skyview(dem, 1, 999999999);
    saveGrid(out, "skyview.dat");
    
    out = prominence(dem, 1, 999999999);
    saveGrid(out, "prominence.dat");
    
    return 0;
}