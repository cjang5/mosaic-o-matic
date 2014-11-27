/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */
	 			
#include <iostream>
#include <map>
#include "maptiles.h"

using namespace std;

MosaicCanvas * mapTiles(SourceImage const & theSource, vector<TileImage> const & theTiles)
{
    //the pointer to a 'MosaicCanvas' we will return
    MosaicCanvas * mosaic = new MosaicCanvas(theSource.getRows(), theSource.getColumns());

    //a std::map of the TileImages in 'theTiles' and their average color values
    map< Point<3>, TileImage> tileImages;

    //a vector containing the average color values of the given TileImages
    vector< Point<3> > tileColors;
    
    /**
     * We use a for loop to iterate through the 'theTiles' vector and 
     * initialize the 'tileImages' map and 'tileColors' vector
     */
    for (size_t i = 0; i < theTiles.size(); i++) {
	RGBAPixel r = theTiles[i].getAverageColor();

	Point<3> p(r.red, r.green, r.blue);

	tileColors.push_back(p);
	tileImages[p] = theTiles[i];
    }

    //construct the kd-tree using 'tileColors'
    KDTree<3> regionColors(tileColors);

    /**
     * here we use a nested for loop to find the nearest neighbor of each region in the source image
     * and set the corresponding TileImage in the MosaicCanvas we are going to return
     */
    for (int i = 0; i < mosaic->getRows(); i++) {
	for (int j = 0; j < mosaic->getColumns(); j++) {
	    //the region color in theSource at coordinates (i, j)
	    RGBAPixel regionColor = theSource.getRegionColor(i, j);

	    //the Point<3> that will represent 'regionColor'
	    Point<3> p(regionColor.red, regionColor.green, regionColor.blue);

	    //the nearest neighbor to 'p'
	    Point<3> nearest = regionColors.findNearestNeighbor(p);

	    //the TileImage that nearest represents
 	    TileImage t = tileImages[nearest];

	    //add that TileImage to the MosaicCanvas
	    mosaic->setTile(i, j, t);
  	}
    }

    return mosaic;
}
