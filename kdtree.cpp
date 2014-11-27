/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <cmath> //for std::abs

template<int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim> & first, const Point<Dim> & second, int curDim) const
{
    //First, check if 'curDim' is an appropriate size
    if (curDim >= Dim || curDim >= Dim) {
	cout << "[smallerDimVal] ERROR: curDim is bigger than dimension" << endl;
	return false;
    }

    //compare the 'curDim' indexes of point 'first' and point 'second'
    if (first[curDim] < second[curDim])
	return true;
    else if (first[curDim] == second[curDim]) //break ties with op<
	return first < second;
    else
	return false;
}


template<int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim> & target, const Point<Dim> & currentBest, const Point<Dim> & potential) const
{
    /**
     * BASE CASE:
     * If we are comparing the same points, then just @return true
     */
    if (currentBest == potential)
	return true; 

    /**
     * Calculate the distance^2 between each point using the helper function 'distance'
     * We're using distance^2 to preserve precision that would be lost if we used 'sqrt'
     */
    double currDistance = distanceSquared(target, currentBest);
    double potDistance = distanceSquared(target, potential);

    //The comparisons
    if (potDistance < currDistance)
	return true;
    else if (potDistance == currDistance) //break ties with op<
	return potential < currentBest;
    else
	return false;
}

/**
 * Helper function to calculate the distance^2 between two points 'a' and 'b'
 * Keep the points 'const' so we don't change anything.
 * @param a - the first point
 * @param b - the second point
 * @return the distance squared between points 'a' and 'b'
 */
template<int Dim>
double KDTree<Dim>::distanceSquared(const Point<Dim> & a, const Point<Dim> & b) const {
    //the distance^2 we will return
    double d = 0;

    /**
     * Loop through each dimension of the points, finding their differences and squaring them
     * We effectively find the distance^2 with this loop
     */
    for (int i = 0; i < Dim; i++) {
	double diff = a[i] - b[i];

     	d += diff * diff;
    }

    return d;
}

/**
 * The constructor for KDTrees
 * We utilize a quick select algorithm using multiple helper functions to partition the vector
 * and place the medians in the correct indices.
 * @param newPoints - the vector of points we will construct the kdtree from
 */
template<int Dim>
KDTree<Dim>::KDTree(const vector< Point<Dim> > & newPoints)
{
    //check for empty 'newPoints' vector
    if (newPoints.size() == 0)
	return;

    //copy 'newPoints' into the tree's member vector so we can edit it
    for (size_t i = 0; i < newPoints.size(); i++)
	points.push_back(newPoints[i]);

    //variables that will be passed into the partition function
    int lo = 0;
    int hi = newPoints.size() - 1;
    int piv = (lo + hi) / 2;

    //call the recursive helper function 'construct' to begin construction of the tree
    construct(points, lo, hi, piv, 0);
}

/** 
 * The recursive helper function for use in the constructor
 * @param list - 	the list of points that we will be partitioning
 * @param low  - 	the lower index
 * @param high - 	the upper index
 * @param n    -	the index of the root (median)
 * @param dim  -	the current dimension we are working with 
 */
template<int Dim>
void KDTree<Dim>::construct(vector< Point<Dim> > & list, int low, int high, int n, int dim) {
    /**
     * BASE CASE:
     * If the lower index is greater than or equal to the higher index, we know that there are no more
     * children to check
     */
    if (low >= high)
	return;
	
    //call the 'select' helper function which will run the QuickSelect algorithm on 'list'
    select(list, low, high, n, dim);

    //recursive calls on the left and right sublists
    construct(list, low, n - 1, (low + n - 1) / 2, (dim + 1) % Dim);
    construct(list, n + 1, high, (n + 1 + high) / 2, (dim + 1) % Dim);
}

/**
 * the helper function 'select' which will recursively partition the 'list' using the QuickSelect algorithm.
 * @param list - 	the list of points that we will be partitioning
 * @param low  - 	the lower index
 * @param high - 	the upper index
 * @param n    -	the index of the root (median)
 * @param dim  -	the current dimension we are working with 
 */
template<int Dim>
void KDTree<Dim>::select(vector< Point<Dim> > & list, int low, int high, int n, int dim) {
    /**
     * BASE CASE:
     * If 'low' is equal to 'high' then we are done and have found the value we want
     */
    if (low == high)
	return;

    //partition the list around 'n'
    int pivotIndex = partition(list, low, high, n, dim);

    //the recursive calls where we determine if we have finished partitioning/selecting
    if (n == pivotIndex)
	return;
    else if (n < pivotIndex)
	select(list, low, pivotIndex - 1, n, dim);
    else
	select(list, pivotIndex + 1, high, n, dim);
}

/**
 * partition
 * This is the helper function that will partition within the given indices in 'list'
 * @param list - 	the list of points that we will be partitioning
 * @param low  - 	the lower index
 * @param high - 	the upper index
 * @param pivotIndex -  the index that we are partitioning around
 * @param dim  -	the current dimension we are working with 
 * @return the index that the element originally in 'pivotIndex' has ended up
 */
template<int Dim>
int KDTree<Dim>::partition(vector< Point<Dim> > & list, int low, int high, int pivotIndex, int dim) {
    //get the Point<Dim> @ pivotIndex
    Point<Dim> pivotValue = list[pivotIndex];
	
    //move pivotValue to the 'high' index, i.e. the end of the list
    std::swap(list[pivotIndex], list[high]);

    //declare 'storeIndex' which will be the next point to swap - initialized to 'low'
    int storeIndex = low;

    //partition the array
    for (int i = low; i < high; i++) {
	if (smallerDimVal(list[i], pivotValue, dim)) {
	    //swap the elements in 'i' and 'storeIndex'
	    std::swap(list[i], list[storeIndex]);
	    
	    //increment 'storeIndex' to the next element that will be swapped
	    storeIndex++;
 	}
    }

    //final swap to move the pivot to the spot it belongs in
    std::swap(list[high], list[storeIndex]);

    //return the index of the pivot
    return storeIndex;
}

/**
 * findNearestNeighbor
 * The function that will use helper functions to find the nearest neighbor to a given point 'query'
 * @param query - the point we want to find the closest one to
 * @return the closest point in the tree to 'query'
 */
template<int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim> & query) const
{
    //call the helper function that will return the index of the closest point to 'query'
    return points[nearestIndex(query, 0, points.size() - 1, 0)];
}

/**
 * nearestIndex
 * the helper function for use in 'findNearestNeighbor' to conduct the NNS algorithm
 * @param query - 	the point we want to find the closest point to
 * @param low   -	the lower index of the list of points we want to compare
 * @param high  - 	the upper index of the list of points we want to compare
 * @param dim   - 	the current splitting dimension we are working with
 * @return the index of the point closest to 'query'
 */
template<int Dim>
int KDTree<Dim>::nearestIndex(const Point<Dim> & query, int low, int high, int dim) const {
    /**
     * BASE CASE:
     * If 'low' is greater than or equal to 'high' we know we are done searching
     */
    if (low > high)
	return low;

    //index to return
    int currentBest;

    //the current median's index based on the 'low' and 'high' range
    int medianIndex = (low + high) / 2;

    /**
     * Here we decide which subtree to traverse based on a comparison of the 
     * values of the current splitting dimension of the point @ medianIndex and query
     */
    if (smallerDimVal(query, points[medianIndex], dim))
	currentBest = nearestIndex(query, low, medianIndex - 1, (dim + 1) % Dim);
    else
	currentBest = nearestIndex(query, medianIndex + 1, high, (dim + 1) % Dim);

    /**
     * Now we traverse back up the tree, comparing the currentBest point to its parents
     * We also decide whether we need to traverse the parent's subtrees
     */
    if (shouldReplace(query, points[currentBest], points[medianIndex]))
	    currentBest = medianIndex;

    int diff = std::abs(points[medianIndex][dim] - query[dim]);
    diff = diff * diff;

    if (diff <= distanceSquared(query, points[currentBest])) {
	//the index that will potentially be better than 'currentBest'
	int potentialBestIndex;

	if (currentBest >= medianIndex)
	    potentialBestIndex = nearestIndex(query, low, medianIndex - 1, (dim + 1) % Dim);
	else 
	    potentialBestIndex = nearestIndex(query, medianIndex + 1, high, (dim + 1) % Dim);

	//check if it is indeed better than 'currentBest' or not
	if (shouldReplace(query, points[currentBest], points[potentialBestIndex]))
	    currentBest = potentialBestIndex;
    }

    //return the best index
    return currentBest;
}
