/*
 * SDependenceCalculator.h
 *
 *  	Created on: January 2010
 *      Author: Ronald A.J. van Elburg
 *      Email:   Ronald A J (at) vanelburg (dot) eu
 */

#ifndef SDEPENDENCECALCULATOR_H_
#define SDEPENDENCECALCULATOR_H_

struct TreeSet
{
	// Variables becoming available at TreeSet construction time
	int NoOfTrees,NoOfLeaves;
	double  *historiesArray,  *multiplicityArray, *SEPArray, *SPAArray;
	int   *child1IndexArray, *child1NTerminalArray, * child2IndexArray, * child2NTerminalArray;
	int  **cofArray;

	// Variables becoming available when calculating SPath Structure
	int *predecessorTreeSetIndexArray, *predecessorCount,**predecessorTreeIndicesArray, **branchedLeafesCOFArray;

	// Variables becoming available when calculating S Dependence Structure
	double  * SFactorArray;
};



class SDependenceCalculator
{
public:
	SDependenceCalculator();
	~SDependenceCalculator();
	SDependenceCalculator(int _NTerminals);
	int updateNTerminals(int _NTerminals);
	int initialize();
	int  updateSPaths();
	int getNMax(){
		return NMax;};
	void updateSDependence(double _SP);

	int  getIndex(int  _n1, int _ti1 ,int  _n2, int _ti2 );
	double logCombinatorialFactor(int _n,int _n1,int _n2);

	TreeSet * TreeSets;
	int ** indicesArray;

private:
	double RallPower;
	int NTerminals;
	int   NMax;
	int NMaxS;
	double S;
};

#endif /* SDEPENDENCECALCULATOR_H_ */
