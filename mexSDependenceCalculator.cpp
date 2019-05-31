/*
 * mexSDependenceCalculator.cpp
 *
 *      Created on: January 2010
 *      Author: Ronald A.J. van Elburg
 *      Email:  Ronald A J (at) vanelburg (dot) eu
 *
 * The calling syntax is:    mexSDependenceCalculator( arg1, arg2 )
 *
 * Usage of the correct libraries can be forced by starting matlab through:
 * LD_PRELOAD="/usr/lib/libstdc++.so.6 /lib/libgcc_s.so.1"  matlab
 *
 */

#include <iostream>
//#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <math.h>

#include "SDependenceCalculator.h"

// Include code for creating persistent object handles.
#include "ObjectHandle.h"

using namespace std;

extern void _main();

/****************************/
class SDCWrapper {

public:
    SDCWrapper();
    ~SDCWrapper();
    void setNMax(int N);
    void setS(double S);
    double * getLogMultiplicity(int _NTerminals);
    double * getLogHistoryCount(int _NTerminals);
    double * getProbability(int _NTerminals);
    double * getSummedElectronicPathlength(int _NTerminals);
    double * getSummedPartitionAsymmetries(int _NTerminals);
    int * getNLargeSubtree(int _NTerminals);
    int  getCount(int _N );
    int NMax;
    double S;
    static bool isCreatedSDCWrapper;
private:
   SDependenceCalculator * mySDC;
};



SDCWrapper::SDCWrapper() {
	//cout << "SDCWrapper::SDCWrapper"<< std::endl;
	mySDC=new SDependenceCalculator();
    cout << "SDCWrapper::SDCWrapper"<< std::endl;
}

SDCWrapper::~SDCWrapper() {

	cout << "SDCWrapper::~SDCWrapper"<< std::endl;
	delete mySDC;
}



void SDCWrapper::setNMax(int _NMax){   
    NMax=mySDC->updateNTerminals(_NMax);
};

void SDCWrapper::setS(double S){
   mySDC->updateSDependence(S);
};

int SDCWrapper::getCount(int  _NTerminals ){
   return mySDC->TreeSets[ _NTerminals-1].NoOfTrees;
};

double * SDCWrapper::getLogMultiplicity(int _NTerminals){
    return  mySDC->TreeSets[_NTerminals-1].multiplicityArray; 
};

double * SDCWrapper::getLogHistoryCount(int _NTerminals){
   return  mySDC->TreeSets[_NTerminals-1].historiesArray;  
};

double * SDCWrapper::getProbability(int _NTerminals){
  return mySDC->TreeSets[_NTerminals-1].SFactorArray;  
};

double *  SDCWrapper::getSummedElectronicPathlength(int _NTerminals){
    return   mySDC->TreeSets[_NTerminals-1].SEPArray; 
};

double *  SDCWrapper::getSummedPartitionAsymmetries(int _NTerminals){
    return   mySDC->TreeSets[_NTerminals-1].SPAArray;
};

bool SDCWrapper::isCreatedSDCWrapper=false;

static void createSDCWrapper( mxArray *plhs[]) {

	SDCWrapper  * mySDCWrapper = new SDCWrapper();      // Create a SDCWrapper object
	plhs[0]= create_handle(mySDCWrapper);              // Create persistent handle to this object

    return;
}


void mexFunction(
        int          nlhs,
        mxArray      *plhs[],
        int          nrhs,
        const mxArray *prhs[]
        ) {

	int  * vin1,* vin4, in1,in4;
	double*  vin2;
	short * vin3;
	double *voutdouble;
	int *voutint;
	char * cstr;
    SDCWrapper  * mySDCWrapper;
    int NoOfTrees;


    if (nrhs == 0) {
    	if(!SDCWrapper::isCreatedSDCWrapper){
    		plhs[0]=mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    		createSDCWrapper(plhs);
    		SDCWrapper::isCreatedSDCWrapper=true;
    	}else{

    		mexWarnMsgIdAndTxt("SDCWrapper:mexFunction", "Only a single creation of the SDCWrapper is allowed");
    	}
    }else if (nrhs >=2   && SDCWrapper::isCreatedSDCWrapper==true) {
            

         if(mxIsUint32(prhs[0])){

            // Get existing SDCWrapper
            mySDCWrapper= &get_object<SDCWrapper>(prhs[0]);
            
            // Get input argument describing function ...
            vin1= (int *) mxGetPr(prhs[1]);
            
            // ... and use it to pick the desired function
            switch (*vin1) {
				case 'A': // Get Tree Asymmetry index

					double * myA;
					vin4 = (int *) mxGetPr(prhs[2]);
					in4 = *vin4;

					// Get data from SDCWrapper
					NoOfTrees = mySDCWrapper->getCount(in4);
					myA = mySDCWrapper->getSummedPartitionAsymmetries(in4);

					plhs[0] = mxCreateNumericMatrix(1, mySDCWrapper->getCount(in4), mxDOUBLE_CLASS, mxREAL);
					voutdouble = (double *) mxGetPr(plhs[0]);

					for (int index = 0; index < NoOfTrees; index++) {

						voutdouble[index] = myA[index] / (in4 - 1);
					}

					break;



                case 'C': // Get index count
                	vin4= (int *) mxGetPr(prhs[2]);
                    in4=*vin4;
                	
                    plhs[0]=mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
                	voutint=(int *) mxGetPr(plhs[0]);
                	voutint[0] =mySDCWrapper->getCount(in4);
                    
                	
                    break;   
                    
                   
                case 'E': // Get Electronic Pathlengths
                	
                    double * myEPs;
                	vin4= (int *) mxGetPr(prhs[2]);
                    in4=*vin4;
                	
                    // Get data from SDCWrapper
                    NoOfTrees=mySDCWrapper->getCount(in4);
                    myEPs=mySDCWrapper->getSummedElectronicPathlength(in4);
                    
                    plhs[0]=mxCreateNumericMatrix(1, mySDCWrapper->getCount(in4), mxDOUBLE_CLASS, mxREAL);
                	voutdouble=(double *) mxGetPr(plhs[0]);
                    
                    for ( int index = 0; index < NoOfTrees; index++ ) {
                       
                        voutdouble[index] = myEPs[index];
                    }
                    
                    break; 
                    




                case 'H': // Get logHistoryCount
                   
                    double * myLogHistoryCount;
                	vin4= (int *) mxGetPr(prhs[2]);
                    in4=*vin4;
                	
                    // Get data from SDCWrapper
                    NoOfTrees=mySDCWrapper->getCount(in4);
                    myLogHistoryCount=mySDCWrapper->getLogHistoryCount(in4);
                    
                    // and copy it to matlab
                    plhs[0]=mxCreateNumericMatrix(1,NoOfTrees , mxDOUBLE_CLASS, mxREAL);
                	voutdouble=(double *) mxGetPr(plhs[0]);
                    
                    for ( int index = 0; index < NoOfTrees; index++ ) {
                       
                        voutdouble[index] = myLogHistoryCount[index];
                    }
                    break;
                    
                case 'M': // Get Multiplicities
                	
                    double * myMultiplicities;
                	vin4= (int *) mxGetPr(prhs[2]);
                    in4=*vin4;
                	
                    // Get data from SDCWrapper
                    NoOfTrees=mySDCWrapper->getCount(in4);
                    myMultiplicities=mySDCWrapper->getLogMultiplicity(in4);
                    
                    plhs[0]=mxCreateNumericMatrix(1, mySDCWrapper->getCount(in4), mxDOUBLE_CLASS, mxREAL);
                	voutdouble=(double *) mxGetPr(plhs[0]);
                    
                    for ( int index = 0; index < NoOfTrees; index++ ) {
                       
                        voutdouble[index] =  myMultiplicities[index];
                    }
                    
                    break;
                    
                case 'N': //Set NMax
                   
                    
                	vin4= (int *) mxGetPr(prhs[2]);
                    in4=*vin4;
                	
                    mySDCWrapper->setNMax(in4);
                	
                    
                    plhs[0]=mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
                	voutint=(int *) mxGetPr(plhs[0]);
                	voutint[0] =mySDCWrapper->NMax;
                   

                    break;
                    
                case 'P': // Get SFactors (Probabilities)
                	
                    double * myProbabilities;
                	vin4= (int *) mxGetPr(prhs[2]);
                    in4=*vin4;
                	
                    // Get data from SDCWrapper
                    NoOfTrees=mySDCWrapper->getCount(in4);
                    myProbabilities=mySDCWrapper->getProbability(in4);
                    
                    
                    //double * myMultiplicities;
                    myMultiplicities=mySDCWrapper->getLogMultiplicity(in4);
                    
                    
                    plhs[0]=mxCreateNumericMatrix(1, mySDCWrapper->getCount(in4), mxDOUBLE_CLASS, mxREAL);
                	voutdouble=(double *) mxGetPr(plhs[0]);
                    
                    for ( int index = 0; index < NoOfTrees; index++ ) {
                       
                        voutdouble[index] = myProbabilities[index]*pow(10,myMultiplicities[index]);
                    }
                    
                    break;
                    
                case 'S': // Set S
                	vin2= (double *) mxGetPr(prhs[2]);
                	mySDCWrapper->setS(*vin2);
                    break;
                    
                case 'D':// Destroy Current Wrapper
                    destroy_object<SDCWrapper>(prhs[0]);
                    SDCWrapper::isCreatedSDCWrapper=false;

                    break;

                    
                    
                    
                default:
                    mexWarnMsgIdAndTxt("SDCWrapper:mexFunction", "Unknown argument.");

            }
        }
     };

    return;
}

