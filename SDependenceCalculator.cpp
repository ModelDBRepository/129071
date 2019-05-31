/*
 * SDependenceCalculator.cpp
 *
 * 	 	Created on: January 2010
 *      Author: Ronald A.J. van Elburg
 *      Email:   Ronald A J (at) vanelburg (dot) eu
 */

#include "SDependenceCalculator.h"
#include <math.h>
#include <iostream>



SDependenceCalculator::~SDependenceCalculator(){
	//std::cout << "... SDependenceCalculator::~SDependenceCalculator"  << std::endl;
	for (int n=NMax-1;n>=0;n--){

		if(n<NMaxS){
			for(int ti=0;ti<TreeSets[n].NoOfTrees;ti++){
				delete[] TreeSets[n].predecessorTreeIndicesArray[ti];
				delete[] TreeSets[n].branchedLeafesCOFArray[ti];

			}


			delete[] TreeSets[n].SFactorArray;
			delete[] TreeSets[n].predecessorCount;
			delete[] TreeSets[n].predecessorTreeSetIndexArray;

			delete[] TreeSets[n].predecessorTreeIndicesArray; //**
			delete[] TreeSets[n].branchedLeafesCOFArray; //**

		}

		for(int ti=0;ti<TreeSets[n].NoOfTrees;ti++){
			delete[] TreeSets[n].cofArray[ti];
		}
		delete[] TreeSets[n].cofArray; //**

		delete[] indicesArray[n];

		delete[] TreeSets[n].historiesArray;

		delete[] TreeSets[n].multiplicityArray;

		delete[] TreeSets[n].child1IndexArray;

		delete[] TreeSets[n].child1NTerminalArray;

		delete[] TreeSets[n].child2IndexArray;

		delete[] TreeSets[n].child2NTerminalArray;

		delete[] TreeSets[n].SEPArray;

		delete[] TreeSets[n].SPAArray;
	}

	delete[] TreeSets; //**
	delete[] indicesArray;//**
}


SDependenceCalculator::SDependenceCalculator(){
    initialize();
};

SDependenceCalculator::SDependenceCalculator(int _NTerminals){
    initialize();
    updateNTerminals( _NTerminals);

}

int SDependenceCalculator::initialize(){

    int n, n1, n2; 				//index enumerating the TreeSets corresponding to the number of terminals-1
    int ti, ti1, ti2;			//index enumerating the trees in a TreeSets corresponding to Harding enumeration-1
    RallPower=1.5;
    NMax=0;
    NMaxS=0;
    S=-1e99;


    indicesArray=new int *[1];
    indicesArray[0]=new int [1];
    indicesArray[0][0]=1;

    TreeSets =new TreeSet[ 1];

    // Construct TreeSet for N=1
    n=0;//
    TreeSets[n].NoOfTrees=1;
    TreeSets[n].NoOfLeaves=1;

    TreeSets[n].historiesArray=new double[1];
    TreeSets[n].historiesArray[0]=0;

    TreeSets[n].multiplicityArray=new double[1];
    TreeSets[n].multiplicityArray[0]=0;

    TreeSets[n].SEPArray=new double[1];
    TreeSets[n].SEPArray[0]=pow(1, -RallPower/2.0);

    TreeSets[n].SPAArray=new double[1];
    TreeSets[n].SPAArray[0]=0;


    TreeSets[n].child1IndexArray=new int[0];
    TreeSets[n].child2IndexArray=new int[0];
    TreeSets[n].child2NTerminalArray=new int[0];
    TreeSets[n].child1NTerminalArray=new int[0];

    TreeSets[n].cofArray=new int*[TreeSets[n].NoOfTrees];
    for(int i=0;i<TreeSets[n].NoOfTrees;i++)  {
        TreeSets[n].cofArray[i]=new int[TreeSets[n].NoOfLeaves];
    }

    TreeSets[n].cofArray[0][0]=1;

    NMaxS=1;
    NMax=1;

    // Construct TreeSet for N=2
    NMax=updateNTerminals(2);

    // Define elementary histories.

    TreeSets[0].predecessorTreeSetIndexArray=new int[1];
    TreeSets[0].predecessorTreeSetIndexArray[0]=-1;

    TreeSets[0].predecessorTreeIndicesArray=new int*[1];
    TreeSets[0].predecessorTreeIndicesArray[0]= new int[0];

    TreeSets[0].branchedLeafesCOFArray=new int*[1];
    TreeSets[0].branchedLeafesCOFArray[0]= new int[0];

    TreeSets[0].predecessorCount=new int[1];
    TreeSets[0].predecessorCount[0]=0;

    TreeSets[1].predecessorTreeSetIndexArray=new int[1];
    TreeSets[1].predecessorTreeSetIndexArray[0]=0;

    TreeSets[1].predecessorCount=new int[1];
    TreeSets[1].predecessorCount[0]=1;

    TreeSets[1].predecessorTreeIndicesArray=new int*[1];
    TreeSets[1].predecessorTreeIndicesArray[0]= new int[1];
    TreeSets[1].predecessorTreeIndicesArray[0][0]=0;

    TreeSets[1].branchedLeafesCOFArray=new int*[1];
    TreeSets[1].branchedLeafesCOFArray[0]= new int[1];
    TreeSets[1].branchedLeafesCOFArray[0][0]= 1;

    TreeSets[1].SFactorArray=new double[1];
    TreeSets[1].SFactorArray[0]=1;

    // Initialize S-dependence arrays
    TreeSets[0].SFactorArray=new double[1];
    TreeSets[0].SFactorArray[0]=1;

    TreeSets[1].SFactorArray=new double[1];
    TreeSets[1].SFactorArray[0]=1;

    NMaxS=2;

    return NMax;

}

inline double 	SDependenceCalculator::logCombinatorialFactor(int _n, int _n1, int _n2){
    // Calculate log( _n!/(_n1! _n2!))
    double result=0;

    for(int j=_n1+1;j<=_n;j++){result+=log10(j);}
    for(int j=1;j<=_n2;j++){result-=log10(j);}

    return result;
}



int	SDependenceCalculator::updateNTerminals(int _NTerminals){
    int n, n1, n2; 				//index enumerating the TreeSets corresponding to the number of terminals-1
    int ti, ti1, ti2;			//index enumerating the trees in a TreeSets corresponding to Harding enumeration-1
    double log2=log10(2);
    double nEP;
    double logCF;

    // If this size is already available do nothing
    if( _NTerminals <= NMax) return NMax;

    // (Re)Allocate arrays varying with _NTerminals
    int ** _indicesArray;
    TreeSet * _TreeSets;

    _indicesArray= new int *[_NTerminals];
    _TreeSets =new TreeSet[_NTerminals];

    for(n=0;n<NMax;n++){
        _indicesArray[n]=indicesArray[n];
        _TreeSets[n] =TreeSets[n];
    }

    delete[] TreeSets;
    TreeSets=_TreeSets;

    delete[] indicesArray;
    indicesArray=_indicesArray;


    // Construct TreeSet for N >=NMax
    for (n=NMax;n< _NTerminals;n++){
        indicesArray[n]=new int[n+1];

        // Calculate number of trees in this TreeSet
        TreeSets[n].NoOfTrees=0;
        for(n1=n-1 ; 2*n1 >= n-1;n1--){
            TreeSets[n].NoOfTrees+=TreeSets[n1].NoOfTrees*TreeSets[n-n1-1].NoOfTrees;
        }
        n1++; // undo last n1-- to get back into the last step to allow correction of this last step
        if(n1==n-n1-1){TreeSets[n].NoOfTrees-=((TreeSets[n1].NoOfTrees-1)*TreeSets[n1].NoOfTrees)/2;}
        // End of Calculation of the number of trees in this TreeSet



        TreeSets[n].NoOfLeaves=n+1;

        TreeSets[n].historiesArray		=new double[TreeSets[n].NoOfTrees];
        TreeSets[n].multiplicityArray	=new double[TreeSets[n].NoOfTrees];
        TreeSets[n].child1IndexArray	=new int[TreeSets[n].NoOfTrees];
        TreeSets[n].child1NTerminalArray=new int[TreeSets[n].NoOfTrees];
        TreeSets[n].child2IndexArray	=new int[TreeSets[n].NoOfTrees];
        TreeSets[n].child2NTerminalArray=new int[TreeSets[n].NoOfTrees];
        TreeSets[n].SEPArray			=new double[TreeSets[n].NoOfTrees];
        TreeSets[n].SPAArray			=new double[TreeSets[n].NoOfTrees];
        TreeSets[n].cofArray			=new int*[TreeSets[n].NoOfTrees];

        for(int i=0;i<TreeSets[n].NoOfTrees;i++){
            TreeSets[n].cofArray[i]=new int[TreeSets[n].NoOfLeaves];
        }


        ti=0;
        nEP=pow(TreeSets[n].NoOfLeaves, -1/(RallPower*2.0));
        for(n1=n-1 ; 2*n1 >= n-1;n1--){
            indicesArray[n][n1]=ti;
            n2=n-n1-1;
            logCF=logCombinatorialFactor(n-1, n1, n2);


            for(ti1=0;ti1<TreeSets[n1].NoOfTrees;ti1++){


                for(n1==n2?ti2=ti1:ti2=0;ti2<TreeSets[n2].NoOfTrees;ti2++){
                    TreeSets[n].child1IndexArray[ti]=ti1;
                    TreeSets[n].child1NTerminalArray[ti]=n1;
                    TreeSets[n].child2IndexArray[ti]=ti2;
                    TreeSets[n].child2NTerminalArray[ti]=n2;

                    TreeSets[n].multiplicityArray[ti]=TreeSets[n1].multiplicityArray[ti1]+TreeSets[n2].multiplicityArray[ti2];
                    if( ti1!=ti2 || n1 != n2 ){TreeSets[n].multiplicityArray[ti]+=log2;}

                    TreeSets[n].historiesArray[ti]=logCF +TreeSets[n1].historiesArray[ti1]+TreeSets[n2].historiesArray[ti2];
                    // The history of a tree with N (n+1) leaves contains N-1 branch events, if we have subtrees with N1 and N2 leaves
                    // the number of ways the N-2 branch events in these subtrees can be combined is  (N-2)!/(N1! n2!). As the branching
                    // of the root segment is always the first event to take place  the number of histories of a tree with N leaves can
                    // simply be calculated from the product of the number of histories of the constituting subtrees multiplied with
                    // this combinatorial factor.



                    TreeSets[n].SEPArray[ti]= nEP + TreeSets[n1].SEPArray[ti1]+	TreeSets[n2].SEPArray[ti2];

                    if(n>1){
                    	TreeSets[n].SPAArray[ti]= fabs(n1-n2)/(n1+n2) + TreeSets[n1].SPAArray[ti1]+	TreeSets[n2].SPAArray[ti2];
                    }else{
                    	TreeSets[n].SPAArray[ti]=0;
                    }

                    for(int j=0;j<=n1;j++){
                        TreeSets[n].cofArray[ti][j]=1+TreeSets[n1].cofArray[ti1][j];
                    }

                    for(int j=0;j<=n2;j++){
                        TreeSets[n].cofArray[ti][j+n1+1]=TreeSets[n2].cofArray[ti2][j]+1;
                    }

                   ti++;

                }
            }
        }

    } // end for n




    NMax=_NTerminals;

    return  NMax;
};

inline int  SDependenceCalculator::getIndex(int  _n1, int _ti1 , int  _n2, int _ti2 ){
    int _nswap, _tiswap, index, N2;


    if(_n2 > _n1 ){
        _nswap=_n1;
        _n1=_n2;
        _n2=_nswap;
        _tiswap=_ti1;
        _ti1=_ti2;
        _ti2=_tiswap;

    }else if(_n2 ==_n1){
        if(_ti1 >_ti2){
            _tiswap=_ti1;
            _ti1=_ti2;
            _ti2=_tiswap;
        }
    }

    //Calculate index
    index = indicesArray[_n1+_n2+1][_n1];
    N2=TreeSets[_n2].NoOfTrees;
    if(_n1 > _n2 ){
        index+=_ti1*N2	+_ti2;
    }else if(_n2 ==_n1){
        index+=(_ti1*(N2*2-_ti1+1))/2	+(_ti2-_ti1);
    }
    return index;

};

int SDependenceCalculator::updateSPaths(){
    int n, n1, n2;
    int ti, ti1, ti2;
    int pre, ipre, ipre2, preCount;

    if (NMaxS ==NMax) return -1;

    for (n=NMaxS;n<NMax;n++){
        //TreeSets[n].SFactorArray=new double[TreeSets[n].NoOfTrees];
        TreeSets[n].predecessorTreeSetIndexArray =new int[TreeSets[n].NoOfTrees];
        TreeSets[n].predecessorCount=new int[TreeSets[n].NoOfTrees];
        TreeSets[n].predecessorTreeIndicesArray=new int*[TreeSets[n].NoOfTrees];
        TreeSets[n].branchedLeafesCOFArray=new int*[TreeSets[n].NoOfTrees];
        TreeSets[n].SFactorArray=new double[TreeSets[n].NoOfTrees];

        for(ti=0;ti<TreeSets[n].NoOfTrees;ti++){

            TreeSets[n].predecessorTreeSetIndexArray[ti]=n-1;

            ti1=TreeSets[n].child1IndexArray[ti];

            n1=TreeSets[n].child1NTerminalArray[ti];
            ti2=TreeSets[n].child2IndexArray[ti];

            n2=TreeSets[n].child2NTerminalArray[ti];
            preCount=TreeSets[n1].predecessorCount[ti1]
                    +TreeSets[n2].predecessorCount[ti2];


            TreeSets[n].predecessorCount[ti]=preCount;

            TreeSets[n].predecessorTreeIndicesArray[ti]=new int[preCount];
            TreeSets[n].branchedLeafesCOFArray[ti]=new int[preCount];

            for(ipre=0;ipre<TreeSets[n1].predecessorCount[ti1];ipre++){
                TreeSets[n].predecessorTreeIndicesArray[ti][ipre]=getIndex(n1-1, TreeSets[n1].predecessorTreeIndicesArray[ti1][ipre] , n2, ti2 );
                TreeSets[n].branchedLeafesCOFArray[ti][ipre]=TreeSets[n1].branchedLeafesCOFArray[ti1][ipre]+1;
            }



            for(ipre2=0;ipre2<TreeSets[n2].predecessorCount[ti2];ipre2++){
                TreeSets[n].predecessorTreeIndicesArray[ti][ipre]=getIndex(n1, ti1 , n2-1, TreeSets[n2].predecessorTreeIndicesArray[ti2][ipre2] );
                TreeSets[n].branchedLeafesCOFArray[ti][ipre]=TreeSets[n2].branchedLeafesCOFArray[ti2][ipre2]+1;
                ipre ++;
            }

        }
    }
    NMaxS=NMax;

}

void 	SDependenceCalculator::updateSDependence(double _SP){
    int n, n1, n2; 				//index enumerating the TreeSets corresponding to the number of terminals-1
    int ti, tipre;			//index enumerating the trees in a TreeSets corresponding to Harding enumeration-1
    int ipre;
    double normalization, pathdependentfactor;

    updateSPaths();


    //std::cout << "_SP "<< _SP << std::endl;

    if(S!=_SP){
        for (n=2;n < NMax;n++){

            for(ti=0;ti<TreeSets[n].NoOfTrees;ti++){
                TreeSets[n].SFactorArray[ti]=0;
                //std::cout << "newtree n: " << n << " ti: "<< ti   << std::endl;
                for(ipre=0;ipre<TreeSets[n].predecessorCount[ti];ipre++){
                    tipre=TreeSets[n].predecessorTreeIndicesArray[ti][ipre];

                    normalization=0;
                    for(int terminal=0;terminal<n;terminal++){
                        normalization+=pow(2,- _SP*TreeSets[n-1].cofArray[tipre][terminal]);
                        //std::cout << "delta normalization" <<pow(2,- _SP*TreeSets[n-1].cofArray[tipre][terminal]) << std::endl;
                    }

                    //std::cout << " tipre "<< tipre << std::endl;

                    pathdependentfactor=pow(2,- _SP*TreeSets[n].branchedLeafesCOFArray[ti][ipre])/normalization;

                    ///std::cout << "pathdependentfactor \t " << pathdependentfactor<<" \t  TreeSets[n-1].SFactorArray[tipre] \t" <<TreeSets[n-1].SFactorArray[tipre] << std::endl;
                    TreeSets[n].SFactorArray[ti]+=pathdependentfactor*TreeSets[n-1].SFactorArray[tipre];

                }
			}
        }
    }

    S=_SP;
    return;
}

