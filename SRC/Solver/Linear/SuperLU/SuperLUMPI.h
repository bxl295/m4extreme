// SuperLU.h: interface for the SuperLU class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details

#if !defined(SOLVER_LINEAR_SUPERLU_MPI_INCLUDED)
#define SOLVER_LINEAR_SUPERLU_MPI_INCLUDED

#pragma once

#include <math.h> 
#include <stdlib.h> 
#include <vector>
#include <iterator>
#include <map>
#include <iostream>
#include "Solver/Linear/Linear.h"
#include "config.h"
#include <cassert>
#include "mpi.h"
#include "SRC/superlu_ddefs.h"
#include <time.h>

using namespace std;


namespace Solver {
  namespace Linear {
    template <class key> class SuperLU_MPI : public System <key> {
    public:
      typedef std::pair<key, key> keypair;

      //      Constructors/destructors:
      SuperLU_MPI();
      SuperLU_MPI(const map<key, double> &, const map<keypair, double> &, int_t nprow, int_t npcol, bool isFullMatrix, const int = 0, const int = 0);
      virtual ~SuperLU_MPI();

      //      General methods:
      void Solve();
      map<key, double> & GetTemperature();

      //      Accessors/mutators:

      double Get(key);
      double Get(key, key);
      void Set(key, double);
      void Set(key, key, double);

      //      General methods:

      void SetToZero1();
      void SetToZero2();
      void Add(key, double);
      void Add(key, key, double);
      double Norm();

    private:

      int_t m, n;
      int_t nnz;
      int many_rhs;

      map<key, unsigned int> KeyMap;
      map<keypair, unsigned int> KeyPairMap;
      map<key, int> indexmap;

      superlu_options_t options;
      SuperLUStat_t stat;
      SuperMatrix A;
      ScalePermstruct_t ScalePermstruct;
      LUstruct_t LUstruct;
      SOLVEstruct_t SOLVEstruct;
      gridinfo_t grid;
      double *berr;
      double *b;
      int_t nprow, npcol;
      int iam, ldb, ldx;
      int_t m_loc, fst_row;
      
      vector<double> LMatrix;
      vector<int_t> ASUB;
      vector<int_t> XA;
      map<key, double> Solution;
      
      

    private:

      SuperLU_MPI(SuperLU_MPI &);
      void operator=(SuperLU_MPI &);
    };

    //

    template <class key>
      SuperLU_MPI<key>::SuperLU_MPI(const map<key, double> &Keys, const map<keypair, double> &Keymaps, 
				    int_t nprow, int_t npcol, bool isFullMatrix, const int SymMode, const int MulRhs)
      : many_rhs(MulRhs) {

      int_t i, j;
      typename map<key, double>::const_iterator pK;
      typename map<keypair, double>::const_iterator pKK;
      double tolerance = 1.0e-14;

    
      int ncount = 0;
      for (pK = Keys.begin(); pK != Keys.end(); pK++, ++ncount ) {
    	indexmap.insert( make_pair(pK->first, ncount) );
      }

      n = (unsigned int) Keys.size();
      nnz = n*n;


      clock_t start,finish;
      double time1, time2, time3,time0,time4;
      start = clock();

      superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
      finish = clock();
      time0 = (double)(finish - start)/CLOCKS_PER_SEC;
      //cout << "1: " << time0 << endl;
      iam = grid.iam;
      assert (iam < nprow*npcol);

      i = 0;
      j = 0;

      double *nzval_loc;         /* local */
      int_t *colind, *rowptr;	 /* local */
      int_t nnz_loc;
      int_t m_loc_fst;
      int_t row, column;


      // Compute the number of rows to be distributed to local process
      m_loc = n / (grid.nprow * grid.npcol);
      m_loc_fst = m_loc;

      if ((m_loc * grid.nprow * grid.npcol) != n) {
      	if (iam == (grid.nprow * grid.npcol - 1)) 
	  {m_loc = n - m_loc * (grid.nprow * grid.npcol - 1);}
      }
    

      /********************************************
       * Read Matrix *
       ********************************************/    


      fst_row = iam * m_loc_fst;
    
      start = clock();

      int row_old;
      row_old = fst_row - 1;
      nnz_loc = 0;
      for (pKK = Keymaps.begin(); pKK != Keymaps.end(); pKK++) {
	keypair pP = pKK->first;
	row = indexmap.find(pP.first)->second;
	if (row >= fst_row && row < fst_row + m_loc) {
	  double value = pKK->second;
	  if (value != 0.0) {
	    LMatrix.push_back(value);
	    column = indexmap.find(pP.second)->second;
	    ASUB.push_back(column);
	    if (row == row_old + 1) {
	      XA.push_back(nnz_loc);
	      row_old = row;
	    }
	    nnz_loc++;
	  }	
	}
      } 		
      XA.push_back(nnz_loc);
      //cout << "nnz_loc: " << nnz_loc << endl;

      vector<double>::iterator p1;
      vector<int_t>::iterator p3, p4;

      p1 = LMatrix.begin();
      nzval_loc = &(*p1);
      p3 = ASUB.begin();
      colind = &(*p3);
      p4 = XA.begin();
      rowptr= &(*p4);

      finish = clock();
      time1 = (double)(finish - start)/CLOCKS_PER_SEC;
      //cout << "Read Matrixes: " << time1 << endl;

      start = clock();
      //Destroy_CompCol_Matrix_dist(&GA);

      dCreate_CompRowLoc_Matrix_dist(&A, n, n, nnz_loc, m_loc, fst_row, nzval_loc,
				     colind, rowptr, SLU_NR_loc, SLU_D, SLU_GE);


      finish = clock();
      time2 = (double)(finish - start)/CLOCKS_PER_SEC;
      //cout << "Construct A: " << time2 << endl;


      start = clock();


      if ( !((b) = doubleMalloc_dist(m_loc)))
      	ABORT("Malloc fails for rhs[]");

      for (pK = Keys.begin(), i = 0; i < fst_row; ++i) {
    	pK++;
      }

      for (i = 0; i < m_loc; ++i) {

      	(b)[i] = pK->second;
      	pK++;
      }
      ldb = m_loc;

      ldx = m_loc;

      if ( !(berr = doubleMalloc_dist(1)) )
      	ABORT("Malloc fails for berr[]");

      set_default_options_dist(&options);

      m = A.nrow;
      n = A.ncol;

      ScalePermstructInit(m, n, &ScalePermstruct);
      LUstructInit(n, &LUstruct);

      PStatInit(&stat);

      finish = clock();
      time3 = (double)(finish - start)/CLOCKS_PER_SEC;
      //cout << "Constructor Functions: " << time3  << endl;

    }

    //


    template <class key>
      SuperLU_MPI<key>::~SuperLU_MPI() {

      PStatFree(&stat);
      //Destroy_CompRowLoc_Matrix_dist(&A);
      ScalePermstructFree(&ScalePermstruct);
      Destroy_LU(n, &grid, &LUstruct);
      LUstructFree(&LUstruct);
      if (options.SolveInitialized) {
	dSolveFinalize(&options, &SOLVEstruct);
      }
      //SUPERLU_FREE(b);
      SUPERLU_FREE(berr);
    } 

    template <class key>
      void SuperLU_MPI<key>::Solve() {
      int info = 0;
   

      pdgssvx(&options, &A, &ScalePermstruct, b, ldb, 1, &grid,
	      &LUstruct, &SOLVEstruct, berr, &stat, &info);


      typename map<key, int>::const_iterator pK;
      int i;

      for (i=0, pK=indexmap.begin(); i < fst_row; ++i) {
	pK++;
      }

      for (i=0; i<m_loc; i++) {
	Solution.insert(make_pair(pK->first, *b));
	b++;
	pK++;
      }

      /*typename map<key, double>::const_iterator pS;
	pS = Solution.begin();
	for (i=0; i<5; ++i) {
	cout << "iam: " << pS->first << " " << pS->second << endl;
	pS++;
	} */ 
    }

    template <class key>
      map<key, double> & SuperLU_MPI<key>::GetTemperature() {
      return Solution;
    }

    template <class key>
      double SuperLU_MPI<key>::Get(key K) {
      assert(false);
    }

    template <class key>
      double SuperLU_MPI<key>::Get(key K1, key K2) {
      assert(false);
    }

    template <class key>
      void SuperLU_MPI<key>::Set(key K, double Input) {
      assert(false);
    }

    template <class key>
      void SuperLU_MPI<key>::Set(key K1, key K2, double Input) {
      assert(false);
    }

    template <class key>
      void SuperLU_MPI<key>::SetToZero1() {
      assert(false);
    }

    template <class key>
      void SuperLU_MPI<key>::SetToZero2() {
      assert(false);      
    }

    template <class key>
      void SuperLU_MPI<key>::Add(key K, double Input) {
      assert(false);
    }

    template <class key>
      double SuperLU_MPI<key>::Norm() {
      assert(false);
    }

    template <class key>
      void SuperLU_MPI<key>::Add(key K1, key K2, double Input) {
      assert(false);
    }
  }
}

#endif 
