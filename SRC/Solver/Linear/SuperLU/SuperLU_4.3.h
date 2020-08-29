// SuperLU.h: interface for the SuperLU class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details

#if !defined(SOLVER_LINEAR_SUPERLU_V4_3__INCLUDED_)
#define SOLVER_LINEAR_SUPERLU_V4_3__INCLUDED_

#pragma once

#include <math.h> 
#include <stdlib.h> 
#include <vector>

#include <iostream>

using namespace std;

#include "SRC/slu_ddefs.h"
#include "SRC/slu_util.h" 
#include "Solver/Linear/Linear.h"


namespace Solver {
  namespace Linear {

    template <class key>
      class SuperLU_V4 : public System <key> {
    public:

      //      Constructors/destructors:

      SuperLU_V4();
      SuperLU_V4(const set<key> &, const set<pair <key, key> > &, const int = 0, const int = 0);
      SuperLU_V4(const vector<key> &, const int = 0, const int = 0);
      virtual ~SuperLU_V4();

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
      void Solve();
      double Norm();

    private:

      unsigned int n;
      unsigned int nnz;

      map<key, unsigned int> KeyMap;
      map<pair<key, key>, unsigned int> KeyPairMap;

      SuperMatrix A;
      SuperMatrix AC;
      SuperMatrix L;
      SuperMatrix U;
      SuperMatrix B;
      int panel_size;
      int relax;
      int many_rhs;
      double *a;
      double *b;
      int *Row;
      int *Column;
      int *PermutationRow;
      int *PermutationColumn;
      int *etree;

      superlu_options_t options;
      SuperLUStat_t stat;

    private:

      SuperLU_V4(SuperLU_V4 &);
      void operator=(SuperLU_V4 &);
    };

    template <class key>
      SuperLU_V4<key>::SuperLU_V4() {
      a = 0;
      b = 0;
      Row = 0;
      Column = 0;
      PermutationRow = 0;
      PermutationColumn = 0;

      dCreate_CompCol_Matrix(&A, 0, 0, 0, a, Row, Column, SLU_NC, SLU_D, SLU_GE);
      dCreate_Dense_Matrix(&B, 0, 1, b, 0, SLU_DN, SLU_D, SLU_GE);
      L.Store = 0;
      U.Store = 0;
    }

    template <class key>
      SuperLU_V4<key>::SuperLU_V4(const set<key> &Keys,
			    const set<pair<key, key> > &KeyPairs, const int SymMode, const int MulRhs) : many_rhs(MulRhs) {
      unsigned int i, j, npp;
      key K1, K2, K3;
      pair <key, key> KK;
      typename std::set <key>::const_iterator iK;
      typename std::set <pair <key, key> >::const_iterator iKK;

      n = (unsigned int) Keys.size();
      npp = n + 1;
      nnz = (unsigned int) KeyPairs.size();

      if (n == 0) throw (0);

      a = (double *) malloc(nnz * sizeof (double));
      b = (double *) malloc(n * sizeof (double));
      Row = (int *) malloc(nnz * sizeof (int));
      Column = (int *) malloc(npp * sizeof (int));

      for (i = 0; i < nnz; i++) {
	a[i] = 0.0;
	Row[i] = 0;
      }

      for (i = 0; i < n; i++) b[i] = 0.0;

      for (i = 0; i < npp; i++) Column[i] = 0;

      for (iK = Keys.begin(), i = 0; iK != Keys.end(); iK++, i++) KeyMap[*iK] = i;

      j = 0;
      Column[0] = 0;
      K3 = *(Keys.begin());

      for (iKK = KeyPairs.begin(), i = 0; iKK != KeyPairs.end(); iKK++, i++) {
	KK = *iKK;
	K2 = KK.first;
	K1 = KK.second;
	Row[i] = KeyMap[K1];
	KeyPairMap[KK] = i;
	if (K2 != K3) {
	  K3 = K2;
	  j++;
	  Column[j] = i;
	}
      }

      Column[n] = nnz;

      PermutationRow = new int [n];
      PermutationColumn = new int [n];

      dCreate_CompCol_Matrix(&A, n, n, nnz, a, Row, Column, SLU_NC, SLU_D, SLU_GE);
      dCreate_Dense_Matrix(&B, n, 1, b, n, SLU_DN, SLU_D, SLU_GE);
      L.Store = 0;
      U.Store = 0;
      options.Fact = DOFACT;
      panel_size = sp_ienv(1);
      relax = sp_ienv(2);
      StatInit(&stat);
      etree = intMalloc(A.ncol);
      get_perm_c(2, &A, PermutationColumn);
      options.IterRefine = NOREFINE;
      options.Trans = NOTRANS;
      if (SymMode != 1) {
	options.SymmetricMode = NO;
	options.DiagPivotThresh = 1.0;
	options.ColPerm = COLAMD;
      } else {
	options.SymmetricMode = YES;
	options.DiagPivotThresh = 0.001;
	options.ColPerm = MMD_AT_PLUS_A;
      }
      sp_preorder(&options, &A, PermutationColumn, etree, &AC);
    }

    //
    //
    //
    template <class key>
      SuperLU_V4<key>::SuperLU_V4(const vector<key> &Keys,
				  const int SymMode, const int MulRhs) 
      : many_rhs(MulRhs) {
      unsigned int i, j, npp;
      key K1, K2, K3;

      n = (unsigned int) Keys.size();
      npp = n + 1;
      nnz = n*n;

      if (n == 0) throw (0);

      a = (double *) malloc(nnz * sizeof (double));
      b = (double *) malloc(n * sizeof (double));
      Row = (int *) malloc(nnz * sizeof (int));
      Column = (int *) malloc(npp * sizeof (int));

      for (i = 0; i < nnz; ++i) {
	a[i] = 0.0;
	Row[i] = 0;
      }

      for (i = 0; i < n; ++i) b[i] = 0.0;

      for (i = 0; i < npp; ++i) Column[i] = 0;

      for (i = 0; i < n; ++i) {
	KeyMap.insert( make_pair(Keys[i], i) );
      }

      int count = 0;
      for ( i = 0; i < n; ++i ) {
	Column[i] = count;
	K1 = Keys[i];
	for ( j = 0; j < n; ++j ) {
	  K2 = Keys[j];
	  Row[count] = j;
	  KeyPairMap[make_pair(K2, K1)] = count;
	  ++count;
	}
      }

      Column[n] = nnz;

      PermutationRow = new int [n];
      PermutationColumn = new int [n];

      dCreate_CompCol_Matrix(&A, n, n, nnz, a, Row, Column, SLU_NC, SLU_D, SLU_GE);
      dCreate_Dense_Matrix(&B, n, 1, b, n, SLU_DN, SLU_D, SLU_GE);
      L.Store = 0;
      U.Store = 0;

      options.Fact = DOFACT;
      panel_size = sp_ienv(1);
      relax = sp_ienv(2);
      StatInit(&stat);
      etree = intMalloc(A.ncol);
      get_perm_c(2, &A, PermutationColumn);
      options.IterRefine = NOREFINE;
      options.Trans = NOTRANS;
      if (SymMode != 1) {
	options.SymmetricMode = NO;
	options.DiagPivotThresh = 1.0;
	options.ColPerm = COLAMD;
      } else {
	options.SymmetricMode = YES;
	options.DiagPivotThresh = 0.001;
	options.ColPerm = MMD_AT_PLUS_A;
      }
      sp_preorder(&options, &A, PermutationColumn, etree, &AC);
    }

    //
    //
    //
    template <class key>
      SuperLU_V4<key>::~SuperLU_V4() {
      SUPERLU_FREE(etree);
      StatFree(&stat);
      SUPERLU_FREE(b);
      if (U.Store != 0) Destroy_CompCol_Matrix(&U);
      if (L.Store != 0) Destroy_SuperNode_Matrix(&L);
      Destroy_SuperMatrix_Store(&B);
      Destroy_CompCol_Matrix(&A);
      Destroy_CompCol_Permuted(&AC);

      delete [] PermutationColumn;
      delete [] PermutationRow;

      KeyPairMap.clear();
      KeyMap.clear();
    }

    template <class key>
      double SuperLU_V4<key>::Get(key K) {
      return b[KeyMap[K]];
    }

    template <class key>
      double SuperLU_V4<key>::Get(key K1, key K2) {
      return a[KeyPairMap[make_pair(K2, K1)]];
    }

    template <class key>
      void SuperLU_V4<key>::Set(key K, double Input) {
      b[KeyMap[K]] = Input;
    }

    template <class key>
      void SuperLU_V4<key>::Set(key K1, key K2, double Input) {
      a[KeyPairMap[make_pair(K2, K1)]] = Input;
    }

    template <class key>
      void SuperLU_V4<key>::SetToZero1() {
      unsigned int i;
      for (i = 0; i < n; i++) b[i] = 0.0;
    }

    template <class key>
      void SuperLU_V4<key>::SetToZero2() {
      unsigned int i;
      for (i = 0; i < nnz; i++) a[i] = 0.0;
    }

    template <class key>
      void SuperLU_V4<key>::Add(key K, double Input) {
      b[KeyMap[K]] += Input;
    }

    template <class key>
      double SuperLU_V4<key>::Norm() {
      unsigned int i;
      double N = 0.0;
      for (i = 0; i < n; i++) N += b[i] * b[i];
      return sqrt(N);
    }

    template <class key>
      void SuperLU_V4<key>::Add(key K1, key K2, double Input) {
      a[KeyPairMap[make_pair(K2, K1)]] += Input;
    }

    template <class key>
      void SuperLU_V4<key>::Solve() {
      int info = 0;
      int drop_tol = 0;
      int lwork = 0;
      trans_t trans = NOTRANS;

      if (U.Store == 0 || L.Store == 0) {// factorization, outputs: L,U 
	dgstrf(&options, &AC, relax, panel_size, etree, NULL,
	       lwork, PermutationColumn, PermutationRow, &L, &U, &stat, &info);
      }

      if (info == 0) // solve
	dgstrs(trans, &L, &U, PermutationColumn, PermutationRow, &B, &stat, &info);

      if (info != 0) {
	cout << "Solver::SuperLU_V4 error in factorization @" << endl;
	dPrint_Dense_Matrix("AC", &AC);
	throw (0);
      }

      if (many_rhs == 0) // want to solve Ax=b1, Ax=b2, ... ?
	{
	  Destroy_CompCol_Matrix(&U);
	  U.Store = 0; // no
	  Destroy_SuperNode_Matrix(&L);
	  L.Store = 0;

	  //next time assume the sparsity pattern would be similar to speed-up factorization
	  options.Fact = SamePattern;
	} else options.Fact = FACTORED; // yes
    }
  }
}
#endif // !defined(SOLVER_LINEAR_SUPERLU_V4_3__INCLUDED_)
