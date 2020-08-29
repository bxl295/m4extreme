// Static.cpp: implementation of the Static class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
/////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <sys/time.h>
#include "./Static.h"
#include "Potential/TwoBody/TwoBody.h"

extern map<Set::Manifold::Point *, Set::VectorSpace::Vector> _gForce;
extern map<Set::Manifold::Point *, size_t> _gIndexSet;

namespace Model {
    namespace Static {
        //////////////////////////////////////////////////////////////////////
        // Class LocalState
        //////////////////////////////////////////////////////////////////////

        LocalState::LocalState() {
        }

        LocalState::LocalState(
                Clock *T_, set<Element::LocalState *> *ELS_,
                map<Set::Manifold::Point *, Set::Manifold::Map *> *Emb_,
                map<Set::Manifold::Point *, Set::Manifold::TMap *> *DEmb_)
        : T(T_), ELS(ELS_), Emb(Emb_), DEmb(DEmb_) {}

        LocalState::~LocalState() {
        }

        Model::LocalState *
        LocalState::Clone() const {
            return new LocalState(*this);
        }

        LocalState::LocalState(const LocalState &rhs)
        : T(rhs.T), ELS(rhs.ELS),
        Emb(rhs.Emb), DEmb(rhs.DEmb) {
        }

        LocalState &
                LocalState::operator =(const LocalState &rhs) {
            if (this == &rhs) return *this;
            T = rhs.T;
            ELS = rhs.ELS;
            Emb = rhs.Emb;
            DEmb = rhs.DEmb;
            return *this;
        }

        void
        LocalState::Embed(const set<Set::Manifold::Point *> &y,
                map <Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> &yemb) {
            typedef map <Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> yemb_type;
            yemb.clear();
            set<Set::Manifold::Point *>::const_iterator py;
            map<Set::Manifold::Point *, Set::Manifold::Map *>::iterator pEmb;
            for (py = y.begin(), pEmb = Emb->begin(); py != y.end(); py++, pEmb++) {
                if (pEmb->second == 0) {
                    yemb.insert(yemb_type::value_type(*py,
                            (Set::Euclidean::Orthonormal::Point &)(*(*py))));
                } else {
                    yemb.insert(yemb_type::value_type(*py,
                            (Set::Euclidean::Orthonormal::Point &)(*pEmb->second)((*(*py)))));
                }
            }
        }
        
        void
        LocalState::Embed(const set<Set::Manifold::Point *> &y,
                          map <Set::Manifold::Point *, LocalState::vector_type> &yemb) {
            yemb.clear();
            set<Set::Manifold::Point *>::const_iterator py;
            map<Set::Manifold::Point *, Set::Manifold::Map *>::iterator pEmb;

            for (py = y.begin(); py != y.end(); py++) {
                Set::Manifold::Point * yloc = *py;
                pEmb = Emb->find(yloc);

                if (pEmb != Emb->end()) {                    

                    if (pEmb->second == 0) {
                        const Set::Euclidean::Orthonormal::Point * ploc = 
                                dynamic_cast<Set::Euclidean::Orthonormal::Point*> (yloc);
                        yemb.insert( make_pair(yloc, vector_type(ploc->size(), ploc->begin())) );
                    } else {
                        const Set::Euclidean::Orthonormal::Point & ploc = 
                                dynamic_cast<const Set::Euclidean::Orthonormal::Point &> ((*pEmb->second)(*yloc));
                        yemb.insert( make_pair(yloc, vector_type(ploc.size(), ploc.begin())) );
                    }                    
                    
                } else {
                    assert(false);
                }
            }
        }

        map <Set::Manifold::Point *, LocalState::vector_type>
        LocalState::Embed(const set<Set::Manifold::Point *> &y) {

            map <Set::Manifold::Point *, vector_type> yemb;

            set<Set::Manifold::Point *>::const_iterator py;
            map<Set::Manifold::Point *, Set::Manifold::Map *>::iterator pEmb;
                        for (py = y.begin(); py != y.end(); py++) {
                Set::Manifold::Point * yloc = *py;
                pEmb = Emb->find(yloc);

                if (pEmb != Emb->end()) {                    

                    if (pEmb->second == 0) {
                        const Set::Euclidean::Orthonormal::Point * ploc = 
                                dynamic_cast<Set::Euclidean::Orthonormal::Point*> (yloc);
                        yemb.insert( make_pair(yloc, vector_type(ploc->size(), ploc->begin())) );
                    } else {
                        const Set::Euclidean::Orthonormal::Point & ploc = 
                                dynamic_cast<const Set::Euclidean::Orthonormal::Point &> ((*pEmb->second)(*yloc));
                        yemb.insert( make_pair(yloc, vector_type(ploc.size(), ploc.begin())) );
                    }                    
                    
                } else {
                    assert(false);
                }
            }
	    
            return yemb;
        }

        void
        LocalState::Embed(
                const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u,
                map<Set::Manifold::Point *, Set::VectorSpace::Vector> &uemb) {
            uemb.clear();
            map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pu;
            map<Set::Manifold::Point *, Set::Manifold::TMap *>::iterator pDEmb;
            for (pu = u.begin(), pDEmb = DEmb->begin(); pu != u.end(); pu++, pDEmb++) {
                if (pDEmb->second == 0) {
                    uemb.insert(make_pair(pu->first, pu->second));
                } else {
                    Set::VectorSpace::Hom A = (*pDEmb->second)((*(pu->first)));
                    uemb.insert(make_pair(pu->first, A(pu->second)));
                }
            }
        }

        void
        LocalState::Submerge(
                const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &femb,
                map<Set::Manifold::Point *, Set::VectorSpace::Vector> &f) {
            f.clear();
            map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pf;
            map<Set::Manifold::Point *, Set::Manifold::TMap *>::iterator pDEmb;
            for (pf = femb.begin(), pDEmb = DEmb->begin(); pf != femb.end(); pf++, pDEmb++) {
                if (pDEmb->second == 0) {
                    f.insert(make_pair(pf->first, pf->second));
                } else {
                    Set::VectorSpace::Hom A = Adjoint((*pDEmb->second)((*(pf->first))));
                    f.insert(make_pair(pf->first, A(pf->second)));
                }
            }
        }

        void
        LocalState::Submerge(
                const map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> &dfemb,
                map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> &df) {
            typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> stiff_type;

            map<Set::Manifold::Point *, Set::Manifold::TMap *> &Dg = *DEmb;
            df.clear();
            stiff_type::const_iterator pdf;
            for (pdf = dfemb.begin(); pdf != dfemb.end(); pdf++) {
                pair<Set::Manifold::Point *, Set::Manifold::Point *> ba = pdf->first;
                Set::Manifold::Point *b = ba.first;
                Set::Manifold::Point *a = ba.second;
                const Set::VectorSpace::Hom & Aemb = pdf->second;

                if ((Dg[a] == 0) && (Dg[b] == 0)) {
                    df.insert(stiff_type::value_type(ba, Aemb));
                } else if ((Dg[a] != 0) && (Dg[b] == 0)) {
                    const Set::Manifold::Point &ya = *a;
                    Set::Manifold::TMap &Dga = *Dg[a];
                    df.insert(stiff_type::value_type(ba, Aemb * Dga(ya)));
                } else if ((Dg[a] == 0) && (Dg[b] != 0)) {
                    const Set::Manifold::Point &yb = *b;
                    Set::Manifold::TMap &Dgb = *Dg[b];
                    df.insert(stiff_type::value_type(ba, Adjoint(Dgb(yb)) * Aemb));
                } else if ((Dg[a] != 0) && (Dg[b] != 0)) {
                    const Set::Manifold::Point &ya = *a;
                    Set::Manifold::TMap &Dga = *Dg[a];
                    const Set::Manifold::Point &yb = *b;
                    Set::Manifold::TMap &Dgb = *Dg[b];
                    df.insert(stiff_type::value_type(ba, Adjoint(Dgb(yb)) * Aemb * Dga(ya)));
                }
            }
        }

        void
        LocalState::operator ++() {
            set<Element::LocalState *>::iterator pELS;
            for (pELS = ELS->begin(); pELS != ELS->end(); pELS++) {
                if ( (*pELS)->isActivated() ) {
                ++(*(*pELS));
                }
            }
        }

        set<DOF>
        LocalState::GetDOFs() const {
            set<DOF> DOFs;
            unsigned int i;
            set<Element::LocalState *>::iterator pELS;
            for (pELS = ELS->begin(); pELS != ELS->end(); pELS++) {
                Element::LocalState &LSLoc = *(*pELS);

                if ( LSLoc.isActivated() ) {
                set<Set::Manifold::Point *> N = LSLoc.GetNodes();
                set<Set::Manifold::Point *>::iterator pN;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
                    for (i = 0; i < a->size(); i++) DOFs.insert(make_pair(a, i));
                }
                }
            }
            return DOFs;
        }

        set<pair<DOF, DOF> >
        LocalState::GetDOFPairs() const {
            set<pair<DOF, DOF> > DOFPairs;
            unsigned int i, j;
            set<Element::LocalState *>::iterator pELS;
            for (pELS = ELS->begin(); pELS != ELS->end(); pELS++) {
                Element::LocalState &LSLoc = *(*pELS);

                if ( LSLoc.isActivated() ) {
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NN = LSLoc.GetNodePairs();
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> >::iterator pNN;
                for (pNN = NN.begin(); pNN != NN.end(); pNN++) {
                    pair<Set::Manifold::Point *, Set::Manifold::Point *> ba = *pNN;
                    Set::Manifold::Point *b = ba.first;
                    Set::Manifold::Point *a = ba.second;
                    unsigned m = a->size();
                    unsigned n = b->size();
                    for (i = 0; i < m; i++) {
                        DOF ai = make_pair(a, i);
                        for (j = 0; j < n; j++) {
                            DOF bj = make_pair(b, j);
                            DOFPairs.insert(make_pair(ai, bj));
                        }
                    }
                }
                }
            }
            return DOFPairs;
        }

        void
        LocalState::operator () (
                Solver::Linear::System<DOF> &S,
                set<Set::Manifold::Point *> &y) const {
            unsigned int i;
            set<Set::Manifold::Point *>::iterator py;
            for (py = y.begin(); py != y.end(); py++) {
                Set::Manifold::Point *a = *py;
                Set::Manifold::Point &ya = *a;
                unsigned int n = ya.size();
                Set::VectorSpace::Vector u(n);
                for (i = 0; i < n; i++)
                    u[i] = S.Get(make_pair(a, i));
                ya -= u;
            }
        }

        void
        LocalState::operator () (
                Solver::Linear::System<DOF> &S,
                map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u) const {
            unsigned int i;
            map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pu;
            for (pu = u.begin(); pu != u.end(); pu++) {
                Set::Manifold::Point *a = pu->first;
                Set::VectorSpace::Vector &ua = pu->second;
                for (i = 0; i < ua.size(); i++) ua[i] = S.Get(make_pair(a, i));
            }
        }

        set<Element::LocalState *> *
        LocalState::GetELS() const {
            return ELS;
        }

        //////////////////////////////////////////////////////////////////////
        // Class Energy<0>
        //////////////////////////////////////////////////////////////////////

        Energy < 0 > ::Energy() {
        }

        Energy < 0 > ::~Energy() {
        }

        Model::Energy < 0 > *
        Energy < 0 > ::Clone() const {
            return new Energy < 0 > (*this);
        }

        Energy < 0 > ::Energy(LocalState *LS_,
                set<Element::Energy < 0 > *> *E_) : LS(LS_), E(E_) {
        }

        Energy < 0 > ::Energy(const Energy < 0 > &rhs) : LS(rhs.LS), E(rhs.E) {
        }

        Energy < 0 > &
                Energy < 0 > ::operator =(const Energy < 0 > &rhs) {
            if (this == &rhs) return *this;
            LS = rhs.LS;
            E = rhs.E;
            return *this;
        }

        double
        Energy < 0 > ::operator () (
                const set<Set::Manifold::Point *> &y) const {            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            double F = 0.0;
            vector_type yemb;
            LS->Embed(y, yemb);
            set<Element::Energy < 0 > *>::iterator pE;
            for (pE = E->begin(); pE != E->end(); pE++) {
                const Element::LocalState &LSLoc = *(*pE)->GetLocalState();
                
                if ( LSLoc.isActivated() ) {
                const Element::Energy < 0 > &ELoc = *(*pE);
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                vector_type yloc;
		vector_type::const_iterator py;
                set <Set::Manifold::Point *>::iterator pN;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (py = yemb.find(a)) != yemb.end() ) {
		      yloc.insert(vector_type::value_type(a, py->second));
		    }
		    else {
		      assert(false);
		    }
                }
                F += ELoc(yloc);
                }
            }
            return F;
        }

        //////////////////////////////////////////////////////////////////////
        // Class Energy<1>
        //////////////////////////////////////////////////////////////////////

        Energy < 1 > ::Energy() {
        }

        Energy < 1 > ::~Energy() {
        }

        Model::Energy < 1 > *
        Energy < 1 > ::Clone() const {
            return new Energy < 1 > (*this);
        }

        Energy < 1 > ::Energy(LocalState *LS_,
                set<Element::Energy < 1 > *> *DE_) : LS(LS_), DE(DE_) {}

        Energy < 1 > ::Energy(const Energy < 1 > &rhs) : LS(rhs.LS), DE(rhs.DE) {
        }

        Energy < 1 > &
                Energy < 1 > ::operator =(const Energy < 1 > &rhs) {
            if (this == &rhs) return *this;
            LS = rhs.LS;
            DE = rhs.DE;
            return *this;
        }

        map<Set::Manifold::Point *, Set::VectorSpace::Vector>
        Energy < 1 > ::operator () (const set<Set::Manifold::Point *> &y) const {            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> force_type;

            double tol = 1.0e-11;
            
            vector_type yemb;
            force_type DF;
            force_type femb;
            LS->Embed(y, yemb);
            set<Set::Manifold::Point *>::const_iterator py;
            for (py = y.begin(); py != y.end(); py++) {
                Set::Manifold::Point *a = *py;
                const Set::Manifold::Point &ya = *a;
                femb.insert(force_type::value_type(
                        a, Set::VectorSpace::Vector(yemb[a].size())));
                DF.insert(force_type::value_type(
                        a, Set::VectorSpace::Vector(ya.size())));
            }
            
            set<Element::Energy < 1 > *>::iterator pDE;
            for (pDE = DE->begin(); pDE != DE->end(); pDE++) {
                const Element::LocalState &LSLoc = *(*pDE)->GetLocalState();

                if ( LSLoc.isActivated() ) {
                const Element::Energy < 1 > &DELoc = *(*pDE);
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                vector_type yloc;
		vector_type::const_iterator py;
                set <Set::Manifold::Point *>::iterator pN;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (py=yemb.find(a)) != yemb.end() ) {
		      yloc.insert(vector_type::value_type(a, py->second));
		    }
		    else {
		      assert(false);
		    }
                }
		
                force_type floc = DELoc(yloc);

#ifdef _M4EXTREME_DEBUG_II    
                std::cout << "\nSingle Element forces:\n";
#endif
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
                    Set::VectorSpace::Vector & force = floc.find(a)->second;
                    
#ifdef _M4EXTREME_DEBUG_II                    
                    for ( int i=0; i<force.size(); i++) {
                        if ( fabs(force[i]) < tol ){
                            force[i] = 0.0;
                        }
                    }
                    
                    std::cout << _gIndexSet.find(*pN)->second << " [" << yemb.find(*pN)->second << "]------>["
                            << force << "]" << std::endl;
#endif
                    
                    
                    //femb.find(a)->second += floc.find(a)->second;
                    femb.find(a)->second += force;
                }
                }
            }            

#ifdef _M4EXTREME_DEBUG_II
            std::cout << "\nAssembled Nodal forces:\n";
            force_type::iterator pf;
            for (pf = femb.begin(); pf != femb.end(); pf++) {
                std::cout << _gIndexSet.find(pf->first)->second << " ["
                        << yemb.find(pf->first)->second << "]------>["
                        << pf->second << "]" << std::endl;
            }

            std::cout << std::endl;
#endif           

            LS->Submerge(femb, DF);
            return DF;
        }

        
        void Energy < 1 > ::operator () (const set<Set::Manifold::Point *> & y,
                map<Set::Manifold::Point *, Set::VectorSpace::Vector> & DF) {
            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> force_type;

            if (!DF.empty()) {
                DF.clear();
            }

            vector_type yemb;
            force_type femb;
            LS->Embed(y, yemb);

            set<Set::Manifold::Point *>::const_iterator py;
            for (py = y.begin(); py != y.end(); py++) {
                femb.insert(force_type::value_type(*py,
                        Set::VectorSpace::Vector(yemb.find(*py)->second.size())));
            }

            set<Element::Energy < 1 > *>::iterator pDE;
	    struct timeval start, stop;
	    int cost = 0;
	    _computational_cost.clear();

            int k = 0;

            for (pDE = DE->begin(); pDE != DE->end(); pDE++) {
                
                ++k;
                cost = 0;
                Element::LocalState * pLSLoc = (*pDE)->GetLocalState();

                if ( pLSLoc->isActivated() ) {
		  gettimeofday(&start, NULL);

		  set <Set::Manifold::Point *> N;
		  pLSLoc->GetNodes(N);

		  vector_type yloc;
		  vector_type::const_iterator py;
		  set <Set::Manifold::Point *>::iterator pN;
		  for (pN = N.begin(); pN != N.end(); pN++) {
		    if ( (py=yemb.find(*pN)) != yemb.end() ) {
		      yloc.insert(vector_type::value_type(*pN, py->second));
		    }
		    else {
		      assert(false);
		    }
		  }

		  int numofNodes = N.size();

		  if ( numofNodes == 2 ) {
		    //cout << "two body potential element" << endl;
		    Potential::TwoBody::LocalState * pTBLS = 
		      dynamic_cast<Potential::TwoBody::LocalState * >(pLSLoc);
		    if ( pTBLS == NULL ) {
		      assert(false);
		    }
		  }

		  force_type floc;
		  (*pDE)->operator()(yloc, floc);

#ifdef _M4EXTREME_DEBUG_II              
		  std::cout << "\nSingle Element forces:\n";
#endif

		  for (pN = N.begin(); pN != N.end(); pN++) {

#ifdef _M4EXTREME_DEBUG_II                    
                    std::cout << _gIndexSet.find(*pN)->second << " [" << yemb.find(*pN)->second << "]------>["
			      << floc.find(*pN)->second << "]" << std::endl;
#endif                    

                    femb.find(*pN)->second += floc.find(*pN)->second;
		  }

		gettimeofday(&stop, NULL);
		cost = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec; 
		}
		_computational_cost.insert( make_pair(pLSLoc, cost) );
            }
            

#ifdef _M4EXTREME_DEBUG_II
            std::cout << "\nAssembled Nodal forces:\n";
            force_type::iterator pf;
            for (pf = femb.begin(); pf != femb.end(); pf++) {
                std::cout << _gIndexSet.find(pf->first)->second << " ["
                        << yemb.find(pf->first)->second << "]------>["
                        << pf->second << "]" << std::endl;
            }

            std::cout << std::endl;
#endif

            LS->Submerge(femb, DF);
            return;
        }

        double
        Energy < 1 > ::operator () (
                const set<Set::Manifold::Point *> &y,
                const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u) const {            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> force_type;

            double DF = 0.0;
            vector_type yemb;
            LS->Embed(y, yemb);
            vector_type uemb;
            LS->Embed(u, uemb);
            set<Element::Energy < 1 > *>::iterator pDE;
            for (pDE = DE->begin(); pDE != DE->end(); pDE++) {
                const Element::LocalState &LSLoc = *(*pDE)->GetLocalState();

                if ( LSLoc.isActivated() ) {
                const Element::Energy < 1 > &DELoc = *(*pDE);
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                vector_type yloc;
                vector_type uloc;
		vector_type::const_iterator py, pu;
                set <Set::Manifold::Point *>::iterator pN;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (py=yemb.find(a)) != yemb.end() &&
			 (pu=uemb.find(a)) != uemb.end() ) {
		      yloc.insert(vector_type::value_type(a, py->second));
		      uloc.insert(vector_type::value_type(a, pu->second));
		    }
		    else {
		      assert(false);
		    }
                }
                force_type floc = DELoc(yloc);
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (pu=uemb.find(a)) != uemb.end() ) { 
		      DF += floc[a](pu->second);
		    }
                }
                }
            }
            return DF;
        }

        void
        Energy < 1 > ::operator () (
                const set<Set::Manifold::Point *> &y,
                Solver::Linear::System<DOF> &S) const {
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> force_type;
            unsigned int i;
            S.SetToZero1();
            force_type DF = this->operator ()(y);
            set<Set::Manifold::Point *>::const_iterator py;
            for (py = y.begin(); py != y.end(); py++) {
                Set::Manifold::Point *a = *py;
                Set::VectorSpace::Vector & DFa = DF[a];
                for (i = 0; i < DFa.size(); i++)
                    S.Add(make_pair(a, i), DFa[i]);
            }
        }

        //////////////////////////////////////////////////////////////////////
        // Class Energy<2>
        //////////////////////////////////////////////////////////////////////

        Energy < 2 > ::Energy() {
        }

        Energy < 2 > ::~Energy() {
        }

        Model::Energy < 2 > *
        Energy < 2 > ::Clone() const {
            return new Energy < 2 > (*this);
        }

        Energy < 2 > ::Energy(LocalState *LS_,
                set<Element::Energy < 2 > *> *DDE_) :
        LS(LS_), DDE(DDE_) {
        }

        Energy < 2 > ::Energy(const Energy < 2 > &rhs) :
        LS(rhs.LS), DDE(rhs.DDE) {
        }

        Energy < 2 > &
                Energy < 2 > ::operator =(const Energy < 2 > &rhs) {
            if (this == &rhs) return *this;
            LS = rhs.LS;
            DDE = rhs.DDE;
            return *this;
        }

        map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom>
        Energy < 2 > ::operator () (const set<Set::Manifold::Point *> &y0) const {            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>,
                    Set::VectorSpace::Hom> stiff_type;

            set<Set::Manifold::Point *> y = y0;
            stiff_type DDF;
            stiff_type dfemb;
            vector_type yemb;
            LS->Embed(y, yemb);
            set<Element::Energy < 2 > *>::iterator pDDE;
            for (pDDE = DDE->begin(); pDDE != DDE->end(); pDDE++) {
                const Element::LocalState &LSLoc = *(*pDDE)->GetLocalState();

                if ( LSLoc.isActivated() ) {
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NN = LSLoc.GetNodePairs();
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> >::iterator pNN;
                for (pNN = NN.begin(); pNN != NN.end(); pNN++) {
                    const pair<Set::Manifold::Point *, Set::Manifold::Point *> & ba = *pNN;
                    Set::Manifold::Point *b = ba.first;
                    Set::Manifold::Point *a = ba.second;
                    dfemb.insert(stiff_type::value_type(
                            ba, Set::VectorSpace::Hom(yemb[b].size(), yemb[a].size())));
                    DDF.insert(stiff_type::value_type(
                            ba, Set::VectorSpace::Hom(b->size(), a->size())));
                }
                }
            }

            for (pDDE = DDE->begin(); pDDE != DDE->end(); pDDE++) {
                const Element::LocalState &LSLoc = *(*pDDE)->GetLocalState();

                if ( LSLoc.isActivated() ) {
                const Element::Energy < 2 > &DDELoc = *(*pDDE);
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                vector_type yloc;
		vector_type::const_iterator py;
                set <Set::Manifold::Point *>::iterator pN;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (py=yemb.find(a)) != yemb.end() ) {
		      yloc.insert(vector_type::value_type(a, py->second));
		    }
		    else {
		      assert(false);
		    }
                }
                stiff_type dfloc = DDELoc(yloc);
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NN = LSLoc.GetNodePairs();
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> >::iterator pNN;
                for (pNN = NN.begin(); pNN != NN.end(); pNN++) {
                    const pair<Set::Manifold::Point *, Set::Manifold::Point *> & ba = *pNN;
                    dfemb[ba] += dfloc[ba];
                }
                }
            }
            LS->Submerge(dfemb, DDF);
            return DDF;
        }

        double
        Energy < 2 > ::operator () (const domain_type &y0,
                const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u) const {            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>,
                    Set::VectorSpace::Hom> stiff_type;

            double DDF = 0.0;
            set<Set::Manifold::Point *> y = y0;
            vector_type yemb;
            LS->Embed(y, yemb);
            vector_type uemb;
            LS->Embed(u, uemb);
            set<Element::Energy < 2 > *>::iterator pDDE;
            for (pDDE = DDE->begin(); pDDE != DDE->end(); pDDE++) {
                const Element::LocalState &LSLoc = *(*pDDE)->GetLocalState();
                if ( LSLoc.isActivated() ) {
                const Element::Energy < 2 > &DDELoc = *(*pDDE);
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                vector_type yloc;
                vector_type uloc;
		vector_type::const_iterator py, pu;
                set <Set::Manifold::Point *>::iterator pN;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (py=yemb.find(a)) != yemb.end() &&
			 (pu=uemb.find(a)) != uemb.end() ) {
		      yloc.insert(vector_type::value_type(a, py->second));
		      uloc.insert(vector_type::value_type(a, pu->second));
		    }
		    else {
		      assert(false);
		    }
                }
		
                stiff_type dfloc = DDELoc(yloc);
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NN = LSLoc.GetNodePairs();
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> >::iterator pNN;
                for (pNN = NN.begin(); pNN != NN.end(); pNN++) {
                    const pair<Set::Manifold::Point *, Set::Manifold::Point *> & ba = *pNN;
                    Set::Manifold::Point *b = ba.first;
                    Set::Manifold::Point *a = ba.second;
                    DDF += dfloc[ba](uloc[b])(uloc[a]);
                }
                }
            }
            return DDF;
        }

        void
        Energy < 2 > ::operator () (
                const set<Set::Manifold::Point *> &y,
                Solver::Linear::System<DOF> &S) const {
            typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> stiff_type;
            unsigned int i, j;
            S.SetToZero2();
            stiff_type DDF = this->operator ()(y);
            stiff_type::iterator pDDF;
            for (pDDF = DDF.begin(); pDDF != DDF.end(); pDDF++) {
                const pair<Set::Manifold::Point *, Set::Manifold::Point *> & ba = pDDF->first;
                Set::Manifold::Point *a = ba.first;
                Set::Manifold::Point *b = ba.second;
                const Set::VectorSpace::Hom & DDFba = pDDF->second;
                for (i = 0; i < DDFba.size1(); i++) {
                    DOF ai = make_pair(a, i);
                    for (j = 0; j < DDFba.size2(); j++) {
                        DOF bj = make_pair(b, j);
                        S.Add(ai, bj, DDFba[j][i]);
                    }
                }
            }
        }

        //////////////////////////////////////////////////////////////////////
        // Class Jet<0>
        //////////////////////////////////////////////////////////////////////

        Jet < 0 > ::Jet() {
        }

        Jet < 0 > ::~Jet() {
        }

        Model::Jet < 0 > *
        Jet < 0 > ::Clone() const {
            return new Jet < 0 > (*this);
        }

        Jet < 0 > ::Jet(LocalState *LS_,
                set<Element::Jet < 0 > *> *J_) : LS(LS_), J(J_) {
        }

        Jet < 0 > ::Jet(const Jet < 0 > &rhs) : LS(rhs.LS), J(rhs.J) {
        }

        Jet < 0 > &
                Jet < 0 > ::operator =(const Jet < 0 > &rhs) {
            if (this == &rhs) return *this;
            LS = rhs.LS;
            J = rhs.J;
            return *this;
        }

        pair<double, map<Set::Manifold::Point *, Set::VectorSpace::Vector> >
        Jet < 0 > ::operator () (const set<Set::Manifold::Point *> &y) const {            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> force_type;

            double F = 0.0;
            force_type femb;
            force_type DF;
            vector_type yemb;
            LS->Embed(y, yemb);
            set<Set::Manifold::Point *>::const_iterator py;
            for (py = y.begin(); py != y.end(); py++) {
                Set::Manifold::Point *a = *py;
                const Set::Manifold::Point &ya = *a;
                femb.insert(force_type::value_type(
                        a, Set::VectorSpace::Vector(yemb[a].size())));
                DF.insert(force_type::value_type(
                        a, Set::VectorSpace::Vector(ya.size())));
            }
            set<Element::Jet < 0 > *>::iterator pJ;
            for (pJ = J->begin(); pJ != J->end(); pJ++) {
                const Element::LocalState &LSLoc = *(*pJ)->GetLocalState();
                const Element::Jet < 0 > &JLoc = *(*pJ);
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                vector_type yloc;
		vector_type::const_iterator py;
                set <Set::Manifold::Point *>::iterator pN;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (py=yemb.find(a)) != yemb.end() ) {
		      yloc.insert(vector_type::value_type(a, py->second));
		    }
		    else {
		      assert(false);
		    }
                }
                pair<double, force_type> K = JLoc(yloc);
                double wloc = K.first;
                const force_type & floc = K.second;
                F += wloc;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
                    femb[a] += floc.find(a)->second;
                }
            }
            LS->Submerge(femb, DF);
            return make_pair(F, DF);
        }

        pair<double, double>
        Jet < 0 > ::operator () (
                const set<Set::Manifold::Point *> &y,
                const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u) const {            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> force_type;

            double F = 0.0;
            double DF = 0.0;
            vector_type yemb;
            LS->Embed(y, yemb);
            vector_type uemb;
            LS->Embed(u, uemb);
            set<Element::Jet < 0 > *>::iterator pJ;
            for (pJ = J->begin(); pJ != J->end(); pJ++) {
                const Element::LocalState &LSLoc = *(*pJ)->GetLocalState();
                const Element::Jet < 0 > &JLoc = *(*pJ);
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                vector_type yloc;
                vector_type uloc;
		vector_type::const_iterator py, pu;		
                set <Set::Manifold::Point *>::iterator pN;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (py=yemb.find(a)) != yemb.end() &&
			 (pu=uemb.find(a)) != uemb.end() ) {
		      yloc.insert(vector_type::value_type(a, py->second));
		      uloc.insert(vector_type::value_type(a, pu->second));
		    }
		    else {
		      assert(false);
		    }
                }
                pair<double, force_type> K = JLoc(yloc);
                double wloc = K.first;
                const force_type & floc = K.second;
                F += wloc;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
                    DF += (floc.find(a)->second)(uemb[a]);
                }
            }
            return make_pair(F, DF);
        }

        //////////////////////////////////////////////////////////////////////
        // Class Jet<1>
        //////////////////////////////////////////////////////////////////////

        Jet < 1 > ::Jet() {
        }

        Jet < 1 > ::~Jet() {
        }

        Model::Jet < 1 > *
        Jet < 1 > ::Clone() const {
            return new Jet < 1 > (*this);
        }

        Jet < 1 > ::Jet(LocalState *LS_,
                set<Element::Jet < 1 > *> *DJ_) : LS(LS_), DJ(DJ_) {
        }

        Jet < 1 > ::Jet(const Jet < 1 > &rhs) : LS(rhs.LS), DJ(rhs.DJ) {
        }

        Jet < 1 > &
                Jet < 1 > ::operator =(const Jet < 1 > &rhs) {
            if (this == &rhs) return *this;
            LS = rhs.LS;
            DJ = rhs.DJ;
            return *this;
        }

        pair<map<Set::Manifold::Point *, Set::VectorSpace::Vector>,
        map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> >
        Jet < 1 > ::operator () (const set<Set::Manifold::Point *> &y0) const {            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> force_type;
            typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>,
                    Set::VectorSpace::Hom> stiff_type;

            set<Set::Manifold::Point *> y = y0;
            vector_type yemb;
            LS->Embed(y, yemb);
            force_type DF;
            force_type femb;
            stiff_type DDF;
            stiff_type dfemb;

            set<Set::Manifold::Point *>::const_iterator py;
            for (py = y.begin(); py != y.end(); py++) {
                Set::Manifold::Point *a = *py;
                const Set::Manifold::Point &ya = *a;
                femb.insert(force_type::value_type(
                        a, Set::VectorSpace::Vector(yemb[a].size())));
                DF.insert(force_type::value_type(
                        a, Set::VectorSpace::Vector(ya.size())));
            }
            set<Element::Jet < 1 > *>::iterator pDJ;
            for (pDJ = DJ->begin(); pDJ != DJ->end(); pDJ++) {
                const Element::LocalState &LSLoc = *(*pDJ)->GetLocalState();
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NN = LSLoc.GetNodePairs();
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> >::iterator pNN;
                for (pNN = NN.begin(); pNN != NN.end(); pNN++) {
                    const pair<Set::Manifold::Point *, Set::Manifold::Point *> & ba = *pNN;
                    Set::Manifold::Point *b = ba.first;
                    Set::Manifold::Point *a = ba.second;
                    dfemb.insert(stiff_type::value_type(
                            ba, Set::VectorSpace::Hom(yemb[b].size(), yemb[a].size())));
                    DDF.insert(stiff_type::value_type(
                            ba, Set::VectorSpace::Hom(b->size(), a->size())));
                }
            }
            for (pDJ = DJ->begin(); pDJ != DJ->end(); pDJ++) {
                const Element::LocalState &LSLoc = *(*pDJ)->GetLocalState();
                const Element::Jet < 1 > &DJLoc = *(*pDJ);
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                set <Set::Manifold::Point *>::iterator pM, pN;
                vector_type yloc;
		vector_type::const_iterator py;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (py = yemb.find(a)) != yemb.end() ) {
		      yloc.insert(vector_type::value_type(a, py->second));
		    }
		    else {
		      assert(false);
		    }
                }
                pair<force_type, stiff_type> K = DJLoc(yloc);
                const force_type & floc = K.first;
                const stiff_type & dfloc = K.second;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
                    femb[a] += floc.find(a)->second;
                }
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NN = LSLoc.GetNodePairs();
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> >::iterator pNN;
                for (pNN = NN.begin(); pNN != NN.end(); pNN++) {
                    const pair<Set::Manifold::Point *, Set::Manifold::Point *> & ba = *pNN;
                    dfemb[ba] += dfloc.find(ba)->second;
                }
            }
            LS->Submerge(femb, DF);
            LS->Submerge(dfemb, DDF);
            return make_pair(DF, DDF);
        }

        pair<double, double>
        Jet < 1 > ::operator () (
                const set<Set::Manifold::Point *> &y,
                const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &u) const {            
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> force_type;
            typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> stiff_type;

            double DF = 0.0;
            double DDF = 0.0;
            vector_type yemb;
            LS->Embed(y, yemb);
            vector_type uemb;
            LS->Embed(u, uemb);

            set<Element::Jet < 1 > *>::iterator pDJ;
            for (pDJ = DJ->begin(); pDJ != DJ->end(); pDJ++) {
                const Element::LocalState &LSLoc = *(*pDJ)->GetLocalState();
                const Element::Jet < 1 > &DJLoc = *(*pDJ);
                set <Set::Manifold::Point *> N = LSLoc.GetNodes();
                set <Set::Manifold::Point *>::iterator pM, pN;
                vector_type yloc;
                vector_type uloc;
		vector_type::const_iterator py, pu;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
		    if ( (py=yemb.find(a)) != yemb.end() &&
			 (pu=uemb.find(a)) != uemb.end() ) {
		      yloc.insert(vector_type::value_type(a, py->second));
		      uloc.insert(vector_type::value_type(a, pu->second));
		    }
		    else {
		      assert(false);
		    }
		}
                pair<force_type, stiff_type> K = DJLoc(yloc);
                const force_type & floc = K.first;
                const stiff_type & dfloc = K.second;
                for (pN = N.begin(); pN != N.end(); pN++) {
                    Set::Manifold::Point *a = *pN;
                    DF += (floc.find(a)->second)(uloc[a]);
                }
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> > NN = LSLoc.GetNodePairs();
                set<pair<Set::Manifold::Point *, Set::Manifold::Point *> >::iterator pNN;
                for (pNN = NN.begin(); pNN != NN.end(); pNN++) {
                    const pair<Set::Manifold::Point *, Set::Manifold::Point *> & ba = *pNN;
                    Set::Manifold::Point *b = ba.first;
                    Set::Manifold::Point *a = ba.second;
                    DDF += (dfloc.find(ba)->second)(uloc[b])(uloc[a]);
                }
            }
            return make_pair(DF, DDF);
        }

        void
        Jet < 1 > ::operator () (
                const set<Set::Manifold::Point *> &y,
                Solver::Linear::System<DOF> &S) const {
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> force_type;
            typedef map<pair<Set::Manifold::Point *, Set::Manifold::Point *>, Set::VectorSpace::Hom> stiff_type;

            unsigned int i, j;
            S.SetToZero1();
            S.SetToZero2();
            pair<force_type, stiff_type> K = this->operator ()(y);
            const force_type & DF = K.first;
            const stiff_type & DDF = K.second;
            set<Set::Manifold::Point *>::const_iterator py;
            for (py = y.begin(); py != y.end(); py++) {
                Set::Manifold::Point *a = *py;
                const Set::VectorSpace::Vector & DFa = DF.find(a)->second;
                for (i = 0; i < DFa.size(); i++)
                    S.Add(make_pair(a, i), DFa[i]);
            }
            stiff_type::const_iterator pDDF;
            for (pDDF = DDF.begin(); pDDF != DDF.end(); pDDF++) {
                const pair<Set::Manifold::Point *, Set::Manifold::Point *> & ba = pDDF->first;
                Set::Manifold::Point *a = ba.first;
                Set::Manifold::Point *b = ba.second;
                const Set::VectorSpace::Hom & DDFba = pDDF->second;
                for (i = 0; i < DDFba.size1(); i++) {
                    DOF ai = make_pair(a, i);
                    for (j = 0; j < DDFba.size2(); j++) {
                        DOF bj = make_pair(b, j);
                        S.Add(ai, bj, DDFba[j][i]);
                    }
                }
            }
        }

        void
        Jet < 1 > ::operator () (
                Solver::Linear::System<DOF> &S,
                set<Set::Manifold::Point *> &y) const {
            unsigned int i;
            set<Set::Manifold::Point *>::iterator py;
            for (py = y.begin(); py != y.end(); py++) {
                Set::Manifold::Point *a = *py;
                Set::Manifold::Point &ya = *a;
                unsigned int n = ya.size();
                Set::VectorSpace::Vector u(n);
                for (i = 0; i < n; i++)
                    u[i] = S.Get(make_pair(a, i));
                ya -= u;
            }
        }

    }

}
