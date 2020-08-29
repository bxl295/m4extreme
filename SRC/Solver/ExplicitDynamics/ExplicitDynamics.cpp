// ExplicitDynamics.cpp: implementation of the ExplicitDynamics class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./ExplicitDynamics.h"

#ifdef _M4EXTREME_DEBUG_II
extern map<Set::Manifold::Point *, Set::VectorSpace::Vector> _gForce;
extern map<Set::Manifold::Point *, size_t> _gIndexSet;
#endif

namespace Solver {

    ExplicitDynamics::ExplicitDynamics() {
    }

    ExplicitDynamics::~ExplicitDynamics() {
    }

  ExplicitDynamics::ExplicitDynamics(
            Clock *T_,
            Model::LocalState *LS_,
            Model::Energy < 1 > *DE_,
            map<Set::Manifold::Point *, double> *m_,
            set<Set::Manifold::Point *> *x_,
            map<Set::Manifold::Point *, Set::VectorSpace::Vector> *v_,
  	    const double GamNew_) :
      Print(0), T(T_), LS(LS_), DE(DE_), m(m_), x(x_), v(v_),
    GamNew(GamNew_), GamOld(1.0 - GamNew_) {
    }

    ExplicitDynamics::ExplicitDynamics(
            Clock *T_,
            Model::LocalState *LS_,
            Model::Energy < 1 > *DE_,
            map<Set::Manifold::Point *, double> *m_,
            set<Set::Manifold::Point *> *x_,
            map<Set::Manifold::Point *, Set::VectorSpace::Vector> *v_,
            const map<Set::Manifold::Point *, Set::VectorSpace::Vector> & a_,
	    const double GamNew_) :
      Print(0), T(T_), LS(LS_), DE(DE_), m(m_), x(x_), v(v_), a(a_),
    GamNew(GamNew_), GamOld(1.0 - GamNew_) {
    }

    bool &
    ExplicitDynamics::SetPrint() {
        return Print;
    }

    void
    ExplicitDynamics::UpdateMass(map<Set::Manifold::Point*, double> * m_) {
        m = m_;
        return;
    }

    void
    ExplicitDynamics::operator ++() {
        if (Print) cout <<"Time = " << T->Time() << ", Explicit dynamics: step begin" << endl;

        map<Set::Manifold::Point *, double>::iterator pm;
        set<Set::Manifold::Point *>::iterator px;
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pf;
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pv;
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pa;

        if (a.size() == 0) {        
         
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> accel_type;
            //map<Set::Manifold::Point *, Set::VectorSpace::Vector> f;	    
            DE->operator()(*x, f);
            for (pm = m->begin(), pf = f.begin(); pm != m->end(); pm++, pf++)
                a.insert(accel_type::value_type(pm->first, -pf->second / pm->second));
            
#ifdef _M4EXTREME_DEBUG_II
            std::cout << "\nInitialization BEGIN >>>>>>>>>>>>>" << std::endl;

            for (px = x->begin(); px != x->end(); px++) {
                std::cout << _gIndexSet[*px] << "th Node(" << *static_cast<Set::Euclidean::Orthonormal::Point*> (*px)
                        << ") v=[" << v->find(*px)->second
                        << "] a=[" << a.find(*px)->second
                        << "]" << std::endl;
            }
            
            std::cout << "<<<<<<<<<<<<<< Initialization DONE\n" << std::endl;
#endif            

        }

        double dt = T->DTime();
        double dtold = GamOld*dt;
        double dtnew = GamNew*dt;
        double dtdthalf = 0.5 * dt*dt;

        for (px = x->begin(), pv = v->begin(), pa = a.begin(); px != x->end(); px++, pv++, pa++) {
            *(*px) += dt * pv->second + dtdthalf * pa->second;
            
#ifdef _M4EXTREME_DEBUG_II
            std::cout << _gIndexSet[*px] << "th Node(" << *static_cast<Set::Euclidean::Orthonormal::Point*>(*px) 
                    << ") v=[" << pv->second
                    << "] a=[" << pa->second 
                    << "]" << std::endl;
#endif
        }

        DE->operator()(*x, f);

        double TOL = 1.0e-16;
        double ma = 0.0;
        for ( pf = f.begin(); pf != f.end(); pf++ ) {
            if ( (pm=m->find(pf->first)) != m->end() ) {
                ma = pm->second;

#ifdef _M4EXTREME_DEBUG_II
            cout << "force for node " << _gIndexSet[pf->first]
                    << " with mass " << ma
                    << "------------->[" << pf->second << "]" << endl; ////
#endif            
	    
                Set::VectorSpace::Vector & vloc = v->find(pf->first)->second;
                Set::VectorSpace::Vector & aloc = a.find(pf->first)->second;

                if ( ma > TOL ) {
                    Set::VectorSpace::Vector anew = -pf->second / ma;
                    vloc += dtold * aloc + dtnew*anew;
                    aloc = anew;
                }
                else {
                    Null(aloc);
                }
            }
        }	
        
        ++(*LS);
        ++(*T);

        if (Print) cout << "Time = " << T->Time() << ", Explicit dynamics: step end" << endl;
    }

    void
    ExplicitDynamics::BallisticUpdate() {
        double dt = T->DTime();

        set<Set::Manifold::Point *>::iterator px;
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pv;
        for (px = x->begin(), pv = v->begin(); px != x->end(); px++, pv++) {
            *(*px) += dt * pv->second;
        }

	return;
    }
  
    /**
     * predict the location of nodes at t+dt
     */
    void 
    ExplicitDynamics::Predictor() {


        set<Set::Manifold::Point *>::iterator px, xend = x->end();
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pv;
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pa;

        if (a.size() == 0) {        
         
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> accel_type;

            DE->operator()(*x, f);
	    map<Set::Manifold::Point *, double>::iterator pm, mend = m->end();
	    map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pf;
            for (pm = m->begin(), pf = f.begin(); pm != mend; pm++, pf++)
                a.insert(accel_type::value_type(pm->first, -pf->second / pm->second));
	}

        double dt = T->DTime();
        double dtdthalf = 0.5 * dt*dt;

        for (px = x->begin(), pv = v->begin(), pa = a.begin(); px != xend; px++, pv++, pa++) {
            *(*px) += dt * pv->second + dtdthalf * pa->second;
            
#ifdef _M4EXTREME_DEBUG_II
            std::cout << _gIndexSet[*px] << "th Node: [" 
		      << *static_cast<Set::Euclidean::Orthonormal::Point*>(*px) 
		      << "]" << std::endl;
#endif
        }

	return;
    }

    /**
     * predict the location of nodes at t+dt
     */
    void 
    ExplicitDynamics::Predictor(const set<Set::Manifold::Point *> & fixedDoF) {


        set<Set::Manifold::Point *>::iterator px, xend = x->end();
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pv;
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pa;

        if (a.size() == 0) {        
         
            typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> accel_type;

            DE->operator()(*x, f);
	    map<Set::Manifold::Point *, double>::iterator pm, mend = m->end();
	    map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pf;
            for (pm = m->begin(), pf = f.begin(); pm != mend; pm++, pf++)
                a.insert(accel_type::value_type(pm->first, -pf->second / pm->second));
	}

        double dt = T->DTime();
        double dtdthalf = 0.5 * dt*dt;

        for (px = x->begin(), pv = v->begin(), pa = a.begin(); px != xend; px++, pv++, pa++) {
	  if ( fixedDoF.find(*px) == fixedDoF.end() ) {
            *(*px) += dt * pv->second + dtdthalf * pa->second;
            
#ifdef _M4EXTREME_DEBUG_II
            std::cout << _gIndexSet[*px] << "th Node: [" 
		      << *static_cast<Set::Euclidean::Orthonormal::Point*>(*px) 
		      << "]" << std::endl;
#endif
	  }
        }

	return;
    }

    /**
     * correct the kinematic info of nodes at t+dt
     */
    void
    ExplicitDynamics::Corrector() {

        DE->operator()(*x, f);

        double TOL = 1.0e-16;
        double dt = T->DTime();
        double dtold = GamOld*dt;
        double dtnew = GamNew*dt;

        double ma = 0.0;
	map<Set::Manifold::Point *, double>::iterator pm, mend = m->end();
	map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pf, fend=f.end();
        for ( pf = f.begin(); pf != fend; pf++ ) {
            if ( (pm=m->find(pf->first)) != mend ) {
                ma = pm->second;

#ifdef _M4EXTREME_DEBUG_II
            cout << "force for node " << _gIndexSet[pf->first]
                    << " with mass " << ma
                    << "------------->[" << pf->second << "]" << endl; ////
#endif            
	    
                Set::VectorSpace::Vector & vloc = v->find(pf->first)->second;
                Set::VectorSpace::Vector & aloc = a.find(pf->first)->second;

                if ( ma > TOL ) {
                    Set::VectorSpace::Vector anew = -pf->second / ma;
                    vloc += dtold * aloc + dtnew*anew;
                    aloc = anew;
                }
                else {
                    Null(aloc);
                }
            }
        }	
        
        ++(*LS);
        ++(*T);

	return;
    }


    /**
     * pull back the location of nodes at from t+dt to t
     */
    void 
    ExplicitDynamics::PullBack() {

        set<Set::Manifold::Point *>::iterator px, xend = x->end();
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pv;
        map<Set::Manifold::Point *, Set::VectorSpace::Vector>::iterator pa;

        if (a.size() == 0) { 
	  assert(false);
	}
         
        double dt = T->DTime();
        double dtdthalf = 0.5 * dt*dt;

        for (px = x->begin(), pv = v->begin(), pa = a.begin(); px != xend; px++, pv++, pa++) {
            *(*px) -= dt * pv->second + dtdthalf * pa->second;
        }

	return;
    }

}
