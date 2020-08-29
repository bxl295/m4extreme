// LumpedMass.cpp: implementation of the LumpedMass class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include <sys/time.h>
#include "./LumpedMass.h"

namespace Model {

#if defined(_M4EXTREME_THREAD_POOL)
    std::vector<int> LumpedMass::_Costs;
    std::vector<int> LumpedMass::_Costs_new;

#if defined(_M4EXTREME_THREAD_EFFICIENT_LOCK_)
    ost::Mutex LumpedMass::_lock;    
#endif

#endif

    LumpedMass::LumpedMass() {
    }

    LumpedMass::~LumpedMass() {
    }

#if defined(_M4EXTREME_THREAD_POOL)
    LumpedMass::LumpedMass(vector<Element::LumpedMass *> *rhs_M) : _LM_mt(rhs_M), M(NULL) {
        _Costs.resize(_LM_mt->size(), 0);
        _Costs_new.resize(_LM_mt->size(), 0);
    }
#endif

    LumpedMass::LumpedMass(set<Element::LumpedMass *> *rhs_M) : 
      M(rhs_M), _computational_cost(0) {}

    LumpedMass::LumpedMass(const LumpedMass &rhs) : M(rhs.M),_computational_cost(0)
#if defined(_M4EXTREME_THREAD_POOL)
    ,_LM_mt(rhs._LM_mt)
#endif
    {}

    LumpedMass &
    LumpedMass::operator =(const LumpedMass &rhs) {
        if (this == &rhs) return *this;
        M = rhs.M;
	_computational_cost = 0;
#if defined(_M4EXTREME_THREAD_POOL)
        _LM_mt = rhs._LM_mt;
#endif       
        return *this;
    }

    double 
    LumpedMass::GetComputationalCost() const{
#if defined(_M4EXTREME_THREAD_POOL)
      double sum = 0.0;
      for ( int i = 0; i < _Costs.size(); ++i ) {
	sum += (double)_Costs[i] / 1.0e6;
      }

      return sum;
#else
      return _computational_cost;
#endif
    }

    map<Set::Manifold::Point *, double>
    LumpedMass::GetMass() const {
        map<Set::Manifold::Point *, double> Mass;

#if defined(_M4EXTREME_THREAD_POOL)
        assert(false);
#else           
        map<Set::Manifold::Point*, double>::iterator pLM, pMass;

        set<Element::LumpedMass *>::iterator pM;
        for (pM = M->begin(); pM != M->end(); pM++) {

            map<Set::Manifold::Point *, double> LocalMass = (*pM)->GetMass();

            for (pLM = LocalMass.begin(); pLM != LocalMass.end(); pLM++) {
                pMass = Mass.find(pLM->first);
                if (pMass != Mass.end()) {
                    pMass->second += pLM->second;
                } else {
                    Mass.insert(make_pair(pLM->first, pLM->second));
                }
            }

        }
#endif

        return Mass;
    }

void LumpedMass::GetMass(range_type & Mass) {

        if (!Mass.empty()) Mass.clear();

#if defined(_M4EXTREME_THREAD_POOL)
        _eureka_thread_arg arg;
        arg._pLM = _LM_mt;

#if defined(_M4EXTREME_THREAD_EFFICIENT_LOCK_)
        arg._pMass = &Mass;
#else
        int numThreads = m4extreme::Utils::GetNumberofThreads();
        std::vector<range_type> MassVec(numThreads);
        arg._pMass = &MassVec;
#endif

	int numMPTs = _LM_mt->size();
	if ( numMPTs > _Costs.size() ) {
	  for ( int i = _Costs.size(); i < numMPTs; ++i ) {
	    _Costs.push_back(0);
	    _Costs_new.push_back(0);
	  }
	}

        m4extreme::Utils::RunThreadMonitor(_getLumpedMass, reinterpret_cast<void*> (&arg));
        
	std::copy(_Costs_new.begin(), _Costs_new.end(), _Costs.begin());

#if !defined(_M4EXTREME_THREAD_EFFICIENT_LOCK_)
        range_type::iterator pMass;
        for (int i = 0; i < numThreads; ++i) {
            range_type::iterator pM = MassVec[i].begin();
            while (pM != MassVec[i].end()) {
                pMass = Mass.find(pM->first);
                if (pMass != Mass.end()) {
                    pMass->second += pM->second;
                } else {
                    Mass.insert(make_pair(pM->first, pM->second));
                }

                pM++;
            }
        }
#endif

#else           
	struct timeval start, stop;
	gettimeofday(&start, NULL);

        map<Set::Manifold::Point*, double>::iterator pLM, pMass;
        map<Set::Manifold::Point *, double> LocalMass;

        set<Element::LumpedMass *>::iterator pM;
        for (pM = M->begin(); pM != M->end(); pM++) {

            (*pM)->GetMass(LocalMass);

            for (pLM = LocalMass.begin(); pLM != LocalMass.end(); pLM++) {
                pMass = Mass.find(pLM->first);
                if (pMass != Mass.end()) {
                    pMass->second += pLM->second;
                } else {
                    Mass.insert(make_pair(pLM->first, pLM->second));
                }
            }

        }

	gettimeofday(&stop, NULL);
	_computational_cost = stop.tv_sec - start.tv_sec
	  + (stop.tv_usec - start.tv_usec) / 1.0e6; //unit in seconds;
#endif

        return;

    }

#if defined(_M4EXTREME_THREAD_POOL)

    void * LumpedMass::_getLumpedMass(void * p) {

        int numThreads = m4extreme::Utils::GetNumberofThreads();
        int myId = m4extreme::Utils::GetMyThreadID();

        _eureka_thread_arg * ploc = reinterpret_cast<_eureka_thread_arg *> (p);
        const std::vector<Element::LumpedMass *> & LM = *ploc->_pLM;

#if defined(_M4EXTREME_THREAD_EFFICIENT_LOCK_)
        range_type mass;
#else
        range_type & mass = (*ploc->_pMass)[myId];
#endif

        int numMPTs = LM.size();
        int i_start, i_end;

#if defined(_M4EXTREME_THREAD_POOL_ADAPTIVE)
        struct timeval start, stop;
        m4extreme::Utils::GetDataShare(myId, numThreads, _Costs, i_start, i_end);
#else
        m4extreme::Utils::GetDataShare(myId, numThreads, numMPTs, i_start, i_end);
#endif

        map<Set::Manifold::Point *, double>::iterator pLM, pMass;
        map<Set::Manifold::Point *, double> LocalMass;

        for (int i = i_start; i <= i_end; i++) {

            gettimeofday(&start, NULL);

            LM[i]->GetMass(LocalMass);
            for (pLM = LocalMass.begin(); pLM != LocalMass.end(); pLM++) {
                pMass = mass.find(pLM->first);
                if (pMass != mass.end()) {
                    pMass->second += pLM->second;
                } else {
                    mass.insert(make_pair(pLM->first, pLM->second));
                }
            }

            gettimeofday(&stop, NULL);
            _Costs_new[i] = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec; //unit is microsecond
        }

#if defined(_M4EXTREME_THREAD_EFFICIENT_LOCK_)
        range_type & TotalMass = *ploc->_pMass;
        _lock.EnterMutex();
        for (pLM = mass.begin(); pLM != mass.end(); pLM++) {
                pMass = TotalMass.find(pLM->first);
                if (pMass != TotalMass.end()) {
                    pMass->second += pLM->second;
                } else {
                    TotalMass.insert(make_pair(pLM->first, pLM->second));
                }
        }
        _lock.LeaveMutex();
#endif

        return NULL;
    }

#endif

}
