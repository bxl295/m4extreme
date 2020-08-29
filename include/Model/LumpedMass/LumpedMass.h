// LumpedMass.h: interface for the LumpedMass class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#if !defined(MODEL_LUMPEDMASS_H__INCLUDED_)
#define MODEL_LUMPEDMASS_H__INCLUDED_

#pragma once

#include "../Model.h"
#include "../../Element/LumpedMass.h"

#if defined(_M4EXTREME_THREAD_POOL)
#include "cc++/thread.h"
#include "Threads/ThreadMonitorIface.h"
#endif

namespace Model
{
  class LumpedMass
  {
  public: 
    typedef map<Set::Manifold::Point *, double> range_type;
    
    LumpedMass();
    virtual ~LumpedMass();

#if defined(_M4EXTREME_THREAD_POOL)
    LumpedMass(vector<Element::LumpedMass *> *);
#endif

    LumpedMass(set<Element::LumpedMass *> *);
    
    LumpedMass(const LumpedMass &);
    LumpedMass & operator = (const LumpedMass &);
    range_type GetMass() const;
        
    void GetMass(range_type &);
    
    double GetComputationalCost() const;

  private:
    set<Element::LumpedMass *> *M;
    double _computational_cost;
    
#if defined(_M4EXTREME_THREAD_POOL)
    static void * _getLumpedMass(void *);
    std::vector<Element::LumpedMass *>* _LM_mt;
    static std::vector<int> _Costs, _Costs_new;

#if defined(_M4EXTREME_THREAD_EFFICIENT_LOCK_)
    typedef struct {
        const std::vector<Element::LumpedMass *>* _pLM;
        range_type * _pMass;
    } _eureka_thread_arg;

    static ost::Mutex _lock;    
#else
    typedef struct {
        const std::vector<Element::LumpedMass *>* _pLM;
        std::vector<range_type> * _pMass;
    } _eureka_thread_arg;
#endif

#endif

  };

}

#endif // !defined(MODEL_LUMPEDMASS_H__INCLUDED_
