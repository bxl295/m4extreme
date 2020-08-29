// -*- C++ -*-

/*!
  \file stochastic/Direct.h
  \brief The direct method for SSA.
*/

#if !defined(__stochastic_Direct_h__)
#define __stochastic_Direct_h__

#include "Solver.h"
#include "Propensities.h"
#include "TimeEpochOffset.h"

#include "../array/StaticArrayOfArrays.h"
#include "../numerical/random/exponential/ExponentialGeneratorZiggurat.h"

namespace stochastic {

//! Perform a stochastic simulation using Gillespie's direct method.
/*!
  \param _DiscreteGenerator Random deviate generator for the discrete,
  finite distribution with reaction propensities as scaled probabilities.
  \param _ExponentialGenerator Random deviate generator for the exponential
  distribution. By default the ziggurat algorithm is used.
  \param _PropensitiesFunctor Can calculate propensities as a function of the
  reaction index and the populations.
*/
template<class _DiscreteGenerator,
         class _ExponentialGenerator =
         numerical::ExponentialGeneratorZiggurat<>,
         class _PropensitiesFunctor = PropensitiesSingle<true> >
class Direct : public Solver {
protected:

   //! The base class.
   typedef Solver Base;

   //
   // Public types.
   //
public:

   //! The propensities functor.
   typedef _PropensitiesFunctor PropensitiesFunctor;
   //! The exponential generator.
   typedef _ExponentialGenerator ExponentialGenerator;
   //! The discrete, finite generator.
   typedef _DiscreteGenerator DiscreteGenerator;
   //! The discrete, uniform generator.
   typedef typename ExponentialGenerator::DiscreteUniformGenerator
   DiscreteUniformGenerator;

   //
   // Enumerations.
   //
private:

   enum {ComputeIndividualPropensities =
            Loki::IsSameType < typename PropensitiesFunctor::result_type,
         double >::value
        };

   //
   // Member data.
   //
protected:

   //! The time is tracked with an epoch and offset.
   TimeEpochOffset _time;
   //! The propensities functor.
   PropensitiesFunctor _propensitiesFunctor;
   //! The reaction influence graph.
   array::StaticArrayOfArrays<std::size_t> _reactionInfluence;
   //! The discrete uniform random number generator.
   DiscreteUniformGenerator _discreteUniformGenerator;
   //! The exponential deviate generator.
   ExponentialGenerator _exponentialGenerator;
   //! The discrete random number generator.
   DiscreteGenerator _discreteGenerator;
   //! The sequence of propensity values.
   std::vector<double> _propensities;
   //! The time step.
   double _tau;

   //
   // Not implemented.
   //
private:

   //! Default constructor not implemented.
   Direct();
   //! Copy constructor not implemented.
   Direct(const Direct&);
   //! Assignment operator not implemented.
   Direct&
   operator=(const Direct&);

   //--------------------------------------------------------------------------
   //! \name Constructors etc.
   //@{
public:

   //! Construct.
   Direct(const State& state,
          const PropensitiesFunctor& propensitiesFunctor,
          const array::StaticArrayOfArrays<std::size_t>& reactionInfluence,
          const double maxSteps) :
      Base(state, maxSteps),
      // Copy.
      _time(),
      _propensitiesFunctor(propensitiesFunctor),
      _reactionInfluence(reactionInfluence),
      // Construct.
      _discreteUniformGenerator(),
      _exponentialGenerator(&_discreteUniformGenerator),
      _discreteGenerator(&_discreteUniformGenerator),
      // Invalid value.
      _tau(-1) {
      // Allocate memory if needed.
      if (! ComputeIndividualPropensities) {
         _propensities.resize(state.getNumberOfReactions());
      }
   }

   // Default destructor is fine.

   //@}
   //--------------------------------------------------------------------------
   //! \name Simulation.
   //@{
public:

   //! Initialize the state with the initial populations and reset the time.
   void
   initialize(const std::vector<double>& populations, const double time) {
      // Initialize the state.
      Base::initialize(populations);
      _time = time;
      // Compute the propensities and initialize the discrete, finite generator.
      computePropensities();
      // Compute the initial time step.
      _tau = computeTau();
   }

   //! Simulate until the end time is reached.
   void
   simulate(const double endTime) {
      // Step until no more reactions can fire or we reach the termination
      // condition.
      while (step(endTime)) {
      }
   }

   //! Simulate until the termination condition is reached.
   /*!
     Record each reaction index and time.
   */
   template<typename _IntOutputIter, typename _NumberOutputIter>
   void
   simulate(const double endTime, _IntOutputIter indices,
            _NumberOutputIter times) {
      // Step until no more reactions can fire or we reach the termination
      // condition.
      while (step(endTime, indices, times)) {
      }
   }

   //
   // Note: Ideally it would be best to have a single step function. But to
   // get the best performance I have to duplicate code.
   //

   //! Try to take a step.  Return true if a step is taken.
   bool
   step(double endTime);

   //! Try to take a step.  Return true if a step is taken.
   /*!
     Record the reaction index and time.
   */
   template<typename _IntOutputIter, typename _NumberOutputIter>
   bool
   step(double endTime, _IntOutputIter indices, _NumberOutputIter times);

   //! Let the process equilibrate.
   void
   equilibrate(const double equilibrationTime) {
      // Let the process equilibrate.
      const double endTime = _time + equilibrationTime;
      simulate(endTime);
      // The next (unfired) reaction will carry the simulation past the end time.
      _time = endTime;
   }

protected:

   //! Record a step count error message. Record the current time.
   /*! Override the same function from the base class. */
   void
   setStepCountError() {
      std::ostringstream out;
      out << "The maximum step count " << _maxSteps << " was reached at time = "
          << _time << ".";
      _error += out.str();
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Accessors.
   //@{
public:

   //! Return the current time.
   double
   getTime() const {
      // Convert the time epoch and offset to a single time.
      return _time;
   }

   //! Return a const reference to the discrete, uniform generator.
   const DiscreteUniformGenerator&
   getDiscreteUniformGenerator() const {
      return _discreteUniformGenerator;
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Manipulators.
   //@{
public:

   //! Return a reference to the discrete, finite generator.
   DiscreteGenerator&
   getDiscreteGenerator() {
      return _discreteGenerator;
   }

   //! Return a reference to the discrete, uniform generator.
   DiscreteUniformGenerator&
   getDiscreteUniformGenerator() {
      return _discreteUniformGenerator;
   }

protected:

   //! Compute the time to the next reaction.
   double
   computeTau() {
      // If any reactions can fire.
      if (_discreteGenerator.isValid()) {
         const double mean = 1.0 / _discreteGenerator.sum();
         _time.updateEpoch(mean);
         // Return the time to the next reaction.
         return mean * _exponentialGenerator();
      }
      // Otherwise return infinity.
      return std::numeric_limits<double>::max();
   }

   void
   computePropensities() {
      _computePropensities(Loki::Int2Type<ComputeIndividualPropensities>());
   }

   void
   _computePropensities(Loki::Int2Type<false> /*dummy*/) {
      // Use the member function propensities array.
      _propensitiesFunctor(_state.getPopulations(), _propensities.begin());
      _discreteGenerator.initialize(_propensities.begin(), _propensities.end());
   }

   void
   _computePropensities(Loki::Int2Type<true> /*dummy*/) {
      // Allocate a propensity array.
      std::vector<double> propensities(_state.getNumberOfReactions());
      // Compute each propensity.
      for (std::size_t i = 0; i != propensities.size(); ++i) {
         propensities[i] = _propensitiesFunctor(i, _state.getPopulations());
      }
      _discreteGenerator.initialize(propensities.begin(), propensities.end());
   }

   void
   updatePropensities(const std::size_t reactionIndex) {
      _updatePropensities
      (Loki::Int2Type<ComputeIndividualPropensities>(), reactionIndex);
   }

   void
   _updatePropensities(Loki::Int2Type<false> /*dummy*/,
                       const std::size_t reactionIndex) {
      // Recompute all of the propensity functions.
      _propensitiesFunctor(_state.getPopulations(), _propensities.begin());
      _discreteGenerator.set(_propensities.begin());
      updateSum(Loki::Int2Type<DiscreteGenerator::AutomaticUpdate>());
   }

   void
   _updatePropensities(Loki::Int2Type<true> /*dummy*/,
                       const std::size_t reactionIndex) {
      for (typename array::StaticArrayOfArrays<std::size_t>::const_iterator
            i = _reactionInfluence.begin(reactionIndex);
            i != _reactionInfluence.end(reactionIndex); ++i) {
         _discreteGenerator.set
         (*i, _propensitiesFunctor(*i, _state.getPopulations()));
      }
      updateSum(Loki::Int2Type<DiscreteGenerator::AutomaticUpdate>());
   }

   //! Do nothing because the discrete generator automatically update the PMF sum.
   void
   updateSum(Loki::Int2Type<true> /*Automatic update*/) {
   }

   //! Tell the discrete generator to update the PMF sum.
   void
   updateSum(Loki::Int2Type<false> /*Automatic update*/) {
      _discreteGenerator.updateSum();
   }

   //@}
};

//@}

} // namespace stochastic

#define __stochastic_Direct_ipp__
#include "Direct.ipp"
#undef __stochastic_Direct_ipp__

#endif
