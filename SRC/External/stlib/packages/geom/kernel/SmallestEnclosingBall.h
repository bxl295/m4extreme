// -*- C++ -*-

#include "../../ext/array.h"

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <list>

namespace geom {

//! The smallest ball with a set of n <= _Dimension+1 points on the boundary.
template<std::size_t _Dimension, typename _T=double>
class SmallestEnclosingBallBoundary {
   //
   // Types
   //
public:

   //! A Cartesian point.
   typedef std::tr1::array<_T, _Dimension> Point;

   //
   // Member data.
   //
private:

   // Size.
   std::size_t m;
   // Number of support points.
   std::size_t s;
   Point q0;

   _T z[_Dimension+1];
   _T f[_Dimension+1];
   Point v[_Dimension+1];
   Point a[_Dimension+1];

   Point c[_Dimension+1];
   _T squaredRadii[_Dimension+1];

   // refers to some c[j]
   Point* currentCenter;
   _T currentSquaredRadius;

public:
   //--------------------------------------------------------------------------
   /*! \name Constructors etc.
     We use the default copy constructor, assignment operator, and destructor.
   */
   // @{

   SmallestEnclosingBallBoundary() {
      reset();
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Accesors.
   // @{

   const Point&
   center() const {
      return *currentCenter;
   }

   _T
   squaredRadius() const {
      return currentSquaredRadius;
   }

   std::size_t
   size() const {
      return m;
   }

   std::size_t
   supportSize() const {
      return s;
   }

   _T
   excess(const Point& p) const {
      return squaredDistance(p, *currentCenter) - currentSquaredRadius;
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Manipulators.
   // @{

   //! Generates an empty sphere with zero size and no support points.
   void
   reset();

   bool
   push(const Point& p);

   void
   pop() {
      --m;
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Validation.
   // @{

   _T
   slack() const;

   // @}
};

template<std::size_t _Dimension, typename _T>
void
SmallestEnclosingBallBoundary<_Dimension,_T>::reset() {
   m = s = 0;
   // we misuse c[0] for the center of the empty sphere
   std::fill(c[0].begin(), c[0].end(), 0);
   currentCenter = &c[0];
   currentSquaredRadius = -1;
}




template<std::size_t _Dimension, typename _T>
bool
SmallestEnclosingBallBoundary<_Dimension,_T>::push(const Point& p) {
   std::size_t i, j;
   _T eps = 1e-32;
   if (m == 0) {
      q0 = p;
      c[0] = q0;
      squaredRadii[0] = 0;
   }
   else {
      // set v_m to Q_m
      v[m] = p - q0;

      // compute the a_{m,i}, i < m
      for (i = 1; i != m; ++i) {
         a[m][i] = 2. * dot(v[i], v[m]) / z[i];
      }

      // update v_m to Q_m-\bar{Q}_m
      for (i = 1; i != m; ++i) {
         for (j = 0; j != _Dimension; ++j) {
            v[m][j] -= a[m][i] * v[i][j];
         }
      }

      // compute z_m
      z[m] = 2. * squaredMagnitude(v[m]);

      // reject push if z_m too small
      if (z[m] < eps * currentSquaredRadius) {
         return false;
      }

      // update c, sqr_r
      _T e = squaredDistance(p, c[m-1]) - squaredRadii[m-1];
      f[m] = e / z[m];

      for (i = 0; i != _Dimension; ++i) {
         c[m][i] = c[m-1][i] + f[m] * v[m][i];
      }
      squaredRadii[m] = squaredRadii[m-1] + 0.5 * e * f[m];
   }
   currentCenter = &c[m];
   currentSquaredRadius = squaredRadii[m];
   s = ++m;
   return true;
}

template<std::size_t _Dimension, typename _T>
_T
SmallestEnclosingBallBoundary<_Dimension,_T>::slack() const {
   _T l[_Dimension+1];
   _T minL = 0;
   l[0] = 1;
   for (std::size_t i = s - 1; i > 0; --i) {
      l[i] = f[i];
      for (std::size_t k = s - 1; k > i; --k) {
         l[i] -= a[k][i] * l[k];
      }
      if (l[i] < minL) {
         minL = l[i];
      }
      l[0] -= l[i];
   }
   if (l[0] < minL) {
      minL = l[0];
   }
   return (minL < 0) ? -minL : 0;
}


//! Smallest enclosing ball of a set of points.
template<std::size_t _Dimension, typename _T=double>
class SmallestEnclosingBall {
   //
   // Types.
   //
public:
   //! A Cartesian point.
   typedef std::tr1::array<_T, _Dimension> Point;
   typedef typename std::list<Point>::iterator iterator;
   typedef typename std::list<Point>::const_iterator const_iterator;

   //
   // Member data.
   //
private:
   //! Internal point set.
   std::list<Point> L;
   //! The current ball.
   SmallestEnclosingBallBoundary<_Dimension,_T> B;
   // The end of the support set.
   iterator supportEnd;

   //--------------------------------------------------------------------------
   /*! \name Constructors etc.
     We use the default destructor. The copy constructor and assignment
     operator are private.
   */
   // @{
public:

   // Create an empty ball.
   SmallestEnclosingBall() {
   }

   // Copy p to the internal point set.
   void
   checkIn(const Point& p) {
      L.push_back(p);
   }

   // Build the smallest enclosing ball of the internal point set.
   void
   build() {
      B.reset();
      supportEnd = L.begin();
      pivot_mb(L.end());
   }

   // Build the smallest enclosing ball with the specified internal point set.
   template<typename _InputIterator>
   void
   build(_InputIterator first, _InputIterator last) {
      L.clear();
      L.insert(L.begin(), first, last);
      build();
   }

private:

   // Not implemented.
   SmallestEnclosingBall(const SmallestEnclosingBall&);

   // Not implemented.
   SmallestEnclosingBall&
   operator=(const SmallestEnclosingBall&);

   // @}
   //--------------------------------------------------------------------------
   //! \name Accessors.
   // @{
public:

   //! Return the center of the ball (undefined if ball is empty).
   const Point&
   center() const {
      return B.center();
   }

   //! Return the squared radius of the ball (-1 if ball is empty).
   _T
   squaredRadius() const {
      return B.squaredRadius();
   }

   //! Return the size of the internal point set.
   std::size_t
   numInternalPoints() const {
      return L.size();
   }

   //! Return the beginning of the internal point set.
   const_iterator
   pointsBegin() const {
      return L.begin();
   }

   //! Return the end of the internal point set.
   const_iterator
   pointsEnd() const {
      return L.end();
   }

   //! Return the size of the support point set.
   /*! The point set has the following properties:
     - There are at most _Dimension + 1 support points.
     - All support points are on the boundary of the computed ball.
     - The smallest enclosing ball of the support point set equals the
     smallest enclosing ball of the internal point set. */
   std::size_t
   numSupportPoints() const {
      return B.supportSize();
   }

   //! Return the beginning of the internal point set.
   const_iterator
   supportPointsBegin() const {
      return L.begin();
   }

   //! Return the end of the internal point set.
   const_iterator
   supportPointsEnd() const {
      return supportEnd;
   }

   //! Assess the quality of the computed ball. 
   /*! The return value is the maximum squared distance of any support
    point or point outside the ball to the boundary of the ball,
    divided by the squared radius of the ball. If everything went
    fine, this will be less than e-15 and says that the computed ball
    approximately contains all the internal points and has all the
    support points on the boundary.
   
    The slack parameter that is set by the method says something about
    whether the computed ball is really the *smallest* enclosing ball
    of the support points; if everything went fine, this value will be
    0; a positive value may indicate that the ball is not smallest
    possible, with the deviation from optimality growing with the
    slack. */
   _T
   accuracy(_T* slack) const;

   // returns true if the accuracy is below the given tolerance and the
   // slack is 0
   bool
   isValid(const _T tolerance = std::numeric_limits<_T>::epsilon())
      const {
      _T slack;
      return ((accuracy(&slack) < tolerance) && (slack == 0));
   }

   // @}
   //
   // Private methods.
   //
private:

   void
   mtf_mb(iterator k);

   void
   pivot_mb(iterator k);

   void
   moveToFront(iterator j);

   _T
   maxExcess(iterator t, iterator i, iterator* pivot) const;
};


template<std::size_t _Dimension, typename _T>
void
SmallestEnclosingBall<_Dimension,_T>::
mtf_mb(iterator i) {
   supportEnd = L.begin();
   if (B.size() == _Dimension + 1) {
      return;
   }
   for (iterator k = L.begin(); k != i;) {
      iterator j = k++;
      if (B.excess(*j) > 0) {
         if (B.push(*j)) {
            mtf_mb(j);
            B.pop();
            moveToFront(j);
         }
      }
   }
}

template<std::size_t _Dimension, typename _T>
void
SmallestEnclosingBall<_Dimension,_T>::
moveToFront(iterator j) {
   if (supportEnd == j) {
      supportEnd++;
   }
   L.splice(L.begin(), L, j);
}


template<std::size_t _Dimension, typename _T>
void
SmallestEnclosingBall<_Dimension,_T>::
pivot_mb(iterator i) {
   iterator t = ++L.begin();
   mtf_mb(t);
   _T maxE;
   _T oldSquaredRadius = -1;
   do {
      iterator pivot;
      maxE = maxExcess(t, i, &pivot);
      if (maxE > 0) {
         t = supportEnd;
         if (t == pivot) {
            ++t;
         }
         oldSquaredRadius = B.squaredRadius();
         B.push(*pivot);
         mtf_mb(supportEnd);
         B.pop();
         moveToFront(pivot);
      }
   }
   while (maxE > 0 && B.squaredRadius() > oldSquaredRadius);
}


template<std::size_t _Dimension, typename _T>
_T
SmallestEnclosingBall<_Dimension,_T>::
maxExcess(iterator t, iterator i, iterator* pivot) const {
   const Point& c = B.center();
   const _T sqr_r = B.squaredRadius();
   _T e, maxE = 0;
   for (iterator k = t; k != i; ++k) {
      const Point& p = *k;
      e = squaredDistance(p, c) - sqr_r;
      if (e > maxE) {
         maxE = e;
         *pivot = k;
      }
   }
   return maxE;
}


template<std::size_t _Dimension, typename _T>
_T
SmallestEnclosingBall<_Dimension,_T>::
accuracy(_T* slack) const {
   _T e, maxE = 0;
   std::size_t numSupp = 0;
   const_iterator i;
   for (i = L.begin(); i != supportEnd; ++i, ++numSupp) {
      if ((e = std::abs(B.excess(*i))) > maxE) {
         maxE = e;
      }
   }

   // you've found a non-numerical problem if the following ever fails
   assert(numSupp == numSupportPoints());

   for (i = supportEnd; i != L.end(); ++i) {
      if ((e = B.excess(*i)) > maxE) {
         maxE = e;
      }
   }

   *slack = B.slack();
   return maxE / squaredRadius();
}

} // namespace geom
