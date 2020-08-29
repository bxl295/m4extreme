// -*- C++ -*-

#if !defined(__numerical_derivative_ridders_ipp__)
#error This file is an implementation detail of ridders.
#endif

namespace numerical {


// Calculate the derivative using Ridders' method.
template<typename _Functor>
inline
double
derivativeRidders(_Functor& f, const double x, const double scale,
                  double* error, ads::SquareMatrix<10, double>* a) {
   const double shrink = 1. / 1.4;
   const double grow2 = 1.4 * 1.4;
   assert(scale != 0);

   double err, fac;
   double h = scale;
   double result = std::numeric_limits<double>::max();
   *error = std::numeric_limits<double>::max();
   (*a)(0, 0) = (f(x + h) - f(x - h)) / (2.0 * h);
   for (std::size_t i = 1; i != 10; ++i) {
      h *= shrink;
      (*a)(i, 0) = (f(x + h) - f(x - h)) / (2.0 * h);
      fac = grow2;
      for (std::size_t j = 1; j != i; ++j) {
         (*a)(i, j) = ((*a)(i, j-1) * fac - (*a)(i-1, j-1)) / (fac - 1.0);
         fac *= grow2;
         err = std::max(std::abs((*a)(i, j) - (*a)(i, j-1)),
                        std::abs((*a)(i, j) - (*a)(i-1, j-1)));
         if (err <= *error) {
            *error = err;
            result = (*a)(i, j);
         }
      }
      // Return when the error is twice as large as the best so far.
      if (std::abs((*a)(i, i) - (*a)(i-1, i-1)) >= 2 * *error) {
         break;
      }
   }
   return result;
}


// Calculate the gradient using Ridders' method.
// Instead of restricting the multi-variable function to single-variable
// functions and using the above derivative implementation, I just duplicate
// the code.
template<typename _Functor, typename _Vector>
inline
double
gradientRidders(_Functor& f, const _Vector& x, _Vector* gradient,
                const _Vector& scales, ads::SquareMatrix<10, double>* a) {
   const double shrink = 1. / 1.4;
   const double grow2 = 1.4 * 1.4;
   assert(x.size() == gradient->size());
   assert(x.size() == scales.size());

   double err, e, fac, fp, fn;
   _Vector y = x;
   double error = 0;
   // For each dimension.
   for (std::size_t n = 0; n != x.size(); ++n) {
      assert(scales[n] != 0);
      double h = scales[n];
      // Calculate the derivative in the n_th dimension.
      (*gradient)[n] = std::numeric_limits<double>::max();
      err = std::numeric_limits<double>::max();
      y[n] = x[n] + h;
      fp = f(y);
      y[n] = x[n] - h;
      fn = f(y);
      (*a)(0, 0) = (fp - fn) / (2.0 * h);
      for (std::size_t i = 1; i != 10; ++i) {
         h *= shrink;
         y[n] = x[n] + h;
         fp = f(y);
         y[n] = x[n] - h;
         fn = f(y);
         (*a)(i, 0) = (fp - fn) / (2.0 * h);
         fac = grow2;
         for (std::size_t j = 1; j <= i; ++j) {
            (*a)(i, j) = ((*a)(i, j-1) * fac - (*a)(i-1, j-1)) / (fac - 1.0);
            fac *= grow2;
            e = std::max(std::abs((*a)(i, j) - (*a)(i, j-1)),
                         std::abs((*a)(i, j) - (*a)(i-1, j-1)));
            if (e <= err) {
               err = e;
               (*gradient)[n] = (*a)(i, j);
            }
         }
         // Return when the error is twice as large as the best so far.
         if (std::abs((*a)(i, i) - (*a)(i-1, i-1)) >= 2 * err) {
            break;
         }
      }
      y[n] = x[n];
      error += err * err;
   }
   error = std::sqrt(error);
   return error;
}


// Calculate the gradient using Ridders' method.
template<typename _Functor, typename _Vector>
inline
double
gradientRidders(_Functor& f, const _Vector& x, _Vector* gradient) {
   // Take a guess at an appropriate length scale over which the function
   // changes appreciably.
   _Vector scales = x;
   for (std::size_t i = 0; i != scales.size(); ++i) {
      scales[i] = std::max(0.001 * std::abs(scales[i]), 0.001);
   }
   return gradientRidders(f, x, gradient, scales);
}


} // namespace numerical

