// -*- C++ -*-

#if !defined(__geom_shape_procrustes_ipp__)
#error This file is an implementation detail of procrustes.
#endif

namespace geom {


template<typename _T, std::size_t _D>
inline
void
procrustes(std::vector<std::tr1::array<_T, _D> >* source,
           std::vector<std::tr1::array<_T, _D> >* target,
           std::tr1::array<_T, _D>* sourceCentroid,
           std::tr1::array<_T, _D>* targetCentroid,
           array::EquilateralArray<_T, 2, _D>* rotation,
           _T* scale) {
   typedef std::tr1::array<_T, _D> Point;
            
   // There must be at least as many points as the space dimension.
   assert(source->size() >= _D);
   // The number of landmark points must be the same in the source and target->
   assert(source->size() == target->size());

   // Center the two sets of points.
   {
      // The source points.
      std::fill(sourceCentroid->begin(), sourceCentroid->end(), 0);
      for (std::size_t i = 0; i != source->size(); ++i) {
         *sourceCentroid += (*source)[i];
      }
      *sourceCentroid /= _T(source->size());
      for (std::size_t i = 0; i != source->size(); ++i) {
         (*source)[i] -= *sourceCentroid;
      }

      // The target points.
      std::fill(targetCentroid->begin(), targetCentroid->end(), 0);
      for (std::size_t i = 0; i != target->size(); ++i) {
         *targetCentroid += (*target)[i];
      }
      *targetCentroid /= _T(target->size());
      for (std::size_t i = 0; i != target->size(); ++i) {
         (*target)[i] -= *targetCentroid;
      }
   }

   // Determine the rotation.
   Eigen::Matrix<_T, Eigen::Dynamic, _D> x1(source->size(), _D),
      x2(target->size(), _D);
   for (std::size_t row = 0; row != source->size(); ++row) {
      for (std::size_t col = 0; col != _D; ++col) {
         x1(row, col) = (*source)[row][col];
         x2(row, col) = (*target)[row][col];
      }
   }
   typedef Eigen::Matrix<_T, _D, _D> SquareMatrix;
   Eigen::SVD<SquareMatrix> svd(_D, _D);
   SquareMatrix product = x2.transpose() * x1;
   svd.compute(product);
   SquareMatrix rot = svd.matrixV() * svd.matrixU().transpose();

   // Determine the scale.
   SquareMatrix num = product * rot;
   SquareMatrix den = x1.transpose() * x1;
   *scale = num.trace() / den.trace();
   assert(*scale > 0);

   // Copy to the output rotation matrix.
   for (std::size_t i = 0; i != _D; ++i) {
      for (std::size_t j = 0; j != _D; ++j) {
         (*rotation)(i, j) = rot(i, j);
      }
   }
}


} // namespace geom
