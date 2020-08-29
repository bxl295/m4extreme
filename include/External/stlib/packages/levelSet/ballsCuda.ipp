/* -*- C++ -*- */

#if !defined(__levelSet_ballsCuda_ipp__)
#error This file is an implementation detail of ballsCuda.
#endif

namespace levelSet {


inline
void
allocateBallIndicesCuda
(const array::StaticArrayOfArrays<std::size_t>& dependencies,
 const std::size_t numRefined,
 std::size_t** ballIndexOffsetsDev,
 std::size_t** packedBallIndicesDev) {
   // Allocate device memory for the ball index offsets.
   CUDA_CHECK(cudaMalloc((void**)ballIndexOffsetsDev,
                         (numRefined + 1) * sizeof(std::size_t)));
   {
      // In "dependencies", we store pointers that determine the beginning
      // and end of each array. We need to convert these to offsets by 
      // subtracting the address of the beginning of the storage.
      std::vector<std::size_t> buffer(numRefined + 1);
      buffer[0] = 0;
      std::size_t n = 0;
      for (std::size_t i = 0; i != dependencies.getNumberOfArrays(); ++i) {
         if (! dependencies.empty(i)) {
            buffer[n+1] = dependencies.end(i) - dependencies.begin();
            ++n;
         }
      }
      assert(n == numRefined);
      // Copy to the device.
      CUDA_CHECK(cudaMemcpy(*ballIndexOffsetsDev, &buffer[0],
                            buffer.size() * sizeof(std::size_t),
                            cudaMemcpyHostToDevice));
   }

   // Allocate device memory for the packed ball indices.
   CUDA_CHECK(cudaMalloc((void**)packedBallIndicesDev,
                         dependencies.size() * sizeof(std::size_t)));
   // Copy to the device.
   CUDA_CHECK(cudaMemcpy(*packedBallIndicesDev, &dependencies[0],
                         dependencies.size() * sizeof(std::size_t),
                         cudaMemcpyHostToDevice));
}


} // namespace levelSet
