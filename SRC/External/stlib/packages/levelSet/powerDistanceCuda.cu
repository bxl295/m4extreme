/* -*- C++ -*- */

#include "powerDistanceCuda.h"
#include "../cuda/check.h"
#include "cuda.h"
#ifndef __CUDA_ARCH__
#include "gridCuda.h"
#include "ballsCuda.h"
#endif

#include <cassert>

namespace levelSet {


__global__
void
powerDistanceKernel(float* patch,
                    const float3 lowerCorner,
                    const float spacing,
                    const std::size_t numBalls,
                    const float4* ballsSquared) {
   // Calculate the Cartesian location of the grid point.
   const float3 p = {lowerCorner.x + spacing * threadIdx.x,
                     lowerCorner.y + spacing * threadIdx.y,
                     lowerCorner.z + spacing * threadIdx.z};
   // Initialize the grid value to infinity.
   float g = 1./0;
   // Find the minimum power distance over the set of balls.
   float d;
   float4 b;
   for (std::size_t i = 0; i != numBalls; ++i) {
      b = ballsSquared[i];
      d = (p.x - b.x) * (p.x - b.x) + (p.y - b.y) * (p.y - b.y) +
         (p.z - b.z) * (p.z - b.z) - b.w;
      if (d < g) {
         g = d;
      }
   }
   // Record the patch value. Convert the multi-index to a single offset.
   patch[threadIdx.x + threadIdx.y * blockDim.x +
         threadIdx.z * blockDim.x * blockDim.y] = g;
}


#ifndef __CUDA_ARCH__
// Translate the data to CUDA format and call powerDistanceKernel().
void
powerDistanceCuda(array::EquilateralArray<float, 3, PatchExtent>* patch,
                  const std::tr1::array<float, 3>& lowerCorner,
                  const float spacing,
                  const std::vector<geom::BallSquared<float, 3> >& balls) {
   // Allocate memory for the patch.
   float* patchDev;
   CUDA_CHECK(cudaMalloc((void**)&patchDev, patch->size() * sizeof(float)));

   // Allocate memory for the balls and copy the memory.
   float4* ballsDev;
   CUDA_CHECK(cudaMalloc((void**)&ballsDev, balls.size() * sizeof(float4)));
   {
      std::vector<float4> buffer(balls.size());
      for (std::size_t i = 0; i != balls.size(); ++i) {
         buffer[i].x = balls[i].center[0];
         buffer[i].y = balls[i].center[1];
         buffer[i].z = balls[i].center[2];
         buffer[i].w = balls[i].squaredRadius;
      }
      CUDA_CHECK(cudaMemcpy(ballsDev, &buffer[0],
                            buffer.size() * sizeof(float4),
                            cudaMemcpyHostToDevice));
   }

   const dim3 ThreadsPerBlock(PatchExtent, PatchExtent, PatchExtent);
   const float3 lower = {lowerCorner[0], lowerCorner[1], lowerCorner[2]};
   powerDistanceKernel<<<1,ThreadsPerBlock>>>
      (patchDev, lower, spacing, balls.size(), ballsDev);

   CUDA_CHECK(cudaMemcpy(&(*patch)[0], patchDev, patch->size() * sizeof(float),
                         cudaMemcpyDeviceToHost));
   CUDA_CHECK(cudaFree(patchDev));
   CUDA_CHECK(cudaFree(ballsDev));

#if 0
   for (std::size_t i = 0; i != patch->size(); ++i) {
      assert((*patch)[i] == 0);
   }
#else
   std::cout << *patch << '\n';
#endif
}
#endif


// Note: Storing ballsSquared in constant memory does not significantly affect
// performance.
__device__
void
powerDistanceKernel(float* patch,
                    const float3 lowerCorner,
                    const float spacing,
                    const float4* ballsSquared,
                    std::size_t numInfluencingBalls,
                    const std::size_t* ballIndices) {
   // Calculate the Cartesian location of the grid point.
   const float3 p = {lowerCorner.x + spacing * threadIdx.x,
                     lowerCorner.y + spacing * threadIdx.y,
                     lowerCorner.z + spacing * threadIdx.z};
   // Initialize the grid value to infinity.
   float g = 1./0;
   // Find the minimum power distance over the set of balls.
   float d;
#if 0
   float4 b;
   for (std::size_t i = 0; i != numInfluencingBalls; ++i) {
      b = ballsSquared[ballIndices[i]];
      d = (p.x - b.x) * (p.x - b.x) + (p.y - b.y) * (p.y - b.y) +
         (p.z - b.z) * (p.z - b.z) - b.w;
      if (d < g) {
         g = d;
      }
   }
   // Record the patch value. Convert the multi-index to a single offset.
   patch[threadIdx.x + threadIdx.y * blockDim.x +
         threadIdx.z * blockDim.x * blockDim.y] = g;
#else
   // Using shared memory improves performance on the 9600M GT, but not on
   // the 9400M.
   __shared__ float4 b[32];
   // Convert the multi-index to a single offset.
   const std::size_t tid = threadIdx.x + threadIdx.y * blockDim.x +
      threadIdx.z * blockDim.x * blockDim.y;
   std::size_t block;
   while (numInfluencingBalls != 0) {
      block = 32;
      if (numInfluencingBalls < 32) {
         block = numInfluencingBalls;
      }
      // Load the balls into shared memory.
      if (tid < block) {
         // Even though these are uncoalesced, this is an insignificant cost.
         b[tid] = ballsSquared[ballIndices[tid]];
      }
      __syncthreads();
      // Process the loaded balls.
      for (std::size_t i = 0; i != block; ++i) {
         // This accounts for about a third of the execution time.
         d = (p.x - b[i].x) * (p.x - b[i].x) + (p.y - b[i].y) * (p.y - b[i].y) +
            (p.z - b[i].z) * (p.z - b[i].z) - b[i].w;
         if (d < g) {
            g = d;
         }
      }
      ballIndices += block;
      numInfluencingBalls -= block;
   }
   // Record the patch value. 
   // Note: This statement accounts for about 35% of execution time.
   patch[tid] = g;
#endif
}


// Select a single patch using the block index. Then call a kernel for that
// patch.
__global__
void
powerDistanceKernel(const std::size_t numRefined,
                    float* patches,
                    const uint3* indices,
                    const float3 lowerCorner,
                    const float spacing,
                    const float4* ballsSquared,
                    const std::size_t* ballIndexOffsets,
                    const std::size_t* packedBallIndices) {
   // Convert the 2-D block index into a single patch index.
   const std::size_t i = blockIdx.x + blockIdx.y * gridDim.x;
   if (i >= numRefined) {
      return;
   }
   const std::size_t begin = ballIndexOffsets[i];
   const std::size_t end = ballIndexOffsets[i+1];
   const std::size_t NumThreads = PatchExtent * PatchExtent * PatchExtent;
   const float3 patchLowerCorner =
      {lowerCorner.x + PatchExtent * indices[i].x * spacing,
       lowerCorner.y + PatchExtent * indices[i].y * spacing,
       lowerCorner.z + PatchExtent * indices[i].z * spacing};
   powerDistanceKernel(patches + i * NumThreads,
                       patchLowerCorner, spacing, ballsSquared,
                       end - begin, &packedBallIndices[begin]);
}


#ifndef __CUDA_ARCH__
// Construct a level set for the power distance to a set of balls.
void
negativePowerDistanceCuda(Grid<float, 3, PatchExtent>* grid,
                          const std::vector<geom::Ball<float, 3> >& balls) {
   // Dispense with the trivial case.
   if (grid->empty()) {
      return;
   }

   // Determine the patch/ball dependencies.
   array::StaticArrayOfArrays<std::size_t> dependencies;
   {
      // 1.1 * (diagonal length of a voxel)
      const float offset = 1.1 * grid->spacing * std::sqrt(float(3));
      // Offset the ball radii to include the volume of calculated distance.
      std::vector<geom::Ball<float, 3> > offsetBalls(balls);
      for (std::size_t i = 0; i != offsetBalls.size(); ++i) {
         offsetBalls[i].radius += offset;
      }
      // Calculate the dependencies.
      patchDependencies(*grid, offsetBalls.begin(), offsetBalls.end(),
                        &dependencies);
   }

   // Refine the appropriate patches and set the rest to have an unknown
   // distance.
   grid->refine(dependencies);
   const std::size_t numRefined = grid->numRefined();

   // Allocate device memory for the refined patches and their indices.
   float* patchesDev;
   uint3* indicesDev;
   allocateGridCuda(*grid, numRefined, &patchesDev, &indicesDev);

   // Allocate device memory for the balls and copy the memory.
   float4* ballsSquaredDev;
   CUDA_CHECK(cudaMalloc((void**)&ballsSquaredDev,
                         balls.size() * sizeof(float4)));
   {
      std::vector<float4> buffer(balls.size());
      for (std::size_t i = 0; i != balls.size(); ++i) {
         buffer[i].x = balls[i].center[0];
         buffer[i].y = balls[i].center[1];
         buffer[i].z = balls[i].center[2];
         buffer[i].w = balls[i].radius * balls[i].radius;
      }
      CUDA_CHECK(cudaMemcpy(ballsSquaredDev, &buffer[0],
                            buffer.size() * sizeof(float4),
                            cudaMemcpyHostToDevice));
   }

   // Allocate device memory for the ball index offsets and packed ball indices.
   // Copy the data to the device.
   std::size_t* ballIndexOffsetsDev;
   std::size_t* packedBallIndicesDev;
   allocateBallIndicesCuda(dependencies, numRefined, &ballIndexOffsetsDev,
                           &packedBallIndicesDev);
   
   const float3 lowerCorner = {grid->lowerCorner[0], grid->lowerCorner[1],
                               grid->lowerCorner[2]};
   // Use a 2-D grid of blocks. Because the number of refined patches may 
   // exceed the maximum allowed single grid dimension.
   cudaDeviceProp prop;
   CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
   const std::size_t Len = prop.maxGridSize[0];
   const dim3 GridDim(std::min(numRefined, Len), (numRefined + Len - 1) / Len);
   // A thread for each patch grid point.
   const dim3 ThreadsPerBlock(PatchExtent, PatchExtent, PatchExtent);
   // Launch the kernel.
   powerDistanceKernel<<<GridDim,ThreadsPerBlock>>>
      (numRefined, patchesDev, indicesDev, lowerCorner, grid->spacing,
       ballsSquaredDev, ballIndexOffsetsDev, packedBallIndicesDev);
   // Copy the patch data back to the host.
   CUDA_CHECK(cudaMemcpy(grid->data(), patchesDev,
                         grid->numVertices() * sizeof(float),
                         cudaMemcpyDeviceToHost));

   // Free the device memory.
   CUDA_CHECK(cudaFree(patchesDev));
   CUDA_CHECK(cudaFree(indicesDev));
   CUDA_CHECK(cudaFree(ballsSquaredDev));
   CUDA_CHECK(cudaFree(ballIndexOffsetsDev));
   CUDA_CHECK(cudaFree(packedBallIndicesDev));
}
#endif

} // namespace levelSet
