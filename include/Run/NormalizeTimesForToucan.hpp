#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

// Specific includes
// #include "Structs/Utility.hpp"
#include "Structs/SRDF.hpp"

namespace Stork {
    namespace Run {

        // For custom sorting struct
        namespace {
            // Custom Struct For Melting
            template <typename FloatType>
            struct sortStruct {
                uint32_t p;
                FloatType t;

                KOKKOS_INLINE_FUNCTION
                bool operator<(const sortStruct& other) const {
                    return t < other.t;
                }
            };
        } // namespace

        // Interpolate entire file
        template <typename FloatType, typename memory_space>
        void Norm_Times_For_Toucan(Structs::SRDF_Dual<FloatType>& SRDF, const FloatType minGap) {

            // Usings
            using idx_deviceView = Kokkos::View<uint32_t*, layout, memory_space>;
            using floatType_deviceView = Kokkos::View<FloatType*, layout, memory_space>;

            // References
            Structs::SRDF_Data<FloatType, memory_space>& data = SRDF.template get_data<memory_space>();
            const uint32_t numSnaps = SRDF.numSnaps;

            // Populate Sort Struct
            Kokkos::View<sortStruct<FloatType>*, layout, memory_space> sortObj(Kokkos::ViewAllocateWithoutInitializing("Stork::NormTimesForToucan::sortObj"), numSnaps);

            // Update times
            Kokkos::parallel_for(
                "Stork::Normalize::Norm_Times_For_Toucan(construct sortObj)",
                Kokkos::RangePolicy<memory_space>(0, numSnaps),
                KOKKOS_LAMBDA(const uint32_t n) {
                    // Fill Object
                    sortObj(n).p = n;
                    sortObj(n).t = data.t_prev(n);
                });

            // Sort in time
            Kokkos::sort(sortObj);

            // Scan to figure out where the differences are (clusters)
            idx_deviceView clusterID(Kokkos::ViewAllocateWithoutInitializing("Stork::NormTimesForToucan::clusterID"), numSnaps);
            uint32_t numClusters;
            Kokkos::parallel_scan(
                "Stork::Normalize::Norm_Times_For_Toucan(findClusters)",
                Kokkos::RangePolicy<memory_space>(0, numSnaps),
                KOKKOS_LAMBDA(const uint32_t n, uint32_t& th_cluster, bool isFinal) {
                    if (isFinal) {
                        clusterID(n) = th_cluster;
                    }
                    if ((n > 0) && ((sortObj(n).t - sortObj(n - 1).t) > minGap)) {
                        th_cluster++;
                    }
                },
                numClusters);
            ++numClusters;

            // TODO::Reduction probably better
            // Make holders for min and max time of each cluster
            floatType_deviceView minTime(Kokkos::ViewAllocateWithoutInitializing("Stork::NormTimesForToucan::minTime"), numClusters);
            floatType_deviceView maxTime(Kokkos::ViewAllocateWithoutInitializing("Stork::NormTimesForToucan::maxTime"), numClusters);
            Kokkos::deep_copy(minTime, std::numeric_limits<FloatType>::max());
            Kokkos::deep_copy(maxTime, -std::numeric_limits<FloatType>::max());
            Kokkos::fence();
            Kokkos::parallel_for(
                "Stork::Normalize::Norm_Times_For_Toucan(Find Min/Max)",
                Kokkos::RangePolicy<memory_space>(0, numSnaps),
                KOKKOS_LAMBDA(const uint32_t n) {
                    // Find cluster ID
                    const uint32_t& cluster = clusterID(n);
                    const FloatType& t = sortObj(n).t;
                    // Reduce specific cluster
                    if (t < minTime(cluster)) {
                        Kokkos::atomic_min(&minTime(cluster), t);
                    }
                    if (t > maxTime(cluster)) {
                        Kokkos::atomic_max(&maxTime(cluster), t);
                    }
                });

            // Now that we have clusters, find cumulative sum of ranges
            floatType_deviceView partial_sum_ranges(Kokkos::ViewAllocateWithoutInitializing("Stork::NormTimesForToucan::sumRanges"), numClusters);
            FloatType total_range_sum;
            Kokkos::parallel_scan(
                "Stork::Normalize::Norm_Times_For_Toucan(Find Range)",
                Kokkos::RangePolicy<memory_space>(0, numSnaps),
                KOKKOS_LAMBDA(const uint32_t n, FloatType& th_range, bool isFinal) {
                    if ((n == 0) || (clusterID(n) != clusterID(n - 1))) {
                        if (isFinal) {
                            partial_sum_ranges(clusterID(n)) = th_range;
                        }
                        th_range += (maxTime(clusterID(n)) - minTime(clusterID(n)));
                    }
                },
                total_range_sum);

            // Where first group should start
            const FloatType baseShift = -(total_range_sum + numClusters * minGap) / static_cast<FloatType>(2.0);
            // const FloatType baseShift = static_cast<FloatType>(0.0);

            // Update times
            Kokkos::parallel_for(
                "Stork::Normalize::Norm_Times_For_Toucan(Shift Times)",
                Kokkos::RangePolicy<memory_space>(0, numSnaps),
                KOKKOS_LAMBDA(const uint32_t p) {
                    // Get base position and clusterID
                    const uint32_t& n = sortObj(p).p;
                    const uint32_t& cluster = clusterID(p);

                    // Shift times
                    data.t_prev(n) = (data.t_prev(n) - minTime(cluster)) + baseShift + (partial_sum_ranges(cluster) - partial_sum_ranges(0)) + minGap * cluster;
                    data.t_cur(n) = (data.t_cur(n) - minTime(cluster)) + baseShift + (partial_sum_ranges(cluster) - partial_sum_ranges(0)) + minGap * cluster;
                });

            // Fence after for
            Kokkos::fence();
        }
    } // namespace Run
} // namespace Stork
