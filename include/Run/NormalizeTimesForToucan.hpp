#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

// Specific includes
#include "Structs/Utility.hpp"
#include "Structs/RDF.hpp"
#include "Structs/SRDF.hpp"

namespace Stork {
    namespace Run {

        // For custom sorting structs (can't put inside class because of CUDA extended lambdas)
        namespace {
            // Custom Struct For RDF Times
            template <typename FloatType>
            struct RDF_SortStruct {
                uint32_t p;
                FloatType t;
                bool tm;

                KOKKOS_INLINE_FUNCTION
                bool operator<(const RDF_SortStruct& other) const {
                    return t < other.t;
                }
            };

            // Custom Struct For SRDF Times
            template <typename FloatType>
            struct SRDF_SortStruct {
                uint32_t p;
                FloatType t;

                KOKKOS_INLINE_FUNCTION
                bool operator<(const SRDF_SortStruct& other) const {
                    return t < other.t;
                }
            };
        } // namespace

        // Normalize times for more accurate grain competition
        template <typename FloatType, typename memory_space>
        void Normalize_RDF_Times(Structs::RDF_Dual<FloatType>& RDF, const FloatType minGap) {

            // References
            Structs::RDF_Data<FloatType, memory_space>& data = RDF.template get_data<memory_space>();
            const uint32_t numEvents = RDF.numEvents;
            const uint32_t numTimes = 2 * numEvents;

            // Return if empty
            if (numEvents == 0) {
                return;
            }

            // Populate Sort Struct
            Kokkos::View<RDF_SortStruct<FloatType>*, layout, memory_space> sortObj(Kokkos::ViewAllocateWithoutInitializing("Stork::NormTimesForToucan::sortObj"), numTimes);

            // Update times
            Kokkos::parallel_for(
                "Stork::Normalize::Norm_Times_For_Toucan(construct sortObj)",
                Kokkos::RangePolicy<memory_space>(0, numEvents),
                KOKKOS_LAMBDA(const uint32_t n) {
                    // Fill Object
                    sortObj(2*n+0).p = n;
                    sortObj(2*n+1).p = n;
                    sortObj(2*n+0).t = data.tm(n);
                    sortObj(2*n+1).t = data.tl(n);
                    sortObj(2*n+0).tm = true;
                    sortObj(2*n+1).tm = false;
                }
            );

            // Sort in time
            Kokkos::sort(sortObj);

            // Scan to adjust the times and remove gaps
            FloatType totalTime;
            Kokkos::parallel_scan(
                "Stork::Normalize::Norm_Times_For_Toucan(adjust Times)",
                Kokkos::RangePolicy<memory_space>(0, numTimes),
                KOKKOS_LAMBDA(const uint32_t n, FloatType& th_time, bool isFinal) {
                    if (n > 0) {
                        const FloatType timeDiff = sortObj(n).t - sortObj(n - 1).t;
                        th_time += (timeDiff > minGap) ? minGap : timeDiff;
                    }

                    if (isFinal) {
                        const uint32_t& p = sortObj(n).p;
                        if (sortObj(n).tm) {
                            data.tm(p) = th_time;
                        }
                        else {
                            data.tl(p) = th_time;
                        }
                    }
                },
                totalTime);

            // Where first group should start
            const FloatType baseShift = -totalTime / static_cast<FloatType>(2.0);
            // const FloatType baseShift = static_cast<FloatType>(0.0);

            // Center times
            Kokkos::parallel_for(
                "Stork::Normalize::Norm_Times_For_Toucan(center Times)",
                Kokkos::RangePolicy<memory_space>(0, numEvents),
                KOKKOS_LAMBDA(const uint32_t n) {
                    data.tm(n) += baseShift;
                    data.tl(n) += baseShift;
                });

            // Fence after for
            Kokkos::fence();
        }

        // Normalize times for more accurate grain competition
        template <typename FloatType, typename memory_space>
        void Normalize_SRDF_Times(Structs::SRDF_Dual<FloatType>& SRDF, const FloatType minGap) {

            // References
            Structs::SRDF_Data<FloatType, memory_space>& data = SRDF.template get_data<memory_space>();
            const uint32_t numSnaps = SRDF.numSnaps;

            // Return if empty
            if (numSnaps == 0) {
                return;
            }

            // Populate Sort Struct
            Kokkos::View<SRDF_SortStruct<FloatType>*, layout, memory_space> sortObj(Kokkos::ViewAllocateWithoutInitializing("Stork::NormTimesForToucan::sortObj"), numSnaps);

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

            // Scan to adjust the times and remove gaps
            FloatType totalTime;
            Kokkos::parallel_scan(
                "Stork::Normalize::Norm_Times_For_Toucan(adjust Times)",
                Kokkos::RangePolicy<memory_space>(0, numSnaps),
                KOKKOS_LAMBDA(const uint32_t n, FloatType& th_time, bool isFinal) {
                    if (n > 0) {
                        const FloatType timeDiff = sortObj(n).t - sortObj(n - 1).t;
                        th_time += (timeDiff > minGap) ? minGap : timeDiff;
                    }

                    if (isFinal) {
                        const uint32_t& p = sortObj(n).p;
                        const FloatType timeShift = th_time - sortObj(n).t;
                        data.t_prev(p) += timeShift;
                        data.t_cur(p) += timeShift;
                    }
                },
                totalTime);

            // Where first group should start
            const FloatType baseShift = -totalTime / static_cast<FloatType>(2.0);
            // const FloatType baseShift = static_cast<FloatType>(0.0);

            // Center times
            Kokkos::parallel_for(
                "Stork::Normalize::Norm_Times_For_Toucan(center Times)",
                Kokkos::RangePolicy<memory_space>(0, numSnaps),
                KOKKOS_LAMBDA(const uint32_t n) {
                    data.t_prev(n) += baseShift;
                    data.t_cur(n) += baseShift;
                });

            // Fence after for
            Kokkos::fence();
        }
    } // namespace Run
} // namespace Stork
