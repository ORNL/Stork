#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

// Specific includes
#include "Structs/Utility.hpp"
#include "Structs/RDF.hpp"

namespace Stork{
    namespace Run{
        template<typename FloatType, typename memory_space>
        void Trim_RDF_in_Z(Structs::RDF_Dual<FloatType>& RDF) {
            // References
            Structs::RegularGrid_Header<FloatType, host_space>& host_header = RDF.template get_header<host_space>();
            Structs::RegularGrid_Header<FloatType, memory_space>& header = RDF.template get_header<memory_space>();
            Structs::RDF_Data<FloatType, memory_space>& data = RDF.template get_data<memory_space>();

            // Find the minimum used K (corresponding to bottom of domain)
            const uint32_t numEvents = RDF.numEvents;
            uint32_t minK = host_header.local_knum() - 1; // Initialize to maximum possible k

            Kokkos::parallel_reduce(
                "Stork::Trim::Trim_RDF_in_Z(findMinK)",
                Kokkos::RangePolicy<memory_space>(0, numEvents),
                KOKKOS_LAMBDA(const uint32_t n, uint32_t& min_k)  // min_k should be passed by reference
                {
                    // Get point number
                    const uint32_t& LOCAL_p = data.p(n);
                    // Get indices (using original header values)
                    uint32_t LOCAL_ijk[3];
                    header.LOCAL_p_to_LOCAL_ijk(LOCAL_ijk, LOCAL_p);
                    // Perform the min reduction correctly:
                    min_k = Kokkos::min(min_k, LOCAL_ijk[2]); // Use Kokkos::min
                }
            , Kokkos::Min<uint32_t>(minK)); // Use Kokkos::Min reducer

            // Keep the pre-trim header available for decoding old point IDs.
            Structs::RegularGrid_Header<FloatType, memory_space> old_header;
            Kokkos::deep_copy(old_header.index_view, header.index_view);
            Kokkos::deep_copy(old_header.floatType_view, header.floatType_view);

            // Update bounds on the host
            host_header.global_z0() = host_header.global_z0() + host_header.gridResolution() * minK;
            host_header.local_knum() = host_header.local_knum() - minK;

            // Create an new device header on the device and copy the updated header
            Structs::RegularGrid_Header<FloatType, memory_space> new_header;
            Kokkos::deep_copy(new_header.index_view, host_header.index_view);
            Kokkos::deep_copy(new_header.floatType_view, host_header.floatType_view);

            // Update grid numbers (using the UPDATED header)
            Kokkos::parallel_for(
                "Stork::Trim::Trim_RDF_in_Z(updateGridNumbers)",
                Kokkos::RangePolicy<memory_space>(0, numEvents),
                KOKKOS_LAMBDA(const uint32_t n)
                {
                    // Get point number
                    const uint32_t& LOCAL_p = data.p(n);

                    // Get indices using the old_header
                    uint32_t LOCAL_ijk[3];
                    old_header.LOCAL_p_to_LOCAL_ijk(LOCAL_ijk, LOCAL_p);
                    const uint32_t& old_k = LOCAL_ijk[2];

                    // Update K
                    const uint32_t new_K = old_k - minK;
                    const uint32_t new_LOCAL_ijk[3] = {LOCAL_ijk[0], LOCAL_ijk[1], new_K};

                    // Get new point number (using UPDATED header values)
                    uint32_t new_LOCAL_p;
                    new_header.LOCAL_ijk_to_LOCAL_p(new_LOCAL_p, new_LOCAL_ijk);
                    data.p(n) = new_LOCAL_p;
                }
            );

            // Now deep copy within device
            Kokkos::deep_copy(header.index_view, new_header.index_view);
            Kokkos::deep_copy(header.floatType_view, new_header.floatType_view);
        }
    }
}           
