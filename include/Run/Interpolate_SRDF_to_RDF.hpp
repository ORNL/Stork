#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

// Specific includes
#include "Structs/Utility.hpp"
#include "Structs/RDF.hpp"
#include "Structs/SRDF.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#include <Kokkos_StdAlgorithms.hpp>

namespace Stork{
    namespace Run{
        
        // For custom structs
        namespace impl{
            // Custom Struct For Melting
            template<typename FloatType>
            struct meltStruct{
                int p;
                FloatType tm;

                KOKKOS_INLINE_FUNCTION
                bool operator<(const meltStruct& other) const {
                    if (p == other.p) {
                        return tm < other.tm;
                    }
                    return p < other.p;
                }
            };

            // Custom Struct For Solidifying
            template<typename FloatType>
            struct solStruct{
                int p;
                FloatType tl, cr;

                KOKKOS_INLINE_FUNCTION
                bool operator<(const solStruct& other) const {
                    if (p == other.p) {
                        return tl < other.tl;
                    }
                    return p < other.p;
                }
            };

            // Kokkos views
            template<typename FloatType, typename memory_space>
            using meltStruct_View = Kokkos::View<meltStruct<FloatType>*, layout, memory_space>;

            template<typename FloatType, typename memory_space>
            using solStruct_View = Kokkos::View<solStruct<FloatType>*, layout, memory_space>;
        }

        namespace impl{
            // Define this structure outside the function, or within an appropriate scope
            struct MeltSolEventCount {
                uint32_t melt_count = 0;
                uint32_t sol_count = 0;

                // Operator += (required for reduction/scan: defines how two values are combined)
                KOKKOS_INLINE_FUNCTION
                MeltSolEventCount& operator+=(const MeltSolEventCount& rhs) {
                    melt_count += rhs.melt_count;
                    sol_count += rhs.sol_count;
                    return *this;
                }
            };

            KOKKOS_INLINE_FUNCTION
            MeltSolEventCount operator+(const MeltSolEventCount& lhs, const MeltSolEventCount& rhs) {
                // Create a copy of the left-hand side
                MeltSolEventCount result = lhs;
                // Use the existing operator+= to add the right-hand side to the copy
                result += rhs;
                // Return the new combined object
                return result;
            }

            template<typename memory_space>
            using eventCount_View = Kokkos::View<MeltSolEventCount*, layout, memory_space>;
        }
        
        // For trilinearInterpolation
        namespace impl{
    
            template<typename FloatType>
            KOKKOS_INLINE_FUNCTION
            FloatType trilinearInterpolation(const FloatType* T_arr, const FloatType xd, const FloatType yd, const FloatType zd) {

                // Interpolate along z axis
                const FloatType T00 = T_arr[0] * (static_cast<FloatType>(1.0) - zd) + T_arr[1] * zd;
                const FloatType T01 = T_arr[2] * (static_cast<FloatType>(1.0) - zd) + T_arr[3] * zd;
                const FloatType T10 = T_arr[4] * (static_cast<FloatType>(1.0) - zd) + T_arr[5] * zd;
                const FloatType T11 = T_arr[6] * (static_cast<FloatType>(1.0) - zd) + T_arr[7] * zd;

                // Interpolate along y axis
                const FloatType T0 = T00 * (static_cast<FloatType>(1.0) - yd) + T01 * yd;
                const FloatType T1 = T10 * (static_cast<FloatType>(1.0) - yd) + T11 * yd;

                // Interpolate along x axis
                const FloatType T = T0 * (static_cast<FloatType>(1.0) - xd) + T1 * xd;

                // Return interpolated temperature
                return T;
            }
        }
        
        // For Triming Layers
        namespace impl{    
            // Trim layers according to spatial bounds
            template<typename FloatType, typename memory_space>
            void Trim_Layer(Structs::SRDF_Dual<FloatType>& SRDF, const Structs::IndexBounds<memory_space>& Bounds){

                // References
                Structs::RegularGrid_Header<FloatType, memory_space>& header = SRDF.template get_header<memory_space>();
                Structs::SRDF_Data<FloatType, memory_space>& data = SRDF.template get_data<memory_space>();

                // Set quick access to views
                using indexType_spaceView = Kokkos::View<uint32_t*, layout, memory_space>; 
                using floatType_spaceView = Kokkos::View<FloatType*, layout, memory_space>;

                // Only include events within the bounds
                const uint32_t numEvents = SRDF.numSnaps;
                uint32_t trimSize;
                Kokkos::parallel_reduce(
                "Stork::Interpolate::TrimLayer(getSize)",
                Kokkos::RangePolicy<host_space>(0, numEvents),
                KOKKOS_LAMBDA(const uint32_t n, uint32_t& th_numTrimEvents)
                    {
                        // Get point number
                        const uint32_t& LOCAL_p = data.p(n);
                        // Get indices
                        uint32_t LOCAL_ijk[3];
                        header.LOCAL_p_to_LOCAL_ijk(LOCAL_ijk,LOCAL_p);
                        // Check bounds
                        if (LOCAL_ijk[0]>=Bounds.imin() && LOCAL_ijk[0]<=Bounds.imax() && LOCAL_ijk[1]>=Bounds.jmin() && LOCAL_ijk[1]<=Bounds.jmin()){
                            // Increment number of events
                            th_numTrimEvents++;
                        }
                    }
                , trimSize);

                // Initialize New Views based on size (in memory space)
                indexType_spaceView trimmed_cellNum_view("trimmed_cellNum_view", trimSize);
                floatType_spaceView trimmed_times_view("trimmed_times_view", 2 * trimSize);
                floatType_spaceView trimmed_thermals_view("trimmed_thermals_view", 16 * trimSize);

                // Place the correct info into the new containers
                Kokkos::parallel_scan(
                "Stork::Interpolate::TrimLayer(placeInfo)",
                Kokkos::RangePolicy<host_space>(0, numEvents),
                KOKKOS_LAMBDA(const uint32_t n, uint32_t& th_numTrimEvents, bool isFinal)
                    {
                        // Get point number
                        const uint32_t& LOCAL_p = data.p(n);
                        // Get indices
                        uint32_t LOCAL_ijk[3];
                        header.LOCAL_p_to_LOCAL_ijk(LOCAL_ijk,LOCAL_p);
                        // Check bounds
                        if (LOCAL_ijk[0]>=Bounds.imin() && LOCAL_ijk[0]<=Bounds.imax() && LOCAL_ijk[1]>=Bounds.jmin() && LOCAL_ijk[1]<=Bounds.jmin()){
                            if (isFinal){
                                // Place cell num info
                                trimmed_cellNum_view(th_numTrimEvents) = data.p(n);
                                // Place times
                                trimmed_times_view(2*th_numTrimEvents+0) = data.times_view(2*n+0);
                                trimmed_times_view(2*th_numTrimEvents+1) = data.times_view(2*n+1);
                                // Place thermals
                                for (uint32_t dn=0; dn<16; dn++){
                                    trimmed_thermals_view(16*th_numTrimEvents+dn) = data.thermals_view(16*n+dn);
                                }   
                            }
                            // Increment number of events
                            th_numTrimEvents++;
                        }
                    }
                , trimSize);

                // Replace the original views with the trimmed views
                SRDF.cellNum_view = trimmed_cellNum_view;
                SRDF.times_view = trimmed_times_view;
                SRDF.thermals_view = trimmed_thermals_view;
            }
        }
        
        // For SRDF_to_meltData
        namespace impl{
            
            template<typename FloatType_In, typename FloatType_Out, typename target_space>
            void init_fineRDF_header(Structs::SRDF_Dual<FloatType_In>& SRDF, Structs::RDF_Dual<FloatType_Out>& RDF, const uint8_t fineFactor){
                // Get host headers
                Structs::RegularGrid_Header<FloatType_In, host_space>& SRDF_header = SRDF.template get_header<host_space>();
                Structs::RegularGrid_Header<FloatType_Out, host_space>& RDF_header = RDF.template get_header<host_space>();
                // Decrement resolution
                RDF_header.gridResolution() = SRDF_header.gridResolution()/fineFactor;
                // Init header to have the same origin point
                RDF_header.global_x0() = SRDF_header.global_x0();
                RDF_header.global_y0() = SRDF_header.global_y0();
                RDF_header.global_z0() = SRDF_header.global_z0();
                // Init headers to have right number of global start indices
                RDF_header.global_i0() = fineFactor*SRDF_header.global_i0();
                RDF_header.global_j0() = fineFactor*SRDF_header.global_j0();
                RDF_header.global_k0() = fineFactor*SRDF_header.global_k0();
                // Init headers to have right number of local extents
                RDF_header.local_inum() = fineFactor*(SRDF_header.local_inum()-1)+1;
                RDF_header.local_jnum() = fineFactor*(SRDF_header.local_jnum()-1)+1;
                RDF_header.local_knum() = fineFactor*(SRDF_header.local_knum()-1)+1;
                // Copy to device (if different spaces)
                RDF.template Copy_Header<host_space, target_space>();    
            }

            template<typename FloatType_In, typename FloatType_Out, typename memory_space>
            void SRDF_to_meltData(Structs::SRDF_Dual<FloatType_In>& SRDF, Structs::RDF_Dual<FloatType_Out>& RDF, impl::meltStruct_View<FloatType_In, memory_space>& meltData, impl::solStruct_View<FloatType_In, memory_space>& solData, const uint8_t fineFactor_uint8){
                
                // References
                const Structs::RegularGrid_Header<FloatType_Out, memory_space>& RDF_header = RDF.template get_header<memory_space>();
                const Structs::RegularGrid_Header<FloatType_In, memory_space>& SRDF_header = SRDF.template get_header<memory_space>();
                const Structs::SRDF_Data<FloatType_In, memory_space>& SRDF_data = SRDF.template get_data<memory_space>();
                
                // Make better fineFactors for numerics
                const uint32_t fineFactor_uint32 = static_cast<uint32_t>(fineFactor_uint8);
                const FloatType_In fineFactor_floatType = static_cast<FloatType_In>(fineFactor_uint8);
                const FloatType_In invFineFactor = static_cast<FloatType_In>(1.0) / fineFactor_floatType;
                const FloatType_In T_liq = SRDF.T_critical;

                // Cache host-known grid constants once so kernels capture plain scalars
                const uint32_t coarse_inum = SRDF.host_header.local_inum();
                const uint32_t coarse_jnum = SRDF.host_header.local_jnum();
                const uint32_t coarse_knum = SRDF.host_header.local_knum();
                const uint32_t coarse_jk_stride = coarse_jnum * coarse_knum;
                const uint32_t fine_jnum = RDF.host_header.local_jnum();
                const uint32_t fine_knum = RDF.host_header.local_knum();
                const uint32_t fine_jk_stride = fine_jnum * fine_knum;
                
                // Make the scan object container and find events per data
                const uint32_t dataSize = SRDF.numSnaps;
                eventCount_View<memory_space> per_element_counts("per_element_counts", dataSize);
                Kokkos::parallel_for(
                    "Stork::Interpolate::countEvents",
                    Kokkos::RangePolicy<memory_space>(0, dataSize),
                    KOKKOS_LAMBDA(const uint32_t n) {
                        const uint32_t local_p = SRDF_data.p(n);
                        const uint32_t i = local_p / coarse_jk_stride;
                        const uint32_t rem = local_p - i * coarse_jk_stride;
                        const uint32_t j = rem / coarse_knum;
                        const uint32_t k = rem - j * coarse_knum;
                        const uint32_t di_limit = fineFactor_uint32 + static_cast<uint32_t>(i == coarse_inum - 2);
                        const uint32_t dj_limit = fineFactor_uint32 + static_cast<uint32_t>(j == coarse_jnum - 2);
                        const uint32_t dk_limit = fineFactor_uint32 + static_cast<uint32_t>(k == coarse_knum - 2);

                        FloatType_In T_prev_arr[8];
                        FloatType_In T_cur_arr[8];
                        const uint32_t thermalOffset = 16 * n;
                        for (uint32_t corner = 0; corner < 8; ++corner) {
                            T_prev_arr[corner] = SRDF_data.thermals_view(thermalOffset + corner);
                            T_cur_arr[corner] = SRDF_data.thermals_view(thermalOffset + 8 + corner);
                        }

                        MeltSolEventCount counts;
                        for (uint32_t di = 0; di < di_limit; ++di) {
                            const FloatType_In dx = static_cast<FloatType_In>(di) * invFineFactor;
                            for (uint32_t dj = 0; dj < dj_limit; ++dj) {
                                const FloatType_In dy = static_cast<FloatType_In>(dj) * invFineFactor;
                                for (uint32_t dk = 0; dk < dk_limit; ++dk) {
                                    const FloatType_In dz = static_cast<FloatType_In>(dk) * invFineFactor;
                                    const FloatType_In T_prev = impl::trilinearInterpolation<FloatType_In>(T_prev_arr,dx,dy,dz);
                                    const FloatType_In T_cur = impl::trilinearInterpolation<FloatType_In>(T_cur_arr,dx,dy,dz);
                                    counts.melt_count += static_cast<uint32_t>((!(T_prev >= T_liq)) && (T_cur >= T_liq));
                                    counts.sol_count += static_cast<uint32_t>((T_prev >= T_liq) && !(T_cur >= T_liq));
                                }
                            }
                        }
                        per_element_counts(n) = counts;
                    }
                );
                
                // Calculate scan offsets for each data (using inclusive scan)
                eventCount_View<memory_space> inclusive_scan_offsets("scan_offsets", dataSize);
                Kokkos::View<MeltSolEventCount, layout, host_space> total_events_h("total_events_h");
                Kokkos::Experimental::inclusive_scan(
                    "Stork::Interpolate::computeInclusiveOffsets", // The label
                    memory_space{},                                // An instance of the execution space
                    per_element_counts,                            // The source view
                    inclusive_scan_offsets                         // The destination view
                );
                auto last_element_idx = dataSize - 1;
                Kokkos::deep_copy(total_events_h, Kokkos::subview(inclusive_scan_offsets, last_element_idx));

                // Throw error if not the same number of events. This means something is wrong with the thermal data.
                const uint32_t numMelt = total_events_h().melt_count;
                const uint32_t numSol = total_events_h().sol_count;
                if (numMelt != numSol){
                    // TODO::DEBUG
                    std::cout << "ERROR" << "\n";
                    std::cout << "\tNum: Melt/Sol\t" << numMelt << " | " << numSol << "\n";
                    std::cout << "\tOrigin: " << RDF.host_header.global_x0() << " | " << RDF.host_header.global_y0() << " | " << RDF.host_header.global_z0() << "\n";
                    std::cout << "\tOffset: " << RDF.host_header.global_i0() << " | " << RDF.host_header.global_j0() << " | " << RDF.host_header.global_k0() << "\n";
                    std::cout << "\tNums: " << RDF.host_header.local_inum() << " | " << RDF.host_header.local_jnum() << " | " << RDF.host_header.local_knum() << "\n";
                    std::cout << "\tResolution: " << RDF.host_header.gridResolution() << std::endl;

                    // Make new views
                    uint32_deviceView p_melt_device("DEBUG::Melt_device", RDF.host_header.local_inum()*RDF.host_header.local_jnum()*RDF.host_header.local_knum());
                    uint32_deviceView p_sol_device("DEBUG::Sol_device", RDF.host_header.local_inum()*RDF.host_header.local_jnum()*RDF.host_header.local_knum());
                    uint32_hostView p_melt_host("DEBUG::Melt_host", RDF.host_header.local_inum()*RDF.host_header.local_jnum()*RDF.host_header.local_knum());
                    uint32_hostView p_sol_host("DEBUG::Sol_host", RDF.host_header.local_inum()*RDF.host_header.local_jnum()*RDF.host_header.local_knum());
                    // Find number of events first
                    Kokkos::parallel_for(
                    "Stork::Interpolate::SRDF_to_meltData(DEBUG)",
                    Kokkos::RangePolicy<memory_space>(0, dataSize),
                    KOKKOS_LAMBDA(const uint32_t n)
                        {
                            // Set references
                            const uint32_t& inum = SRDF_header.local_inum();
                            const uint32_t& jnum = SRDF_header.local_jnum();
                            const uint32_t& knum = SRDF_header.local_knum();
                            const FloatType_In& T_liq = SRDF.T_critical;
                            const uint32_t& inum_fine = RDF_header.local_inum();
                            const uint32_t& jnum_fine = RDF_header.local_jnum();
                            const uint32_t& knum_fine = RDF_header.local_knum();
                            // Get i,j,k
                            uint32_t ijk[3];
                            SRDF_header.LOCAL_p_to_LOCAL_ijk(ijk,SRDF_data.p(n));
                            const uint32_t& i = ijk[0];
                            const uint32_t& j = ijk[1];
                            const uint32_t& k = ijk[2];
                            // Determine if global edge (xnum-1 is maximum index of cells, so that's why it's xnum-2)
                            const bool di_edge = (i==inum-2);
                            const bool dj_edge = (j==jnum-2);
                            const bool dk_edge = (k==knum-2);
                            // Loop over possible points
                            const uint32_t steps = (fineFactor_uint32+1);
                            const uint32_t range = steps*steps*steps;
                            for (uint32_t dp=0; dp<range; dp++)
                            {
                                // Get local shifts
                                const uint32_t di = (dp/(steps*steps));
                                const uint32_t dj = (dp/steps)%(steps);
                                const uint32_t dk = (dp%steps); 
                                // Make normalized shifts
                                const FloatType_In dx = di/fineFactor_floatType;
                                const FloatType_In dy = dj/fineFactor_floatType;
                                const FloatType_In dz = dk/fineFactor_floatType;
                                // Determine if not local and global edge
                                const bool cond1 = (di==fineFactor_uint32 && !di_edge);
                                const bool cond2 = (dj==fineFactor_uint32 && !dj_edge);
                                const bool cond3 = (dk==fineFactor_uint32 && !dk_edge);
                                // If not global edge, but local edge -> mask results
                                if (cond1 || cond2 || cond3){continue;}
                                // Interpolate temperatres
                                const FloatType_In T_prev = impl::trilinearInterpolation<FloatType_In>(&SRDF_data.T_prev(n),dx,dy,dz);
                                const FloatType_In T_cur = impl::trilinearInterpolation<FloatType_In>(&SRDF_data.T_cur(n),dx,dy,dz);
                                // Is melting or solidifying?
                                const bool isMelt = (!(T_prev>=T_liq) && (T_cur>=T_liq));
                                const bool isSol = ((T_prev>=T_liq) && !(T_cur>=T_liq));
                                // Calculate new local indices
                                const uint32_t interp_i = i*fineFactor_uint32+di;
                                const uint32_t interp_j = j*fineFactor_uint32+dj;
                                const uint32_t interp_k = k*fineFactor_uint32+dk;
                                const uint32_t interp_p2D = interp_i*jnum_fine*knum_fine+interp_j*knum_fine+interp_k;
                                // If melting or solidifying, increment size
                                if (isMelt){
                                    Kokkos::atomic_add(&p_melt_device(interp_p2D),1);
                                }
                                if (isSol){
                                    Kokkos::atomic_add(&p_sol_device(interp_p2D),1);
                                }
                            }
                        }
                    );

                    // Copy data to host
                    Kokkos::deep_copy(p_melt_host, p_melt_device);
                    Kokkos::deep_copy(p_sol_host, p_sol_device);

                    // Output to file (i,j,numSol,numMelt) [HELP]
                    const Common::string filename = "DEBUG-STORK-NUMEVENTS.csv";
                    std::ofstream outfile(filename);
                    if (!outfile.is_open()) {
                        throw std::runtime_error("Could not open output file: " + filename);
                    }

                    outfile << "i,j,k,numMelt,numSol\n";
                    for (uint32_t p = 0; p < p_melt_host.extent(0); ++p) {
                        // Compute indices
                        const uint32_t i = p / (RDF.host_header.local_knum()*RDF.host_header.local_jnum());
                        const uint32_t j = (p / RDF.host_header.local_knum()) % RDF.host_header.local_jnum();
                        const uint32_t k = (p % RDF.host_header.local_knum());

                        outfile << i << "," << j << "," << k << "," << p_melt_host(p) << "," << p_sol_host(p) << "\n";
                    }
                    outfile.close();

                    throw std::runtime_error("Error with Thermal Input. Number of melting and solidification events are not equal. Interpolation Failed.");
                }

                // Make new views
                meltData = impl::meltStruct_View<FloatType_In, memory_space>(Kokkos::ViewAllocateWithoutInitializing("meltData"), numMelt);
                solData = impl::solStruct_View<FloatType_In, memory_space>(Kokkos::ViewAllocateWithoutInitializing("solStruct"), numSol);
                
                // Now fill with a for
                Kokkos::parallel_for(
                "Stork::Interpolate::SRDF_to_meltData(fillEvents)",
                Kokkos::RangePolicy<memory_space>(0, dataSize),
                KOKKOS_LAMBDA(const uint32_t n)
                    {   
                        const MeltSolEventCount counts = per_element_counts(n);
                        if (counts.melt_count == 0 && counts.sol_count == 0) {
                            return;
                        }

                        const uint32_t local_p = SRDF_data.p(n);
                        const uint32_t i = local_p / coarse_jk_stride;
                        const uint32_t rem = local_p - i * coarse_jk_stride;
                        const uint32_t j = rem / coarse_knum;
                        const uint32_t k = rem - j * coarse_knum;
                        const uint32_t di_limit = fineFactor_uint32 + static_cast<uint32_t>(i == coarse_inum - 2);
                        const uint32_t dj_limit = fineFactor_uint32 + static_cast<uint32_t>(j == coarse_jnum - 2);
                        const uint32_t dk_limit = fineFactor_uint32 + static_cast<uint32_t>(k == coarse_knum - 2);

                        FloatType_In T_prev_arr[8];
                        FloatType_In T_cur_arr[8];
                        const uint32_t thermalOffset = 16 * n;
                        for (uint32_t corner = 0; corner < 8; ++corner) {
                            T_prev_arr[corner] = SRDF_data.thermals_view(thermalOffset + corner);
                            T_cur_arr[corner] = SRDF_data.thermals_view(thermalOffset + 8 + corner);
                        }

                        const FloatType_In t_prev = SRDF_data.t_prev(n);
                        const FloatType_In t_cur = SRDF_data.t_cur(n);
                        const FloatType_In dt = t_cur - t_prev;
                        const FloatType_In inv_dt = static_cast<FloatType_In>(1.0) / dt;
                        const uint32_t fine_base_p =
                            (i * fineFactor_uint32) * fine_jk_stride +
                            (j * fineFactor_uint32) * fine_knum +
                            (k * fineFactor_uint32);

                        uint32_t melt_idx = inclusive_scan_offsets(n).melt_count - per_element_counts(n).melt_count;
                        uint32_t sol_idx  = inclusive_scan_offsets(n).sol_count - per_element_counts(n).sol_count;
                        for (uint32_t di = 0; di < di_limit; ++di) {
                            const FloatType_In dx = static_cast<FloatType_In>(di) * invFineFactor;
                            const uint32_t fine_p_i = fine_base_p + di * fine_jk_stride;
                            for (uint32_t dj = 0; dj < dj_limit; ++dj) {
                                const FloatType_In dy = static_cast<FloatType_In>(dj) * invFineFactor;
                                const uint32_t fine_p_ij = fine_p_i + dj * fine_knum;
                                for (uint32_t dk = 0; dk < dk_limit; ++dk) {
                                    const FloatType_In dz = static_cast<FloatType_In>(dk) * invFineFactor;
                                    const FloatType_In T_prev = impl::trilinearInterpolation<FloatType_In>(T_prev_arr,dx,dy,dz);
                                    const FloatType_In T_cur = impl::trilinearInterpolation<FloatType_In>(T_cur_arr,dx,dy,dz);
                                    const bool isMelt = (!(T_prev >= T_liq)) && (T_cur >= T_liq);
                                    const bool isSol = (T_prev >= T_liq) && !(T_cur >= T_liq);

                                    if (!(isMelt || isSol)) {
                                        continue;
                                    }

                                    const FloatType_In inv_dT = static_cast<FloatType_In>(1.0) / (T_cur - T_prev);
                                    const FloatType_In crossTime = t_prev + (T_liq - T_prev) * dt * inv_dT;
                                    const uint32_t p_interp = fine_p_ij + dk;
                                    if (isMelt) {
                                        meltData(melt_idx).p = p_interp;
                                        meltData(melt_idx).tm = crossTime;
                                        ++melt_idx;
                                    }
                                    if (isSol) {
                                        solData(sol_idx).p = p_interp;
                                        solData(sol_idx).tl = crossTime;
                                        solData(sol_idx).cr = (T_prev - T_cur) * inv_dt;
                                        ++sol_idx;
                                    }
                                }
                            }
                        }
                    }
                );

                // Sort by p and then (tm,tl)  
                Kokkos::sort(meltData);
                Kokkos::sort(solData);
            }
        }
            
        // Interpolate entire file
        template<typename FloatType_In, typename memory_space, typename FloatType_Out, typename target_space>
        Structs::RDF_Dual<FloatType_Out> Interpolate_SRDF_to_RDF(Structs::SRDF_Dual<FloatType_In>& SRDF, const uint8_t fineFactor){
            
            // Make structure
            Structs::RDF_Dual<FloatType_Out> RDF;

            // Initialize RDF header (copy to device if necessary)
            impl::init_fineRDF_header<FloatType_In, FloatType_Out, target_space>(SRDF, RDF, fineFactor);

            // Make SRDF mirrors and copy to device (if necessary)
            SRDF.template Make_Data_Mirrors<memory_space, target_space>();
            SRDF.template Copy_All<memory_space, target_space>();

            // Make containers for returned data
            impl::meltStruct_View<FloatType_In, target_space> meltData;
            impl::solStruct_View<FloatType_In, target_space> solData;

            // Find the number of events, and the sorted melt and sol data
            impl::SRDF_to_meltData<FloatType_In, FloatType_Out, target_space>(SRDF, RDF, meltData, solData, fineFactor);

            // Construct data containers and get view reference
            const uint32_t numEvents = meltData.extent(0);
            RDF.template Make_Data_Views<target_space>(numEvents);
            RDF.numEvents = numEvents;
            Structs::RegularGrid_Header<FloatType_Out, target_space>& RDF_header = RDF.template get_header<target_space>();
            Structs::RDF_Data<FloatType_Out, target_space>& RDF_data = RDF.template get_data<target_space>();

            // Fill data containers
            Kokkos::parallel_for(
            "Stork::Interpolate::SRDF_to_meltData(fillDataContainer)",
            Kokkos::RangePolicy<target_space>(0, numEvents),
            KOKKOS_LAMBDA(const uint32_t n)
                {
                    RDF_data.p(n) = meltData(n).p;
                    RDF_data.tm(n) = meltData(n).tm;
                    RDF_data.tl(n) = solData(n).tl;
                    RDF_data.cr(n) = solData(n).cr;
                }
            );

            // Return data
            return RDF;
        }

        // Interpolate entire file
        template<typename FloatType_In, typename memory_space, typename FloatType_Out, typename target_space>
        Structs::RDF_Dual<FloatType_Out> Interpolate_And_Trim_SRDF_to_RDF(Structs::SRDF_Dual<FloatType_In>& SRDF, const uint8_t fineFactor, const uint32_t GLOBAL_iMin, const uint32_t GLOBAL_iMax, const uint32_t GLOBAL_jMin, const uint32_t GLOBAL_jMax){

            // Make structure
            Structs::RDF_Dual<FloatType_Out> RDF;
            
            // Initialize RDF header (copy to device if necessary)
            impl::init_fineRDF_header<FloatType_In, FloatType_Out, target_space>(SRDF, RDF, fineFactor);

            // Make SRDF mirrors and copy to device (if necessary)
            SRDF.template Make_Data_Mirrors<memory_space, target_space>();
            SRDF.template Copy_All<memory_space, target_space>();

            // Make containers for returned data
            impl::meltStruct_View<FloatType_In, target_space> meltData;
            impl::solStruct_View<FloatType_In, target_space> solData;

            // Find the number of events, and the sorted melt and sol data
            impl::SRDF_to_meltData<FloatType_In, FloatType_Out, target_space>(SRDF, RDF, meltData, solData, fineFactor);

            // Construct data containers and get view reference
            const uint32_t numEvents = meltData.extent(0);
            RDF.template Make_Data_Views<target_space>(numEvents);
            RDF.numEvents = numEvents;
            Structs::RegularGrid_Header<FloatType_Out, target_space>& RDF_header = RDF.template get_header<target_space>();
            Structs::RDF_Data<FloatType_Out, target_space>& RDF_data = RDF.template get_data<target_space>();

            // Keep the untrimmed header for decoding the event list below.
            Structs::RegularGrid_Header<FloatType_Out, target_space> old_header;
            Kokkos::deep_copy(old_header.index_view, RDF_header.index_view);
            Kokkos::deep_copy(old_header.floatType_view, RDF_header.floatType_view);

            // Update header with snipped information
            RDF.host_header.global_i0() = GLOBAL_iMin;
            RDF.host_header.local_inum() = (GLOBAL_iMax-GLOBAL_iMin)+1;
            RDF.host_header.global_j0() = GLOBAL_jMin;
            RDF.host_header.local_jnum() = (GLOBAL_jMax-GLOBAL_jMin)+1; 

            // Create an new device header on the device and copy the updated header
            Structs::RegularGrid_Header<FloatType_Out, target_space> new_header;
            Kokkos::deep_copy(new_header.index_view, RDF.host_header.index_view);
            Kokkos::deep_copy(new_header.floatType_view, RDF.host_header.floatType_view);

            // Construct and fill final data containers (using a scan)
            uint32_t trimEvents;
            Kokkos::parallel_scan(
            "Stork::Interpolate::SRDF_to_meltData(fillDataContainer-withTrim)",
            Kokkos::RangePolicy<target_space>(0, numEvents),
            KOKKOS_LAMBDA(const uint32_t n, uint32_t& th_trimEvents, bool isFinal)
                {
                    // Get global ijk from local p (global)
                    const uint32_t p = meltData(n).p;
                    uint32_t GLOBAL_ijk[3];
                    // USE OLD HEADER
                    old_header.LOCAL_p_to_GLOBAL_ijk(GLOBAL_ijk,p);
                    const uint32_t& GLOBAL_i = GLOBAL_ijk[0];
                    const uint32_t& GLOBAL_j = GLOBAL_ijk[1];
                    const uint32_t& GLOBAL_k = GLOBAL_ijk[2];
                    // Check bounds
                    if (GLOBAL_i>= GLOBAL_iMin && GLOBAL_i<=GLOBAL_iMax && GLOBAL_j>=GLOBAL_jMin && GLOBAL_j<=GLOBAL_jMax){
                        if (isFinal){
                            // Make new p value
                            uint32_t new_LOCAL_p;
                            // USE NEW HEADER
                            new_header.GLOBAL_ijk_to_LOCAL_p(new_LOCAL_p, GLOBAL_ijk);
                            // Store value
                            RDF_data.p(th_trimEvents) = new_LOCAL_p;
                            RDF_data.tm(th_trimEvents) = meltData(n).tm;
                            RDF_data.tl(th_trimEvents) = solData(n).tl;
                            RDF_data.cr(th_trimEvents) = solData(n).cr;
                        }
                        // Increment number of events
                        th_trimEvents++;
                    }
                }
            , trimEvents);  

            // Now deep copy within device to set new header information
            Kokkos::deep_copy(RDF_header.index_view, new_header.index_view);
            Kokkos::deep_copy(RDF_header.floatType_view, new_header.floatType_view);

            //Set new size (trimmed)
            RDF.numEvents = trimEvents;

            // Return data
            return RDF;
        }
    }
}
