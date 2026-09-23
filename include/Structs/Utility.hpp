#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

namespace Stork {
    namespace Structs {

        // Class for Reduced Data Format header (on a single memory space)
        template <typename FloatType, typename memory_space>
        class RegularGrid_Header {
          private:
            /////////////////////////////////////////
            //// Set quick access to views types ////
            /////////////////////////////////////////
            using indexType_view = Kokkos::View<uint32_t[6], layout, memory_space>;
            using floatingType_view = Kokkos::View<FloatType[5], layout, memory_space>;

          public:
            /////////////////////////////////////////
            //// Views and easy access functions ////
            /////////////////////////////////////////
            indexType_view index_view;

            KOKKOS_INLINE_FUNCTION
            uint32_t& global_i0() const {
                return index_view(0);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& global_j0() const {
                return index_view(1);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& global_k0() const {
                return index_view(2);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& local_inum() const {
                return index_view(3);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& local_jnum() const {
                return index_view(4);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& local_knum() const {
                return index_view(5);
            }

            floatingType_view floatType_view;

            KOKKOS_INLINE_FUNCTION
            FloatType& global_x0() const {
                return floatType_view(0);
            }

            KOKKOS_INLINE_FUNCTION
            FloatType& global_y0() const {
                return floatType_view(1);
            }

            KOKKOS_INLINE_FUNCTION
            FloatType& global_z0() const {
                return floatType_view(2);
            }

            KOKKOS_INLINE_FUNCTION
            FloatType& gridResolution() const {
                return floatType_view(3);
            }

            // Default constructor with view allocation
            RegularGrid_Header() : index_view("RegularGrid_Header.index_view"), floatType_view("RegularGrid_Header.floatType_view") {}

            /////////////////////////
            /// UTILITY FUNCTIONS ///
            /////////////////////////

            // Get p from <i,j,k> (z->y->x indexing)
            KOKKOS_INLINE_FUNCTION
            void LOCAL_ijk_to_LOCAL_p(uint32_t& LOCAL_p, const uint32_t (&LOCAL_ijk)[3]) const {
                // Get references
                // const uint32_t& i_num = local_inum();
                const uint32_t& j_num = local_jnum();
                const uint32_t& k_num = local_knum();

                // Make and return p
                LOCAL_p = LOCAL_ijk[0] * (j_num * k_num) + LOCAL_ijk[1] * (k_num) + LOCAL_ijk[2];
            }

            // Get ijk from p (z->y->x indexing)
            KOKKOS_INLINE_FUNCTION
            void LOCAL_p_to_LOCAL_ijk(uint32_t (&LOCAL_ijk)[3], const uint32_t& LOCAL_p) const {
                // Get references
                // const uint32_t& i_num = local_inum();
                const uint32_t& j_num = local_jnum();
                const uint32_t& k_num = local_knum();

                // Compute indices
                LOCAL_ijk[0] = (LOCAL_p / (k_num * j_num));
                LOCAL_ijk[1] = ((LOCAL_p / k_num) % j_num);
                LOCAL_ijk[2] = (LOCAL_p % k_num);
            }

            // Get ijk(global) from ijk(local)
            KOKKOS_INLINE_FUNCTION
            void LOCAL_ijk_to_GLOBAL_ijk(uint32_t (&GLOBAL_ijk)[3], const uint32_t (&LOCAL_ijk)[3]) const {
                // Shift LOCAL values into GLOBAL grid
                GLOBAL_ijk[0] = LOCAL_ijk[0] + global_i0();
                GLOBAL_ijk[1] = LOCAL_ijk[1] + global_j0();
                GLOBAL_ijk[2] = LOCAL_ijk[2] + global_k0();
            }

            // Get ijk(global) from ijk(local)
            KOKKOS_INLINE_FUNCTION
            void GLOBAL_ijk_to_LOCAL_ijk(uint32_t (&LOCAL_ijk)[3], const uint32_t (&GLOBAL_ijk)[3]) const {
                // Shift GLOBAL values into LOCAL grid
                LOCAL_ijk[0] = GLOBAL_ijk[0] - global_i0();
                LOCAL_ijk[1] = GLOBAL_ijk[1] - global_j0();
                LOCAL_ijk[2] = GLOBAL_ijk[2] - global_k0();
            }

            // Get ijk(global) from p (z->y->x indexing)
            KOKKOS_INLINE_FUNCTION
            void LOCAL_p_to_GLOBAL_ijk(uint32_t (&GLOBAL_ijk)[3], const uint32_t& LOCAL_p) const {
                // Get LOCAL ijk
                uint32_t LOCAL_ijk[3];
                LOCAL_p_to_LOCAL_ijk(LOCAL_ijk, LOCAL_p);

                // Now get global indices
                LOCAL_ijk_to_GLOBAL_ijk(GLOBAL_ijk, LOCAL_ijk);
            }

            // Get p(local) from ijk(global)
            KOKKOS_INLINE_FUNCTION
            void GLOBAL_ijk_to_LOCAL_p(uint32_t& LOCAL_p, const uint32_t (&GLOBAL_ijk)[3]) const {
                // Shift GLOBAL values into LOCAL grid
                uint32_t LOCAL_ijk[3];
                GLOBAL_ijk_to_LOCAL_ijk(LOCAL_ijk, GLOBAL_ijk);

                // Set LOCAL p
                LOCAL_ijk_to_LOCAL_p(LOCAL_p, LOCAL_ijk);
            }

            // Get xyz(global) from ijk(global)
            template <typename ReturnFloatType>
            KOKKOS_INLINE_FUNCTION void GLOBAL_ijk_to_GLOBAL_xyz(ReturnFloatType (&GLOBAL_xyz)[3], const uint32_t (&GLOBAL_ijk)[3]) const {
                GLOBAL_xyz[0] = global_x0() + GLOBAL_ijk[0] * gridResolution();
                GLOBAL_xyz[1] = global_y0() + GLOBAL_ijk[1] * gridResolution();
                GLOBAL_xyz[2] = global_z0() + GLOBAL_ijk[2] * gridResolution();
            }

            // Get ijk(global) from p (z->y->x indexing)
            template <typename ReturnFloatType>
            KOKKOS_INLINE_FUNCTION void LOCAL_p_to_GLOBAL_xyz(ReturnFloatType (&GLOBAL_xyz)[3], const uint32_t& LOCAL_p) const {
                // Get global ijk
                uint32_t GLOBAL_ijk[3];
                LOCAL_p_to_GLOBAL_ijk(GLOBAL_ijk, LOCAL_p);

                // Compute position
                GLOBAL_ijk_to_GLOBAL_xyz(GLOBAL_xyz, GLOBAL_ijk);
            }

            // Get ijk from p (z->y->x indexing)
            template <typename ReturnFloatType>
            KOKKOS_INLINE_FUNCTION void LOCAL_ijk_to_GLOBAL_xyz(ReturnFloatType (&GLOBAL_xyz)[3], const uint32_t (&LOCAL_ijk)[3]) const {
                // Shift local ijk to global ijk
                uint32_t GLOBAL_ijk[3] = {
                    LOCAL_ijk[0] + global_i0(),
                    LOCAL_ijk[1] + global_j0(),
                    LOCAL_ijk[2] + global_k0()};

                // Compute position
                GLOBAL_ijk_to_GLOBAL_xyz(GLOBAL_xyz, GLOBAL_ijk);
            }
        };

        // Class for Reduced Data Format header (on a single memory space)
        template <typename memory_space>
        class IndexBounds {
          private:
            /////////////////////////////////////////
            //// Set quick access to views types ////
            /////////////////////////////////////////
            using indexType_view = Kokkos::View<uint32_t[6], layout, memory_space>;

          public:
            /////////////////////////////////////////
            //// Views and easy access functions ////
            /////////////////////////////////////////
            indexType_view index_view;

            KOKKOS_INLINE_FUNCTION
            uint32_t& imin() const {
                return index_view(0);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& imax() const {
                return index_view(1);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& jmin() const {
                return index_view(2);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& jmax() const {
                return index_view(3);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& kmin() const {
                return index_view(4);
            }

            KOKKOS_INLINE_FUNCTION
            uint32_t& kmax() const {
                return index_view(5);
            }
        };
    } // namespace Structs
} // namespace Stork
