#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

namespace Stork { 
    namespace Structs {

        // Class for Space/Super Reduced Data Format (on a single memory space)
        template<typename FloatType, typename memory_space>
        class SRDF_Data {
            private:
                /////////////////////////////////////////
                //// Set quick access to views types ////
                /////////////////////////////////////////
                // using indexType_view = Kokkos::View<uint32_t*, layout, memory_space>;
                // using floatingType_view = Kokkos::View<FloatType*, layout, memory_space>;
            
            public:
                /////////////////////////////////////////
                //// Views and easy access functions ////
                /////////////////////////////////////////    

                // TODO::DEBUG
                using indexType_view = Kokkos::View<uint32_t*, layout, memory_space>;
                using floatingType_view = Kokkos::View<FloatType*, layout, memory_space>;

                indexType_view cellNum_view;
                
                    KOKKOS_INLINE_FUNCTION
                    uint32_t& p(const uint32_t n) const {
                        return cellNum_view(n);
                    }

                floatingType_view times_view;

                    KOKKOS_INLINE_FUNCTION
                    FloatType& t_prev(const uint32_t n) const {
                        return times_view(2* n);
                    }

                    KOKKOS_INLINE_FUNCTION
                    FloatType& t_cur(const uint32_t n) const {
                        return times_view(2* n + 1);
                    }

                floatingType_view thermals_view;

                    KOKKOS_INLINE_FUNCTION
                    FloatType& T_prev(const uint32_t n) const {
                        return thermals_view(16* n);
                    }

                    KOKKOS_INLINE_FUNCTION
                    FloatType& T_cur(const uint32_t n) const {
                        return thermals_view(16* n + 8);
                    }
        };

        // Class for Space/Super Reduced Data Format (on all memory spaces)
        template <typename FloatType>
        class SRDF_Dual{
            public:
                SRDF_Dual() {
                    Mirror_Header_When_Same_Space();
                    Mirror_Data_When_Same_Space();
                } 

                ///////////////
                /// OBJECTS ///
                ///////////////
                
                // Number Of Temperature Snapshots
                uint32_t numSnaps = 0;

                // Critical Temperature
                FloatType T_critical;

                // Header views on host and device
                RegularGrid_Header<FloatType, host_space> host_header;
                RegularGrid_Header<FloatType, device_space> device_header;

                    // Get specific RDF based on memory space
                    template<typename memory_space>
                    KOKKOS_INLINE_FUNCTION
                    RegularGrid_Header<FloatType, memory_space>& get_header() {
                        if constexpr (std::is_same_v<memory_space, host_space>) {
                            return host_header;
                        } 
                        else {
                            return device_header;
                        }
                    }

                    // Const-qualified version
                    template<typename memory_space>
                    KOKKOS_INLINE_FUNCTION
                    const RegularGrid_Header<FloatType, memory_space>& get_header() const {
                        if constexpr (std::is_same_v<memory_space, host_space>) {
                            return host_header;
                        } 
                        else {
                            return device_header;
                        }
                    }

                // Data views on host and device
                SRDF_Data<FloatType, host_space> host_data;
                SRDF_Data<FloatType, device_space> device_data;

                    // Get specific RDF based on memory space
                    template<typename memory_space>
                    KOKKOS_INLINE_FUNCTION
                    SRDF_Data<FloatType, memory_space>& get_data() {
                        if constexpr (std::is_same_v<memory_space, host_space>) {
                            return host_data;
                        } 
                        else {
                            return device_data;
                        }
                    }

                    // Const-qualified version
                    template<typename memory_space>
                    KOKKOS_INLINE_FUNCTION
                    const SRDF_Data<FloatType, memory_space>& get_data() const {
                        if constexpr (std::is_same_v<memory_space, host_space>) {
                            return host_data;
                        } 
                        else {
                            return device_data;
                        }
                    }

                /////////////////
                /// FUNCTIONS ///
                /////////////////

                void Mirror_Header_When_Same_Space(){
                    if constexpr (std::is_same<host_space, device_space>::value){
                        device_header.index_view = host_header.index_view;
                        device_header.floatType_view = host_header.floatType_view;
                    }
                }

                void Mirror_Data_When_Same_Space(){
                    if constexpr (std::is_same<host_space, device_space>::value){
                        device_data.cellNum_view = host_data.cellNum_view;
                        device_data.times_view = host_data.times_view;
                        device_data.thermals_view = host_data.thermals_view;
                    }
                }

                // Initialze views
                template<typename memory_space>
                void Make_Data_Views(const uint32_t size){
                    // Get references
                    SRDF_Data<FloatType, memory_space>& data = this->template get_data<memory_space>();
                    // Quick access to views
                    using indexType_view = Kokkos::View<uint32_t*, layout, memory_space>;
                    using floatingType_view = Kokkos::View<FloatType*, layout, memory_space>;
                    // Make views
                    data.cellNum_view = indexType_view(Kokkos::ViewAllocateWithoutInitializing("srdf_indices"), size);
                    data.times_view = floatingType_view(Kokkos::ViewAllocateWithoutInitializing("srdf_indices"), 2 * size);
                    data.thermals_view = floatingType_view(Kokkos::ViewAllocateWithoutInitializing("srdf_indices"), 16 * size);
                    Mirror_Data_When_Same_Space();
                }

                // Make Mirrors
                template<typename source_space, typename destination_space>
                void Make_Data_Mirrors(){
                    if constexpr (std::is_same<source_space, destination_space>::value){
                        Mirror_Data_When_Same_Space();
                    }
                    else{
                        // Get references
                        SRDF_Data<FloatType, source_space>& srcData = this->template get_data<source_space>();
                        SRDF_Data<FloatType, destination_space>& dstData = this->template get_data<destination_space>();
                        // Make Mirrors
                        dstData.cellNum_view = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, destination_space(), srcData.cellNum_view);
                        dstData.times_view = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, destination_space(), srcData.times_view);
                        dstData.thermals_view = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, destination_space(), srcData.thermals_view);
                    }
                }

                // Copy header information
                template<typename source_space, typename destination_space>
                void Copy_Header(){
                    if constexpr (std::is_same<source_space, destination_space>::value){
                        Mirror_Header_When_Same_Space();
                    }
                    else{
                        // Get references
                        RegularGrid_Header<FloatType, source_space>& srcHeader = this->template get_header<source_space>();
                        RegularGrid_Header<FloatType, destination_space>& dstHeader = this->template get_header<destination_space>();
                        // Copy Data between spaces
                        Kokkos::deep_copy(dstHeader.index_view, srcHeader.index_view);
                        Kokkos::deep_copy(dstHeader.floatType_view, srcHeader.floatType_view);
                    }   
                }

                // Copy views
                template<typename source_space, typename destination_space>
                void Copy_Data(){
                    if constexpr (std::is_same<source_space, destination_space>::value){
                        Mirror_Data_When_Same_Space();
                    }
                    else{
                        // Get references
                        SRDF_Data<FloatType, source_space>& srcData = this->template get_data<source_space>();
                        SRDF_Data<FloatType, destination_space>& dstData = this->template get_data<destination_space>();
                        // Copy Data between spaces
                        Kokkos::deep_copy(dstData.cellNum_view, srcData.cellNum_view);
                        Kokkos::deep_copy(dstData.times_view, srcData.times_view);
                        Kokkos::deep_copy(dstData.thermals_view, srcData.thermals_view);
                    }   
                }

                // Copy views
                template<typename source_space, typename destination_space>
                void Copy_All(){
                    Copy_Header<source_space, destination_space>();
                    Copy_Data<source_space, destination_space>();
                }

                // Free memory
                template<typename memory_space>
                void Free_Data_Memory(){
                    // Do nothing if the same space (as it will free when out of scope anyways)
                    if constexpr (!std::is_same<host_space, device_space>::value){
                        // Get references
                        SRDF_Data<FloatType, memory_space>& data = this->template get_data<memory_space>();
                        // Free host memory by resizing or reinitializing host views to zero
                        Kokkos::resize(data.cellNum_view, 0);
                        Kokkos::resize(data.times_view, 0);
                        Kokkos::resize(data.thermals_view, 0);
                    }
                }
        };   
    }
}
