#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

namespace Stork { 
    namespace Structs {

        // Class for Reduced Data Format (on a single memory space)
        template<typename FloatType, typename memory_space>
        class RDF_Data {
            private:
                /////////////////////////////////////////
                //// Set quick access to views types ////
                /////////////////////////////////////////
                using indexType_view = Kokkos::View<uint32_t*, layout, memory_space>;
                using floatingType_view = Kokkos::View<FloatType*, layout, memory_space>;
            public:        
                /////////////////////////////////////////
                //// Views and easy access functions ////
                /////////////////////////////////////////  
                indexType_view cellNum_view;
                floatingType_view solInfo_view;
                
                KOKKOS_INLINE_FUNCTION
                uint32_t& p(const uint32_t n) const {
                    return cellNum_view(n);
                }

                KOKKOS_INLINE_FUNCTION
                FloatType& tm(const uint32_t n) const {
                    return solInfo_view(3* n);
                }

                KOKKOS_INLINE_FUNCTION
                FloatType& tl(const uint32_t n) const {
                    return solInfo_view(3* n + 1);
                }

                KOKKOS_INLINE_FUNCTION
                FloatType& cr(const uint32_t n) const {
                    return solInfo_view(3* n + 2);
                }     
        };

        // Class for Reduced Data Format (on all memory spaces)
        template <typename FloatType>
        class RDF_Dual{
            public:
                RDF_Dual() {
                    Mirror_Header_When_Same_Space();
                    Mirror_Data_When_Same_Space();
                } 
                
                ///////////////
                /// OBJECTS ///
                ///////////////

                // Number Of Events
                uint32_t numEvents = 0;
                
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
                RDF_Data<FloatType, host_space> host_data;
                RDF_Data<FloatType, device_space> device_data;

                    // Get specific RDF based on memory space
                    template<typename memory_space>
                    KOKKOS_INLINE_FUNCTION
                    RDF_Data<FloatType, memory_space>& get_data() {
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
                    const RDF_Data<FloatType, memory_space>& get_data() const {
                        if constexpr (std::is_same_v<memory_space, host_space>) {
                            return host_data;
                        } 
                        else {
                            return device_data;
                        }
                    }

                ////////////////////////
                /// MEMORY FUNCTIONS ///
                ////////////////////////

                void Mirror_Header_When_Same_Space(){
                    if constexpr (std::is_same<host_space, device_space>::value){
                        device_header.index_view = host_header.index_view;
                        device_header.floatType_view = host_header.floatType_view;
                    }
                }

                void Mirror_Data_When_Same_Space(){
                    if constexpr (std::is_same<host_space, device_space>::value){
                        device_data.cellNum_view = host_data.cellNum_view;
                        device_data.solInfo_view = host_data.solInfo_view;
                    }
                }
                
                // Initialze views
                template<typename memory_space>
                void Make_Data_Views(const uint32_t size){
                    // Get references
                    RDF_Data<FloatType, memory_space>& data = this->template get_data<memory_space>();
                    // Quick access to views
                    using indexType_view = Kokkos::View<uint32_t*, layout, memory_space>;
                    using floatingType_view = Kokkos::View<FloatType*, layout, memory_space>;
                    // Make views
                    data.cellNum_view = indexType_view(Kokkos::ViewAllocateWithoutInitializing("RDF_cellNum"), size);
                    data.solInfo_view = floatingType_view(Kokkos::ViewAllocateWithoutInitializing("RDF_solInfo"), 3*size);
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
                        RDF_Data<FloatType, source_space>& srcData = this->template get_data<source_space>();
                        RDF_Data<FloatType, destination_space>& dstData = this->template get_data<destination_space>();
                        // Make Mirrors
                        dstData.cellNum_view = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, destination_space(), srcData.cellNum_view);
                        dstData.solInfo_view = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, destination_space(), srcData.solInfo_view);
                    }
                }

                // Copy views
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
                        RDF_Data<FloatType, source_space>& srcData = this->template get_data<source_space>();
                        RDF_Data<FloatType, destination_space>& dstData = this->template get_data<destination_space>();
                        // Copy Data between spaces
                        Kokkos::deep_copy(dstData.cellNum_view, srcData.cellNum_view);
                        Kokkos::deep_copy(dstData.solInfo_view, srcData.solInfo_view);
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
                        RDF_Data<FloatType, memory_space>& data = this->template get_data<memory_space>();
                        // Free host memory by resizing or reinitializing host views to zero
                        Kokkos::resize(data.cellNum_view, 0);
                        Kokkos::resize(data.solInfo_view, 0);
                    }
                }
        };
    }
}
