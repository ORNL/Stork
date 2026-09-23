#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

// Utility includes
#include "IO/Utility.hpp"

// Specific includes
#include "Structs/SRDF.hpp"

namespace Stork::IO {

    template <typename FloatType>
    void Output_SRDF_binary(const Structs::SRDF_Dual<FloatType>& SRDF, const std::string& name) {
        // Make file name
        const std::string binFile = name + ".stork";

        // Open binary file with error checking and exceptions
        std::ofstream os(binFile, std::ios::binary);
        os.exceptions(std::ofstream::failbit | std::ofstream::badbit);

        try {
            // Get sizes
            const uint32_t size = SRDF.numSnaps;
            const FloatType T_crit = SRDF.T_critical;
            const auto& H_index = SRDF.host_header.index_view;
            const auto& H_float = SRDF.host_header.floatType_view;
            const auto& D_cellNum = SRDF.host_data.cellNum_view;
            const auto& D_t = SRDF.host_data.times_view;
            const auto& D_T = SRDF.host_data.thermals_view;

            // Compute total buffer size more precisely
            const size_t bufferSize =
                sizeof(FileType) +                    // File type
                sizeof(DataType) +                    // Type specifier
                sizeof(uint32_t) +                    // Size
                sizeof(FloatType) +                   // Critical Temperature
                H_index.size() * sizeof(uint32_t) +   // Header index information
                H_float.size() * sizeof(FloatType) +  // Header float information
                D_cellNum.size() * sizeof(uint32_t) + // Data cellNum
                D_t.size() * sizeof(FloatType) +      // Data timesView
                D_T.size() * sizeof(FloatType);       // Data thermalsView

            // Preallocate buffer with exact size
            std::vector<char> buffer;
            buffer.reserve(bufferSize);

            // Make Lambdas
            auto insertPrimitive = impl::makeInsertPrimitive(buffer);
            auto writeViewData = impl::makeWriteViewData(buffer);

            // Write file type
            const FileType fileType = FileType::SRDF;
            insertPrimitive(fileType);

            // Write type specifier
            const DataType dataType = impl::getExpectedDataType<FloatType>();
            insertPrimitive(dataType);

            // Write size
            insertPrimitive(size);

            // Write temperature
            insertPrimitive(T_crit);

            // Append all header sections
            writeViewData(H_index);
            writeViewData(H_float);

            // Make "subviews" to write (since the capacity may not equal the number used)
            using SRDF_DataType = Structs::SRDF_Data<FloatType, host_space>;
            using IndexViewType = typename SRDF_DataType::indexType_view;    // Access type alias
            using FloatViewType = typename SRDF_DataType::floatingType_view; // Access type alias
            const auto D_cellNum_subview = IndexViewType(D_cellNum.data(), size);
            const auto D_t_subview = FloatViewType(D_t.data(), 2 * size);
            const auto D_T_subview = FloatViewType(D_T.data(), 16 * size);
            writeViewData(D_cellNum_subview);
            writeViewData(D_t_subview);
            writeViewData(D_T_subview);

            // Single write operation for the entire buffer
            os.write(buffer.data(), buffer.size());
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to write SRDF binary file: " + binFile + ". Error: " + e.what());
        }
    }

    template <typename FloatType>
    void Input_SRDF_binary(Structs::SRDF_Dual<FloatType>& SRDF, const std::string& name) {

        // Make file name
        const std::string binFile = name + ".stork";

        // Open binary file
        std::ifstream is(binFile, std::ios::binary);
        if (!is) {
            throw std::runtime_error("Failed to open file: " + binFile);
        }

        // Make lambdas
        auto readPrimitive = impl::makeReadPrimitive(is);
        auto readViewData = impl::makeReadViewData(is);

        // Read and verify FileType
        FileType fileType;
        readPrimitive(fileType);
        if (fileType != FileType::SRDF) {
            throw std::runtime_error("Invalid file type: expected SRDF, got " + impl::fileTypeToString(fileType) + "\n");
        }

        // Read and verify DataType
        DataType dataType;
        readPrimitive(dataType);
        if (dataType != impl::getExpectedDataType<FloatType>()) {
            throw std::runtime_error("Invalid data type: expected " + impl::getExpectedDataTypeString<FloatType>() + ", got " + impl::dataTypeToString(dataType) + "\n");
        }

        // Read size
        readPrimitive(SRDF.numSnaps);
        const uint32_t& size = SRDF.numSnaps;

        // Read critical temperature
        readPrimitive(SRDF.T_critical);

        // Resize SRDF views
        Kokkos::realloc(Kokkos::WithoutInitializing, SRDF.host_data.cellNum_view, size);
        Kokkos::realloc(Kokkos::WithoutInitializing, SRDF.host_data.times_view, 2 * size);
        Kokkos::realloc(Kokkos::WithoutInitializing, SRDF.host_data.thermals_view, 16 * size);

        // Read header info
        readViewData(SRDF.host_header.index_view);
        readViewData(SRDF.host_header.floatType_view);

        // Read cellNum, times, and thermals data
        readViewData(SRDF.host_data.cellNum_view);
        readViewData(SRDF.host_data.times_view);
        readViewData(SRDF.host_data.thermals_view);

        // Check for read errors
        if (!is) {
            throw std::runtime_error("Error reading from file: " + binFile);
        }
    }
} // namespace Stork::IO
