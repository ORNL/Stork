#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

// Utility includes
#include "IO/Utility.hpp"

// Specific includes
#include "Structs/RDF.hpp"

namespace Stork::IO {

    // Outputting Data to a CSV file
    template <typename FloatType>
    void Output_RDF_csv(Structs::RDF_Dual<FloatType>& RDF, const Common::string name) {
        // Set reference
        Structs::RegularGrid_Header<FloatType, host_space>& header = RDF.host_header;
        Structs::RDF_Data<FloatType, host_space>& data = RDF.host_data;

        // Get total number of events
        const size_t numEvents = RDF.numEvents;

        // Output CSV
        std::ofstream datafile;
        datafile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
        Common::string out_file = name + ".csv";
        try {
            datafile.open(out_file.c_str());
            datafile << "x,y,z,tm,tl,cr\n";
            for (uint32_t n = 0; n < numEvents; n++) {
                // Get global xyz from local p
                const uint32_t& LOCAL_p = data.p(n);
                FloatType GLOBAL_xyz[3];
                header.LOCAL_p_to_GLOBAL_xyz(GLOBAL_xyz, LOCAL_p);
                // Output to csv file
                datafile << GLOBAL_xyz[0] << "," << GLOBAL_xyz[1] << "," << GLOBAL_xyz[2] << ",";
                datafile << data.tm(n) << "," << data.tl(n) << "," << data.cr(n) << "\n";
            }
        } catch (const std::ofstream::failure& e) {
            std::cout << "Exception writing data file, check that Data directory exists\n";
        }
        datafile.close();
    }

    template <typename FloatType>
    void Output_RDF_binary(const Structs::RDF_Dual<FloatType>& RDF, const std::string& name) {
        // Make file name
        const std::string binFile = name + ".stork";

        // Open binary file
        std::ofstream os(binFile, std::ios::binary);
        if (!os) {
            throw std::runtime_error("Failed to open file: " + binFile);
        }

        // Get sizes
        const uint32_t size = RDF.host_data.cellNum_view.size();
        const auto& indexData = RDF.host_header.index_view;
        const auto& floatData = RDF.host_header.floatType_view;
        const auto& cellNumData = RDF.host_data.cellNum_view;
        const auto& solInfoData = RDF.host_data.solInfo_view;

        // Reserve buffer
        size_t bufferSize = sizeof(FileType) +                      // File type
                            sizeof(DataType) +                      // Type specifier
                            sizeof(uint32_t) +                      // Size
                            indexData.size() * sizeof(uint32_t) +   // Header index information
                            floatData.size() * sizeof(FloatType) +  // Header float information
                            cellNumData.size() * sizeof(uint32_t) + // Data cellNum
                            solInfoData.size() * sizeof(FloatType); // Data solInfo
        std::vector<char> buffer;
        buffer.reserve(bufferSize);

        // Make Lambdas
        auto insertPrimitive = impl::makeInsertPrimitive(buffer);
        auto writeViewData = impl::makeWriteViewData(buffer);

        // Write file type
        const FileType fileType = FileType::RDF;
        insertPrimitive(fileType);

        // Write type specifier
        const DataType dataType = impl::getExpectedDataType<FloatType>();
        insertPrimitive(dataType);

        // Write size
        insertPrimitive(size);

        // Append all data sections
        writeViewData(indexData);
        writeViewData(floatData);
        writeViewData(cellNumData);
        writeViewData(solInfoData);

        // Write buffer to file
        os.write(buffer.data(), buffer.size());
        if (!os) {
            throw std::runtime_error("Failed to write data to file: " + binFile);
        }
    }

    template <typename FloatType>
    void Input_RDF_binary(Structs::RDF_Dual<FloatType>& RDF, const std::string& name) {

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
        if (fileType != FileType::RDF) {
            throw std::runtime_error("Invalid file type: expected RDF, got " + impl::fileTypeToString(fileType) + "\n");
        }

        // Read and verify DataType
        DataType dataType;
        readPrimitive(dataType);
        if (dataType != impl::getExpectedDataType<FloatType>()) {
            throw std::runtime_error("Invalid data type: expected " + impl::getExpectedDataTypeString<FloatType>() + ", got " + impl::dataTypeToString(dataType) + "\n");
        }

        // Read size
        readPrimitive(RDF.numEvents);
        const uint32_t& size = RDF.numEvents;

        // Resize RDF views
        Kokkos::realloc(Kokkos::WithoutInitializing, RDF.host_data.cellNum_view, size);
        Kokkos::realloc(Kokkos::WithoutInitializing, RDF.host_data.solInfo_view, 3 * size);

        // Read header info
        readViewData(RDF.host_header.index_view);
        readViewData(RDF.host_header.floatType_view);

        // Read cellNum, times, and thermals data
        readViewData(RDF.host_data.cellNum_view);
        readViewData(RDF.host_data.solInfo_view);

        // Check for read errors
        if (!is) {
            throw std::runtime_error("Error reading from file: " + binFile);
        }
    }
} // namespace Stork::IO
