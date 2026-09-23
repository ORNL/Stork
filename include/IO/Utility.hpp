#pragma once

// Common includes
#include "Definitions.hpp"
#include "Common.hpp"

namespace Stork::IO {

    enum class FileType : uint8_t {
        RDF = 0,
        SRDF = 1
    };

    enum class DataType : uint8_t {
        FLOAT32 = 0,
        FLOAT64 = 1
    };

    // For "stringifying" the enums
    namespace impl {
        using namespace Common;

        // For fileTypes
        constexpr array<string_view, 2> FileTypeNames = {"RDF", "SRDF"};
        inline constexpr string_view fileTypeToStringView(FileType fileType) {
            uint8_t index = static_cast<uint8_t>(fileType);
            return index < FileTypeNames.size() ? FileTypeNames[index] : "Unknown";
        }
        inline string fileTypeToString(FileType fileType) {
            return string(fileTypeToStringView(fileType));
        }

        // For dataTypes
        constexpr array<string_view, 2> DataTypeNames = {"FLOAT32", "FLOAT64"};
        inline constexpr string_view dataTypeToStringView(DataType dataType) {
            uint8_t index = static_cast<uint8_t>(dataType);
            return index < DataTypeNames.size() ? DataTypeNames[index] : "Unknown";
        }
        inline string dataTypeToString(DataType dataType) {
            return string(dataTypeToStringView(dataType));
        }

        // DataType Validation Helper
        template <typename FloatType>
        constexpr DataType getExpectedDataType() {
            if constexpr (std::is_same_v<FloatType, float>) {
                return DataType::FLOAT32;
            } else if constexpr (std::is_same_v<FloatType, double>) {
                return DataType::FLOAT64;
            } else {
                throw std::runtime_error("Unsupported FloatType provided!");
            }
        }

        // DataType Validation Helper
        template <typename FloatType>
        string getExpectedDataTypeString() {
            if constexpr (std::is_same_v<FloatType, float>) {
                return dataTypeToString(DataType::FLOAT32);
            } else if constexpr (std::is_same_v<FloatType, double>) {
                return dataTypeToString(DataType::FLOAT64);
            } else {
                throw std::runtime_error("Unsupported FloatType provided!");
            }
        }
    } // namespace impl

    // Get Filename Extension
    namespace impl {
        using namespace Common;

        inline string getFileExtension(const string& filename) {
            // Find the last period in the filename
            size_t dotPosition = filename.find_last_of('.');
            if (dotPosition == string::npos) {
                return ""; // No extension found
            }
            return filename.substr(dotPosition);
        }
    } // namespace impl

    // For making binary write/read lambdas
    namespace impl {
        // Write-side lambdas (templated to be usable with different buffer and view types)
        template <typename BufferType>
        auto makeInsertPrimitive(BufferType& buffer) {
            return [&buffer](const auto& value) {
                buffer.insert(buffer.end(),
                              reinterpret_cast<const char*>(&value),
                              reinterpret_cast<const char*>(&value + 1));
            };
        }

        template <typename BufferType>
        auto makeWriteViewData(BufferType& buffer) {
            return [&buffer](const auto& data) {
                buffer.insert(buffer.end(),
                              reinterpret_cast<const char*>(data.data()),
                              reinterpret_cast<const char*>(data.data() + data.size()));
            };
        }

        // Read-side lambdas (templated to be usable with different stream and view types)
        template <typename StreamType>
        auto makeReadPrimitive(StreamType& is) {
            return [&is](auto& value) {
                is.read(reinterpret_cast<char*>(&value), sizeof(decltype(value)));
            };
        }

        template <typename StreamType>
        auto makeReadViewData(StreamType& is) {
            return [&is](auto& view) {
                using ValueType = typename std::decay_t<decltype(view)>::value_type;
                is.read(
                    reinterpret_cast<char*>(view.data()),
                    view.size() * sizeof(ValueType));
            };
        }
    } // namespace impl
} // namespace Stork::IO
