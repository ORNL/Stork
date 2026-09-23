// Includes
#include <Stork_Core.hpp>

// Templated Helper Function
template<typename FloatType_In, typename FloatType_Out>
void Run(const std::string file, const int fineFactor){
    // Make SRDF Object
    Stork::Structs::SRDF_Dual<FloatType_In> SRDF;
    
    // Read file into SRDF Object 
    Stork::IO::Input_SRDF_binary<FloatType_In>(SRDF, file);
    
    // Fine Interpolation
    Stork::Structs::RDF_Dual<FloatType_Out> RDF = Stork::Run::Interpolate_SRDF_to_RDF<FloatType_In, Stork::host_space, FloatType_Out, Stork::device_space>(SRDF, fineFactor);
    
    // Copy Data from Device to Host and Free Memory on Device
    RDF.template Make_Data_Mirrors<Stork::device_space, Stork::host_space>();
    RDF.template Copy_All<Stork::device_space, Stork::host_space>();
    RDF.template Free_Data_Memory<Stork::device_space>();

    // Output Data to File
    Stork::IO::Output_RDF_csv<FloatType_Out>(RDF, file+".interp." + std::to_string(fineFactor));

    // How Big?
    std::cout << "NumEvents: " << RDF.numEvents << std::endl;
}

// Main Function
int main(int argc, char* argv[]) {

    // Initialize Kokkos
    Kokkos::initialize(argc, argv);
	{
        
        // Tell user whats up if they get it wrong
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " <filename> <fineFactor> <precisionIn> <precisionOut>" << std::endl;
            return 1;
        }

        // Get inputs (filename and interpFactor)
        const std::string file = argv[1];
        const int fineFactor = std::stoi(argv[2]);

        // Precisions
        const std::string precisionIn = argv[3];
        const std::string precisionOut = argv[4];

        const std::string floatStr = "float";
        const std::string doubleStr = "double";

        if (precisionIn == floatStr && precisionOut == floatStr) {
            Run<float,float>(file, fineFactor);
        } 
        else if (precisionIn == floatStr && precisionOut == doubleStr) {
            Run<float,double>(file, fineFactor);
        } 
        else if (precisionIn == doubleStr && precisionOut == floatStr) {
            Run<double,float>(file, fineFactor);
        } 
        else if (precisionIn == doubleStr && precisionOut == doubleStr) {
            Run<double,double>(file, fineFactor);
        } 
        else {
            std::cerr << "Invalid precision(s). Use 'float' or 'double'." << std::endl;
            return 1;
        }
        
    }
    // Finalize Kokkos
    Kokkos::finalize();

    return 0;
}
