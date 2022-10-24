

#include <TRestWimpSensitivity.h>

void REST_WIMP_Sensitivity(const std::string& rmlFile, const double wimpStart = 0.15,
                           const double wimpEnd = 20, const double wimpStep = 0.1) {
    TRestWimpSensitivity WS(rmlFile.c_str());
    WS.PrintMetadata();

    std::string outputFile = WS.BuildOutputFileName(".dat");

    std::ofstream ofs(outputFile, std::ofstream::out);
    if (!ofs.is_open()) {
        std::cout << "It was not possible to create file: " << outputFile << std::endl;
        return;
    }

    for (double wimpMass = wimpStart; wimpMass <= wimpEnd; wimpMass += wimpStep) {
        const double sens = WS.GetSensitivity(wimpMass);
        if (sens > 0) {
            std::cout << "WIMP mass " << wimpMass << " " << sens << std::endl;
            ofs << wimpMass << "\t" << sens << "\n";
        }
    }
}
