#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>  // for system()

void aPlotRunner() {
    std::ifstream fin("../mycustom.txt");
    if (!fin.is_open()) {
        std::cerr << "txt missing\n";
        return;
    }

    std::string line;
    int idx = 0;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        double SEF, sieie, toe, isoC, ecal, hcal;
        if (!(iss >> SEF >> sieie >> toe >> isoC >> ecal >> hcal)) continue;

        std::ofstream config("input_params.h");
        config << "double medium_sieie = " << sieie << ";\n";
        config << "double medium_toe = " << toe << ";\n";
        config << "double medium_isoC = " << isoC << ";\n";
        config << "double medium_ecal = " << ecal << ";\n";
        config << "double medium_hcal = " << hcal << ";\n";
        config << "const char* output_filename = \"" << int(SEF * 100 + 0.5) << "_MediumEffBck_TruePV.png\";\n";
        config << "const char* testoutput = \"" << "SEF" <<int(SEF * 100 + 0.5) << "_Barrel_Plots_updat_SR_TruePVID.root\";\n";
        config.close();

        std::cout << "Running for SEF = " << SEF << std::endl;
        gSystem->Exec("root -l -b -q PlotBuilder.C");
    }

    fin.close();
}

