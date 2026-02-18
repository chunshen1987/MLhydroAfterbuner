/////////////////////////////////////////////////////////////////////////
//                      MLhydroAfterburn
//
//              author: Chun Shen <chunshen@wayne.edu>
//              Copyright: Chun Shen 2026
//
//  This program load hydrodynamic evolution files and find hyper-surface
//  at a given energy density
//
//
/////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <memory>

#include "Hydroinfo_MUSIC.h"
#include "Stopwatch.h"
#include "ParameterReader.h"
#include "SurfaceFinder.h"


int main(int argc, char *argv[]) {
    ParameterReader *paraRdr = new ParameterReader();
    paraRdr->readFromFile("parameters.input");
    paraRdr->readFromArguments(argc, argv);
    paraRdr->echo();

    int load_viscous = paraRdr->getVal("load_viscous_info");
    int hydro_type   = paraRdr->getVal("hydro_type");

    std::shared_ptr<Hydroinfo_MUSIC> hydroinfo_ptr = (
            std::make_shared<Hydroinfo_MUSIC>());

    Stopwatch sw;
    sw.tic();
    // hydro data file pointer
    // 3: (2+1)D hydro MUSIC (no grid Chun's format)
    // 4: (3+1)D hydro MUSIC (no grid Chun's format)
    hydroinfo_ptr->readHydroData(hydro_type);

    //FluidcellStatistic fluidcellanalysis(hydroinfo_ptr_in, paraRdr);
    //fluidcellanalysis.outputTempasTauvsX();
    //fluidcellanalysis.outputKnudersonNumberasTauvsX();
    //fluidcellanalysis.outputinverseReynoldsNumberasTauvsX();
    ////double T_cut = paraRdr->getVal("T_cut");
    ////fluidcellanalysis.analysis_hydro_volume_for_photon(T_cut);
    //fluidcellanalysis.output_temperature_vs_avg_utau();
    //fluidcellanalysis.output_flowvelocity_vs_tau();

    // construct freeze-out hyper-surface
    // SurfaceFinder* surface_ptr = new SurfaceFinder(hydroinfo_ptr, paraRdr);
    // surface_ptr->Find_full_hypersurface();

    sw.toc();
    std::cout << "totally takes : " << sw.takeTime() << " seconds."
              << std::endl;

    return(0);
}


