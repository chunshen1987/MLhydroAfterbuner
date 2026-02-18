// Hydroinfo_MUSIC.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009-2010 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in hydro data from files and functions
// that return interpolated data at a given space-time point

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>

#include "data_struct.h"
#include "Hydroinfo_MUSIC.h"

using namespace std;


Hydroinfo_MUSIC::Hydroinfo_MUSIC() {
    boost_invariant_ = false;
}


Hydroinfo_MUSIC::~Hydroinfo_MUSIC() {
    lattice_3D_ideal.clear();
}


void Hydroinfo_MUSIC::readHydroData(int whichHydro) {
    hydroType_ = whichHydro;

    // all hydro data is stored in tau steps (not t)
    // evolution is converted to tau when accessing the hydro data
    lattice_3D_ideal.clear();

    // read in setups of the hydro simulation
    cout << "Using new MUSIC hydro format (no grid) reading data ..."
         << endl;
    if (whichHydro == 4) {
        boost_invariant_ = false;
    } else {
        boost_invariant_ = true;
    }

    // read in temperature and flow velocity
    // The name of the evolution file: evolution_name
    string evolution_name = "results/evolution_all_xyeta.dat";
    cout << "Evolution file name = " << evolution_name << endl;
    std::FILE *fin;
    fin = std::fopen(evolution_name.c_str(), "rb");
    if (fin == NULL) {
        cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
             << "Unable to open file: " << evolution_name << endl;
        exit(1);
    }

    float header[16];
    int status = std::fread(&header, sizeof(float), 16, fin);
    if (status == 0) {
        cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
             << "Can not read the evolution file header" << endl;
        exit(1);
    }

    hydroTau0 = header[0];
    hydroDtau = header[1];
    ixmax = static_cast<int>(header[2]);
    hydroDx = header[3];
    hydroXmax = std::abs(header[4]);
    iymax = static_cast<int>(header[5]);
    hydroDy = header[6];
    hydroYmax = std::abs(header[7]);
    ietamax = static_cast<int>(header[8]);
    hydroDeta = header[9];
    hydro_eta_max = std::abs(header[10]);
    turn_on_rhob = static_cast<int>(header[11]);
    turn_on_shear = static_cast<int>(header[12]);
    turn_on_bulk = static_cast<int>(header[13]);
    turn_on_diff = static_cast<int>(header[14]);
    const int nVar_per_cell = static_cast<int>(header[15]);

    float* cell_info = new float [nVar_per_cell];

    int itau_max = 0;
    fluidCell_3D_ideal zeroCell;
    zeroCell.itau = 0;
    zeroCell.ix = 0;
    zeroCell.iy = 0;
    zeroCell.ieta = 0;
    zeroCell.temperature = 0.;
    zeroCell.ed = 0.;
    zeroCell.pressure = 0.;
    zeroCell.ux = 0.;
    zeroCell.uy = 0.;
    zeroCell.uz = 0.;
    lattice_3D_ideal.push_back(zeroCell);
    int ik = 0;
    while (true) {
        status = 0;
        status = std::fread(cell_info, sizeof(float), nVar_per_cell, fin);
        if (status == 0) break;
        if (status != nVar_per_cell) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "the evolution file format is not correct" << endl;
            exit(1);
        }

        if (itau_max < static_cast<int>(cell_info[0]))
            itau_max = static_cast<int>(cell_info[0]);
        fluidCell_3D_ideal newCell;
        newCell.itau = static_cast<int>(cell_info[0]);
        newCell.ix   = static_cast<int>(cell_info[1]);
        newCell.iy   = static_cast<int>(cell_info[2]);
        newCell.ieta = static_cast<int>(cell_info[3]);
        newCell.temperature = cell_info[6];
        newCell.ed = cell_info[4];
        newCell.pressure = cell_info[5];
        newCell.ux = cell_info[8];
        newCell.uy = cell_info[9];
        newCell.uz = cell_info[10];
        lattice_3D_ideal.push_back(newCell);
        ik++;
        if (ik%50000 == 0)
            cout << "o" << flush;
    }
    cout << endl;
    std::fclose(fin);
    itaumax = itau_max;
    // create the index map
    long long ncells = (itaumax + 1)*ixmax*iymax*ietamax;
    idx_map_.resize(ncells, 0);
    for (unsigned int i = 0; i < lattice_3D_ideal.size(); i++) {
        const auto cell_i = lattice_3D_ideal[i];
        int cell_idx = (
            (  (cell_i.itau*ietamax + cell_i.ieta)*iymax
             + cell_i.iy)*ixmax + cell_i.ix);
        idx_map_[cell_idx] = i;
    }
    hydroTauMax = hydroTau0 + hydroDtau*itaumax;

    // One final step for easy automation of MARTINI:
    // hydroTauMax is reset for the case where writing to evolution.dat
    // ended early (due to all cells freezing out):
    cout << "hydro_tau0 = " << hydroTau0 << " fm"<< endl;
    cout << "hydro_tau_max = " << hydroTauMax << " fm" << endl;
    cout << "hydry_dtau = " << hydroDtau << " fm" << endl;
    cout << "hydro_Xmax = " << hydroXmax << " fm" << endl;
    cout << "hydro_dx = " << hydroDx << " fm" << endl;
    cout << "hydro_Ymax = " << hydroYmax << " fm" << endl;
    cout << "hydro_dy = " << hydroDy << " fm" << endl;
    cout << "hydro_eta_max = " << hydro_eta_max << endl;
    cout << "hydro_deta = " << hydroDeta << endl;

    delete[] cell_info;
}


void Hydroinfo_MUSIC::getHydroValues(float x, float y, float eta,
                                     float tau, fluidCell* info) {
// For interpolation of evolution files in tau-eta coordinates. Only the
// reading of MUSIC's evolution_xyeta.dat file is implemented here.

    int itau = static_cast<int>((tau-hydroTau0)/hydroDtau + 0.0001);
    int ix   = static_cast<int>((hydroXmax+x)/hydroDx + 0.0001);
    int iy   = static_cast<int>((hydroYmax+y)/hydroDy + 0.0001);
    int ieta = static_cast<int>((hydro_eta_max+eta)/hydroDeta + 0.0001);

    float taufrac = (tau - hydroTau0)/hydroDtau - static_cast<float>(itau);
    float xfrac   = (x - (static_cast<float>(ix)*hydroDx - hydroXmax))/hydroDx;
    float yfrac   = (y - (static_cast<float>(iy)*hydroDy - hydroYmax))/hydroDy;
    float etafrac = (eta/hydroDeta - static_cast<float>(ieta)
                     + 0.5*static_cast<float>(ietamax));

    if (boost_invariant_) {
        ieta = 0;
        etafrac = 0.;
    }

    if (ix < 0 || ix >= ixmax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - x out of range x=" << x
             << ", ix=" << ix << ", ixmax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info->ed = 0.0;
        info->temperature = 0.0;
        info->utau = 1.0;
        info->ux = 0.0;
        info->uy = 0.0;
        info->ueta = 0.0;
        return;
    }

    if (iy < 0 || iy >= iymax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - y out of range, y=" << y << ", iy="  << iy
             << ", iymax=" << iymax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info->ed = 0.0;
        info->temperature = 0.0;
        info->utau = 1.0;
        info->ux = 0.0;
        info->uy = 0.0;
        info->ueta = 0.0;
        return;
    }

    if (itau < 0 || itau > itaumax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "tau out of range, itau=" << itau << ", itaumax=" << itaumax
             << endl;
        cout << "[Hydroinfo_MUSIC::getHydroValues]: tau= " << tau
             << ", hydroTauMax = " << hydroTauMax
             << ", hydroDtau = " << hydroDtau << endl;

        info->ed = 0.0;
        info->temperature = 0.0;
        info->utau = 1.0;
        info->ux = 0.0;
        info->uy = 0.0;
        info->ueta = 0.0;
        return;
    }

    if (ieta < 0 || ieta >= ietamax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "eta out of range, ieta=" << ieta << ", ietamax=" << ietamax
             << endl;
        info->ed = 0.0;
        info->temperature = 0.0;
        info->utau = 1.0;
        info->ux = 0.0;
        info->uy = 0.0;
        info->ueta = 0.0;
        return;
    }

    // The array of positions on the 4-dimensional rectangle:
    int position[2][2][2][2];
    for (int ipx = 0; ipx < 2; ipx++) {
        int px;
        if (ipx == 0 || ix == ixmax-1) {
            px = ix;
        } else {
            px = ix + 1;
        }
        for (int ipy = 0; ipy < 2; ipy++) {
            int py;
            if (ipy == 0 || iy == iymax-1) {
                py = iy;
            } else {
                py = iy + 1;
            }
            for (int ipeta = 0; ipeta < 2; ipeta++) {
                int peta;
                if (ipeta == 0 || ieta == ietamax-1) {
                    peta = ieta;
                } else {
                    peta = ieta + 1;
                }
                for (int iptau = 0; iptau < 2; iptau++) {
                    int ptau;
                    if (iptau == 0 || itau == itaumax) {
                        ptau = itau;
                    } else {
                        ptau = itau + 1;
                    }
                    position[ipx][ipy][ipeta][iptau] = (
                                px + ixmax*(py + iymax*(peta + ietamax*ptau)));
                }
            }
        }
    }

    // And now, the interpolation:
    float T = 0.0;
    float ed = 0.;
    float p = 0.;
    float utau = 1.0;
    float ux = 0.0;
    float uy = 0.0;
    float uz = 0.0;
    float ueta = 0.0;
    float pi00 = 0.0;
    float pi01 = 0.0;
    float pi02 = 0.0;
    float pi03 = 0.0;
    float pi11 = 0.0;
    float pi12 = 0.0;
    float pi13 = 0.0;
    float pi22 = 0.0;
    float pi23 = 0.0;
    float pi33 = 0.0;
    float bulkPi = 0.0;

    fluidCell_3D_ideal *HydroCell_3D_ideal_ptr1, *HydroCell_3D_ideal_ptr2;
    for (int iptau = 0; iptau < 2; iptau++) {
        float taufactor;
        if (iptau == 0)
            taufactor = 1. - taufrac;
        else
            taufactor = taufrac;
        for (int ipeta = 0; ipeta < 2; ipeta++) {
            float etafactor;
            if (ipeta == 0)
                etafactor = 1. - etafrac;
            else
                etafactor = etafrac;
            for (int ipy = 0; ipy < 2; ipy++) {
                float yfactor;
                if (ipy == 0)
                    yfactor = 1. - yfrac;
                else
                    yfactor = yfrac;

                float prefrac = yfactor*etafactor*taufactor;

                HydroCell_3D_ideal_ptr1 = (
                    &lattice_3D_ideal[idx_map_[position[0][ipy][ipeta][iptau]]]);
                HydroCell_3D_ideal_ptr2 = (
                    &lattice_3D_ideal[idx_map_[position[1][ipy][ipeta][iptau]]]);
                T += prefrac*(
                    (1. - xfrac)*HydroCell_3D_ideal_ptr1->temperature
                    + xfrac*HydroCell_3D_ideal_ptr2->temperature);
                ed += prefrac*(
                    (1. - xfrac)*HydroCell_3D_ideal_ptr1->ed
                    + xfrac*HydroCell_3D_ideal_ptr2->ed);
                p += prefrac*(
                    (1. - xfrac)*HydroCell_3D_ideal_ptr1->pressure
                    + xfrac*HydroCell_3D_ideal_ptr2->pressure);
                ux += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->ux
                                + xfrac*HydroCell_3D_ideal_ptr2->ux);
                uy += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->uy
                                + xfrac*HydroCell_3D_ideal_ptr2->uy);
                uz += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->uz
                               + xfrac*HydroCell_3D_ideal_ptr2->uz);
            }
        }
    }

    float sinh_eta = sinh(eta);
    float cosh_eta = cosh(eta);
    if (hydroType_ == 4) {
        // full (3+1)D case
        float ut = sqrt(1. + ux*ux + uy*uy + uz*uz);
        utau = ut*cosh_eta - uz*sinh_eta;
        ueta = - ut*sinh_eta + uz*cosh_eta;
    } else {
        // boost invariant
        utau = sqrt(1. + ux*ux + uy*uy);
        ueta = 0.;
    }

    info->temperature = T;
    info->utau = utau;
    info->ux = ux;
    info->uy = uy;
    info->ueta = ueta;

    info->ed = ed;
    info->sd = (ed + p)/(T + 1e-16);
    info->pressure = p;

    info->pi[0][0] = pi00;
    info->pi[0][1] = pi01;
    info->pi[0][2] = pi02;
    info->pi[0][3] = pi03;
    info->pi[1][0] = pi01;
    info->pi[1][1] = pi11;
    info->pi[1][2] = pi12;
    info->pi[1][3] = pi13;
    info->pi[2][0] = pi02;
    info->pi[2][1] = pi12;
    info->pi[2][2] = pi22;
    info->pi[2][3] = pi23;
    info->pi[3][0] = pi03;
    info->pi[3][1] = pi13;
    info->pi[3][2] = pi23;
    info->pi[3][3] = pi33;

    info->bulkPi = bulkPi;
    return;
}


void Hydroinfo_MUSIC::output_temperature_evolution(string filename_base) {
    fluidCell *hydroInfo = new fluidCell;
    for (int i = 0; i < itaumax; i++) {
        float tau = hydroTau0 + i*hydroDtau;
        ostringstream filename;
        filename << filename_base << "_tau_" << tau << ".dat";
        ofstream temp_evo(filename.str().c_str());
        for (int ix = 0; ix < ixmax; ix++) {
            float x_local = -hydroXmax + ix*hydroDx;
            for (int iy = 0; iy < iymax; iy++) {
                float y_local = -hydroXmax + iy*hydroDx;
                getHydroValues(x_local, y_local, 0.0, tau, hydroInfo);
                float temp_local = hydroInfo->temperature;
                temp_evo << scientific << setw(16) << setprecision(8)
                         << temp_local << "   ";
            }
            temp_evo << endl;
        }
        temp_evo.close();
    }
    delete hydroInfo;
}


void Hydroinfo_MUSIC::update_grid_info(
    float tau0, float tau_max, float dtau,
    float x_max, float dx, float eta_max, float deta) {
    hydroTau0 = tau0;
    hydroTauMax = tau_max;
    hydroDtau = dtau;
    hydroXmax = x_max;
    hydroDx = dx;
    hydro_eta_max = eta_max;
    hydroDeta = deta;
}
