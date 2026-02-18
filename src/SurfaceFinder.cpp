// Copyright @ Chun Shen 2014
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>

#include "SurfaceFinder.h"
#include "cornelius.h"

SurfaceFinder::SurfaceFinder(std::shared_ptr<Hydroinfo_MUSIC> hydroinfo_ptr_in,
                             ParameterReader* paraRdr_in) {
    paraRdr = paraRdr_in;
    hydroinfo_MUSIC_ptr = hydroinfo_ptr_in;
}

SurfaceFinder::~SurfaceFinder() {}

bool SurfaceFinder::check_intersect(double e_sw, double tau, double x,
                                    double y, double dtau, double dx, double dy,
                                    double ***cube) {
    fluidCell *fluidCellptr = new fluidCell();
    bool intersect = true;

    double tau_low = tau - dtau/2.;
    double tau_high = tau + dtau/2.;
    double x_left = x - dx/2.;
    double x_right = x + dx/2.;
    double y_left = y - dy/2.;
    double y_right = y + dy/2.;

    hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_left, 0.0, tau_low,
                                        fluidCellptr);
    cube[0][0][0] = fluidCellptr->ed;
    hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_right, 0.0, tau_low,
                                        fluidCellptr);
    cube[0][0][1] = fluidCellptr->ed;
    hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_left, 0.0, tau_low,
                                        fluidCellptr);
    cube[0][1][0] = fluidCellptr->ed;
    hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_right, 0.0, tau_low,
                                        fluidCellptr);
    cube[0][1][1] = fluidCellptr->ed;
    hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_left, 0.0, tau_high,
                                        fluidCellptr);
    cube[1][0][0] = fluidCellptr->ed;
    hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_right, 0.0, tau_high,
                                        fluidCellptr);
    cube[1][0][1] = fluidCellptr->ed;
    hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_left, 0.0, tau_high,
                                        fluidCellptr);
    cube[1][1][0] = fluidCellptr->ed;
    hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_right, 0.0, tau_high,
                                        fluidCellptr);
    cube[1][1][1] = fluidCellptr->ed;

    if ((e_sw - cube[0][0][0])*(cube[1][1][1] - e_sw) < 0.0)
        if ((e_sw - cube[0][1][0])*(cube[1][0][1] - e_sw) < 0.0)
            if ((e_sw - cube[0][1][1])*(cube[1][0][0] - e_sw) < 0.0)
                if ((e_sw - cube[0][0][1])*(cube[1][1][0] - e_sw) < 0.0)
                    intersect = false;

    delete fluidCellptr;
    return(intersect);
}


bool SurfaceFinder::check_intersect_3D(
        double e_sw, double tau, double x, double y, double eta_s,
        double dtau, double dx, double dy, double deta_s, double ****cube) {
    fluidCell *fluidCellptr = new fluidCell();
    bool intersect = true;

    for (int itau = 0; itau < 2; itau++) {
        for (int ix = 0; ix < 2; ix++) {
            for(int iy = 0; iy < 2; iy++) {
                for (int ieta = 0; ieta < 2; ieta++) {
                    hydroinfo_MUSIC_ptr->getHydroValues(
                            x + (ix - 0.5)*dx,
                            y + (iy - 0.5)*dy,
                            eta_s + (ieta - 0.5)*deta_s,
                            tau + itau*dtau,
                            fluidCellptr);
                    cube[itau][ix][iy][ieta] = fluidCellptr->ed;
                }
            }
        }
    }

    bool cond1 = (e_sw - cube[0][0][0][0])*(cube[1][1][1][1] - e_sw) < 0.0;
    bool cond2 = (e_sw - cube[0][1][0][0])*(cube[1][0][1][1] - e_sw) < 0.0;
    bool cond3 = (e_sw - cube[0][1][1][0])*(cube[1][0][0][1] - e_sw) < 0.0;
    bool cond4 = (e_sw - cube[0][1][0][1])*(cube[1][0][1][0] - e_sw) < 0.0;
    bool cond5 = (e_sw - cube[0][1][1][1])*(cube[1][0][0][0] - e_sw) < 0.0;
    bool cond6 = (e_sw - cube[0][0][1][0])*(cube[1][1][0][1] - e_sw) < 0.0;
    bool cond7 = (e_sw - cube[0][0][1][1])*(cube[1][1][0][0] - e_sw) < 0.0;
    bool cond8 = (e_sw - cube[0][0][1][0])*(cube[1][1][0][1] - e_sw) < 0.0;

    if (cond1 && cond2 && cond3 && cond4 && cond5 && cond6 && cond7 && cond8)
        intersect = false;

    delete fluidCellptr;
    return(intersect);
}


int SurfaceFinder::Find_full_hypersurface(double e_sw) {
    std::ofstream output;
    output.open("hyper_surface_2+1d.dat", std::ios::binary);

    double grid_tau0, grid_tauf, grid_x0, grid_y0;
    grid_tau0 = hydroinfo_MUSIC_ptr->get_hydro_tau0();
    grid_tauf = hydroinfo_MUSIC_ptr->get_hydro_tau_max();
    grid_x0 = (- hydroinfo_MUSIC_ptr->get_hydro_x_max()
               + hydroinfo_MUSIC_ptr->get_hydro_dx());
    grid_y0 = grid_x0;

    double grid_dt = paraRdr->getVal("grid_dt");
    double grid_dx = paraRdr->getVal("grid_dx");
    double grid_dy = paraRdr->getVal("grid_dy");

    int dim = 3;
    double *lattice_spacing = new double [dim];
    lattice_spacing[0] = grid_dt;
    lattice_spacing[1] = grid_dx;
    lattice_spacing[2] = grid_dy;

    Cornelius cornelius;
    cornelius.init(dim, e_sw, lattice_spacing);

    int ntime = static_cast<int>((grid_tauf - grid_tau0)/grid_dt);
    int nx = static_cast<int>(fabs(2.*grid_x0)/grid_dx);
    int ny = static_cast<int>(fabs(2.*grid_y0)/grid_dy);

    double ***cube = new double** [2];
    for (int i = 0; i < 2; i++) {
        cube[i] = new double* [2];
        for (int j = 0; j < 2; j++) {
            cube[i][j] = new double [2];
            for (int k = 0; k < 2; k++)
                cube[i][j][k] = 0.0;
        }
    }

    fluidCell *fluidCellptr = new fluidCell();

    const int FOsize = 34;
    const float etas = 0.0;
    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        double tau_local = grid_tau0 + (itime + 0.5)*grid_dt;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + (i + 0.5)*grid_dx;
            for (int j = 0; j < ny; j++) {
                double y_local = grid_y0 + (j + 0.5)*grid_dy;
                bool intersect = check_intersect(e_sw, tau_local, x_local,
                                                 y_local, grid_dt, grid_dx,
                                                 grid_dy, cube);
                if (intersect) {
                    cornelius.find_surface_3d(cube);
                    for (int isurf = 0; isurf < cornelius.get_Nelements();
                         isurf++) {
                        double tau_center = (
                                cornelius.get_centroid_elem(isurf, 0)
                                + tau_local - grid_dt/2.);
                        double x_center = (
                                cornelius.get_centroid_elem(isurf, 1)
                                + x_local - grid_dx/2.);
                        double y_center = (
                                cornelius.get_centroid_elem(isurf, 2)
                                + y_local - grid_dy/2.);

                        double da_tau = cornelius.get_normal_elem(isurf, 0);
                        double da_x = cornelius.get_normal_elem(isurf, 1);
                        double da_y = cornelius.get_normal_elem(isurf, 2);
                        hydroinfo_MUSIC_ptr->getHydroValues(
                            x_center, y_center, etas, tau_center,
                            fluidCellptr);

                        double eps_plus_p_over_T = (
                                fluidCellptr->ed + fluidCellptr->pressure)
                                / fluidCellptr->temperature;

                        std::vector<float> array(FOsize);
                        array[0] = static_cast<float>(tau_center);
                        array[1] = static_cast<float>(x_center);
                        array[2] = static_cast<float>(y_center);
                        array[3] = static_cast<float>(etas);
                        array[4] = static_cast<float>(da_tau);
                        array[5] = static_cast<float>(da_x);
                        array[6] = static_cast<float>(da_y);
                        array[7] = static_cast<float>(0.0);
                        array[8] = static_cast<float>(fluidCellptr->utau);
                        array[9] = static_cast<float>(fluidCellptr->ux);
                        array[10] = static_cast<float>(fluidCellptr->uy);
                        array[11] = static_cast<float>(fluidCellptr->ueta);
                        array[12] = static_cast<float>(fluidCellptr->ed);
                        array[13] = static_cast<float>(fluidCellptr->temperature);
                        array[14] = static_cast<float>(0.0);    // mu_B
                        array[15] = static_cast<float>(0.0);    // mu_S
                        array[16] = static_cast<float>(0.0);    // mu_Q
                        array[17] = static_cast<float>(eps_plus_p_over_T);
                        array[18] = static_cast<float>(fluidCellptr->pi[0][0]);
                        array[19] = static_cast<float>(fluidCellptr->pi[0][1]);
                        array[20] = static_cast<float>(fluidCellptr->pi[0][2]);
                        array[21] = static_cast<float>(fluidCellptr->pi[0][3]);
                        array[22] = static_cast<float>(fluidCellptr->pi[1][1]);
                        array[23] = static_cast<float>(fluidCellptr->pi[1][2]);
                        array[24] = static_cast<float>(fluidCellptr->pi[1][3]);
                        array[25] = static_cast<float>(fluidCellptr->pi[2][2]);
                        array[26] = static_cast<float>(fluidCellptr->pi[2][3]);
                        array[27] = static_cast<float>(fluidCellptr->pi[3][3]);
                        array[28] = static_cast<float>(fluidCellptr->bulkPi);
                        array[29] = static_cast<float>(0.0);    // rhob
                        array[30] = static_cast<float>(0.0);    // VB^0
                        array[31] = static_cast<float>(0.0);    // VB^1
                        array[32] = static_cast<float>(0.0);    // VB^2
                        array[33] = static_cast<float>(0.0);    // VB^3 
                        for (int i = 0; i < FOsize; i++)
                            output.write((char *)&(array[i]), sizeof(float));
                    }
                }
            }
        }
    }
    output.close();

    delete fluidCellptr;
    delete [] lattice_spacing;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++)
            delete [] cube[i][j];
        delete [] cube[i];
    }
    delete [] cube;
    return 0;
}


int SurfaceFinder::Find_full_hypersurface_3D(double e_sw) {
    std::ofstream output;
    output.open("hyper_surface_3+1d.dat", std::ios::binary);

    double grid_tau0, grid_tauf;
    grid_tau0 = hydroinfo_MUSIC_ptr->get_hydro_tau0();
    grid_tauf = hydroinfo_MUSIC_ptr->get_hydro_tau_max();
    double grid_x0, grid_y0, grid_etas0;
    grid_x0 = (- hydroinfo_MUSIC_ptr->get_hydro_x_max()
               + hydroinfo_MUSIC_ptr->get_hydro_dx());
    grid_y0 = (- hydroinfo_MUSIC_ptr->get_hydro_y_max()
               + hydroinfo_MUSIC_ptr->get_hydro_dy());
    grid_etas0 = (- hydroinfo_MUSIC_ptr->get_hydro_eta_max()
                 + hydroinfo_MUSIC_ptr->get_hydro_deta());

    double grid_dtau = paraRdr->getVal("grid_dt");
    double grid_dx = paraRdr->getVal("grid_dx");
    double grid_dy = paraRdr->getVal("grid_dy");
    double grid_deta = paraRdr->getVal("grid_deta");

    int dim = 4;
    double *lattice_spacing = new double [dim];
    lattice_spacing[0] = grid_dtau;
    lattice_spacing[1] = grid_dx;
    lattice_spacing[2] = grid_dy;
    lattice_spacing[3] = grid_deta;

    Cornelius cornelius;
    cornelius.init(dim, e_sw, lattice_spacing);

    const int ntime = static_cast<int>((grid_tauf - grid_tau0)/grid_dtau);
    const int nx = static_cast<int>(std::abs(2.*grid_x0)/grid_dx);
    const int ny = static_cast<int>(std::abs(2.*grid_y0)/grid_dy);
    const int neta = static_cast<int>(std::abs(2.*grid_etas0)/grid_deta);

    double ****cube = new double*** [2];
    for (int i = 0; i < 2; i++) {
        cube[i] = new double** [2];
        for (int j = 0; j < 2; j++) {
            cube[i][j] = new double* [2];
            for (int k = 0; k < 2; k++) {
                cube[i][j][k] = new double [2];
                for (int l = 0; l < 2; l++) {
                    cube[i][j][k][l] = 0.0;
                }
            }
        }
    }

    fluidCell *fluidCellptr = new fluidCell();

    const int FOsize = 34;
    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        double tau_local = grid_tau0 + (itime + 0.5)*grid_dtau;
        for (int k = 0; k < neta; k++) {
            double eta_local = grid_etas0 + (k + 0.5)*grid_deta;
            for (int i = 0; i < nx; i++) {
                // loops over the transverse plane
                double x_local = grid_x0 + (i + 0.5)*grid_dx;
                for (int j = 0; j < ny; j++) {
                    double y_local = grid_y0 + (j + 0.5)*grid_dy;
                    bool intersect = check_intersect_3D(
                            e_sw, tau_local, x_local, y_local, eta_local,
                            grid_dtau, grid_dx, grid_dy, grid_deta, cube);

                    if (!intersect) continue;

                    cornelius.find_surface_4d(cube);
                    for (int isurf = 0; isurf < cornelius.get_Nelements();
                         isurf++) {
                        double tau_center = (
                                cornelius.get_centroid_elem(isurf, 0)
                                + tau_local - grid_dtau/2.);
                        double x_center = (
                                cornelius.get_centroid_elem(isurf, 1)
                                + x_local - grid_dx/2.);
                        double y_center = (
                                cornelius.get_centroid_elem(isurf, 2)
                                + y_local - grid_dy/2.);
                        double eta_center = (
                                cornelius.get_centroid_elem(isurf, 3)
                                + eta_local - grid_deta/2.);

                        double da_tau = cornelius.get_normal_elem(isurf, 0);
                        double da_x = cornelius.get_normal_elem(isurf, 1);
                        double da_y = cornelius.get_normal_elem(isurf, 2);
                        double da_eta = cornelius.get_normal_elem(isurf, 3);

                        hydroinfo_MUSIC_ptr->getHydroValues(
                            x_center, y_center, eta_center, tau_center,
                            fluidCellptr);

                        double eps_plus_p_over_T = (
                                fluidCellptr->ed + fluidCellptr->pressure)
                                / fluidCellptr->temperature;

                        std::vector<float> array(FOsize);
                        array[0] = static_cast<float>(tau_center);
                        array[1] = static_cast<float>(x_center);
                        array[2] = static_cast<float>(y_center);
                        array[3] = static_cast<float>(eta_center);
                        array[4] = static_cast<float>(da_tau);
                        array[5] = static_cast<float>(da_x);
                        array[6] = static_cast<float>(da_y);
                        array[7] = static_cast<float>(da_eta);
                        array[8] = static_cast<float>(fluidCellptr->utau);
                        array[9] = static_cast<float>(fluidCellptr->ux);
                        array[10] = static_cast<float>(fluidCellptr->uy);
                        array[11] = static_cast<float>(fluidCellptr->ueta);
                        array[12] = static_cast<float>(fluidCellptr->ed);
                        array[13] = static_cast<float>(fluidCellptr->temperature);
                        array[14] = static_cast<float>(0.0);    // mu_B
                        array[15] = static_cast<float>(0.0);    // mu_S
                        array[16] = static_cast<float>(0.0);    // mu_Q
                        array[17] = static_cast<float>(eps_plus_p_over_T);
                        array[18] = static_cast<float>(fluidCellptr->pi[0][0]);
                        array[19] = static_cast<float>(fluidCellptr->pi[0][1]);
                        array[20] = static_cast<float>(fluidCellptr->pi[0][2]);
                        array[21] = static_cast<float>(fluidCellptr->pi[0][3]);
                        array[22] = static_cast<float>(fluidCellptr->pi[1][1]);
                        array[23] = static_cast<float>(fluidCellptr->pi[1][2]);
                        array[24] = static_cast<float>(fluidCellptr->pi[1][3]);
                        array[25] = static_cast<float>(fluidCellptr->pi[2][2]);
                        array[26] = static_cast<float>(fluidCellptr->pi[2][3]);
                        array[27] = static_cast<float>(fluidCellptr->pi[3][3]);
                        array[28] = static_cast<float>(fluidCellptr->bulkPi);
                        array[29] = static_cast<float>(0.0);    // rhob
                        array[30] = static_cast<float>(0.0);    // VB^0
                        array[31] = static_cast<float>(0.0);    // VB^1
                        array[32] = static_cast<float>(0.0);    // VB^2
                        array[33] = static_cast<float>(0.0);    // VB^3 
                        for (int i = 0; i < FOsize; i++)
                            output.write((char *)&(array[i]), sizeof(float));
                    }
                }
            }
        }
    }
    output.close();

    delete fluidCellptr;
    delete [] lattice_spacing;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++)
                delete [] cube[i][j][k];
            delete [] cube[i][j];
        }
        delete [] cube[i];
    }
    delete [] cube;
    return 0;
}
