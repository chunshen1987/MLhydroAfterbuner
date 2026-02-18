// Copyright Chun Shen @ 2015
#ifndef SRC_SurfaceFinder_H_
#define SRC_SurfaceFinder_H_

#include<memory>

#include "Hydroinfo_MUSIC.h"
#include "ParameterReader.h"

class SurfaceFinder {
 private:
    //int hydro_type;
    std::shared_ptr<Hydroinfo_MUSIC> hydroinfo_MUSIC_ptr;
    ParameterReader *paraRdr;
    //double e_sw_;

 public:
    SurfaceFinder() = delete;

    SurfaceFinder(std::shared_ptr<Hydroinfo_MUSIC> hydroinfo_ptr_in,
                  ParameterReader* paraRdr_in);
    ~SurfaceFinder();

    bool check_intersect(double e_sw, double tau, double x, double y,
                         double dt, double dx, double dy, double ***cube);
    int Find_full_hypersurface(double e_sw);
};

#endif  // SRC_SurfaceFinder_H_
