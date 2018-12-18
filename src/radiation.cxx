/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "radiation.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "input.h"
#include "data_block.h"
#include "stats.h"
#include "cross.h"
#include "dump.h"
#include "column.h"
#include "field3d_operators.h"
#include "timeloop.h"
#include <time.h>
namespace
{
    /*
    // Wrapper functions to the RRTMG long wave kernels.
    extern "C"
    {
        void c_rrtmg_lw_init(double *cpdair);
        void c_rrtmg_lw (
                int *ncol    ,int *nlay    ,int *icld    ,int *idrv    ,
                double *play    ,double *plev    ,double *tlay    ,double *tlev    ,double *tsfc    ,
                double *h2ovmr  ,double *o3vmr   ,double *co2vmr  ,double *ch4vmr  ,double *n2ovmr  ,double *o2vmr,
                double *cfc11vmr,double *cfc12vmr,double *cfc22vmr,double *ccl4vmr ,double *emis    ,
                int *inflglw ,int *iceflglw,int *liqflglw,double *cldfr   ,
                double *taucld  ,double *cicewp  ,double *cliqwp  ,double *reice   ,double *reliq   ,
                double *tauaer  ,
                double *uflx    ,double *dflx    ,double *hr      ,double *uflxc   ,double *dflxc,  double *hrc,
                double *duflx_dt,double *duflxc_dt );
    }
    */

    template<typename TF>
    void calc_gcss_rad_SW(const TF* const restrict ql, const TF* const restrict qt, const TF* const restrict rhoref,
        const TF* const z, const TF* const dzi,
        const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
        const int icells, const int ijcells
    )
    {
        const int jj = icells;
        const int kk = ijcells;
        const TF rho_l = 1000.;
        const TF reff = 1.E-5;
        TF tauc;
        TF fact;
        int ki; //PBLH index
        std::vector<TF> tau;
        tau.resize(kend); //kcells so that it takes care of the ghost cells
        const TF mu = 0.05;//zenith(32.5,time); //zenith
        for (int j=jstart; j<jend; ++j)
        {
            for (int i=istart; i<iend; ++i)
            {
                if (mu>0.035)
                {
                    tauc = TF(0.0);
                    for (int k=kstart;k<kend;++k)
                    {
                        const int ij   = i + j*jj;
                        const int ijk  = i + j*jj + k*kk;
                        const int km1 = std::max(1,k-1);
                        tau[k] = TF(0.0);
                        if (ql[ijk]>1.E-5)
                        {
                            tau[k] = std::max(TF(0.0) , TF(1.5) * ql[ijk] * rhoref[k] * (z[k]-z[km1]) / reff / rho_l);
                            tauc = tauc + tau[k];
                        }
                    }
                    //sunray
                    // swn = 1.0; //no SW for now
                } //end if mu
            }
        }
    }

    template<typename TF>
    void calc_gcss_rad_LW(const TF* const restrict ql, const TF* const restrict qt,
    TF* const restrict lwp, TF* const restrict flx, const TF* const restrict rhoref,
    const TF* const z, const TF* const dzi,
    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
    const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;
        const TF xka = 85.;
        const TF fr0 = 70.;
        const TF fr1 = 22.;
        const TF cp = 1005; //can read this from constant.h
        const TF div = 3.75E-6; //fix divergence for now
        int ki; //pblh index
        TF fact;
        for (int j=jstart; j<jend; ++j)
        {
            for (int i=istart; i<iend; ++i)
            {
                lwp[i+j*jj] = TF(0.0);
                ki = kend; //set to top of domain
                for (int k=kstart; k<kend; ++k)
                {
                    const int ij   = i + j*jj;
                    const int ijk  = i + j*jj + k*kk;
                    const int km1 = std::max(1,k-1);
                    lwp[ij] = lwp[ij] + std::max( TF(0.0) , ql[ijk] * rhoref[k] * (z[k]-z[km1]));
                    flx[ijk] = fr1 * std::exp(TF(-1.0) * xka * lwp[ij]);
                    if ( (ql[ijk] > TF(0.01E-3) ) && ( qt[ijk] >= TF(0.008) ) ) ki = k; //this is the PBLH index
                }
                fact = div * cp * rhoref[ki];
                const int ij   = i + j*jj;
                flx[ij + kstart*kk] = flx[ij + kstart*kk] + fr0 * std::exp(TF(-1.0) * xka *lwp[ij]);
                for (int k=kstart+1;k<kend;++k)
                {
                    const int ij   = i + j*jj;
                    const int ijk  = i + j*jj + k*kk;
                    const int km1 = std::max(kstart+1,k-1);
                    const int ijkm = i + j*jj + km1*kk;
                    lwp[ij] = lwp[ij] - std::max( TF(0.0) , ql[ijk] * rhoref[k] * (z[k]-z[k-1]));
                    flx[ijk] = flx[ijk] + fr0 * std::exp(-1.0 * xka * lwp[ij]);
                    if ((k>ki) && (ki>1) && (fact>0.))
                    { //above PBLH
                        flx[ijk] = flx[ijk] + fact * ( TF(0.25) * std::pow(z[k]-z[ki],TF(1.333)) + z[ki] * std::pow(z[k]-z[ki],TF(0.33333)) );
                    }
                }
            } // end of i
        } // end of j
    }
    template<typename TF> //EW: simplified radiative parameterization for LW and SW fluxes for DYCOMS
    void exec_gcss_rad(
            TF* const restrict tt, const TF* const restrict ql, const TF* const restrict qt,
            TF* const restrict lwp, TF* const restrict flx, const TF* const restrict rhoref,
            const TF* const z, const TF* const dzi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;
        const TF mu = 0.05;//zenith(32.5,time); //zenith
        const TF cp = 1005; //can read this from constant.h

        //call LW
        calc_gcss_rad_LW<TF>(ql,qt,
        lwp,flx,rhoref,
        z,dzi,
        istart,iend,jstart,jend,kstart,kend,
        icells,ijcells);
        for (int j=jstart; j<jend; ++j)
        {
            for (int i=istart; i<iend; ++i)
            {
                // if (mu>0.035)
                // {
                //     //call shortwave
                // } //end if mu
                for (int k=kstart+1;k<kend;++k)
                {
                    const int ij   = i + j*jj;
                    const int ijk  = i + j*jj + k*kk;
                    const int km1 = std::max(kstart+1,k-1);
                    const int ijkm = i + j*jj + km1*kk;
                    tt[ijk] = tt[ijk] - (flx[ijk] - flx[ijkm]) * dzi[k] / (rhoref[k] * cp);
                    // tt[ijk] = tt[ijk]+(swn[ijk]-swn[ijkm])*dzh[k]/(fields.rhoref[k]*cp); //no SW for now
                }
            } // end of i
        } // end of j
    } // end of calc_gcss_rad
}

template<typename TF>
Radiation<TF>::Radiation(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(master, grid, fields)
{
    // Read the switches from the input
    std::string swradiation_in = inputin.get_item<std::string>("radiation", "swradiation", "", "0");

    if (swradiation_in == "0")
        swradiation = Radiation_type::Disabled;
    else if (swradiation_in == "1")
        swradiation = Radiation_type::Enabled;
    else if (swradiation_in == "2")
        swradiation = Radiation_type::Gcss;
    else
        throw std::runtime_error("Invalid option for \"swradiation\"");

    double cp = 1004.;
    // c_rrtmg_lw_init(&cp);
    ncol = 1;
    nlay = 60;
    nbndlw = 16;
}

template<typename TF>
Radiation<TF>::~Radiation()
{
}

template<typename TF>
void Radiation<TF>::init()
{
    auto& gd = grid.get_grid_data();
    if (swradiation == Radiation_type::Disabled)
        return;
    if (swradiation == Radiation_type::Enabled)
    {
    play.resize(ncol*nlay);     // (ncol, nlay)
    plev.resize(ncol*(nlay+1)); // (ncol, nlay+1)
    tlay.resize(ncol*nlay);     // (ncol, nlay)
    tlev.resize(ncol*(nlay+1)); // (ncol, nlay+1)

    tsfc.resize(ncol); // (ncol)

    h2ovmr.resize(ncol*nlay); // (ncol, nlay)
    o3vmr.resize(ncol*nlay);  // (ncol, nlay)
    co2vmr.resize(ncol*nlay); // (ncol, nlay)
    ch4vmr.resize(ncol*nlay); // (ncol, nlay)
    n2ovmr.resize(ncol*nlay); // (ncol, nlay)
    o2vmr.resize(ncol*nlay);  // (ncol, nlay)

    cfc11vmr.resize(ncol*nlay); // (ncol, nlay)
    cfc12vmr.resize(ncol*nlay); // (ncol, nlay)
    cfc22vmr.resize(ncol*nlay); // (ncol, nlay)
    ccl4vmr.resize(ncol*nlay);  // (ncol, nlay)
    emis.resize(ncol*nbndlw);   // (ncol, nbndlw)

    cldfr.resize(ncol*nlay);  // (ncol, nlay)
    cicewp.resize(ncol*nlay); // (ncol, nlay)
    cliqwp.resize(ncol*nlay); // (ncol, nlay)
    reice.resize(ncol*nlay);  // (ncol, nlay)
    reliq.resize(ncol*nlay);  // (ncol, nlay)
    taucld.resize(nbndlw*ncol*nlay); // (nbndlw, ncol, nlay)
    tauaer.resize(ncol*nlay*nbndlw); // (ncol, nlay, nbndlw)

    // OUTPUT
    uflx.resize(ncol*(nlay+1));      // (ncol, nlay+1)
    dflx.resize(ncol*(nlay+1));      // (ncol, nlay+1)
    hr.resize(ncol*nlay);            // (ncol, nlay)
    uflxc.resize(ncol*(nlay+1));     // (ncol, nlay+1)
    dflxc.resize(ncol*(nlay+1));     // (ncol, nlay+1)
    hrc.resize(ncol*nlay);           // (ncol, nlay)
    duflx_dt.resize(ncol*(nlay+1));  // (ncol, nlay+1)
    duflxc_dt.resize(ncol*(nlay+1)); // (ncol, nlay+1)
    }
}

template<typename TF>
void Radiation<TF>::create(Thermo<TF>& thermo,Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross, Dump<TF>& dump)
{
    if (swradiation == Radiation_type::Disabled)
        return;
    if (swradiation == Radiation_type::Enabled)
    {
    std::string block_name = "radiation.prof";
    Data_block data_block(master, block_name);

    std::vector<TF> p_rad(nlay);
    data_block.get_vector(p_rad, "pavel", nlay, 0, 0);
    std::vector<TF> ph_rad(nlay+1);
    data_block.get_vector(ph_rad, "pz", nlay+1, 0, 0);

    // Convert the data to
    for (auto& d : p_rad)
        d *= 100.;

    for (auto& d : ph_rad)
        d *= 100.;

    // Get copies of pressure profiles.
    std::vector<TF> p = thermo.get_p_vector();
    std::vector<TF> ph = thermo.get_ph_vector();

    auto& gd = grid.get_grid_data();

    // Remove the ghost cells from the vector.
    p.erase(p.begin() + gd.kend, p.end());
    ph.erase(ph.begin() + gd.kend+1, ph.end());
    p.erase(p.begin(), p.begin() + gd.kstart);
    ph.erase(ph.begin(), ph.begin() + gd.kstart);

    // Remove all elements from radiation profile vector that are in the domain.
    auto it = p_rad.begin();
    int counter = 0;
    while (*it > ph.back())
    {
        ++counter;
        ++it;
    }

    // Delete all items until first index, remove one level extra for ph
    p_rad.erase(p_rad.begin(), p_rad.begin()+counter);
    ph_rad.erase(ph_rad.begin(), ph_rad.begin()+counter+1);
    p_rad.insert(p_rad.begin(), p.begin(), p.end());
    ph_rad.insert(ph_rad.begin(), ph.begin(), ph.end());

    // CvH if the first pressure level above the domain top is too close to the half level, it could be moved...

    // Temporary...
    data_block.get_vector(play, "pavel", nlay, 0, 0);
    data_block.get_vector(plev, "pz", nlay+1, 0, 0);
    data_block.get_vector(tlay, "tavel", nlay, 0, 0);
    data_block.get_vector(tlev, "tz", nlay+1, 0, 0);

    // The specific humidity comes from climatology for now.
    data_block.get_vector(h2ovmr, "wkl1", nlay, 0, 0);

    // These profiles come from climatology
    data_block.get_vector(co2vmr, "wkl2", nlay, 0, 0);
    data_block.get_vector(o3vmr, "wkl3", nlay, 0, 0);
    data_block.get_vector(n2ovmr, "wkl4", nlay, 0, 0);
    data_block.get_vector(ch4vmr, "wkl6", nlay, 0, 0);
    data_block.get_vector(o2vmr, "wkl7", nlay, 0, 0);
    }
    if (swradiation == Radiation_type::Gcss) //only write stats if it's gcss for now
    {
        // Set up output classes
        create_stats(stats);
        create_column(column);
        create_dump(dump);
        create_cross(cross);
    }
}

template<typename TF>
void Radiation<TF>::exec(Thermo<TF>& thermo, double time, Timeloop<TF>& timeloop)
{
    if (swradiation == Radiation_type::Disabled)
        return;
    auto& gd = grid.get_grid_data();
    if (swradiation == Radiation_type::Enabled)
    {
        // For now...
        inflglw = 0;
        iceflglw = 0;
        liqflglw = 0;

        std::fill(emis.begin(), emis.end(), 1.);

        // Step 1. Get the mean absolute atmospheric temperature.
        auto T = fields.get_tmp();
        thermo.get_thermo_field(*T, "T", false, false);

        // Calculate radiative cooling only for single column.
        field3d_operators.calc_mean_profile(T->fld_mean.data(), T->fld.data());

        thermo.get_thermo_field(*T, "T_h", false, false);

        // Step 2.

        // Get absolute atmospheric and surface temperature.
        tsfc[0] = 300.;

        /*
        c_rrtmg_lw(
                &ncol    ,&nlay    ,&icld    ,&idrv    ,
                play.data()    ,plev.data()    ,tlay.data()    ,tlev.data()    ,tsfc.data()    ,
                h2ovmr.data()  ,o3vmr.data()   ,co2vmr.data()  ,ch4vmr.data()  ,n2ovmr.data()  ,o2vmr.data(),
                cfc11vmr.data(),cfc12vmr.data(),cfc22vmr.data(),ccl4vmr.data() ,emis.data()    ,
                &inflglw ,&iceflglw,&liqflglw,cldfr.data() ,
                taucld.data()  ,cicewp.data()  ,cliqwp.data()  ,reice.data()   ,reliq.data() ,
                tauaer.data()  ,
                uflx.data()    ,dflx.data()    ,hr.data()      ,uflxc.data()   ,dflxc.data(),  hrc.data(),
                duflx_dt.data(),duflxc_dt.data());
                */

        fields.release_tmp(T);

        std::cout << "Heating rate" << std::endl;
        for (int i=0; i<nlay; ++i)
            std::cout << i       << ", "
                << play[i] << ", "
                << hr[i] << ", "
                << std::endl;

        std::cout << "Upflux/downflux rate" << std::endl;
        for (int i=0; i<nlay+1; ++i)
            std::cout << i       << ", "
                << plev[i] << ", "
                << uflx[i] << ", "
                << dflx[i] << ", "
                << std::endl;
    }


    if (swradiation == Radiation_type::Gcss)
    {
        auto lwp = fields.get_tmp();
        auto flx = fields.get_tmp();
        auto ql  = fields.get_tmp();
        struct tm current_datetime;
        current_datetime = timeloop.get_realtime();
        mktime ( &current_datetime ); //refresh time
        std::cout << "Current month: " << current_datetime.tm_mon
                  << ", Day: " << current_datetime.tm_mday
                  << ", Hour: " << current_datetime.tm_hour
                  << ", Min: " << current_datetime.tm_min <<
                  << ", Sec: " << current_datetime.tm_sec << std::endl;
        exec_gcss_rad<TF>(
            fields.st.at("thl")->fld.data(), ql->fld.data(), fields.sp.at("qt")->fld.data(),
            lwp->fld.data(), flx->fld.data(), fields.rhoref.data(),
            gd.z.data(), gd.dzhi.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
        fields.release_tmp(lwp);
        fields.release_tmp(flx);
        fields.release_tmp(ql);
    }
}
template<typename TF>
bool Radiation<TF>::check_field_exists(const std::string name)
{
    if (name == "rflx")
        return true;
    else
        return false;
}
template<typename TF>
void Radiation<TF>::get_radiation_field(Field3d<TF>& fld, std::string name, Thermo<TF>& thermo)
{
    if (name == "rflx")
    {
        auto& gd = grid.get_grid_data();

        const TF no_offset = 0.;
        const TF no_threshold = 0.;

        auto lwp = fields.get_tmp();
        auto ql  = fields.get_tmp();
        thermo.get_thermo_field(*ql,"ql",false,false);
        calc_gcss_rad_LW(ql->fld.data(), fields.ap.at("qt")->fld.data(),
        lwp->fld.data(), fld.fld.data(), fields.rhoref.data(),
        gd.z.data(), gd.dzi.data(),
        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
        fields.release_tmp(lwp);
        fields.release_tmp(ql);
    }
}

template<typename TF>
void Radiation<TF>::create_stats(Stats<TF>& stats)
{
    if (stats.get_switch())
    {
        stats.add_prof("rflx", "Total radiative flux", "W m-2", "z");
        //when full radiation is available, add the following:
        // stats.add_prof("sflxd", "Downward shortwave radiative flux", "W m-2", "z");
        // stats.add_prof("sflxu", "Upward shortwave radiative flux", "W m-2", "z");
        // stats.add_prof("lflxd", "Downward longwave radiative flux", "W m-2", "z");
        // stats.add_prof("lflxu", "Upward longwave radiative flux", "W m-2", "z");
    }
}

template<typename TF>
void Radiation<TF>::create_column(Column<TF>& column)
{
    // add the profiles to the columns
    if (column.get_switch())
    {
        column.add_prof("rflx", "Total radiative flux", "W m-2", "z");
    }
}

template<typename TF>
void Radiation<TF>::create_cross(Cross<TF>& cross)
{
    if (cross.get_switch())
    {
        swcross_rflx = false;

        // Vectors with allowed cross variables for radiative flux
        std::vector<std::string> allowed_crossvars_rflx = {"rflx"};

        std::vector<std::string> rflxvars  = cross.get_enabled_variables(allowed_crossvars_rflx);
        if (rflxvars.size() > 0)
            swcross_rflx  = true;
        crosslist = rflxvars;
    }
}

template<typename TF>
void Radiation<TF>::create_dump(Dump<TF>& dump)
{
    if (dump.get_switch())
    {
        // Get global cross-list from cross.cxx
        std::vector<std::string> *dumplist_global = dump.get_dumplist();

        // Check if fields in dumplist are retrievable thermo fields
        std::vector<std::string>::iterator dumpvar=dumplist_global->begin();
        while (dumpvar != dumplist_global->end())
        {
            if (check_field_exists(*dumpvar))
            {
                // Remove variable from global list, put in local list
                dumplist.push_back(*dumpvar);
                dumplist_global->erase(dumpvar); // erase() returns iterator of next element..
            }
            else
                ++dumpvar;
        }
    }
}

template<typename TF>
void Radiation<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    auto flx = fields.get_tmp();
    flx->loc = gd.sloc;
    get_radiation_field(*flx,"rflx",thermo);

    //if daytime, rflx = LW + SW
    // calculate the mean
    std::vector<std::string> operators = {"mean"}; //add 2nd moment, if needed
    stats.calc_stats("rflx", *flx, no_offset, no_threshold, operators);

    fields.release_tmp(flx);
}

template<typename TF>
void Radiation<TF>::exec_column(Column<TF>& column, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();
    const TF no_offset = 0.;

    auto flx = fields.get_tmp();
    flx->loc = gd.sloc;
    get_radiation_field(*flx,"rflx",thermo);

    column.calc_column("rflx", flx->fld.data(), no_offset);

    fields.release_tmp(flx);
}

template<typename TF>
void Radiation<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();
    auto flx = fields.get_tmp();
    flx->loc = gd.sloc;

    if(swcross_rflx)
    {
        get_radiation_field(*flx,"rflx",thermo);
    }
    for (auto& it : crosslist)
    {
        if (it == "rflx")
            cross.cross_simple(flx->fld.data(), "rflx", iotime);
    }
    fields.release_tmp(flx);
}

template<typename TF>
void Radiation<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();
    auto flx = fields.get_tmp();

    for (auto& it : dumplist)
    {
        if (it == "rflx")
            {
                get_radiation_field(*flx,"rflx",thermo);
            }
        else
        {
            master.print_error("Radiation dump of field \"%s\" not supported\n", it.c_str());
            throw std::runtime_error("Error in Radiation Dump");
        }
        dump.save_dump(flx->fld.data(), it, iotime);
    }
    fields.release_tmp(flx);
}

template class Radiation<double>;
template class Radiation<float>;
