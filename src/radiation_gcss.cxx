/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
 * Copyright (c) 2018-2019 Elynn Wu
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

 #include "radiation_gcss.h"

 namespace
 {
     template<typename TF>
	 TF calc_zenith(struct tm datetime, TF lat, TF lon)
     {
         const TF pi        = M_PI;
         const TF year2days = 365.;
         const TF piAngle   = 180.;
         const TF day2secs  = 86400.;
         const TF z1        = 279.934;
         const TF z2        = 1.914827;
         const TF z3        = 0.7952;
         const TF z4        = 0.019938;
         const TF z5        = 0.00162;
         const TF z6        = 23.4439;
         TF time2sec = (datetime.tm_yday + 1) +
                            lon / 360. +
                           (datetime.tm_hour * 3600. +
                            datetime.tm_min * 60. +
                            datetime.tm_sec) / day2secs;
         TF day    = floor(time2sec);
         TF lamda  = lat * pi / piAngle;
         TF d      = 2. * pi * int(time2sec) / year2days;
         TF sig    = d + pi/piAngle * (z1 + z2*std::sin(d)
                                                   - z3*std::cos(d)
                                                   + z4*std::sin(2.*d)
                                                   - z5*std::cos(2.*d));
         TF del     = std::asin(std::sin(z6*pi / piAngle)*std::sin(sig));
         TF h       = 2. * pi * ((time2sec - day) - 0.5);
         TF mu      = std::sin(lamda) * std::sin(del) + std::cos(lamda) * std::cos(del) * std::cos(h);
         return mu;
     }

     template<typename TF>
     void sunray(const TF mu, const int i, const int j,
         const int kstart, const int kend, const int icells, const int ijcells,
         std::vector<TF> tau, const TF tauc,
         TF* const restrict swn)
     {
         const int jj = icells;
         const int kk = ijcells;
         TF o_c1 = TF(0.9);
         TF o_c2 = TF(2.75);
         TF o_c3 = TF(0.09);
         TF sw0 = TF(1100.);
         TF gc  = TF(0.85);
         TF sfc_albedo = TF(0.05);
         TF taucde = TF(0.);
         TF taupath = TF(0.);
         std::vector<TF> taude(kend , TF(0.));
         TF omega  = TF(1.) - TF(1.e-3) * (o_c1 + o_c2 * (mu+TF(1.)) * std::exp(-o_c3 * tauc)); //fouquart
         TF ff     = gc * gc;
         TF gcde   = gc / (TF(1.) + gc);
         taucde = ( TF(1.0) - omega*ff) * tauc;
         for (int k=kstart; k<kend; ++k)
         {
             taude[k] = ( TF(1.) - omega*ff ) * tau[k];
         }
         TF omegade = (TF(1.)-ff) * omega/(TF(1.) - omega*ff);
         TF x1  = TF(1.) - omegade * gcde;
         TF x2  = TF(1.) - omegade;
         TF rk  = std::sqrt(TF(3.) * x2 * x1);
         TF mu2 = mu * mu;
         TF x3  = TF(4.) * (TF(1.) - rk*rk*mu2);
         TF rp  = std::sqrt(TF(3.) * x2/x1);
         TF alpha = TF(3.) * omegade * mu2 * (TF(1.) + gcde*x2) / x3;
         TF beta  = TF(3.) * omegade * mu * (TF(1.) + TF(3.)*gcde*mu2*x2) / x3;

         TF rtt = TF(2.0/3.0);
         TF exmu0 = std::exp(-taucde / mu);
         TF expk  = std::exp(rk * taucde);
         TF exmk  = TF(1.) / expk;
         TF xp23p = TF(1.) + rtt*rp;
         TF xm23p = TF(1.) - rtt*rp;
         TF ap23b = alpha + rtt*beta;

         TF t1 = TF(1.) - sfc_albedo - rtt * (TF(1.) + sfc_albedo) * rp;
         TF t2 = TF(1.) - sfc_albedo + rtt * (TF(1.) + sfc_albedo) * rp;
         TF t3 = (TF(1.) - sfc_albedo) * alpha - rtt * (TF(1.) + sfc_albedo) * beta + sfc_albedo*mu;
         TF c2 = (xp23p*t3*exmu0 - t1*ap23b*exmk) / (xp23p*t2*expk - xm23p*t1*exmk);
         TF c1 = (ap23b - c2*xm23p)/xp23p;

         for (int k=kend-1;k>=kstart;--k)
         {
             const int ijk  = i + j*jj + k*kk;
             taupath = taupath + taude[k];
             swn[ijk] = sw0 * TF(4./3.) * (rp * (c1*std::exp(-rk*taupath)
             - c2 * std::exp(rk*taupath)) - beta * std::exp(-taupath/mu))
             + mu * sw0 * std::exp(-taupath / mu);
         }
     }

     template<typename TF>
     void calc_gcss_rad_SW(TF* const restrict swn, const TF* const restrict ql, const TF* const restrict qt,
         const TF* const restrict rhoref, const TF* const z, const TF* const dzi,
         const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
         const int icells, const int ijcells, const int ncells, TF mu)
     {
         const int jj = icells;
         const int kk = ijcells;
         const TF rho_l = 1000.;
         const TF reff = 1.E-5;
         TF tauc;
         TF fact;
         int ki; //PBLH index
         std::vector<TF> tau(kend,TF(0.));
         for (int n=0; n<ncells; ++n)
             swn[n] = TF(0.); //initialize as 0 otherwise weird things might be stored
         for (int j=jstart; j<jend; ++j)
         {
             for (int i=istart; i<iend; ++i)
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
                 sunray<TF>(TF(mu), i, j,
                     kstart, kend, icells, ijcells,
                     tau, tauc, swn);
             }
         }
     }

     template<typename TF>
     void calc_gcss_rad_LW(const TF* const restrict ql, const TF* const restrict qt,
     TF* const restrict lwp, TF* const restrict flx, const TF* const restrict rhoref,
     const TF fr0, const TF fr1, const TF xka, const TF div,
     const TF* const z, const TF* const dzi,
     const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
     const int icells, const int ijcells)
     {
         const int jj = icells;
         const int kk = ijcells;
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
                 fact = div * Constants::cp<TF> * rhoref[ki];
                 const int ij   = i + j*jj;
                 flx[ij + kstart*kk] = flx[ij + kstart*kk] + fr0 * std::exp(TF(-1.0) * xka *lwp[ij]);
                 for (int k=kstart+1; k<kend; ++k)
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
             TF* const restrict lwp, TF* const restrict flx, TF* const restrict swn, const TF* const restrict rhoref,
             const TF mu, const TF mu_min, const TF fr0, const TF fr1, const TF xka, const TF div,
             const TF* const z, const TF* const dzi,
             const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
             const int icells, const int ijcells, const int ncells)
     {
         const int jj = icells;
         const int kk = ijcells;
         //call LW
         calc_gcss_rad_LW<TF>(ql,qt,
         lwp,flx,rhoref, fr0, fr1, xka, div,
         z,dzi,
         istart,iend,jstart,jend,kstart,kend,
         icells,ijcells);
         for (int j=jstart; j<jend; ++j)
         {
             for (int i=istart; i<iend; ++i)
             {
                 for (int k=kstart+1;k<kend;++k)
                 {
                     const int ij   = i + j*jj;
                     const int ijk  = i + j*jj + k*kk;
                     const int km1 = std::max(kstart+1,k-1);
                     const int ijkm = i + j*jj + km1*kk;
                     tt[ijk] = tt[ijk] - (flx[ijk] - flx[ijkm]) * dzi[k] / (rhoref[k] * Constants::cp<TF>);
                 }
             } // end of i
         } // end of j
         if (mu>mu_min) //if daytime, call SW
         {
             calc_gcss_rad_SW<TF>(swn, ql, qt,
                 rhoref, z, dzi,
                 istart, iend, jstart, jend, kstart, kend,
                 icells, ijcells, ncells, mu);
             for (int j=jstart; j<jend; ++j)
             {
                 for (int i=istart; i<iend; ++i)
                 {
                     for (int k=kstart+1;k<kend;++k)
                     {
                         const int ij   = i + j*jj;
                         const int ijk  = i + j*jj + k*kk;
                         const int km1 = std::max(kstart+1,k-1);
                         const int ijkm = i + j*jj + km1*kk;
                         tt[ijk] = tt[ijk] + (swn[ijk] - swn[ijkm]) * dzi[k] / (rhoref[k] * Constants::cp<TF>);
                     }
                 }
             }
         } //end if mu
     } // end of calc_gcss_rad
 }

template<typename TF>
Radiation_gcss<TF>::Radiation_gcss(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
	Radiation<TF>(masterin, gridin, fieldsin, inputin)
{
	swradiation = "gcss"; //2 for gcss
    lat = inputin.get_item<TF>("radiation", "lat", "");
    lon = inputin.get_item<TF>("radiation", "lon", "");
    xka = inputin.get_item<TF>("radiation", "xka", "");
    fr0 = inputin.get_item<TF>("radiation", "fr0", "");
    fr1 = inputin.get_item<TF>("radiation", "fr1", "");
    div = inputin.get_item<TF>("radiation", "div", "");
}

template<typename TF>
Radiation_gcss<TF>::~Radiation_gcss()
{

}

template<typename TF>
void Radiation_gcss<TF>::init()
{
	auto& gd = grid.get_grid_data();
}

template<typename TF>
void Radiation_gcss<TF>::create(Thermo<TF>& thermo,Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross, Dump<TF>& dump)
{
    // Set up output classes
	create_stats(stats);
	create_column(column);
	create_dump(dump);
	create_cross(cross);
}

#ifndef USECUDA
template<typename TF>
void Radiation_gcss<TF>::exec(Thermo<TF>& thermo, double time, Timeloop<TF>& timeloop)
{
	auto& gd = grid.get_grid_data();
	auto lwp = fields.get_tmp();
	auto flx = fields.get_tmp();
	auto swn = fields.get_tmp();
	auto ql  = fields.get_tmp();
	thermo.get_thermo_field(*ql,"ql",false,false);
	struct tm current_datetime;
	current_datetime = timeloop.get_phytime();
	TF mu = calc_zenith(current_datetime, lat, lon);

	exec_gcss_rad<TF>(
		fields.st.at("thl")->fld.data(), ql->fld.data(), fields.sp.at("qt")->fld.data(),
		lwp->fld.data(), flx->fld.data(), swn->fld.data(), fields.rhoref.data(),
        mu, mu_min, fr0, fr1, xka, div,
		gd.z.data(), gd.dzhi.data(),
		gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
		gd.icells, gd.ijcells, gd.ncells);
	fields.release_tmp(lwp);
	fields.release_tmp(flx);
	fields.release_tmp(swn);
	fields.release_tmp(ql);
}

#endif
template<typename TF>
bool Radiation_gcss<TF>::check_field_exists(const std::string name)
{
    if (name == "rflx" || name == "sflx")
        return true;
    else
        return false;
}

template<typename TF>
void Radiation_gcss<TF>::get_radiation_field(Field3d<TF>& fld, std::string name, Thermo<TF>& thermo, Timeloop<TF>& timeloop)
{
    if (name == "lflx")
    {
        auto& gd = grid.get_grid_data();
        auto lwp = fields.get_tmp();
        auto ql  = fields.get_tmp();
        thermo.get_thermo_field(*ql,"ql",false,false);
        calc_gcss_rad_LW(ql->fld.data(), fields.ap.at("qt")->fld.data(),
        lwp->fld.data(), fld.fld.data(), fields.rhoref.data(), fr0, fr1, xka, div,
        gd.z.data(), gd.dzi.data(),
        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
        fields.release_tmp(lwp);
        fields.release_tmp(ql);
    }

    else if (name == "sflx")
    {
        struct tm current_datetime;
        current_datetime = timeloop.get_phytime();
        TF mu = calc_zenith(current_datetime, lat, lon);
        auto& gd = grid.get_grid_data();
        if (mu > mu_min) //if daytime, call SW (make a function for day/night determination)
        {
            auto ql  = fields.get_tmp();
            thermo.get_thermo_field(*ql,"ql",false,false);
            calc_gcss_rad_SW(fld.fld.data(), ql->fld.data(), fields.ap.at("qt")->fld.data(),
                fields.rhoref.data(), gd.z.data(), gd.dzi.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells, gd.ncells, mu);
            fields.release_tmp(ql);
        }

        else //night time, set SW to 0
        {
            for (int n=0; n<gd.ncells; ++n)
            {
                fld.fld[n] = TF(0.);
            }

        }

    }
}

template<typename TF>
void Radiation_gcss<TF>::create_stats(Stats<TF>& stats)
{
    if (stats.get_switch())
    {
        stats.add_prof("sflx", "Total shortwave radiative flux", "W m-2", "z");
        stats.add_prof("lflx", "Total longwave radiative flux", "W m-2", "z");
    }
}

template<typename TF>
void Radiation_gcss<TF>::create_column(Column<TF>& column)
{
    // add the profiles to the columns
    if (column.get_switch())
    {
        column.add_prof("sflx", "Total shortwave radiative flux", "W m-2", "z");
        column.add_prof("lflx", "Total longwave radiative flux", "W m-2", "z");
    }
}

template<typename TF>
void Radiation_gcss<TF>::create_cross(Cross<TF>& cross)
{
    if (cross.get_switch())
    {
        // Vectors with allowed cross variables for radiative flux
        std::vector<std::string> allowed_crossvars_rflx = {"sflx","lflx"};

        crosslist  = cross.get_enabled_variables(allowed_crossvars_rflx);
    }
}

template<typename TF>
void Radiation_gcss<TF>::create_dump(Dump<TF>& dump)
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
void Radiation_gcss<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

    const TF no_offset = 0.;
    const TF no_threshold = 0.;
    // calculate the mean
    std::vector<std::string> operators = {"mean"}; //add 2nd moment, if needed

    auto tmp = fields.get_tmp();

    get_radiation_field(*tmp,"lflx",thermo, timeloop);
    stats.calc_stats("lflx", *tmp, no_offset, no_threshold, operators);

    get_radiation_field(*tmp,"sflx",thermo, timeloop);
    stats.calc_stats("sflx", *tmp, no_offset, no_threshold, operators);

    fields.release_tmp(tmp);
}

#ifndef USECUDA
template<typename TF>
void Radiation_gcss<TF>::exec_column(Column<TF>& column, Thermo<TF>& thermo, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();
    const TF no_offset = 0.;

    auto flx = fields.get_tmp();
    get_radiation_field(*flx,"lflx",thermo,timeloop);
    column.calc_column("lflx", flx->fld.data(), no_offset);

    get_radiation_field(*flx,"sflx",thermo,timeloop);
    column.calc_column("sflx", flx->fld.data(), no_offset);

    fields.release_tmp(flx);
}
#endif

template<typename TF>
void Radiation_gcss<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime, Thermo<TF>& thermo, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();
    auto tmp = fields.get_tmp();

    for (auto& it : crosslist)
    {
        get_radiation_field(*tmp, it, thermo, timeloop);
        cross.cross_simple(tmp->fld.data(), it, iotime);
    }
    fields.release_tmp(tmp);
}

template<typename TF>
void Radiation_gcss<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime, Thermo<TF>& thermo, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();
    auto output = fields.get_tmp();

    for (auto& it : dumplist)
    {
        get_radiation_field(*output, it, thermo, timeloop);
        dump.save_dump(output->fld.data(), it, iotime);
    }
    fields.release_tmp(output);
}

template class Radiation_gcss<double>;
template class Radiation_gcss<float>;