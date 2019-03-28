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

 #include "master.h"
 #include "grid.h"
 #include "fields.h"
 #include "boundary_non_cyclic.h"
 #include "netcdf_interface.h"
 #include "timedep.h"
 #include "timeloop.h"
 #include "field3d_operators.h"

 namespace // code from boundary_outflow.cxx in the outflow branch
 {

     template<typename TF>
     void set_inflow_BC_Dirichlet(
         TF* const restrict data, const TF* const restrict value,
         const int igc, const int jj, const int kk,
         const int jcells, const int kcells)
     {
         for (int k=0; k<kcells; ++k)
             for (int j=0; j<jcells; ++j)
                 #pragma ivdep
                 for (int i=0; i<igc; ++i)
                 {
                     const int ijk0 = i + j*jj + k*kk;
                     data[ijk0] = value;
                 }
     }

     template<typename TF>
     void compute_outflow_2nd(
             TF* const restrict a,
             const int iend,
             const int icells, const int jcells, const int kcells,
             const int ijcells)
     {
         const int ii = 1;

         // Set the ghost cells using extrapolation.
         for (int k=0; k<kcells; ++k)
             for (int j=0; j<jcells; ++j)
             {
                 const int ijk = (iend-1) + j*icells + k*ijcells;
                 a[ijk+ii] = a[ijk];
             }
     }

     template<typename TF>
     void compute_inflow_2nd(
             TF* const restrict a, const TF value,
             const int istart,
             const int icells, const int jcells, const int kcells,
             const int ijcells)
     {
         const int ii = 1;

         // Set the ghost cells using extrapolation.
         for (int k=0; k<kcells; ++k)
             for (int j=0; j<jcells; ++j)
             {
                 const int ijk = istart + j*icells + k*ijcells;
                 a[ijk-ii] = value - a[ijk];
             }
     }

     template<typename TF>
     void compute_outflow_4th(
             TF* const restrict a,
             const int iend,
             const int icells, const int jcells, const int kcells,
             const int ijcells)
     {
         const int ii1 = 1;
         const int ii2 = 2;
         const int ii3 = 3;

         // Set the ghost cells using extrapolation.
         for (int k=0; k<kcells; ++k)
             for (int j=0; j<jcells; ++j)
             {
                 const int ijk = (iend-1) + j*icells + k*ijcells;
                 a[ijk+ii1] = TF(2.)*a[ijk] - TF( 3./2.)*a[ijk-ii1] + TF(1./2.)*a[ijk-ii2];
                 a[ijk+ii2] = TF(3.)*a[ijk] - TF( 7./2.)*a[ijk-ii1] + TF(3./2.)*a[ijk-ii2];
                 a[ijk+ii3] = TF(5.)*a[ijk] - TF(15./2.)*a[ijk-ii1] + TF(7./2.)*a[ijk-ii2];
             }
     }

     template<typename TF>
     void compute_inflow_4th(
             TF* const restrict a, const TF value,
             const int istart,
             const int icells, const int jcells, const int kcells,
             const int ijcells)
     {
         const int ii1 = 1;
         const int ii2 = 2;
         const int ii3 = 3;

         // Set the ghost cells using extrapolation.
         for (int k=0; k<kcells; ++k)
             for (int j=0; j<jcells; ++j)
             {
                 const int ijk = istart + j*icells + k*ijcells;
                 a[ijk-ii1] = value + TF( 9./8.)*a[ijk] - TF( 14./8.)*a[ijk+ii1] + TF( 5./8.)*a[ijk+ii2];
                 a[ijk-ii2] = value + TF(33./8.)*a[ijk] - TF( 54./8.)*a[ijk+ii1] + TF(21./8.)*a[ijk+ii2];
                 a[ijk-ii3] = value + TF(65./8.)*a[ijk] - TF(110./8.)*a[ijk+ii1] + TF(45./8.)*a[ijk+ii2];
             }
     }
 }

template<typename TF>
Boundary_non_cyclic<TF>::Boundary_non_cyclic(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(masterin, gridin, fieldsin)
{
    std::string swinflowBC_in = inputin.get_item<std::string>("boundary", "inflowBC"    , "", "0");
    if (swinflowBC_in == "0")
        swinflowBC = Impose_inflow_BC_type::disabled;
    else if (swinflowBC_in == "1")
    {
        swinflowBC = Impose_inflow_BC_type::enabled;
        inflowBClist = inputin.get_list<std::string>("boundary", "inflowBClist", "", std::vector<std::string>());

        if (inputin.get_item<bool>("boundary", "swtimedep_inflowBC", "", false))
        {
            std::vector<std::string> tdepvars = inputin.get_list<std::string>("boundary", "timedep_inflowBC", "", std::vector<std::string>());
            for (auto& it : tdepvars)
                tdep_inflowBC.emplace(it, new Timedep<TF>(master, grid, it+"_inflowBC", true));
        }
    }
    else
    {
        throw std::runtime_error("Invalid option for \"swinflowbc\"");
    }

}

template<typename TF>
Boundary_non_cyclic<TF>::~Boundary_non_cyclic()
{
}

template<typename TF>
void Boundary_non_cyclic<TF>::init()
{
    auto& gd = grid.get_grid_data();
    if (swinflowBC == Impose_inflow_BC_type::enabled)
    {
        for (auto& it : inflowBClist)
            inflowBCprofs[it] = std::vector<TF>(gd.kcells);
    }
}

template <typename TF>
void Boundary_non_cyclic<TF>::create(Input& inputin, Netcdf_handle& input_nc)
{
    auto& gd = grid.get_grid_data();
    Netcdf_group group_nc = input_nc.get_group("init");
    if (swinflowBC == Impose_inflow_BC_type::enabled)
    {
        // Check whether the fields in the list exist in the prognostic fields.
        for (std::string& it : inflowBClist)
            if (!fields.ap.count(it))
            {
                std::string msg = "field " + it + " in [boundary][inflowBClist] is illegal";
                throw std::runtime_error(msg);
            }

        for (std::string& it : inflowBClist)
        {
            group_nc.get_variable(inflowBCprofs[it], it+"_inflowBC", {0}, {gd.ktot});
            std::rotate(inflowBCprofs[it].rbegin(), inflowBCprofs[it].rbegin() + gd.kstart, inflowBCprofs[it].rend());
        }

        // Process the time dependent data
        const TF offset = 0;
        for (auto& it : tdep_inflowBC)
            it.second->create_timedep_prof(input_nc, offset);
    }
}

template<typename TF>
void Boundary_non_cyclic<TF>::exec()
{
    auto& gd = grid.get_grid_data();

    const int jj = gd.icells;
    const int kk = gd.icells*gd.jcells;

    for (auto& it : inflowBClist)
    set_inflow_BC_Dirichlet<TF>(
        fields.sp.at(it)->fld.data(), inflowBCprofs.at(it).data(),
        gd.igc, jj, kk,
        gd.jcells, gd.kcells);
}

template class Boundary_non_cyclic<double>;
template class Boundary_non_cyclic<float>;
