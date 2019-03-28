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

#ifndef BOUNDARY_NON_CYCLIC
#define BOUNDARY_NON_CYCLIC

#include "timedep.h"

#ifdef USEMPI
#include <mpi.h>
#endif

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Field3d_operators;
template<typename> class Timedep;

enum class Edge {East_west_edge, North_south_edge, Both_edges};
enum class Impose_inflow_BC_type {disabled, enabled};

template<typename TF>
class Boundary_non_cyclic
{
    public:
        Boundary_non_cyclic(Master&, Grid<TF>&, Fields<TF>&, Input&); // Constuctor of the boundary class.
        ~Boundary_non_cyclic();                  // Destructor of the boundary class.

        void init();           ///< Initialize the arrays that contain the profiles.
        void create(Input&, Netcdf_handle&);   ///< Read the profiles of the forces from the input.
        void exec();     ///< Add the tendencies belonging to the large-scale processes.

        void update_time_dependent(Timeloop<TF>&); ///< Update the time dependent parameters.

        std::vector<std::string> inflowBClist;        ///< List of variables that have large-scale forcings.
        std::map<std::string, std::vector<TF>> inflowBCprofs; ///< Map of profiles with inflow BC stored by its name.

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_operators<TF> field3d_operators;

        Impose_inflow_BC_type swinflowBC;

        std::map<std::string, Timedep<TF>*> tdep_inflowBC;


};
#endif
