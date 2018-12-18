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

#ifndef RADIATION_H
#define RADIATION_H

#include <vector>
#include "field3d_operators.h"

enum class Radiation_type {Enabled, Disabled, Gcss};

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Thermo;
template<typename> class Stats;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;

template<typename TF>
class Radiation
{
    public:
        Radiation(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Radiation();
        void init();
        void create(Thermo<TF>&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(Thermo<TF>&, double, Timeloop<TF>&);

        bool check_field_exists(std::string name);
        void get_radiation_field(Field3d<TF>&, std::string, Thermo<TF>&);

        void exec_stats(Stats<TF>&, Thermo<TF>&);
        void exec_cross(Cross<TF>&, unsigned long, Thermo<TF>&);
        void exec_dump(Dump<TF>&, unsigned long, Thermo<TF>&);
        void exec_column(Column<TF>&, Thermo<TF>&);
    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_operators<TF> field3d_operators;

        Radiation_type swradiation;

        // cross sections
        std::vector<std::string> crosslist;        ///< List with all crosses from ini file
        bool swcross_rflx;
        std::vector<std::string> dumplist;         ///< List with all 3d dumps from the ini file.

        void create_stats(Stats<TF>&);   ///< Initialization of the statistics.
        void create_column(Column<TF>&); ///< Initialization of the single column output.
        void create_dump(Dump<TF>&);     ///< Initialization of the single column output.
        void create_cross(Cross<TF>&);   ///< Initialization of the single column output.
        std::vector<std::string> available_masks;   // Vector with the masks that fields can provide


        int ncol;
        int nlay;
        int nbndlw;

        int icld;
        int idrv;

        std::vector<double> play; // (ncol, nlay)
        std::vector<double> plev; // (ncol, nlay+1)
        std::vector<double> tlay; // (ncol, nlay)
        std::vector<double> tlev; // (ncol, nlay+1)

        std::vector<double> tsfc; // (ncol)

        std::vector<double> h2ovmr; // (ncol, nlay)
        std::vector<double> o3vmr;  // (ncol, nlay)
        std::vector<double> co2vmr; // (ncol, nlay)
        std::vector<double> ch4vmr; // (ncol, nlay)
        std::vector<double> n2ovmr; // (ncol, nlay)
        std::vector<double> o2vmr;  // (ncol, nlay)

        std::vector<double> cfc11vmr; // (ncol, nlay)
        std::vector<double> cfc12vmr; // (ncol, nlay)
        std::vector<double> cfc22vmr; // (ncol, nlay)
        std::vector<double> ccl4vmr;  // (ncol, nlay)
        std::vector<double> emis;     // (ncol, nbndlw)

        int inflglw;
        int iceflglw;
        int liqflglw;

        std::vector<double> cldfr;  // (ncol, nlay)
        std::vector<double> cicewp; // (ncol, nlay)
        std::vector<double> cliqwp; // (ncol, nlay)
        std::vector<double> reice;  // (ncol, nlay)
        std::vector<double> reliq;  // (ncol, nlay)
        std::vector<double> taucld; // (nbndlw, ncol, nlay)
        std::vector<double> tauaer; // (nbndlw, ncol, nlay)

        // OUTPUT
        std::vector<double> uflx;      // (ncol, nlay)
        std::vector<double> dflx;      // (ncol, nlay)
        std::vector<double> hr;        // (ncol, nlay)
        std::vector<double> uflxc;     // (ncol, nlay)
        std::vector<double> dflxc;     // (ncol, nlay)
        std::vector<double> hrc;       // (ncol, nlay)
        std::vector<double> duflx_dt;  // (ncol, nlay)
        std::vector<double> duflxc_dt; // (ncol, nlay)
};
#endif
