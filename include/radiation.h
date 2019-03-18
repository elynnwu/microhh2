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

#ifndef RADIATION_H
#define RADIATION_H

#include <memory>
#include <vector>
#include "field3d_operators.h"

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
template<typename> class Timeloop;

template<typename TF>
class Radiation
{
    public:
        Radiation(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Radiation();
		static std::shared_ptr<Radiation> factory(Master&, Grid<TF>&, Fields<TF>&, Input&);
		std::string get_switch();

		//functions that the derived class has to implement
        virtual void init() = 0;
        virtual void create(Thermo<TF>&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&) = 0;
        virtual void exec(Thermo<TF>&, double, Timeloop<TF>&) = 0;

        virtual bool check_field_exists(std::string name) = 0;
        virtual void get_radiation_field(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&) = 0;

        virtual void exec_stats(Stats<TF>&, Thermo<TF>&, Timeloop<TF>&) = 0;
        virtual void exec_cross(Cross<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&) = 0;
        virtual void exec_dump(Dump<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&) = 0;
        virtual void exec_column(Column<TF>&, Thermo<TF>&, Timeloop<TF>&) = 0;
    protected:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_operators<TF> field3d_operators;
		std::string swradiation;
};
#endif
