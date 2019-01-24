/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
 * Copyright (c) 2019 Elynn Wu
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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"

#include "radiation.h"
#include "radiation_disabled.h"
#include "radiation_gcss.h"
#include "radiation_rrtmg.h"

template<typename TF>
Radiation<TF>::Radiation(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(master, grid, fields)
{
}

template<typename TF>
Radiation<TF>::~Radiation()
{
}

template<typename TF>
std::string Radiation<TF>::get_switch()
{
	return swradiation;
}

template<typename TF>
std::shared_ptr<Radiation<TF>> Radiation<TF>::factory(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin)
{
	std::string swradiation = inputin.get_item<std::string>("radiation", "swradiation", "", "0");
    if (swradiation == "0")
        return std::make_shared<Radiation_disabled<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swradiation == "1") //rrtmg - call fortran
        return std::make_shared<Radiation_rrtmg<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swradiation == "2") //gcss - for Sc clouds
        return std::make_shared<Radiation_gcss<TF>>(masterin, gridin, fieldsin, inputin);
    else
    {
        masterin.print_error("\"%s\" is an illegal value for swradiation\n", swradiation.c_str());
        throw std::runtime_error("Illegal options swradiation");
    }
}

template class Radiation<double>;
template class Radiation<float>;
