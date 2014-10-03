/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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

#ifndef THERMO
#define THERMO

// forward declarations to speed up build time
class Master;
class Input;
class Grid;
class Fields;
struct mask;

class Thermo
{
  public:
    Thermo(Model *, Input *);
    virtual ~Thermo();
    static Thermo* factory(Master *, Input *, Model *); ///< Factory function for thermo class generation.

    virtual void init();
    virtual void create(Input *);
    virtual int exec();
    virtual int execstats(mask *);

    virtual void execcross();

    virtual int getmask(Field3d *, Field3d *, mask *);

    // interfacing functions to get buoyancy properties from other classes
    virtual int checkthermofield(std::string name);
    virtual int getthermofield(Field3d *, Field3d *, std::string name);
    virtual int getbuoyancysurf(Field3d *);
    virtual int getbuoyancyfluxbot(Field3d *);
    virtual int getprogvars(std::vector<std::string> *);

    std::string getsw();

    // GPU functions and variables
    virtual int prepareDevice();
    virtual int clearDevice();

  protected:
    Grid   *grid;
    Fields *fields;
    Master *master;
    Model  *model;

    std::string swthermo;
};
#endif
