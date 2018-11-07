# Ubuntu 16.04
if(USEMPI) 
  set(ENV{CC}  mpicc ) # C compiler for parallel build
  set(ENV{CXX} mpicxx) # C++ compiler for parallel build
else()
  set(ENV{CC}  gcc) # C compiler for serial build
  set(ENV{CXX} g++) # C++ compiler for serial build
endif()

set(USER_CXX_FLAGS "-std=c++14")
set(USER_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -mtune=native -march=native")
set(USER_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wno-unknown-pragmas")

set(FFTW_INCLUDE_DIR   "/home/mzamoraz/local/include")
set(FFTW_LIB           "/home/mzamoraz/local/lib/libfftw3.so")
set(FFTWF_LIB          "/home/mzamoraz/local/lib/libfftw3f.a")
set(NETCDF_INCLUDE_DIR "/home/mzamoraz/local/include")
set(NETCDF_LIB_C       "/home/mzamoraz/local/lib/libnetcdf.so")
set(NETCDF_LIB_CPP     "/home/mzamoraz/local/lib/libnetcdf_c++4.so")
set(HDF5_LIB_1         "/home/mzamoraz/local/lib/libhdf5.so")
set(HDF5_LIB_2         "/home/mzamoraz/local/lib/libhdf5_hl.so")
set(SZIP_LIB           "")
set(LIBS ${FFTW_LIB} ${FFTWF_LIB} ${NETCDF_LIB_CPP} ${NETCDF_LIB_C} ${HDF5_LIB_2} ${HDF5_LIB_1} ${SZIP_LIB} m z curl)
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR} ${NETCDF_INCLUDE_DIR})

add_definitions(-DRESTRICTKEYWORD=__restrict__)
