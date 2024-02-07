#include "output.hpp"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

/*
 * We currently support the HDF5 1.6 API (and when using 1.8 the
 * compatibility mode introduced by H5_USE_16_API).  Several machines
 * in SimFactory use HDF5 1.6, so we cannot drop support for it.  It
 * seems it is hard to support both the 1.6 and 1.8 API
 * simultaneously; for example H5Fopen takes a different number of
 * arguments in the two versions.
 */

#ifdef HAVE_CAPABILITY_HDF5
#define H5_USE_16_API
#include <hdf5.h>
#endif

#include <cstdio>
#include <string>

// TODO: Change this to a static inline function
// check return code of HDF5 call abort with an error message if there was an
// error. adapted from CarpetIOHDF5/src/CarpetIOHDF5.hh.
#define HDF5_ERROR(fn_call)                                                    \
  do {                                                                         \
    hid_t _error_code = fn_call;                                               \
                                                                               \
    if (_error_code < 0) {                                                     \
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,                      \
                 "HDF5 call '%s' returned error code %d", #fn_call,            \
                 (int)_error_code);                                            \
    }                                                                          \
  } while (0)

static FILE *OpenOutputFile(CCTK_ARGUMENTS, const std::string &name) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  bool first_time = cctk_iteration == 0;
  const char *mode = first_time ? "w" : "a";
  const char *my_out_dir = strcmp(out_dir, "") ? out_dir : io_out_dir;
  const int err = CCTK_CreateDirectory(0755, my_out_dir);
  if (err < 0)
    CCTK_VWarn(
        CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
        "Multipole output directory %s could not be created (error code %d)",
        my_out_dir, err);

  std::string output_name(std::string(my_out_dir) + "/" + name);

  FILE *fp = fopen(output_name.c_str(), mode);

  if (fp == 0) {
    CCTK_VWarn(
        1, __LINE__, __FILE__, CCTK_THORNSTRING, "%s",
        (std::string("Could not open output file ") + output_name).c_str());
  }

  return fp;
}

static void OutputArray(CCTK_ARGUMENTS, FILE *f, const real_vec &th,
                        const real_vec &ph, const real_vec &xs,
                        const real_vec &ys, const real_vec &zs,
                        const real_vec &data) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_REAL last_ph = ph[0];

  for (size_t i = 0; i < data.size(); i++) {
    if (ph[i] != last_ph) // Separate blocks for gnuplot
      fprintf(f, "\n");
    fprintf(f, "%f %f %f %f %f %f %.19g\n", cctk_time, th[i], ph[i], xs[i],
            ys[i], zs[i], data[i]);
    last_ph = ph[i];
  }
}

void MultipoleX::OutputArrayToFile(CCTK_ARGUMENTS, const std::string &name,
                                   const real_vec &th, const real_vec &ph,
                                   const real_vec &xs, const real_vec &ys,
                                   const real_vec &zs, const real_vec &data) {
  DECLARE_CCTK_ARGUMENTS;

  if (FILE *fp = OpenOutputFile(CCTK_PASS_CTOC, name)) {
    OutputArray(CCTK_PASS_CTOC, fp, th, ph, xs, ys, zs, data);
    fclose(fp);
  }
}

void MultipoleX::Output1D(CCTK_ARGUMENTS, const std::string &name,
                          const real_vec &th, const real_vec &ph,
                          mp_coord coord, const real_vec &data) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (FILE *f = OpenOutputFile(CCTK_PASS_CTOC, name)) {
    fprintf(f, "\"Time = %.19g\n", cctk_time);

    if (coord == mp_theta) {
      for (int i = 0; i <= ntheta; i++) {
        int idx = Multipole_Index(i, 0, ntheta);
        fprintf(f, "%f %.19g\n", th[idx], data[idx]);
      }
    } else if (coord == mp_phi) {
      for (int i = 0; i <= nphi; i++) {
        int idx = Multipole_Index(ntheta / 4, i, ntheta);
        fprintf(f, "%f %.19g\n", ph[idx], data[idx]);
      }
    }
    fprintf(f, "\n\n");
    fclose(f);
  }
}

void Multipole_OutputComplex(CCTK_ARGUMENTS, FILE *fp, CCTK_REAL redata,
                             CCTK_REAL imdata) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  fprintf(fp, "%f %.19g %.19g\n", cctk_time, redata, imdata);
}

void MultipoleX::OutputComplexToFile(CCTK_ARGUMENTS, const std::string &name,
                                     CCTK_REAL redata, CCTK_REAL imdata) {
  DECLARE_CCTK_ARGUMENTS;

  if (FILE *fp = OpenOutputFile(CCTK_PASS_CTOC, name)) {
    Multipole_OutputComplex(CCTK_PASS_CTOC, fp, redata, imdata);
    fclose(fp);
  }
}

#ifdef HAVE_CAPABILITY_HDF5

static bool file_exists(const std::string &name) {
  struct stat sts;
  return !(stat(name.c_str(), &sts) == -1 && errno == ENOENT);
}

static bool dataset_exists(hid_t file, const std::string &dataset_name) {
  // To test whether a dataset exists, the recommended way in API 1.6
  // is to use H5Gget_objinfo, but this prints an error to stderr if
  // the dataset does not exist.  We explicitly avoid this by wrapping
  // the call in H5E_BEGIN_TRY/H5E_END_TRY statements.  In 1.8,
  // H5Gget_objinfo is deprecated, and H5Lexists does the job.  See
  // http://www.mail-archive.com/hdf-forum@hdfgroup.org/msg00125.html

#if 1
  bool exists;
  H5E_BEGIN_TRY {
    exists = H5Gget_objinfo(file, dataset_name.c_str(), 1, NULL) >= 0;
  }
  H5E_END_TRY;
  return exists;
#else
  return H5Lexists(file, dataset_name.c_str(), H5P_DEFAULT);
#endif
}

void MultipoleX::OutputComplexToH5File(
    CCTK_ARGUMENTS, const vector<Multipole::variable_desc> &vars,
    const CCTK_REAL radii[], const Multipole::mode_array &modes) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const char *my_out_dir = strcmp(out_dir, "") ? out_dir : io_out_dir;
  if (CCTK_CreateDirectory(0755, my_out_dir) < 0)
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Multipole output directory %s could not be created",
               my_out_dir);

  static map<std::string, bool> checked; // Has the given file been checked
                                         // for truncation? map<*,bool>
                                         // defaults to false
  for (int v = 0; v < modes.get_nvars(); v++) {
    std::string basename = "mp_" + vars[v].name + ".h5";
    std::string output_name = my_out_dir + std::string("/") + basename;

    hid_t file;

    if (!file_exists(output_name) ||
        (!checked[output_name] && IO_TruncateOutputFiles(cctkGH))) {
      file = H5Fcreate(output_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                       H5P_DEFAULT);
    } else {
      file = H5Fopen(output_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

    checked[output_name] = true;

    for (int i = 0; i < modes.get_nradii(); i++) {
      const CCTK_REAL rad = radii[i];
      for (int l = 0; l <= modes.get_lmax(); l++) {
        for (int m = -l; m <= l; m++) {
          ostringstream datasetname;
          datasetname << "l" << l << "_m" << m << "_r"
                      << setiosflags(ios::fixed) << setprecision(2) << rad;

          hid_t dataset = -1;

          if (dataset_exists(file, datasetname.str())) {
            dataset = H5Dopen(file, datasetname.str().c_str());
          } else {
            hsize_t dims[2] = {0, 3};
            hsize_t maxdims[2] = {H5S_UNLIMITED, 3};
            hid_t dataspace = H5Screate_simple(2, dims, maxdims);

            hid_t cparms = -1;
            hsize_t chunk_dims[2] = {hsize_t(hdf5_chunk_size), 3};
            cparms = H5Pcreate(H5P_DATASET_CREATE);
            HDF5_ERROR(H5Pset_chunk(cparms, 2, chunk_dims));

            dataset = H5Dcreate(file, datasetname.str().c_str(),
                                H5T_NATIVE_DOUBLE, dataspace, cparms);
            H5Pclose(cparms);
          }

          hid_t filespace = H5Dget_space(dataset);
          hsize_t filedims[2];
          hsize_t maxdims[2];
          HDF5_ERROR(H5Sget_simple_extent_dims(filespace, filedims, maxdims));

          filedims[0] += 1;
          hsize_t size[2] = {filedims[0], filedims[1]};
          HDF5_ERROR(H5Dextend(dataset, size));
          HDF5_ERROR(H5Sclose(filespace));

          /* Select a hyperslab  */
          hsize_t offset[2] = {filedims[0] - 1, 0};
          hsize_t dims2[2] = {1, 3};
          filespace = H5Dget_space(dataset);
          HDF5_ERROR(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset,
                                         NULL, dims2, NULL));

          CCTK_REAL data[] = {cctk_time, modes(v, i, l, m, 0),
                              modes(v, i, l, m, 1)};

          hid_t memdataspace = H5Screate_simple(2, dims2, NULL);

          /* Write the data to the hyperslab  */
          HDF5_ERROR(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memdataspace,
                              filespace, H5P_DEFAULT, data));

          HDF5_ERROR(H5Dclose(dataset));
          HDF5_ERROR(H5Sclose(filespace));
          HDF5_ERROR(H5Sclose(memdataspace));
        }
      }
    }
    HDF5_ERROR(H5Fclose(file));
  }
}

#else

void MultipoleX::OutputComplexToH5File(
    CCTK_ARGUMENTS, const vector<Multipole::variable_desc> &vars,
    const CCTK_REAL radii[], const Multipole::mode_array &modes) {
  CCTK_WARN(0, "HDF5 output has been requested but Cactus has been compiled "
               "without HDF5 support");
}

#endif

static void output_modes(CCTK_ARGUMENTS,
                         const vector<Multipole::variable_desc> &vars,
                         const CCTK_REAL radii[],
                         const Multipole::mode_array &modes) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (output_ascii) {
    if (CCTK_MyProc(cctkGH) == 0) {
      for (int v = 0; v < modes.get_nvars(); v++) {
        for (int i = 0; i < modes.get_nradii(); i++) {
          const CCTK_REAL rad = radii[i];
          for (int l = 0; l <= modes.get_lmax(); l++) {
            for (int m = -l; m <= l; m++) {
              ostringstream name;
              name << "mp_" << vars[v].name << "_l" << l << "_m" << m << "_r"
                   << setiosflags(ios::fixed) << setprecision(2) << rad
                   << ".asc";
              Multipole_OutputComplexToFile(CCTK_PASS_CTOC, name.str(),
                                            modes(v, i, l, m, 0),
                                            modes(v, i, l, m, 1));
            }
          }
        }
      }
    }
  }
  if (output_hdf5) {
    if (CCTK_MyProc(cctkGH) == 0) {
      Multipole_OutputComplexToH5File(CCTK_PASS_CTOC, vars, radii, modes);
    }
  }
}

static void output_1D(CCTK_ARGUMENTS, const Multipole::variable_desc &v,
                      CCTK_REAL rad, const vector<CCTK_REAL> &th,
                      const vector<CCTK_REAL> &ph,
                      const vector<CCTK_REAL> &real,
                      const vector<CCTK_REAL> &imag) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (CCTK_MyProc(cctkGH) == 0 && output_ascii) {
    if (out_1d_every != 0 && cctk_iteration % out_1d_every == 0) {
      ostringstream real_base;
      real_base << "mp_" << string(CCTK_VarName(v.index)) << "_r"
                << setiosflags(ios::fixed) << setprecision(2) << rad;
      Multipole_Output1D(CCTK_PASS_CTOC, real_base.str() + string(".th.asc"),
                         th, ph, mp_theta, real);
      Multipole_Output1D(CCTK_PASS_CTOC, real_base.str() + string(".ph.asc"),
                         th, ph, mp_phi, real);

      if (v.imag_index != -1) {
        ostringstream imag_base;
        imag_base << "mp_" << string(CCTK_VarName(v.imag_index)) << "_r"
                  << setiosflags(ios::fixed) << setprecision(2) << rad;
        Multipole_Output1D(CCTK_PASS_CTOC, imag_base.str() + string(".th.asc"),
                           th, ph, mp_theta, imag);
        Multipole_Output1D(CCTK_PASS_CTOC, imag_base.str() + string(".ph.asc"),
                           th, ph, mp_phi, imag);
      }
    }
  }
}