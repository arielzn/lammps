/***************************************************************************
                                base_three.h
                             -------------------
                            W. Michael Brown (ORNL)

  Base class for 3-body potentials

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Tue April 2, 2013
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_BASE_THREE_H
#define LAL_BASE_THREE_H

#include "lal_device.h"
#include "lal_balance.h"
#include "mpi.h"

#if defined(USE_OPENCL)
#include "geryon/ocl_texture.h"
#elif defined(USE_CUDART)
#include "geryon/nvc_texture.h"
#else
#include "geryon/nvd_texture.h"
#endif

//#define THREE_CONCURRENT

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class BaseThree {
 public:
  BaseThree();
  virtual ~BaseThree();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device
    * \param k_two name for the kernel for 2-body force calculation
    * \param k_three name for the kernel for 3-body force calculation
    * 
    * Returns:
    * -  0 if successfull
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card
    * - -10 if invalid thread_per_atom setting **/
  int init_three(const int nlocal, const int nall, const int max_nbors,
                 const int maxspecial, const double cell_size, 
                 const double gpu_split, FILE *screen, 
                 const void *pair_program, const char *k_two,
                 const char *k_three_center, const char *k_three_end);

  /// Estimate the overhead for GPU context changes and CPU driver
  void estimate_gpu_overhead();

  /// Check if there is enough storage for atom arrays and realloc if not
  /** \param success set to false if insufficient memory **/
  inline void resize_atom(const int inum, const int nall, bool &success) {
    if (atom->resize(nall, success))
      pos_tex.bind_float(atom->x,4);
    ans->resize(inum,success);
    #ifdef THREE_CONCURRENT
    ans2->resize(inum,success);
    #endif
  }

  /// Check if there is enough storage for neighbors and realloc if not
  /** \param nlocal number of particles whose nbors must be stored on device
    * \param host_inum number of particles whose nbors need to copied to host
    * \param current maximum number of neighbors
    * \note olist_size=total number of local particles **/
  inline void resize_local(const int inum, const int max_nbors, bool &success) {
    nbor->resize(inum,max_nbors,success);
  }

  /// Check if there is enough storage for neighbors and realloc if not
  /** \param nlocal number of particles whose nbors must be stored on device
    * \param host_inum number of particles whose nbors need to copied to host
    * \param current maximum number of neighbors
    * \note host_inum is 0 if the host is performing neighboring
    * \note nlocal+host_inum=total number local particles
    * \note olist_size=0 **/
  inline void resize_local(const int inum, const int host_inum, 
                           const int max_nbors, bool &success) {
    nbor->resize(inum,host_inum,max_nbors,success);
  }

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear_atomic();

  /// Returns memory usage on device per atom
  int bytes_per_atom_atomic(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage_atomic() const;

  /// Accumulate timers
  inline void acc_timers() {
    if (device->time_device()) {
      nbor->acc_timers();
      time_pair.add_to_total();
      atom->acc_timers();
      ans->acc_timers();
      #ifdef THREE_CONCURRENT
      ans2->acc_timers();
      #endif
    }
  }

  /// Zero timers
  inline void zero_timers() {
    time_pair.zero();
    atom->zero_timers();
    ans->zero_timers();
    #ifdef THREE_CONCURRENT
    ans2->zero_timers();
    #endif
  }

  /// Copy neighbor list from host
  int * reset_nbors(const int nall, const int inum, const int nlist, int *ilist,
                    int *numj, int **firstneigh, bool &success);

  /// Build neighbor list on device
  int build_nbor_list(const int inum, const int host_inum,
                       const int nall, double **host_x, int *host_type,
                       double *sublo, double *subhi, tagint *tag, int **nspecial, 
                       tagint **special, bool &success);

  /// Pair loop with host neighboring
  void compute(const int f_ago, const int inum_full, const int nall, 
               const int nlist, double **host_x, int *host_type,
               int *ilist, int *numj, int **firstneigh, const bool eflag,
               const bool vflag, const bool eatom, const bool vatom,
               int &host_start, const double cpu_time, bool &success);

  /// Pair loop with device neighboring
  int * compute(const int ago, const int inum_full, const int nall, 
                double **host_x, int *host_type, double *sublo,
                double *subhi, tagint *tag, int **nspecial,
                tagint **special, const bool eflag, const bool vflag, 
                const bool eatom, const bool vatom, int &host_start, 
                const double cpu_time, bool &success);

  /// Pair loop with device neighboring
  int ** compute(const int ago, const int inum_full,
                 const int nall, double **host_x, int *host_type, double *sublo,
                 double *subhi, tagint *tag, int **nspecial,
                 tagint **special, const bool eflag, const bool vflag, 
                 const bool eatom, const bool vatom, int &host_start, 
                 int **ilist, int **numj, const double cpu_time, bool &success);

  // -------------------------- DEVICE DATA ------------------------- 

  /// Device Properties and Atom and Neighbor storage
  Device<numtyp,acctyp> *device;

  /// Geryon device
  UCL_Device *ucl_device;

  /// Device Timers
  UCL_Timer time_pair;

  /// Host device load balancer
  Balance<numtyp,acctyp> hd_balancer;

  /// LAMMPS pointer for screen output
  FILE *screen;

  // --------------------------- ATOM DATA --------------------------

  /// Atom Data
  Atom<numtyp,acctyp> *atom;

  // ------------------------ FORCE/ENERGY DATA -----------------------

  Answer<numtyp,acctyp> *ans;
  #ifdef THREE_CONCURRENT
  Answer<numtyp,acctyp> *ans2;
  #endif  

  // --------------------------- NBOR DATA ----------------------------

  /// Neighbor data
  Neighbor *nbor;

  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Program *pair_program;
  UCL_Kernel k_pair, k_three_center, k_three_end, k_three_end_vatom;
  inline int block_pair() { return _block_pair; }
  inline int block_size() { return _block_size; }

  // --------------------------- TEXTURES -----------------------------
  UCL_Texture pos_tex;

 protected:
  bool _compiled;
  int _block_pair, _block_size, _threads_per_atom, _end_command_queue;
  double _max_bytes, _max_an_bytes;
  double _gpu_overhead, _driver_overhead;
  UCL_D_Vec<int> *_nbor_data;

  void compile_kernels(UCL_Device &dev, const void *pair_string, 
                       const char *k_two, const char *k_three_center,
                       const char *k_three_end);

  virtual void loop(const bool _eflag, const bool _vflag, 
                    const int evatom) = 0;
};

}

#endif

