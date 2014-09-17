/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2014 the CP2K developers group                      *
 *****************************************************************************/

#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include "acc_cuda_error.h"
#include "../include/acc.h"

/****************************************************************************/
extern "C" int acc_get_ndevices(int *n_devices){
  cudaError_t cErr;

  cErr = cudaGetDeviceCount (n_devices);
  if (cuda_error_check (cErr))
    return -1;
  return 0;
}


/****************************************************************************/
extern "C" int acc_set_active_device(int device_id){
  cudaError_t cErr;
  int myDevice;

  cErr = cudaSetDevice (device_id);
  if (cuda_error_check (cErr))
    return -1;

  cErr = cudaGetDevice (&myDevice);
  if (cuda_error_check (cErr))
    return -1;

  if (myDevice != device_id)
    return -1;

  return 0;
}

//EOF
