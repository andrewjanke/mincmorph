/* kernel_ops.h */

#ifndef KERNEL_OPS
#define KERNEL_OPS

#include <float.h>
#include <volume_io.h>
#include "kernel_io.h"

/* kernel functions */
Volume  *pad_volume(Kernel * K, Volume * vol, double pad_value);
Volume  *convolve_kernel(Kernel * K, Volume * vol);
Volume  *erosion_kernel(Kernel * K, Volume * vol);
Volume  *dilation_kernel(Kernel * K, Volume * vol);
Volume  *distance_kernel(Kernel * K, Volume * vol);
Volume  *group_kernel(Kernel * K, Volume * vol, char *group_file, 
                      double gfloor, double gceil, int max_groups);
Volume  *binarise(Volume * vol, double floor, double ceil);
Volume  *pad_volume(Kernel * K, Volume * vol, double pad_value);

#endif
