/* kernel_ops.h */

#ifndef KERNEL_OPS
#define KERNEL_OPS

#include <float.h>
#include <volume_io.h>
#include "kernel_io.h"

/* kernel functions */
Volume  *binarise(Volume * vol, nc_type dtype, double floor, double ceil, double foreground, double background);
Volume  *clamp(Volume * vol, double floor, double ceil, double background);
Volume  *pad(Kernel * K, Volume * vol, double background);
Volume  *dilation_kernel(Kernel * K, Volume * vol);
Volume  *erosion_kernel(Kernel * K, Volume * vol);
Volume  *convolve_kernel(Kernel * K, Volume * vol);
Volume  *distance_kernel(Kernel * K, Volume * vol, double background);
Volume  *group_kernel(Kernel * K, Volume * vol, nc_type dtype, char *group_file);

#endif
