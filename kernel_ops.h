/* kernel_ops.h */

#ifndef KERNEL_OPS
#define KERNEL_OPS

#include <volume_io.h>
#include "kernel_io.h"

/* kernel functions */
Volume  *binarise(Volume * vol, double floor, double ceil, double fg, double bg);
Volume  *clamp(Volume * vol, double floor, double ceil, double bg);
Volume  *pad(Kernel * K, Volume * vol, double bg);
Volume  *erosion_kernel(Kernel * K, Volume * vol);
Volume  *dilation_kernel(Kernel * K, Volume * vol);
Volume  *median_dilation_kernel(Kernel * K, Volume * vol);
Volume  *convolve_kernel(Kernel * K, Volume * vol);
Volume  *distance_kernel(Kernel * K, Volume * vol, double bg);
Volume  *group_kernel(Kernel * K, Volume * vol, double bg);

#endif
