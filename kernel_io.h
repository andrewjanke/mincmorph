/* kernel_io.h */

#ifndef KERNEL_IO
#define KERNEL_IO

#define KERNEL_DIMS 5

#include <volume_io.h>

/* Structure for Kernel information */
typedef struct {
   int      nelems;
   int      pre_pad[KERNEL_DIMS];
   int      post_pad[KERNEL_DIMS];
   Real   **K;
   } Kernel;

/* returns a new B_Matrix struct (pointer) */
Kernel  *new_kernel(int nelems);

/* reads in a B_Matrix from a file (pointer) */
Status   input_kernel(const char *kernel_file, Kernel * kernel);

/* pretty print a kernel */
int      print_kernel(Kernel * kernel);

/* calculate start and step offsets for this kernel */
int      setup_pad_values(Kernel * kernel);

/* return the default kernel */
void     setup_def_kernel(Kernel * K);

#endif
