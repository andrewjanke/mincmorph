/* kernel_ops.c */

#include <limits.h>
#include "kernel_ops.h"

extern int verbose;

/* function prototypes */
void     split_kernel(Kernel * K, Kernel * k1, Kernel * k2);
int      compare_ints(const void *a, const void *b);
int      compare_groups(const void *a, const void *b);


/* structure for group information */
typedef struct {
   unsigned int orig_label;
   unsigned int count;
} Group_struct;

int compare_ints(const void *a, const void *b)
{
   return (*(int *)a - *(int *)b);
   }

int compare_groups(const void *a, const void *b)
{
   return ((Group_struct *) b)->count - ((Group_struct *) a)->count;
   }

void split_kernel(Kernel * K, Kernel * k1, Kernel * k2)
{
   int      c, k1c, k2c;

   /* fill the two sub kernels */
   k1c = k2c = 0;
   for(c = 0; c < K->nelems; c++){
      if((K->K[c][2] < 0) ||
          (K->K[c][1] < 0 && K->K[c][2] <= 0) ||
          (K->K[c][0] < 0 && K->K[c][1] <= 0 && K->K[c][2] <= 0)){

         k1->K[k1c] = K->K[c];
         k1c++;
         }
      else{
         k2->K[k2c] = K->K[c];
         k2c++;
         }
      }
   k1->nelems = k1c;
   k2->nelems = k2c;
   }

/* binarise a volume between a range */
Volume  *binarise(Volume * vol, nc_type dtype, double floor, double ceil,
                  double foreground, double background)
{
   int      x, y, z;
   int      sizes[MAX_VAR_DIMS];
   double   value;
   progress_struct progress;
   Volume  *tmp_vol;
   Real     min, max;

   if(verbose){
      fprintf(stdout, "Binarising, range: [%g:%g] fg/bg: [%g:%g]\n", floor, ceil,
              foreground, background);
      }

   get_volume_sizes(*vol, sizes);

   /* copy the volume setting the range to ensure 1:1 between voxel and real */
   tmp_vol = malloc(sizeof(*tmp_vol));
   *tmp_vol = copy_volume_definition(*vol, dtype, FALSE, 0.0, 0.0);
   get_volume_voxel_range(*tmp_vol, &min, &max);
   set_volume_real_range(*tmp_vol, min, max);

   initialize_progress_report(&progress, FALSE, sizes[2], "Binarise");
   for(z = sizes[0]; z--;){
      for(y = sizes[1]; y--;){
         for(x = sizes[2]; x--;){

            value = get_volume_real_value(*vol, z, y, x, 0, 0);
            if((value >= floor) && (value <= ceil)){
               set_volume_voxel_value(*tmp_vol, z, y, x, 0, 0, foreground);
               }
            else{
               set_volume_voxel_value(*tmp_vol, z, y, x, 0, 0, background);
               }

            }
         }
      update_progress_report(&progress, z + 1);
      }

   delete_volume(*vol);
   terminate_progress_report(&progress);

   return (tmp_vol);
   }

/* clamp a volume between a range */
Volume  *clamp(Volume * vol, double floor, double ceil, double background)
{
   int      x, y, z;
   int      sizes[MAX_VAR_DIMS];
   double   value;
   progress_struct progress;

   if(verbose){
      fprintf(stdout, "Clamping, range: [%g:%g] bg: %g\n", floor, ceil, background);
      }

   get_volume_sizes(*vol, sizes);

   initialize_progress_report(&progress, FALSE, sizes[2], "Binarise");
   for(z = sizes[0]; z--;){
      for(y = sizes[1]; y--;){
         for(x = sizes[2]; x--;){
            value = get_volume_real_value(*vol, z, y, x, 0, 0);
            if((value < floor) || (value > ceil)){
               set_volume_real_value(*vol, z, y, x, 0, 0, background);
               }
            }
         }
      update_progress_report(&progress, z + 1);
      }

   terminate_progress_report(&progress);
   return (vol);
   }

/* pad a volume using the background value */
Volume  *pad(Kernel * K, Volume * vol, double background)
{
   int      x, y, z;
   int      sizes[MAX_VAR_DIMS];

   get_volume_sizes(*vol, sizes);

   /* z */
   for(y = 0; y < sizes[1]; y++){
      for(x = 0; x < sizes[2]; x++){
         for(z = 0; z < -K->pre_pad[2]; z++){
            set_volume_real_value(*vol, z, y, x, 0, 0, background);
            }
         for(z = sizes[0] - K->post_pad[2]; z < sizes[0]; z++){
            set_volume_real_value(*vol, z, y, x, 0, 0, background);
            }
         }
      }

   /* y */
   for(z = 0; z < sizes[0]; z++){
      for(x = 0; x < sizes[2]; x++){
         for(y = 0; y < -K->pre_pad[1]; y++){
            set_volume_real_value(*vol, z, y, x, 0, 0, background);
            }
         for(y = sizes[1] - K->post_pad[1]; y < sizes[1]; y++){
            set_volume_real_value(*vol, z, y, x, 0, 0, background);
            }
         }
      }

   /* x */
   for(z = 0; z < sizes[0]; z++){
      for(y = 0; y < sizes[1]; y++){
         for(x = 0; x < -K->pre_pad[0]; x++){
            set_volume_real_value(*vol, z, y, x, 0, 0, background);
            }
         for(x = sizes[2] - K->post_pad[0]; x < sizes[2]; x++){
            set_volume_real_value(*vol, z, y, x, 0, 0, background);
            }
         }
      }

   return (vol);
   }

/* perform a dilation on a volume */
Volume  *dilation_kernel(Kernel * K, Volume * vol)
{
   int      x, y, z, c;
   double   value;
   int      sizes[MAX_VAR_DIMS];
   progress_struct progress;
   Volume   tmp_vol;

   if(verbose){
      fprintf(stdout, "Dilation kernel\n");
      }
   get_volume_sizes(*vol, sizes);
   initialize_progress_report(&progress, FALSE, sizes[2], "Dilation");

   /* copy the volume */
   tmp_vol = copy_volume(*vol);

   for(z = -K->pre_pad[2]; z < sizes[0] - K->post_pad[2]; z++){
      for(y = -K->pre_pad[1]; y < sizes[1] - K->post_pad[1]; y++){
         for(x = -K->pre_pad[0]; x < sizes[2] - K->post_pad[0]; x++){

            value = get_volume_real_value(tmp_vol, z, y, x, 0, 0);
            for(c = 0; c < K->nelems; c++){
               if(get_volume_real_value(*vol,
                                         z + K->K[c][2],
                                         y + K->K[c][1], x + K->K[c][0], 0 + K->K[c][3],
                                         0 + K->K[c][4]) < value){
                  set_volume_real_value(*vol, z + K->K[c][2], y + K->K[c][1],
                                        x + K->K[c][0], 0 + K->K[c][3], 0 + K->K[c][4],
                                        value * K->K[c][5]);
                  }
               }
            }
         }

      update_progress_report(&progress, z + 1);
      }

   delete_volume(tmp_vol);
   terminate_progress_report(&progress);
   return (vol);
   }

/* perform an erosion on a volume */
Volume  *erosion_kernel(Kernel * K, Volume * vol)
{
   int      x, y, z, c, count;
   double   value;
   int      sizes[MAX_VAR_DIMS];
   progress_struct progress;
   Volume   tmp_vol;

   if(verbose){
      fprintf(stdout, "Erosion kernel\n");
      }
   get_volume_sizes(*vol, sizes);
   initialize_progress_report(&progress, FALSE, sizes[2], "Erosion");

   /* copy the volume */
   tmp_vol = copy_volume(*vol);

   for(z = -K->pre_pad[2]; z < sizes[0] - K->post_pad[2]; z++){
      for(y = -K->pre_pad[1]; y < sizes[1] - K->post_pad[1]; y++){
         for(x = -K->pre_pad[0]; x < sizes[2] - K->post_pad[0]; x++){

            value = get_volume_real_value(tmp_vol, z, y, x, 0, 0);

            count = 0;
            for(c = 0; c < K->nelems; c++){
               if(get_volume_real_value(tmp_vol,
                                         z + K->K[c][2],
                                         y + K->K[c][1], x + K->K[c][0], 0 + K->K[c][3],
                                         0 + K->K[c][4]) >= value){
                  count++;
                  }
               }

            if(count == K->nelems){
               set_volume_real_value(*vol, z, y, x, 0, 0, value);
               }
            else{
               set_volume_real_value(*vol, z, y, x, 0, 0, 0.0);
               }
            }
         }
      update_progress_report(&progress, z + 1);
      }

   delete_volume(tmp_vol);
   terminate_progress_report(&progress);
   return (vol);
   }

/* convolve a volume with a input kernel */
Volume  *convolve_kernel(Kernel * K, Volume * vol)
{
   int      x, y, z, c;
   double   value;
   int      sizes[MAX_VAR_DIMS];
   progress_struct progress;
   Volume   tmp_vol;

   if(verbose){
      fprintf(stdout, "Convolve kernel\n");
      }
   get_volume_sizes(*vol, sizes);
   initialize_progress_report(&progress, FALSE, sizes[2], "Convolve");

   /* copy the volume */
   tmp_vol = copy_volume(*vol);

   for(z = -K->pre_pad[2]; z < sizes[0] - K->post_pad[2]; z++){
      for(y = -K->pre_pad[1]; y < sizes[1] - K->post_pad[1]; y++){
         for(x = -K->pre_pad[0]; x < sizes[2] - K->post_pad[0]; x++){

            value = 0;
            for(c = 0; c < K->nelems; c++){
               value += get_volume_real_value(tmp_vol,
                                              z + K->K[c][2],
                                              y + K->K[c][1],
                                              x + K->K[c][0], 0 + K->K[c][3],
                                              0 + K->K[c][4]) * K->K[c][5];
               }
            set_volume_real_value(*vol, z, y, x, 0, 0, value);
            }
         }

      update_progress_report(&progress, z + 1);
      }

   delete_volume(tmp_vol);
   terminate_progress_report(&progress);
   return (vol);
   }


/* should really only work on binary images    */
/* from the original 2 pass Borgefors alg      */
Volume  *distance_kernel(Kernel * K, Volume * vol, double background)
{
   int      x, y, z, c;
   double   value, min;
   int      sizes[MAX_VAR_DIMS];
   progress_struct progress;
   Kernel  *k1, *k2;

   /* split the Kernel */
   k1 = new_kernel(K->nelems);
   k2 = new_kernel(K->nelems);
   split_kernel(K, k1, k2);

   if(verbose){
      fprintf(stdout, "Distance kernel\n");
      fprintf(stdout, "forward direction kernel:\n");
      print_kernel(k1);
      fprintf(stdout, "\nreverse direction kernel:\n");
      print_kernel(k2);
      }

   get_volume_sizes(*vol, sizes);
   initialize_progress_report(&progress, FALSE, sizes[2] * 2, "Distance");

   /* forward raster direction */
   for(z = -K->pre_pad[2]; z < sizes[0] - K->post_pad[2]; z++){
      for(y = -K->pre_pad[1]; y < sizes[1] - K->post_pad[1]; y++){
         for(x = -K->pre_pad[0]; x < sizes[2] - K->post_pad[0]; x++){

            if(get_volume_real_value(*vol, z, y, x, 0, 0) > background){

               /* find the minimum */
               min = DBL_MAX;
               for(c = 0; c < k1->nelems; c++){
                  value = get_volume_real_value(*vol,
                                                z + k1->K[c][2],
                                                y + k1->K[c][1],
                                                x + k1->K[c][0],
                                                0 + k1->K[c][3], 0 + k1->K[c][4]) + 1;
                  if(value < min){
                     min = value;
                     }
                  }

               set_volume_real_value(*vol, z, y, x, 0, 0, min);
               }
            }
         }
      update_progress_report(&progress, z + 1);
      }

   /* reverse raster direction */
   for(z = sizes[0] - K->post_pad[0] - 1; z > -K->pre_pad[0]; z--){
      for(y = sizes[1] - K->post_pad[1] - 1; y > -K->pre_pad[1]; y--){
         for(x = sizes[2] - K->post_pad[2] - 1; x > -K->pre_pad[2]; x--){

            min = get_volume_real_value(*vol, z, y, x, 0, 0);
            if(min > background){

               /* find the minimum distance to 0 in the neighbouring vectors */
               for(c = 0; c < k2->nelems; c++){
                  value = get_volume_real_value(*vol,
                                                z + k2->K[c][2],
                                                y + k2->K[c][1],
                                                x + k2->K[c][0],
                                                0 + k2->K[c][3], 0 + k2->K[c][4]) + 1;
                  if(value < min){
                     min = value;
                     }
                  }

               set_volume_real_value(*vol, z, y, x, 0, 0, min);
               }
            }
         }
      update_progress_report(&progress, sizes[2] + z + 1);
      }

   free(k1);
   free(k2);
   terminate_progress_report(&progress);
   return (vol);
   }

/* do connected components labelling on a volume */
/* resulting groups are sorted WRT size          */
Volume  *group_kernel(Kernel * K, Volume * vol, nc_type dtype, char *group_file)
{
   int      x, y, z;
   int      sizes[MAX_VAR_DIMS];
   progress_struct progress;
   Volume  *tmp_vol;
   Real     min, max;
   Kernel  *k1, *k2;

   unsigned int *equiv;
   unsigned int *counts;
   unsigned int *trans;
   unsigned int neighbours[K->nelems];

   /* counters */
   unsigned int c;
   unsigned int value;
   unsigned int group_idx;             /* label for the next group     */
   unsigned int num_groups;

   unsigned int min_label;
   unsigned int curr_label;
   unsigned int prev_label;
   unsigned int num_matches;


   /* structure for group data */
   Group_struct *group_data;


   /* split the Kernel into forward and backwards kernels */
   k1 = new_kernel(K->nelems);
   k2 = new_kernel(K->nelems);
   split_kernel(K, k1, k2);

   setup_pad_values(k1);
   setup_pad_values(k2);

   if(verbose){
      fprintf(stdout, "Group kernel\n");
      fprintf(stdout, "forward direction kernel:\n");
      print_kernel(k1);
      fprintf(stdout, "\nreverse direction kernel:\n");
      print_kernel(k2);
      }

   get_volume_sizes(*vol, sizes);
   initialize_progress_report(&progress, FALSE, sizes[2], "Groups");

   /* copy the volume setting the range to ensure 1:1 between voxel and real */
   tmp_vol = malloc(sizeof(*tmp_vol));
   *tmp_vol = copy_volume_definition(*vol, dtype, FALSE, 0.0, 0.0);
   get_volume_voxel_range(*tmp_vol, &min, &max);
   set_volume_real_range(*tmp_vol, min, max);
   for(z = sizes[0]; z--;){
      for(y = sizes[1]; y--;){
         for(x = sizes[2]; x--;){
            set_volume_voxel_value(*tmp_vol, z, y, x, 0, 0, 0);
            }
         }
      }

   /* pass 1 - forward direction (we assume a symmetric kernel) */
   group_idx = 1;

   /* initialise the equiv and counts arrays */
   SET_ARRAY_SIZE(equiv, 0, group_idx, 500);
   equiv[0] = 0;

   SET_ARRAY_SIZE(counts, 0, group_idx, 500);
   counts[0] = 0;

   for(z = -k1->pre_pad[2]; z < sizes[0] - k1->post_pad[2]; z++){
      for(y = -k1->pre_pad[1]; y < sizes[1] - k1->post_pad[1]; y++){
         for(x = -k1->pre_pad[0]; x < sizes[2] - k1->post_pad[0]; x++){

            if(get_volume_real_value(*vol, z, y, x, 0, 0) != 0){

               /* get this voxels neighbours */
               num_matches = 0;
               min_label = INT_MAX;

               for(c = 0; c < k1->nelems; c++){
                  value = (unsigned int)get_volume_voxel_value(*tmp_vol,
                                                               z + k1->K[c][2],
                                                               y + k1->K[c][1],
                                                               x + k1->K[c][0],
                                                               0 + k1->K[c][3],
                                                               0 + k1->K[c][4]);
                  if(value != 0){
                     if(value < min_label){
                        min_label = value;
                        }
                     neighbours[num_matches] = value;
                     num_matches++;
                     }
                  }

               switch (num_matches){
               case 0:
                  /* no neighbours, make a new label and increment */
                  set_volume_voxel_value(*tmp_vol, z, y, x, 0, 0, (Real)group_idx);

                  SET_ARRAY_SIZE(equiv, group_idx, group_idx + 1, 500);
                  equiv[group_idx] = group_idx;

                  SET_ARRAY_SIZE(counts, group_idx, group_idx + 1, 500);
                  counts[group_idx] = 1;

                  group_idx++;
                  break;

               case 1:
                  /* only one neighbour, no equivalences needed */
                  set_volume_voxel_value(*tmp_vol, z, y, x, 0, 0, (Real)min_label);
                  counts[min_label]++;
                  break;

               default:
                  /* more than one neighbour */

                  /* first sort the neighbours array */
                  qsort(&neighbours[0], (size_t) num_matches, sizeof(unsigned int),
                        &compare_ints);

                  /* find the minimum possible label for this voxel,    */
                  /* this is done by descending through each neighbours */
                  /* equivalences untill an equivalence equal to itself */
                  /* is found                                           */
                  prev_label = -1;
                  for(c = 0; c < num_matches; c++){
                     curr_label = neighbours[c];

                     /* recurse this label if we haven't yet */
                     if(curr_label != prev_label){
                        while (equiv[curr_label] != equiv[equiv[curr_label]]){
                           curr_label = equiv[curr_label];
                           }

                        /* check against the current minimum value */
                        if(equiv[curr_label] < min_label){
                           min_label = equiv[curr_label];
                           }
                        }

                     prev_label = neighbours[c];
                     }

                  /* repeat, setting equivalences to the min_label */
                  prev_label = -1;
                  for(c = 0; c < num_matches; c++){
                     curr_label = neighbours[c];

                     if(curr_label != prev_label){
                        while (equiv[curr_label] != equiv[equiv[curr_label]]){
                           curr_label = equiv[curr_label];

                           equiv[curr_label] = min_label;
                           }

                        /* set the label itself */
                        if(equiv[neighbours[c]] != min_label){
                           equiv[neighbours[c]] = min_label;
                           }
                        }

                     prev_label = neighbours[c];
                     }

                  /* finally set the voxel in question to the minimum value */
                  set_volume_voxel_value(*tmp_vol, z, y, x, 0, 0, (Real)min_label);
                  counts[min_label]++;
                  break;
               }                       /* end case */

               }
            }
         }
      update_progress_report(&progress, z + 1);
      }


   /* reduce the equiv and counts array */
   num_groups = 0;
   for(c = 0; c < group_idx; c++){

      /* if this equivalence is not resolved yet */
      if(c != equiv[c]){

         /* find the min label value */
         min_label = equiv[c];
         while (min_label != equiv[min_label]){
            min_label = equiv[min_label];
            }

         /* update the label and its counters */
         equiv[c] = min_label;
         counts[min_label] += counts[c];
         counts[c] = 0;
         }
      else{
         num_groups++;
         }
      }

   /* set up the group structure */
   group_data = malloc(num_groups * sizeof(*group_data));

   num_groups = 0;
   for(c = 0; c < group_idx; c++){
      if(counts[c] > 0){
         group_data[num_groups].orig_label = equiv[c];
         group_data[num_groups].count = counts[c];
         num_groups++;
         }
      }

   /* sort the groups by the count size */
   if(verbose){
      fprintf(stdout, "Found %d groups, sorting\n", num_groups);
      }
   qsort(group_data, num_groups, sizeof(*group_data), &compare_groups);

   /* set up the transpose array */
   trans = (unsigned int *)malloc(sizeof(*trans) * num_groups);
   for(c = 0; c < num_groups; c++){
      trans[group_data[c].orig_label] = c + 1; /* +1 to bump past 0 */
      }

   /* pass 2 - resolve equivalences in the output data */
   for(z = 0; z < sizes[0]; z++){
      for(y = 0; y < sizes[1]; y++){
         for(x = 0; x < sizes[2]; x++){

            value = (unsigned int)get_volume_voxel_value(*tmp_vol, z, y, x, 0, 0);
            if(value != 0){
               value = trans[equiv[value]];
               if(value > max){
                  value = 0;
                  }
               }
            set_volume_voxel_value(*tmp_vol, z, y, x, 0, 0, (Real)value);
            }
         }
      }

   if(verbose){
      fprintf(stdout, "# Group kernel counts\n");
      fprintf(stdout, "#\n");
      fprintf(stdout, "#  group      position         count\n");
      // insert stuff in here
      fprintf(stdout, "#\n");
      fprintf(stdout, "#   Total voxels: %14d\n", (sizes[0] * sizes[1] * sizes[2]));
      fprintf(stdout, "#       # groups: %14d\n", group_idx);
      }

   delete_volume(*vol);
   free(k1);
   free(k2);
   terminate_progress_report(&progress);
   return (tmp_vol);
   }
