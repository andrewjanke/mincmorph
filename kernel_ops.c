/* kernel_ops.c */

#include "kernel_ops.h"

extern int verbose;

/* function prototypes */
void     split_kernel(Kernel * K, Kernel * k1, Kernel * k2);

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

Volume  *pad_volume(Kernel * K, Volume * vol, double pad_value)
{
   int      x, y, z;
   int      sizes[MAX_VAR_DIMS];

   get_volume_sizes(*vol, sizes);

   /* z */
   for(y = 0; y < sizes[1]; y++){
      for(x = 0; x < sizes[2]; x++){
         for(z = 0; z < -K->pre_pad[2]; z++){
            set_volume_real_value(*vol, z, y, x, 0, 0, pad_value);
            }
         for(z = sizes[0] - K->post_pad[2]; z < sizes[0]; z++){
            set_volume_real_value(*vol, z, y, x, 0, 0, pad_value);
            }
         }
      }

   /* y */
   for(z = 0; z < sizes[0]; z++){
      for(x = 0; x < sizes[2]; x++){
         for(y = 0; y < -K->pre_pad[1]; y++){
            set_volume_real_value(*vol, z, y, x, 0, 0, pad_value);
            }
         for(y = sizes[1] - K->post_pad[1]; y < sizes[1]; y++){
            set_volume_real_value(*vol, z, y, x, 0, 0, pad_value);
            }
         }
      }

   /* x */
   for(z = 0; z < sizes[0]; z++){
      for(y = 0; y < sizes[1]; y++){
         for(x = 0; x < -K->pre_pad[0]; x++){
            set_volume_real_value(*vol, z, y, x, 0, 0, pad_value);
            }
         for(x = sizes[2] - K->post_pad[0]; x < sizes[2]; x++){
            set_volume_real_value(*vol, z, y, x, 0, 0, pad_value);
            }
         }
      }

   return (vol);
   }

/* morphological functions */
Volume  *binarise(Volume * vol, double floor, double ceil)
{
   int      x, y, z;
   int      sizes[MAX_VAR_DIMS];
   double   value;
   progress_struct progress;
   Volume   tmp_vol;

   if(verbose){
      fprintf(stdout, "Binarising, range: [%g:%g]\n", floor, ceil);
      }

   get_volume_sizes(*vol, sizes);

   /* copy the volume */
   tmp_vol = copy_volume(*vol);

   /* fiddle with the original volume */
   set_volume_type(*vol, NC_BYTE, FALSE, 0.0, 255);
   set_volume_real_range(*vol, 0, 255);

   initialize_progress_report(&progress, FALSE, sizes[2], "Binarise");
   for(z = 0; z < sizes[0]; z++){
      for(y = 0; y < sizes[1]; y++){
         for(x = 0; x < sizes[2]; x++){

            value = get_volume_real_value(tmp_vol, z, y, x, 0, 0);
            if((value >= floor) && (value <= ceil)){
               set_volume_real_value(*vol, z, y, x, 0, 0, 1);
               }
            else{
               set_volume_real_value(*vol, z, y, x, 0, 0, 0);
               }

            }
         }
      update_progress_report(&progress, z + 1);
      }

   delete_volume(tmp_vol);
   terminate_progress_report(&progress);
   return (vol);
   }

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
                  set_volume_real_value(*vol, z + K->K[c][2], y + K->K[c][1], x + K->K[c][0],
                                        0 + K->K[c][3], 0 + K->K[c][4], value * K->K[c][5]);
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

/* should really only work on binary images    */
/* from the original 2 pass Borgefors alg      */
Volume  *distance_kernel(Kernel * K, Volume * vol)
{
   int      x, y, z, c;
   double   value, min;
   int      sizes[MAX_VAR_DIMS];
   progress_struct progress;
   Kernel  *k1, *k2;

   if(verbose){
      fprintf(stdout, "Distance kernel\n");
      }
   get_volume_sizes(*vol, sizes);
   initialize_progress_report(&progress, FALSE, sizes[2] * 2, "Distance");

   k1 = new_kernel(K->nelems);
   k2 = new_kernel(K->nelems);

   /* split the Kernel */
   split_kernel(K, k1, k2);

   if(verbose){
      fprintf(stdout, "forward direction kernel:\n");
      print_kernel(k1);
      fprintf(stdout, "\nreverse direction kernel:\n");
      print_kernel(k2);
      }

   /* forward raster direction */
   for(z = -K->pre_pad[2]; z < sizes[0] - K->post_pad[2]; z++){
      for(y = -K->pre_pad[1]; y < sizes[1] - K->post_pad[1]; y++){
         for(x = -K->pre_pad[0]; x < sizes[2] - K->post_pad[0]; x++){

            if(get_volume_real_value(*vol, z, y, x, 0, 0) > 0){

               /* find the minimum */
               min = DBL_MAX;
               for(c = 0; c < k1->nelems; c++){
                  value = get_volume_real_value(*vol,
                                                z + k1->K[c][2],
                                                y + k1->K[c][1], 
                                                x + k1->K[c][0], 
                                                0 + k1->K[c][3],
                                                0 + k1->K[c][4]) + 1;
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
            if(min > 0){

               /* find the minimum distance to 0 in the neighbouring vectors */
               for(c = 0; c < k2->nelems; c++){
                  value = get_volume_real_value(*vol,
                                                z + k2->K[c][2],
                                                y + k2->K[c][1], 
                                                x + k2->K[c][0], 
                                                0 + k2->K[c][3],
                                                0 + k2->K[c][4]) + 1;
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

Volume  *group_kernel(Kernel * K, Volume * vol, char *group_file, double gfloor, double gceil,
                      int max_groups)
{
   int      x, y, z, c;
   double   value;
   int      sizes[MAX_VAR_DIMS];
   progress_struct progress;
   Volume   tmp_vol;
   Kernel  *k1, *k2;

   double  *equiv;
   double  *counts;
   double   neighbours[K->nelems];

   /* counters */
   long int group_idx = 1;             /* label for the next group     */
   double   group_total = 0.0;         /* sum of all labeled groups    */

   double   min_label;
   int      num_matches;

   /* split the Kernel into forward and backwards kernels */
   k1 = new_kernel(K->nelems);
   k2 = new_kernel(K->nelems);
   split_kernel(K, k1, k2);

   setup_pad_values(k1);
   setup_pad_values(k2);

   fprintf(stdout, "Kernel 1\n");
   print_kernel(k1);
   fprintf(stdout, "Kernel 2\n");
   print_kernel(k2);

   get_volume_sizes(*vol, sizes);
   initialize_progress_report(&progress, FALSE, sizes[2], "Groups");

   /* copy the volume */
   tmp_vol = copy_volume_definition(*vol, NC_DOUBLE, FALSE, 0.0, 0.0);

   /* Count groups outputting as we go */
   if(verbose){
      fprintf(stdout, "# Group kernel counts\n");
      fprintf(stdout, "#\n");
      fprintf(stdout, "#     Group Range: %14g -> %14g\n", gfloor, gceil);
      fprintf(stdout, "#     Max Groups:  %d\n", max_groups);
      fprintf(stdout, "#\n");
      fprintf(stdout, "#  group      position         count\n");
      }

   /* initialise the equiv and counts arrays */
   SET_ARRAY_SIZE(equiv, 0, group_idx + 1, 500);
   equiv[0] = 0.0;
   equiv[1] = 1.0;

   SET_ARRAY_SIZE(counts, 0, group_idx + 1, 500);
   counts[0] = 0.0;
   counts[1] = 0.0;

   for(c = 0; c < max_groups; c++){
      counts[c] = 0;
      }

   /* zero the sucker */
   for(z = 0; z < sizes[0]; z++){
      for(y = 0; y < sizes[1]; y++){
         for(x = 0; x < sizes[2]; x++){
            set_volume_real_value(tmp_vol, z, y, x, 0, 0, 0.0);
            }
         }
      }

   /* pass 1 - forward direction (we assume a symmetric kernel) */
   for(z = -k1->pre_pad[2]; z < sizes[0] - k1->post_pad[2]; z++){
      for(y = -k1->pre_pad[1]; y < sizes[1] - k1->post_pad[1]; y++){
         for(x = -k1->pre_pad[0]; x < sizes[2] - k1->post_pad[0]; x++){

            if(get_volume_real_value(*vol, z, y, x, 0, 0) > 0){

               /* check this voxels neighbours */
               num_matches = 0;
               min_label = 10000000.0;

//               fprintf(stdout, "Value: [%d][%d][%d] ", x, y, z);
               for(c = 0; c < k1->nelems; c++){
                  value = get_volume_real_value(tmp_vol,
                                                z + k1->K[c][2],
                                                y + k1->K[c][1],
                                                x + k1->K[c][0], 0 + k1->K[c][3], 0 + k1->K[c][4]);

                  if(value > 0.0){
                     if(value < min_label){
                        min_label = value;
                        }
                     neighbours[num_matches] = value;
                     num_matches++;
//                     fprintf(stdout, "%g-", value);
                     }
                  }
//               fprintf(stdout, " [%g]", min_label);

               switch (num_matches){

               case 0:
                  /* no neighbours, make a new label and increment */
                  set_volume_real_value(tmp_vol, z, y, x, 0, 0, (double)group_idx);

                  SET_ARRAY_SIZE(equiv, group_idx + 1, group_idx + 2, 500);
                  equiv[group_idx] = (double)group_idx;

                  SET_ARRAY_SIZE(counts, group_idx + 1, group_idx + 2, 500);
                  counts[group_idx] = 0.0;

//                     fprintf(stdout, "  =====  %d", group_idx);

                  counts[group_idx]++;
                  group_idx++;
                  break;

               case 1:
                  /* only one neighbour, no equivalences needed */
                  set_volume_real_value(tmp_vol, z, y, x, 0, 0, min_label);
                  counts[(int)min_label]++;
//                     fprintf(stdout, "  +++++  %g", min_label);
                  break;

               default:
                  /* more than one neighbour, set some equivalences */
                  set_volume_real_value(tmp_vol, z, y, x, 0, 0, min_label);
                  counts[(int)min_label]++;
//                     fprintf(stdout, "  +++++  %g", min_label);

                  for(c = 0; c < num_matches; c++){

                     if(min_label < equiv[(int)neighbours[c]]){
                        /* set an equivalence */
//                           fprintf(stdout, "\n      E [%d]... %d => %g\n", c, (int)neighbours[c], min_label);
                        equiv[(int)neighbours[c]] = min_label;

                        }
                     }
                  break;
               }                       /* end case */

//               fprintf(stdout, "\n");
               }
            }
         }
      update_progress_report(&progress, z + 1);
      }


   /* reduce the equiv and counts array */
   for(c = 0; c < max_groups; c++){
      int      current;

      fprintf(stdout, "%d => %g\t\t%g", c, equiv[c], counts[c]);

      /* find the lowest label for this label */
      if(c != (int)equiv[c]){

         current = (int)equiv[c];
         while (current != equiv[current]){
            current = (int)equiv[current];
            }

         equiv[c] = current;
         }

      fprintf(stdout, " ===> %g\n", equiv[c]);
      }

   /* pass 2 - resolve equivalences in the output data */
   for(z = 0; z < sizes[0]; z++){
      for(y = 0; y < sizes[1]; y++){
         for(x = 0; x < sizes[2]; x++){

            value = get_volume_real_value(tmp_vol, z, y, x, 0, 0);
            set_volume_real_value(*vol, z, y, x, 0, 0, equiv[(int)value]);
            }
         }
      }

   if(verbose){
      fprintf(stdout, "#\n");
      fprintf(stdout, "#   Total voxels: %14g\n", (double)(sizes[0] * sizes[1] * sizes[2]));
      fprintf(stdout, "#       # voxels: %14g\n", (double)group_total);
      fprintf(stdout, "#       # groups: %14g\n", (double)group_idx);
      fprintf(stdout, "#     %% of total: %14g\n",
              (double)(group_total / (sizes[0] * sizes[1] * sizes[2]) * 100));
      }

   delete_volume(tmp_vol);
   free(k1);
   free(k2);
   terminate_progress_report(&progress);
   return (vol);
   }
