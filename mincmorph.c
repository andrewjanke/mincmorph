/* mincmorph.c                                                               */
/*                                                                           */
/* Andrew Janke - rotor@cmr.uq.edu.au                                        */
/* Center for Magnetic Resonance                                             */
/* University of Queensland                                                  */
/*                                                                           */
/* Copyright Andrew Janke, The University of Queensland.                     */
/* Permission to use, copy, modify, and distribute this software and its     */
/* documentation for any purpose and without fee is hereby granted,          */
/* provided that the above copyright notice appear in all copies.  The       */
/* author and the University of Queensland make no representations about the */
/* suitability of this software for any purpose.  It is provided "as is"     */
/* without express or implied warranty.                                      */
/*                                                                           */
/* Morphology on a minc volume...                                            */
/*    o Erosions                                                             */
/*    o Dilations                                                            */
/*    o Group counts (mm and voxel)                                          */
/*    o Arbitrary kernels                                                    */
/*                                                                           */
/* Fri Jan 18 11:41:08 EST 2002 - initial version from mincgroup             */
/* Mon Jan 28 11:41:37 EST 2002 - first version that works                   */


#include <stdlib.h>
#include <sys/param.h>
#include <unistd.h>
#include <float.h>
#include <volume_io.h>
#include <ParseArgv.h>
#include <time_stamp.h>
#include "kernel_io.h"
#include "kernel_ops.h"

/* typedefs */
typedef enum {
   UNDEF = 0,
   BINARISE, GROUP, CONVOLVE, ERODE, DILATE,
   OPEN, CLOSE, LPASS, HPASS, PAD, KEEP,
   DISTANCE, READ_KERNEL, WRITE
} op_types;

typedef struct {
   op_types type;
   Kernel  *kern;
   char    *outfile;
   double   binarise_range[2];
   double   group_range[2];
   int      max_groups;
   int      num_groups;
   double   pad_value;
} Operation;

/* Argument variables */
int      verbose = FALSE;
int      clobber = FALSE;
double   vol_range[2] = { -DBL_MAX, DBL_MAX };
double   group_range[2] = { 0.0, DBL_MAX };
char    *kernel_fn = NULL;
int      max_groups = 250;
char    *succ_txt = "B";

char     successive_help[] = "Successive operations (Maximum: 100) \
\n\tB[floor:ceil] - binarise \
\n\tG[gfloor:gceil:max_groups] - group \
\n\tX - convolve \
\n\tE - erosion \
\n\tD - dilation \
\n\tO - open \
\n\tC - close \
\n\tL - lowpass filter \
\n\tH - highpass filter \
\n\tP[n] - pad volume with respect to current kernel using n as the fill-value (default: 0)\
\n\tK[n] - keep the largest n objects (default: 1) \
\n\tF - distance transform (binary input only - not checked) \
\n\tR[file.kern] - read in a kernel file \
\n\tW[file.mnc]  - write out current results \
\n\tDefault: ";

/* Argument table */
ArgvInfo argTable[] = {
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "be verbose"},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "clobber existing files"},
   {"-kernel", ARGV_STRING, (char *)1, (char *)&kernel_fn,
    "<kernel.kern> read in a kernel file"},
   {"-floor", ARGV_FLOAT, (char *)1, (char *)&vol_range[0],
    "ignore voxels below this value"},
   {"-ceil", ARGV_FLOAT, (char *)1, (char *)&vol_range[1],
    "ignore voxels above this value (incl)"},
   {"-range", ARGV_FLOAT, (char *)2, (char *)vol_range,
    "ignore voxels outside the range (incl)"},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nGroup size options:"},
   {"-group_floor", ARGV_FLOAT, (char *)1, (char *)&group_range[0],
    "ignore groups below this size"},
   {"-group_ceil", ARGV_FLOAT, (char *)1, (char *)&group_range[1],
    "ignore groups above this size  (incl)"},
   {"-group_range", ARGV_FLOAT, (char *)2, (char *)group_range,
    "ignore groups outside the range (incl)"},

   {"-max_groups", ARGV_INT, (char *)1, (char *)&max_groups,
    "the maximum number of groups to segment"},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nSingle morphological operations:"},
   {"-binarise", ARGV_CONSTANT, (char *)"B", (char *)&succ_txt,
    "binarise volume (use -range for range)"},
   {"-group", ARGV_CONSTANT, (char *)"G", (char *)&succ_txt,
    "find groups (use -group_range for size)"},
   {"-convolve", ARGV_CONSTANT, (char *)"X", (char *)&succ_txt,
    "convolve file with kernel"},
   {"-erosion", ARGV_CONSTANT, (char *)"E", (char *)&succ_txt,
    "erosion"},
   {"-dilation", ARGV_CONSTANT, (char *)"D", (char *)&succ_txt,
    "dilation"},
   {"-open", ARGV_CONSTANT, (char *)"O", (char *)&succ_txt,
    "open:            dilation(erosion(X))"},
   {"-close", ARGV_CONSTANT, (char *)"C", (char *)&succ_txt,
    "close:           erosion(dilation(X))"},
   {"-lowpass", ARGV_CONSTANT, (char *)"L", (char *)&succ_txt,
    "lowpass filter:  close(open(X))"},
   {"-highpass", ARGV_CONSTANT, (char *)"H", (char *)&succ_txt,
    "highpass filter: X - lowpass(X)"},
   {"-pad", ARGV_CONSTANT, (char *)"P", (char *)&succ_txt,
    "pad volume with respect to current kernel"},
   {"-keep", ARGV_CONSTANT, (char *)"K", (char *)&succ_txt,
    "keep the largest cluster"},
   {"-distance", ARGV_CONSTANT, (char *)"F", (char *)&succ_txt,
    "distance transform (implicit pad before)"},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nSuccessive morphological operations:"},
   {"-successive", ARGV_STRING, (char *)1, (char *)&succ_txt, successive_help},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nXXXXXXX DON'T USE THESE! XXXXXXX"},
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "XXX compatibility with binop XXX"},
   {"-objlimit", ARGV_INT, (char *)1, (char *)&max_groups,
    "maps to -max_groups"},
   {"-thres", ARGV_FLOAT, (char *)2, (char *)vol_range,
    "maps to -range"},
   {"-obinary", ARGV_CONSTANT, (char *)"B", (char *)&succ_txt,
    "maps to -binarise"},

   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

int main(int argc, char *argv[])
{
   int      offset, op_c, c, num_matched;
   char    *arg_string;
   char    *infile;
   char    *outfile;

   Volume  *volume;
   Kernel  *kern;
   int      num_ops;
   Operation operation[100];
   char     tmp_str[256];
   char     ext_txt[256];
   char     tmp_filename[MAXPATHLEN];

   char    *axis_order[3] = { MIzspace, MIyspace, MIxspace };

   /* Save time stamp and args */
   arg_string = time_stamp(argc, argv);

   /* Get arguments */
   if(ParseArgv(&argc, argv, argTable, 0) || (argc < 2)){
      (void)fprintf(stderr, "\nUsage: %s [options] <in.mnc> <out.mnc>\n", argv[0]);
      (void)fprintf(stderr, "       %s -help\n\n", argv[0]);
      exit(EXIT_FAILURE);
      }
   infile = argv[1];
   outfile = argv[2];

   /* check for the infile */
   if(access(infile, F_OK) != 0){
      fprintf(stderr, "%s: Couldn't find %s\n", argv[0], infile);
      exit(EXIT_FAILURE);
      }

   /* check for the outfile */
   if(access(outfile, F_OK) != -1 && !clobber){
      fprintf(stderr, "%s: %s exists! (use -clobber to overwrite)\n", argv[0], outfile);
      exit(EXIT_FAILURE);
      }

   /* setup the default kernel or read in a kernel file */
   kern = new_kernel(6);
   if(kernel_fn == NULL){
      setup_def_kernel(kern);
      }
   else{
      if(access(kernel_fn, F_OK) != 0){
         fprintf(stderr, "%s: Couldn't find kernel file: %s\n", argv[0], kernel_fn);
         exit(EXIT_FAILURE);
         }

      if(input_kernel(kernel_fn, kern) != OK){
         fprintf(stderr, "%s: Died whilst reading in kernel file: %s\n", argv[0], kernel_fn);
         exit(EXIT_FAILURE);
         }
      if(verbose){
         fprintf(stdout, "Input kernel:\n");
         print_kernel(kern);
         }
      }
   setup_pad_values(kern);

   /* setup operations and check them... */
   num_ops = strlen(succ_txt);
   op_c = 0;
   if(verbose){
      fprintf(stdout, "---Checking Operation(s): %s---\n", succ_txt);
      }
   for(c = 0; c < num_ops; c++){

      /* set up counters and extra text */
      offset = 0;
      strcpy(ext_txt, "");

      switch (succ_txt[c]){
      case 'B':
         operation[op_c].type = BINARISE;

         /* look for a range else use the default C/L ones */
         if(succ_txt[c + 1] == '['){
            while (succ_txt[c + offset] != ']'){

               /* check that we aren't off the edge of the world */
               if(succ_txt[c + offset] == '\0'){
                  fprintf(stderr,
                          "%s: I think you missed the closing ']' on the B[floor:ceil] ? : %s\n",
                          argv[0], tmp_str);
                  exit(EXIT_FAILURE);
                  }

               sscanf(&succ_txt[c + offset], "%1s", &tmp_str[offset]);
               offset++;
               }

            if(strlen(&tmp_str[2]) == 0){
               fprintf(stderr, "%s: missing numbers in the B[floor:ceil] argument\n", argv[0]);
               exit(EXIT_FAILURE);
               }

            sscanf(&tmp_str[2], "%lf:%lf", &operation[op_c].binarise_range[0],
                   &operation[op_c].binarise_range[1]);

            /* figger out the offset to the next operator */
            while (succ_txt[c + offset] != ']'){
               offset++;
               }
            }
         else{
            operation[op_c].binarise_range[0] = vol_range[0];
            operation[op_c].binarise_range[1] = vol_range[1];
            }

         sprintf(ext_txt, "range: [%g:%g]", operation[op_c].binarise_range[0],
                 operation[op_c].binarise_range[1]);
         break;

      case 'G':
         operation[op_c].type = GROUP;

         /* look for a range else use the default C/L ones */
         if(succ_txt[c + 1] == '['){
            while (succ_txt[c + offset] != ']'){

               /* check that we aren't off the edge of the world */
               if(succ_txt[c + offset] == '\0'){
                  fprintf(stderr,
                          "%s: I think you missed the closing ']' on the G[gfloor:gceil:max_groups] ? : %s\n",
                          argv[0], tmp_str);
                  exit(EXIT_FAILURE);
                  }

               sscanf(&succ_txt[c + offset], "%1s", &tmp_str[offset]);
               offset++;
               }

            if(strlen(&tmp_str[2]) == 0){
               fprintf(stderr, "%s: missing numbers in the G[gfloor:gceil:max_groups] argument\n",
                       argv[0]);
               exit(EXIT_FAILURE);
               }


            num_matched = sscanf(&tmp_str[2], "%lf:%lf:%d", &operation[op_c].group_range[0],
                                 &operation[op_c].group_range[1], &operation[op_c].max_groups);

            if(operation[op_c].max_groups > 250){
               fprintf(stderr, "%s: Input number of groups (%d) exceeds 250\n",
                       argv[0], operation[op_c].max_groups);
               exit(EXIT_FAILURE);
               }

            /* check for a missing number of groups */
            if(num_matched == 2){
               operation[op_c].max_groups = max_groups;
               }

            }
         else{
            operation[op_c].group_range[0] = group_range[0];
            operation[op_c].group_range[1] = group_range[1];
            operation[op_c].max_groups = max_groups;
            }

         sprintf(ext_txt, "group range: [%g:%g]  max # groups: %d",
                 operation[op_c].group_range[0], operation[op_c].group_range[1],
                 operation[op_c].max_groups);
         break;

      case 'X':
         operation[op_c].type = CONVOLVE;
         break;

      case 'E':
         operation[op_c].type = ERODE;
         break;

      case 'D':
         operation[op_c].type = DILATE;
         break;

      case 'O':
         operation[op_c].type = OPEN;
         break;

      case 'C':
         operation[op_c].type = CLOSE;
         break;

      case 'L':
         operation[op_c].type = LPASS;
         break;

      case 'H':
         operation[op_c].type = HPASS;
         break;

      case 'P':
         operation[op_c].type = PAD;

         /* look for a range else use the default C/L ones */
         if(succ_txt[c + 1] == '['){
            while (succ_txt[c + offset] != ']'){

               /* check that we aren't off the edge of the world */
               if(succ_txt[c + offset] == '\0'){
                  fprintf(stderr, "%s: I think you missed the closing ']' on the P[n] ? : %s\n",
                          argv[0], tmp_str);
                  exit(EXIT_FAILURE);
                  }

               sscanf(&succ_txt[c + offset], "%1s", &tmp_str[offset]);
               offset++;
               }

            if(strlen(&tmp_str[2]) == 0){
               fprintf(stderr, "%s: missing number in P[n] argument\n", argv[0]);
               exit(EXIT_FAILURE);
               }

            sscanf(&tmp_str[2], "%lf", &operation[op_c].pad_value);

            /* figger out the offset to the next operator */
            while (succ_txt[c + offset] != ']'){
               offset++;
               }
            }
         else{
            operation[op_c].pad_value = 0;
            }

         sprintf(ext_txt, "fill value: %g", operation[op_c].pad_value);
         break;

      case 'K':
         operation[op_c].type = KEEP;

         /* look for a range else use the default C/L ones */
         if(succ_txt[c + 1] == '['){

            /* find the end of the string */
            while (succ_txt[c + offset] != ']'){

               /* check that we aren't off the edge of the world */
               if(succ_txt[c + offset] == '\0'){
                  fprintf(stderr, "%s: I think you missed the closing ']' on the K[n] ? : %s\n",
                          argv[0], tmp_str);
                  exit(EXIT_FAILURE);
                  }

               sscanf(&succ_txt[c + offset], "%1s", &tmp_str[offset]);
               offset++;
               }

            if(strlen(&tmp_str[2]) == 0){
               fprintf(stderr, "%s: missing number in K[n] argument\n", argv[0]);
               exit(EXIT_FAILURE);
               }

            sscanf(&tmp_str[2], "%d", &operation[op_c].num_groups);
            }
         else{
            operation[op_c].num_groups = 1;
            }

         sprintf(ext_txt, "keep %d largest group(s)", operation[op_c].num_groups);
         break;

      case 'F':
         operation[op_c].type = DISTANCE;
         break;

      case 'R':
         operation[op_c].type = READ_KERNEL;

         /* get the filename */
         if(succ_txt[c + 1] == '['){

            /* find the end of the string */
            while (succ_txt[c + offset] != ']'){
               if(succ_txt[c + offset] == '\0'){
                  fprintf(stderr,
                          "%s: I think you missed the closing ']' on the R[file.kern] ? : %s\n",
                          argv[0], tmp_str);
                  exit(EXIT_FAILURE);
                  }

               sscanf(&succ_txt[c + offset], "%1s", &tmp_str[offset]);
               offset++;
               }

            operation[op_c].kern = new_kernel(0);

            /* set up the real filename */
            realpath(&tmp_str[2], tmp_filename);

            fprintf(stderr, "\nXXXX\nKernel File: %s\nXXXX\n\n", tmp_filename);
            if(access(tmp_filename, F_OK) != 0){
               fprintf(stderr, "%s: Couldn't find kernel file: %s\n", argv[0], tmp_filename);
               exit(EXIT_FAILURE);
               }

            if(input_kernel(&tmp_str[2], operation[op_c].kern) != OK){
               fprintf(stderr, "%s: Died whilst reading in kernel file: %s\n", argv[0],
                       &tmp_str[2]);
               exit(EXIT_FAILURE);
               }
            }
         else{
            fprintf(stderr, "%s: R[file.kern] _requires_ a kernel filename\n", argv[0]);
            exit(EXIT_FAILURE);
            }

         sprintf(ext_txt, "kernel: %s", &tmp_str[2]);
         break;

      case 'W':
         operation[op_c].type = WRITE;

         /* get the filename */
         if(succ_txt[c + 1] == '['){
            while (succ_txt[c + offset] != ']'){

               /* check that we aren't off the edge of the world */
               if(succ_txt[c + offset] == '\0'){
                  fprintf(stderr,
                          "%s: I think you missed the closing ']' on the W[file.mnc] ? : %s\n",
                          argv[0], tmp_str);
                  exit(EXIT_FAILURE);
                  }

               sscanf(&succ_txt[c + offset], "%1s", &tmp_str[offset]);
               offset++;
               }

            operation[op_c].outfile = (char *)malloc((offset - 2) * sizeof(char));
            strcpy(operation[op_c].outfile, &tmp_str[2]);

            if(access(operation[op_c].outfile, F_OK) != -1 && !clobber){
               fprintf(stderr, "%s: %s exists! (use -clobber to overwrite)\n", argv[0],
                       operation[op_c].outfile);
               exit(EXIT_FAILURE);
               }

            }
         else{
            fprintf(stderr, "%s: W[file.mnc] _requires_ a filename\n", argv[0]);
            exit(EXIT_FAILURE);
            }

         sprintf(ext_txt, "filename: %s", operation[op_c].outfile);
         break;

      default:
         fprintf(stderr, "\nUnknown op: %c\n\n  %s -help  for operations\n\n", succ_txt[c],
                 argv[0]);
         exit(EXIT_FAILURE);
         }

      if(verbose){
         fprintf(stdout, "  [%02d+%02d] Op[%02d] %c = %d\t\t%s\n", c, offset, op_c,
                 succ_txt[c], operation[op_c].type, ext_txt);
         }

      op_c++;
      c += offset;
      }
   num_ops = op_c;

   /* add an implicit write statment to the end if needed */
   if(operation[num_ops - 1].type != WRITE){
      operation[num_ops].type = WRITE;
      operation[num_ops].outfile = outfile;
      num_ops++;
      }

   /* malloc space for volume structure and read in infile */
   volume = (Volume *) malloc(sizeof(Volume));
   input_volume(infile, MAX_VAR_DIMS, axis_order, NC_UNSPECIFIED, TRUE, 0.0, 0.0, TRUE, volume,
                NULL);

   /* do some operations.... */
   if(verbose){
      fprintf(stdout, "\n---Doing %d Operation(s)---\n", num_ops);
      }
   for(c = 0; c < num_ops; c++){
      switch (operation[c].type){
      case BINARISE:
         volume = binarise(volume, operation[c].binarise_range[0], operation[c].binarise_range[1]);
         break;

      case GROUP:
         volume = group_kernel(kern, volume, outfile,
                               operation[c].group_range[0], operation[c].group_range[1],
                               operation[c].max_groups);
         break;

      case CONVOLVE:
         volume = convolve_kernel(kern, volume);
         break;

      case ERODE:
         volume = erosion_kernel(kern, volume);
         break;

      case DILATE:
         volume = dilation_kernel(kern, volume);
         break;

      case OPEN:
         volume = erosion_kernel(kern, volume);
         volume = dilation_kernel(kern, volume);
         break;

      case CLOSE:
         volume = dilation_kernel(kern, volume);
         volume = erosion_kernel(kern, volume);
         break;

      case LPASS:
         volume = erosion_kernel(kern, volume);
         volume = dilation_kernel(kern, volume);
         volume = dilation_kernel(kern, volume);
         volume = erosion_kernel(kern, volume);
         break;

      case HPASS:
         fprintf(stderr, "GNFARK! Highpass Not implemented yet..\n");
         break;

      case PAD:
         volume = pad_volume(kern, volume, operation[c].pad_value);
         break;

      case KEEP:
         fprintf(stderr, "GNFARK! Keep not implemented yet..\n");
         break;

      case DISTANCE:
         volume = distance_kernel(kern, volume);
         break;

      case READ_KERNEL:
         /* free the existing kernel then set the pointer to the new one */
         free(kern);
         kern = operation[c].kern;
         setup_pad_values(kern);
         if(verbose){
            fprintf(stdout, "Input kernel:\n");
            print_kernel(kern);
            }

         break;

      case WRITE:
         if(operation[c].outfile != NULL){
            if(verbose){
               fprintf(stdout, "Outputting to %s\n", operation[c].outfile);
               }
            output_modified_volume(operation[c].outfile,
                                   NC_UNSPECIFIED, FALSE,
                                   0.0, 0.0, *volume, infile, arg_string, NULL);
            }
         else{
            fprintf(stdout, "WRITE passed a NULL pointer! - this is bad\n");
            exit(EXIT_FAILURE);
            }
         break;

      default:
         fprintf(stderr, "\nUnknown op: %c\n\n  %s -help  for ops\n\n", succ_txt[c], argv[0]);
         exit(EXIT_FAILURE);
         }
      }

   delete_volume(*volume);
   exit(EXIT_SUCCESS);
   }
