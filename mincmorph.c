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
/* Mon Nov 25 16:32:54 EST 2002 - rewrote most of the successive text parser */
/*                              - simplified Group option for flexibility    */
/*                              - replaced K (Keep) option with Clamp.       */
/*                              - introduced notion of fore and back colors  */

#include <config.h>
#include <stdlib.h>
#include <sys/param.h>
#include <unistd.h>
#include <float.h>
#include <volume_io.h>
#include <ParseArgv.h>
#include <time_stamp.h>
#include "kernel_io.h"
#include "kernel_ops.h"

#define INTERNAL_PREC NC_FLOAT         /* should be NC_FLOAT or NC_DOUBLE */
#define DEF_DOUBLE -1

/* function prototypes */
char    *get_real_from_string(char *string, double *value);
char    *get_string_from_string(char *string, char **value);
void     calc_volume_range(Volume * vol, double *min, double *max);
void     print_version_info(void);

/* typedefs */
typedef enum {
   UNDEF = 0,
   BINARISE, CLAMP, PAD, ERODE, DILATE,
   OPEN, CLOSE, LPASS, HPASS, CONVOLVE,
   DISTANCE, GROUP, READ_KERNEL, WRITE
   } op_types;

typedef struct {
   op_types type;
   char     op_c;
   Kernel  *kern;
   char    *outfile;
   double   range[2];
   double   foreground;
   double   background;
   } Operation;

/* Argument variables */
int      verbose = FALSE;
int      clobber = FALSE;
nc_type  dtype = NC_SHORT;
double   range[2] = { -DBL_MAX, DBL_MAX };
double   foreground = 1.0;
double   background = 0.0;
char    *kernel_fn = NULL;
char    *succ_txt = "B";

char     successive_help[] = "Successive operations (Maximum: 100) \
\n\tB[floor:ceil:fg:bg] - binarise in the range, using foreground and background \
\n\tK[floor:ceil:bg] - clamp betwen the specified range. Set other voxels to 'bg' (default: 0) \
\n\tP[bg] - pad volume with respect to the current kernel using 'bg' (default: 0)\
\n\tE - erosion \
\n\tD - dilation \
\n\tO - open \
\n\tC - close \
\n\tL - lowpass filter \
\n\tH - highpass filter \
\n\tX - convolve \
\n\tF - distance transform (binary input only - not checked) \
\n\tG - Label the groups in the volume in ascending order \
\n\tR[file.kern] - read in a kernel file \
\n\tW[file.mnc]  - write out current results \
\n\tDefault: ";

/* Argument table */
ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "General options:"},
   {"-version", ARGV_FUNC, (char *)print_version_info, (char *)NULL,
    "print version info and exit"},
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "be verbose"},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "clobber existing files"},

   {NULL, ARGV_HELP, NULL, NULL,
    "\nOutfile Options (NB: this effects internal precision, int or short suggested)"},
   {"-filetype", ARGV_CONSTANT, (char *)NC_UNSPECIFIED, (char *)&dtype,
    "Use data type of the input file."},
   {"-byte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&dtype,
    "Write out byte data."},
   {"-short", ARGV_CONSTANT, (char *)NC_SHORT, (char *)&dtype,
    "Write out short integer data. (Default)"},
   {"-int", ARGV_CONSTANT, (char *)NC_INT, (char *)&dtype,
    "Write out long integer data."},
   {"-float", ARGV_CONSTANT, (char *)NC_FLOAT, (char *)&dtype,
    "Write out single-precision data."},
   {"-double", ARGV_CONSTANT, (char *)NC_DOUBLE, (char *)&dtype,
    "Write out double-precision data."},

   {NULL, ARGV_HELP, NULL, NULL, "\nMorphology Options"},
   {"-kernel", ARGV_STRING, (char *)1, (char *)&kernel_fn,
    "<kernel.kern> read in a kernel file"},
   {"-floor", ARGV_FLOAT, (char *)1, (char *)&range[0],
    "lowwer value for binarising or clamping"},
   {"-ceil", ARGV_FLOAT, (char *)1, (char *)&range[1],
    "upper value for binarising or clamping (incl)"},
   {"-range", ARGV_FLOAT, (char *)2, (char *)range,
    "range for binarising or clamping (incl)"},
   {"-foreground", ARGV_FLOAT, (char *)1, (char *)&foreground,
    "foreground value"},
   {"-background", ARGV_FLOAT, (char *)1, (char *)&background,
    "background value"},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, "\nSingle morphological operations:"},
   {"-binarise", ARGV_CONSTANT, (char *)"B", (char *)&succ_txt,
    "binarise volume using the input range"},
   {"-clamp", ARGV_CONSTANT, (char *)"K", (char *)&succ_txt,
    "clamp volume using the input range"},
   {"-pad", ARGV_CONSTANT, (char *)"P", (char *)&succ_txt,
    "pad volume with respect to the current kernel (-background to specify value)"},
   {"-erosion", ARGV_CONSTANT, (char *)"E", (char *)&succ_txt,
    "do a single erosion"},
   {"-dilation", ARGV_CONSTANT, (char *)"D", (char *)&succ_txt,
    "do a single dilation"},
   {"-open", ARGV_CONSTANT, (char *)"O", (char *)&succ_txt,
    "open:            dilation(erosion(X))"},
   {"-close", ARGV_CONSTANT, (char *)"C", (char *)&succ_txt,
    "close:           erosion(dilation(X))"},
   {"-lowpass", ARGV_CONSTANT, (char *)"L", (char *)&succ_txt,
    "lowpass filter:  close(open(X))"},
   {"-highpass", ARGV_CONSTANT, (char *)"H", (char *)&succ_txt,
    "highpass filter: X - lowpass(X)"},
   {"-convolve", ARGV_CONSTANT, (char *)"X", (char *)&succ_txt,
    "convolve file with kernel"},
   {"-distance", ARGV_CONSTANT, (char *)"F", (char *)&succ_txt,
    "distance transform"},
   {"-group", ARGV_CONSTANT, (char *)"G", (char *)&succ_txt,
    "label groups in ascending order"},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "\nSuccessive morphological operations:"},
   {"-successive", ARGV_STRING, (char *)1, (char *)&succ_txt, successive_help},

   {NULL, ARGV_HELP, NULL, NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
   };

int main(int argc, char *argv[])
{
   int      c;
   char    *arg_string;
   char    *infile;
   char    *outfile;

   Volume  *volume;
   Kernel  *kern;
   int      num_ops;
   Operation operation[100];
   Operation *op;
   char    *tmp_str;
   char     ext_txt[256];
   char     tmp_filename[MAXPATHLEN];
   double   tmp_double[4];
   double   min, max;
   char    *ptr;

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
      fprintf(stderr, "%s: Couldn't find %s\n\n", argv[0], infile);
      exit(EXIT_FAILURE);
      }

   /* check for the outfile */
   if(access(outfile, F_OK) == 0 && !clobber){
      fprintf(stderr, "%s: %s exists! (use -clobber to overwrite)\n\n", argv[0], outfile);
      exit(EXIT_FAILURE);
      }

   /* setup the default kernel or read in a kernel file */
   kern = new_kernel(6);
   if(kernel_fn == NULL){
      setup_def_kernel(kern);
      }
   else {
      if(access(kernel_fn, F_OK) != 0){
         fprintf(stderr, "%s: Couldn't find kernel file: %s\n\n", argv[0], kernel_fn);
         exit(EXIT_FAILURE);
         }

      if(input_kernel(kernel_fn, kern) != OK){
         fprintf(stderr, "%s: Died reading in kernel file: %s\n\n", argv[0], kernel_fn);
         exit(EXIT_FAILURE);
         }
      if(verbose){
         fprintf(stdout, "Input kernel:\n");
         print_kernel(kern);
         }
      }
   setup_pad_values(kern);

   /* setup operations and check them... */
   if(verbose){
      fprintf(stdout, "---Checking Operation(s): %s---\n", succ_txt);
      }

   num_ops = 0;
   ptr = succ_txt;
   while(ptr[0] != '\0'){

      /* set up counters and extra text */
      strcpy(ext_txt, "");
      op = &operation[num_ops++];

      /* get the operation type */
      op->op_c = ptr[0];
      ptr++;

      switch (op->op_c){
      case 'B':
         op->type = BINARISE;

         /* get 4 possible values */
         ptr = get_real_from_string(ptr, &tmp_double[0]);
         ptr = get_real_from_string(ptr, &tmp_double[1]);
         ptr = get_real_from_string(ptr, &tmp_double[2]);
         ptr = get_real_from_string(ptr, &tmp_double[3]);
         op->range[0] = (tmp_double[0] == DEF_DOUBLE) ? range[0] : tmp_double[0];
         op->range[1] = (tmp_double[1] == DEF_DOUBLE) ? range[1] : tmp_double[1];
         op->foreground = (tmp_double[2] == DEF_DOUBLE) ? foreground : tmp_double[2];
         op->background = (tmp_double[3] == DEF_DOUBLE) ? background : tmp_double[3];

         sprintf(ext_txt, "range: [%g:%g] fg/bg: [%g:%g]", op->range[0],
                 op->range[1], op->foreground, op->background);
         break;

      case 'K':
         op->type = CLAMP;

         /* get 3 possible values */
         ptr = get_real_from_string(ptr, &tmp_double[0]);
         ptr = get_real_from_string(ptr, &tmp_double[1]);
         ptr = get_real_from_string(ptr, &tmp_double[2]);
         op->range[0] = (tmp_double[0] == DEF_DOUBLE) ? range[0] : tmp_double[0];
         op->range[1] = (tmp_double[1] == DEF_DOUBLE) ? range[1] : tmp_double[1];
         op->background = (tmp_double[2] == DEF_DOUBLE) ? background : tmp_double[2];

         sprintf(ext_txt, "range: [%g:%g] fg/bg: [%g:%g]", op->range[0],
                 op->range[1], op->foreground, op->background);
         break;

      case 'P':
         op->type = PAD;

         /* get 1 possible value */
         ptr = get_real_from_string(ptr, &tmp_double[0]);
         op->background = (tmp_double[0] == DEF_DOUBLE) ? background : tmp_double[0];

         sprintf(ext_txt, "fill value: %g", op->background);
         break;

      case 'E':
         op->type = ERODE;
         break;

      case 'D':
         op->type = DILATE;
         break;

      case 'O':
         op->type = OPEN;
         break;

      case 'C':
         op->type = CLOSE;
         break;

      case 'L':
         op->type = LPASS;
         break;

      case 'H':
         op->type = HPASS;
         break;

      case 'X':
         op->type = CONVOLVE;
         break;

      case 'F':
         op->type = DISTANCE;
         break;

      case 'G':
         op->type = GROUP;
         break;

      case 'R':
         op->type = READ_KERNEL;
         op->kern = new_kernel(0);

         /* get the filename */
         ptr = get_string_from_string(ptr, &tmp_str);
         if(tmp_str == NULL){
            fprintf(stderr, "%s: R[file.kern] _requires_ a filename\n\n", argv[0]);
            exit(EXIT_FAILURE);
            }

         /* set up and check for the real filename */
         realpath(tmp_str, tmp_filename);
         if(access(tmp_filename, F_OK) != 0){
            fprintf(stderr, "%s: Couldn't find kernel file: %s\n\n", argv[0],
                    tmp_filename);
            exit(EXIT_FAILURE);
            }

         if(input_kernel(tmp_str, op->kern) != OK){
            fprintf(stderr, "%s: Died reading in kernel file: %s\n\n", argv[0], tmp_str);
            exit(EXIT_FAILURE);
            }

         sprintf(ext_txt, "kernel: %s", tmp_str);
         break;

      case 'W':
         op->type = WRITE;

         /* get the filename */
         ptr = get_string_from_string(ptr, &op->outfile);

         if(op->outfile == NULL){
            fprintf(stderr, "%s: W[file.mnc] _requires_ a filename\n\n", argv[0]);
            exit(EXIT_FAILURE);
            }

         /* check for the outfile */
         if(access(op->outfile, F_OK) == 0 && !clobber){
            fprintf(stderr, "%s: %s exists! (use -clobber to overwrite)\n\n", argv[0],
                    op->outfile);
            exit(EXIT_FAILURE);
            }

         sprintf(ext_txt, "filename: %s", op->outfile);
         break;

      default:
         fprintf(stderr, "\nUnknown op: %c\n\n  %s -help  for operations\n\n",
                 op->op_c, argv[0]);
         exit(EXIT_FAILURE);
         }

      if(verbose){
         fprintf(stdout, "  Op[%02d] %c = %d\t\t%s\n", num_ops, op->op_c, op->type,
                 ext_txt);
         }
      }

   /* add an implicit write statment to the end if needed */
   if(operation[num_ops - 1].type != WRITE){
      operation[num_ops].type = WRITE;
      operation[num_ops].outfile = outfile;
      num_ops++;
      }

   /* malloc space for volume structure and read in infile */
   volume = (Volume *) malloc(sizeof(Volume));
   input_volume(infile, MAX_VAR_DIMS, axis_order,
                INTERNAL_PREC, TRUE, 0.0, 0.0, TRUE, volume, NULL);
   get_type_range(get_volume_data_type(*volume), &min, &max);
   set_volume_real_range(*volume, min, max);

   /* do some operations */
   if(verbose){
      fprintf(stdout, "\n---Doing %d Operation(s)---\n", num_ops);
      }
   for(c = 0; c < num_ops; c++){
      op = &operation[c];

      switch (op->type){
      case BINARISE:
         volume = binarise(volume, op->range[0], op->range[1],
                           op->foreground, op->background);
         break;

      case CLAMP:
         volume = clamp(volume, op->range[0], op->range[1], op->background);
         break;

      case PAD:
         volume = pad(kern, volume, op->background);
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
         fprintf(stderr, "%s: GNFARK! Highpass Not implemented yet..\n\n", argv[0]);
         break;

      case CONVOLVE:
         volume = convolve_kernel(kern, volume);
         break;

      case DISTANCE:
         volume = distance_kernel(kern, volume, background);
         break;

      case GROUP:
         volume = group_kernel(kern, volume, background);
         break;

      case READ_KERNEL:
         /* free the existing kernel then set the pointer to the new one */
         free(kern);
         kern = op->kern;
         setup_pad_values(kern);
         if(verbose){
            fprintf(stdout, "Input kernel:\n");
            print_kernel(kern);
            }
         break;

      case WRITE:
         if(op->outfile != NULL){
            if(verbose){
               fprintf(stdout, "Outputting to %s\n", op->outfile);
               }

            /* set the range to something sensible (if possible) */
            calc_volume_range(volume, &min, &max);
            if(dtype == NC_BYTE){
               if(min >= 0 && max < 255){
                  fprintf(stdout, "BYTE data, setting 1:1 mapping (0-256)\n");
                  min = 0;
                  max = 255;
                  }
               }
            set_volume_real_range(*volume, min, max);

            output_modified_volume(op->outfile,
                                   dtype, FALSE,
                                   0.0, 0.0, *volume, infile, arg_string, NULL);
            }
         else {
            fprintf(stdout, "%s: WRITE passed a NULL pointer! - this is bad\n\n",
                    argv[0]);
            exit(EXIT_FAILURE);
            }
         break;

      default:
         fprintf(stderr, "\n%s: This is vary bad, call Houston\n\n", argv[0]);
         exit(EXIT_FAILURE);
         }
      }

   delete_volume(*volume);
   return (EXIT_SUCCESS);
   }

/* get a real from a char* stream                   */
/* with possible trailing or leading square bracket */
/* and possible leading ':'                         */
/* return the string advanced to the next token or  */
/* as it was input if nothing found                 */
char    *get_real_from_string(char *string, double *value)
{
   char    *ptr;

   /* skip a [ or : else we probably don't belong here */
   if(string[0] == '[' || string[0] == ':'){
      string++;
      }

   /* get a double */
   *value = strtod(&string[0], &ptr);

   /* if nothing found return a default value */
   if(&string[0] == ptr){
      *value = DEF_DOUBLE;
      }

   /* skip over a possible ']' */
   if(ptr[0] == ']'){
      ptr++;
      }

   return ptr;
   }

/* get the string between [ and ] from a char*      */
/* return the string advanced to the next token or  */
/* as it was input if nothing found                 */
char    *get_string_from_string(char *string, char **value)
{
   int      offset;

   /* initialise the return values first */
   *value = NULL;

   /* get the string if there is one */
   if(string[0] == '['){
      string++;

      /* get the string */
      offset = strcspn(string, "]");

      *value = (char *)malloc(offset * sizeof(char));
      strncpy(*value, string, offset);
      string += offset;

      /* skip over a possible ']' */
      if(string[0] == ']'){
         string++;
         }
      }

   return string;
   }

void calc_volume_range(Volume * vol, double *min, double *max)
{

   int      x, y, z;
   int      sizes[MAX_VAR_DIMS];
   double   value;
   progress_struct progress;

   *min = DBL_MAX;
   *max = -DBL_MIN;

   get_volume_sizes(*vol, sizes);

   initialize_progress_report(&progress, FALSE, sizes[2], "Finding Range");
   for(z = sizes[0]; z--;){
      for(y = sizes[1]; y--;){
         for(x = sizes[2]; x--;){

            value = get_volume_voxel_value(*vol, z, y, x, 0, 0);
            if(value < *min){
               *min = value;
               }
            else if(value > *max){
               *max = value;
               }
            }
         }
      update_progress_report(&progress, z + 1);
      }
   terminate_progress_report(&progress);

   if(verbose){
      fprintf(stdout, "Found range of [%g:%g]\n", *min, *max);
      }
   }

void print_version_info(void)
{
   fprintf(stdout, "%s version %s\n", PACKAGE, VERSION);
   fprintf(stdout, "Comments to %s\n", PACKAGE_BUGREPORT);
   fprintf(stdout, "\n");
   exit(EXIT_SUCCESS);
   }
