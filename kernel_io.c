/* kernel_io.c - reads kernel files */

#include "kernel_io.h"
#define MAX_KERNEL_ELEMS 1000

extern int verbose;

static const STRING KERNEL_FILE_HEADER = "MNI Morphology Kernel File";
static const STRING KERNEL_TYPE = "Kernel_Type";
static const STRING NORMAL_KERNEL = "Normal_Kernel";
static const STRING KERNEL = "Kernel";

/* returns a new Kernel struct (pointer)                */
Kernel  *new_kernel(int nelems)
{
   int      i;
   Kernel  *tmp;

   ALLOC_VAR_SIZED_STRUCT(tmp, Real, 10);

   /* allocate for the Kernel Array */
   SET_ARRAY_SIZE(tmp->K, 0, nelems, 10);
   for(i = 0; i < nelems; i++){
      ALLOC(tmp->K[i], KERNEL_DIMS + 1);
      }
   tmp->nelems = nelems;

   return tmp;
   }

/* reads in a Kernel from a file                        */
Status input_kernel(const char *kernel_file, Kernel * kernel)
{
   int      i, j;

   STRING   line;
   STRING   type_name;
   STRING   str;
   Real     tmp_real;
   FILE    *file;

   /* parameter checking */
   if(kernel_file == NULL){
      print_error("input_kernel(): passed NULL FILE.\n");
      return (ERROR);
      }

   file = fopen(kernel_file, "r");
   if(file == NULL){
      print_error("input_kernel(): error opening Kernel file.\n");
      return (ERROR);
      }

   /* okay read the header */
   if(mni_input_string(file, &line, (char)0, (char)0) != OK){
      delete_string(line);
      print_error("input_kernel(): could not read header in file.\n");
      return (ERROR);
      }

   if(!equal_strings(line, KERNEL_FILE_HEADER)){
      delete_string(line);
      print_error("input_kernel(): invalid header in file.\n");
      return (ERROR);
      }

   /* --- read the type of Kernel */
   if(mni_input_keyword_and_equal_sign(file, KERNEL_TYPE, FALSE) != OK){
      return (ERROR);
      }

   if(mni_input_string(file, &type_name, (char)';', (char)0) != OK){
      print_error("input_kernel(): missing kernel type.\n");
      return (ERROR);
      }

   if(mni_skip_expected_character(file, (char)';') != OK){
      return (ERROR);
      }

   if(!equal_strings(type_name, NORMAL_KERNEL)){
      print_error("input_kernel(): invalid kernel type.\n");
      delete_string(type_name);
      return (ERROR);
      }
   delete_string(type_name);

   /* --- read the next string */
   if(mni_input_string(file, &str, (char)'=', (char)0) != OK)
      return (ERROR);

   if(!equal_strings(str, KERNEL)){
      print_error("Expected %s =\n", KERNEL);
      delete_string(str);
      return (ERROR);
      }
   delete_string(str);

   if(mni_skip_expected_character(file, (char)'=') != OK){
      return (ERROR);
      }

   /* now read the elements (lines) of the kernel */
   if(verbose){
      fprintf(stderr, "Reading [%s]", kernel_file);
      }
   for(i = 0; i < MAX_KERNEL_ELEMS; i++){

      /* allocate a bit of memory */
      SET_ARRAY_SIZE(kernel->K, kernel->nelems, kernel->nelems + 1, 10);
      ALLOC(kernel->K[i], KERNEL_DIMS + 1);

      /* get the 5 dimension vectors and the coefficient */
      for(j = 0; j < 6; j++){
         if(mni_input_real(file, &tmp_real) == OK){
            kernel->K[i][j] = tmp_real;
            }
         else {
            /* check for end */
            if(mni_skip_expected_character(file, (char)';') == OK){
               kernel->nelems = i;
               if(verbose){
                  fprintf(stderr, " %dx%d Kernel elements read\n", i, kernel->nelems);
                  }
               return (OK);
               }
            else {
               print_error("input_kernel(): error reading kernel [%d,%d]\n", i + 1,
                           j + 1);
               return (ERROR);
               }
            }
         }
      kernel->nelems++;

      if(verbose){
         fprintf(stderr, ".");
         fflush(stderr);
         }
      }

   /* SHOLDN'T BE REACHED */
   print_error("input_kernel(): Glark! Something is amiss in the State of Kansas\n");
   return (ERROR);
   }

/* pretty print a kernel */
int print_kernel(Kernel * kernel)
{
   int      i, j;

   fprintf(stderr, "           x       y       z       t       v   coeff\n");
   fprintf(stderr, "     -----------------------------------------------\n");
   for(i = 0; i < kernel->nelems; i++){
      fprintf(stderr, "[%02d]", i);
      for(j = 0; j < KERNEL_DIMS + 1; j++){
         fprintf(stderr, "%8.02f", kernel->K[i][j]);
         }
      fprintf(stderr, "\n");
      }

   return (TRUE);
   }

int setup_pad_values(Kernel * kernel)
{
   int      c, n;

   /* init padding values */
   for(n = 0; n < KERNEL_DIMS; n++){
      kernel->pre_pad[n] = 0;
      kernel->post_pad[n] = 0;
      }

   /* find the padding sizes */
   for(c = 0; c < kernel->nelems; c++){
      for(n = 0; n < KERNEL_DIMS; n++){
         if(kernel->K[c][n] < kernel->pre_pad[n]){
            kernel->pre_pad[n] = kernel->K[c][n];
            }

         if(kernel->K[c][n] > kernel->post_pad[n]){
            kernel->post_pad[n] = kernel->K[c][n];
            }
         }
      }

   return (TRUE);
   }

/* create the default 3D 8 connectivity kernel           */
/*            x       y       z       t       v   coeff  */
/*      -----------------------------------------------  */
/* [00]    1.00    0.00    0.00    0.00    0.00    1.00  */
/* [01]   -1.00    0.00    0.00    0.00    0.00    1.00  */
/* [02]    0.00    1.00    0.00    0.00    0.00    1.00  */
/* [03]    0.00   -1.00    0.00    0.00    0.00    1.00  */
/* [04]    0.00    0.00    1.00    0.00    0.00    1.00  */
/* [05]    0.00    0.00   -1.00    0.00    0.00    1.00  */
void setup_def_kernel(Kernel * K)
{
   int      i, j;

   /* first setup the "null" kernel */
   for(i = 0; i < 6; i++){
      for(j = 0; j < 6; j++){
         if(j == 5){
            K->K[i][j] = 1.0;
            }
         else {
            K->K[i][j] = 0.0;
            }
         }
      }

   /* then fool with it */
   K->K[0][0] = 1.0;
   K->K[1][0] = -1.0;
   K->K[2][1] = 1.0;
   K->K[3][1] = -1.0;
   K->K[4][2] = 1.0;
   K->K[5][2] = -1.0;
   }
