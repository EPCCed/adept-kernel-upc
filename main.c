/* Copyright (c) 2015 The University of Edinburgh. */

/* 
* This software was developed as part of the                       
* EC FP7 funded project Adept (Project ID: 610490)                 
* www.adept-project.eu                                            
*/

/* Licensed under the Apache License, Version 2.0 (the "License"); */
/* you may not use this file except in compliance with the License. */
/* You may obtain a copy of the License at */

/*     http://www.apache.org/licenses/LICENSE-2.0 */

/* Unless required by applicable law or agreed to in writing, software */
/* distributed under the License is distributed on an "AS IS" BASIS, */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/* See the License for the specific language governing permissions and */
/* limitations under the License. */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <limits.h>
#include <upc.h>

#include "level1.h"

void usage();

void settable(UL i1,UL i2,UL i3,UL i4,UL i5, UL i6){
  int i; z=i1;w=i2,jsr=i3; jcong=i4; a=i5; b=i6;
  for(i=0;i<256;i=i+1) t[i]=KISS;
}


int main(int argc, char **argv){
  
  int c;

  /* Init RNG here */
  settable(12345,65435,34221,12345,9983651,95746118);

    char *bench = "blas_op";
    unsigned int size = 200;
    unsigned long rep = ULONG_MAX;
    char *op  = "dot_product";
    char *dt = "int";
    
    static struct option option_list[] =
    { {"bench", required_argument, NULL, 'b'},
      {"size", required_argument, NULL, 's'},
      {"reps", required_argument, NULL, 'r'},
      {"op", required_argument, NULL, 'o'},
      {"dtype", required_argument, NULL, 'd'},
      {"help", no_argument, NULL, 'h'},
      {0, 0, 0, 0}
    };

    if (MYTHREAD == 0) printf("Executing benchmark on %d UPC threads.\n", THREADS);
    
    while((c = getopt_long(argc, argv, "b:s:r:o:d:h", option_list, NULL)) != -1){
      switch(c){
        case 'b':
          bench = optarg;
          if (MYTHREAD==0) printf("Benchmark is %s.\n", bench);
          break;
        case 's':
          size = atoi(optarg);
          if (MYTHREAD==0) printf("Size is %d.\n", size);
          break;
        case 'r':
          rep = atol(optarg);
          if (MYTHREAD==0) printf("Number of repetitions %lu.\n", rep);
          break;
        case 'o':
          op = optarg;
          if (MYTHREAD==0) printf("Operation %s\n", op);
          break;
        case 'd':
          dt = optarg;
          if (MYTHREAD==0) printf("Data type is %s\n", dt);
          break;
        case 'h':
          if (MYTHREAD==0) usage();
          return 0;
        default:
          if (MYTHREAD==0) printf("Undefined.\n");
          return 0;
      }
    }
    
    bench_level1(bench, size, rep, op, dt);
    
  return 0;
  
}

void usage(){
  printf("Usage for UPC KERNEL benchmarks:\n\n");
  printf("\t -b, --bench NAME \t name of the benchmark - possible values are blas_op and stencil.\n");
  printf("\t -s, --size N \t\t vector length. Default is 200.\n");
  printf("\t -r, --reps N \t\t number of repetitions. Default value is ULONG_MAX.\n");
  printf("\t -o, --op TYPE \t\t TYPE of operation.\n");
  printf("\t\t\t\t --> for blas_op benchmark: \"dot_product\", \"scalar_mult\", \"dmatvec_product\", \"norm\", \"axpy\", \"spmv\" and \"spgemm\". Default is \"dot_product\".\n");
  printf("\t\t\t\t --> for stencil benchmark: \"27\", \"19\", \"9\" and \"5\". Default is \"27\".\n");
  printf("\t -d, --dtype DATATYPE \t DATATYPE to be used - possible values are int, long, float, double. Default is int.\n");
  printf("\t -h, --help \t\t Displays this help.\n");
  printf("\n\n");
}
