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

/*
 * UPC Stencil benchmark
 *
 * Based on that of Alan Gray, EPCC
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include <upc.h>

#include "level1.h"
#include "utils.h"

#define REPS 100

void stencil27(unsigned int size){

  int i, j, k, iter;
  int n = size-2;
  int size3d=size*size*size;
  int size2d=size*size;
  double fac = 1.0/26;
  
  /* Work buffers, with halos */
  shared double *a0 = (shared double*)upc_all_alloc(THREADS, sizeof(double)*size3d);
  shared double *a1 = (shared double*)upc_all_alloc(THREADS, sizeof(double)*size3d);
  
  if(a0==NULL||a1==NULL){
    /* Something went wrong in the memory allocation here, fail gracefully */
    if (MYTHREAD == 0) printf("27-point Stencil Error: Unable to allocate memory\n");
  }
  
  if (MYTHREAD == 0){
    printf("Running with %d UPC thread(s):\n\n",THREADS);
  }
  
  struct timespec start, end;
  
  /* zero all of array (including halos) */
  upc_forall (i = 0; i < size3d; i++; &a0[i]) {
    a0[i] = 0.0;
  }
  
  /* use random numbers to fill interior */
  for (i = 1; i < n+1; i++) {
    for (j = 1; j < n+1; j++) {
      for (k = 1; k < n+1; k++) {
	a0[i*size2d+j*size+k] = (double) rand()/ (double)(1.0 + RAND_MAX);
      }
    }
  }
  
  /* run main computation */
  clock_gettime(CLOCK, &start);

  for (iter = 0; iter < REPS; iter++) {
    
    upc_forall (i = 1; i < n+1; i++; i) {
      for (j = 1; j < n+1; j++) {
	for (k = 1; k < n+1; k++) {
	  a1[i*size2d+j*size+k] = (
				   a0[i*size2d+(j-1)*size+k] + a0[i*size2d+(j+1)*size+k] +
				   a0[(i-1)*size2d+j*size+k] + a0[(i+1)*size2d+j*size+k] +
				   a0[(i-1)*size2d+(j-1)*size+k] + a0[(i-1)*size2d+(j+1)*size+k] +
				   a0[(i+1)*size2d+(j-1)*size+k] + a0[(i+1)*size2d+(j+1)*size+k] +
				      
				   a0[i*size2d+(j-1)*size+(k-1)] + a0[i*size2d+(j+1)*size+(k-1)] +
				   a0[(i-1)*size2d+j*size+(k-1)] + a0[(i+1)*size2d+j*size+(k-1)] +
				   a0[(i-1)*size2d+(j-1)*size+(k-1)] + a0[(i-1)*size2d+(j+1)*size+(k-1)] +
				   a0[(i+1)*size2d+(j-1)*size+(k-1)] + a0[(i+1)*size2d+(j+1)*size+(k-1)] +
				      
				   a0[i*size2d+(j-1)*size+(k+1)] + a0[i*size2d+(j+1)*size+(k+1)] +
				   a0[(i-1)*size2d+j*size+(k+1)] + a0[(i+1)*size2d+j*size+(k+1)] +
				   a0[(i-1)*size2d+(j-1)*size+(k+1)] + a0[(i-1)*size2d+(j+1)*size+(k+1)] +
				   a0[(i+1)*size2d+(j-1)*size+(k+1)] + a0[(i+1)*size2d+(j+1)*size+(k+1)] +
				      
				   a0[i*size2d+j*size+(k-1)] + a0[i*size2d+j*size+(k+1)]
				   ) * fac;
        }
      }
    }
    
    upc_forall (i = 1; i < n+1; i++; i) {
      for (j = 1; j < n+1; j++) {
	for (k = 1; k < n+1; k++) {
	  a0[i*size2d+j*size+k] = a1[i*size2d+j*size+k];
	}
      }
    }
    
    upc_barrier;

  } /* end iteration loop */
  clock_gettime(CLOCK, &end);

  if (MYTHREAD == 0){
    elapsed_time_hr(start, end, "Stencil - 27 point");
    /* Free malloc'd memory to prevent leaks */
    upc_free(a0);
    upc_free(a1);
  }
    
}

void stencil19(unsigned int size){

  int i, j, k, iter;
  int n = size-2;
  int size3d = size*size*size;
  int size2d = size*size;
  double fac = 1.0/18;
  
  /* Work buffers, with halos */
  shared double *a0 = (shared double*)upc_all_alloc(THREADS, sizeof(double)*size3d);
  shared double *a1 = (shared double*)upc_all_alloc(THREADS, sizeof(double)*size3d);

  if(a0==NULL||a1==NULL){
    /* Something went wrong in the memory allocation here, fail gracefully */
    if (MYTHREAD == 0)printf("19-point Stencil Error: Unable to allocate memory\n");
  }

  if (MYTHREAD == 0){
    printf("Running with %d UPC thread(s):\n\n",THREADS);
  }

  struct timespec start,end;
  /* zero all of array (including halos) */

  upc_forall (i = 0; i < size3d; i++; &a0[i]) {
    a0[i] = 0.0;
  }

  /* use random numbers to fill interior */
  upc_forall (i = 1; i < n+1; i++; &a0[i]) {
    for (j = 1; j < n+1; j++) {
      for (k = 1; k < n+1; k++) {
	a0[i*size2d+j*size+k] = (double) rand()/ (double)(1.0 + RAND_MAX);
      }
    }
  }
  
  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(n+1);
  upc_barrier_timer();

  upc_barrier;

  /* run main computation */

  clock_gettime(CLOCK, &start);
  for (iter = 0; iter < REPS; iter++) {

    upc_forall (i = 1; i < n+1; i++; i) {
      for (j = 1; j < n+1; j++) {
	for (k = 1; k < n+1; k++) {
	  a1[i*size2d+j*size+k] = (
				   a0[i*size2d+(j-1)*size+k] + a0[i*size2d+(j+1)*size+k] +
				   a0[(i-1)*size2d+j*size+k] + a0[(i+1)*size2d+j*size+k] +
				   a0[(i-1)*size2d+(j-1)*size+k] + a0[(i-1)*size2d+(j+1)*size+k] +
				   a0[(i+1)*size2d+(j-1)*size+k] + a0[(i+1)*size2d+(j+1)*size+k] +
					
				   a0[i*size2d+(j-1)*size+(k-1)] + a0[i*size2d+(j+1)*size+(k-1)] +
				   a0[(i-1)*size2d+j*size+(k-1)] + a0[(i+1)*size2d+j*size+(k-1)] +
					
				   a0[i*size2d+(j-1)*size+(k+1)] + a0[i*size2d+(j+1)*size+(k+1)] +
				   a0[(i-1)*size2d+j*size+(k+1)] + a0[(i+1)*size2d+j*size+(k+1)] +
					
				   a0[i*size2d+j*size+(k-1)] + a0[i*size2d+j*size+(k+1)]
				   ) * fac;
	}
      }
    }
      
    upc_forall (i = 1; i < n+1; i++; i) {
      for (j = 1; j < n+1; j++) {
	for (k = 1; k < n+1; k++) {
	  a0[i*size2d+j*size+k] = a1[i*size2d+j*size+k];
	}
      }
    }

    upc_barrier;
      
  } /* end iteration loop */

  clock_gettime(CLOCK, &end);
  
  if(MYTHREAD==0){
    elapsed_time_hr(start, end, "Stencil - 19 point");
  
    /* Free malloc'd memory to prevent leaks */
    upc_free(a0);
    upc_free(a1);
  }
}

void stencil9(unsigned int size){

  int i, j, iter;
  int n = size-2;
  int array_size = size*size;
  double fac = 1.0/8;
  
  /* Work buffers, with halos */
  shared double *a0 = (shared double*)upc_all_alloc(THREADS, sizeof(double)*array_size);
  shared double *a1 = (shared double*)upc_all_alloc(THREADS, sizeof(double)*array_size);
  
  if(a0==NULL||a1==NULL){
    /* Something went wrong in the memory allocation here, fail gracefully */
    if(MYTHREAD==0) printf("9-point Stencil Error: Unable to allocate memory\n");
  }

  if (MYTHREAD == 0){
    printf("Running with %d UPC thread(s):\n\n",THREADS);
  }

  struct timespec start,end;
  /* zero all of array (including halos) */
  upc_forall (i = 0; i < array_size; i++; &a0[i]) {
    a0[i] = 0.0;
  }
  
  /* use random numbers to fill interior */
  upc_forall (i = 1; i < n+1; i++; &a0[i]) {
    for (j = 1; j < n+1; j++) {
      a0[i*size+j] = (double) rand()/ (double)(1.0 + RAND_MAX);
    }
  }
  
  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(n+1);
  upc_barrier_timer();

  upc_barrier;

  /* run main computation */
  clock_gettime(CLOCK, &start);
  for (iter = 0; iter < REPS; iter++) {

    upc_forall (i = 1; i < n+1; i++; &a0[i*size+i]) {
      for (j = 1; j < n+1; j++) {
	a1[i*size+j] = (a0[i*size+(j-1)] + a0[i*size+(j+1)] +
			a0[(i-1)*size+j] + a0[(i+1)*size+j] +
			a0[(i-1)*size+(j-1)] + a0[(i-1)*size+(j+1)] +
			a0[(i+1)*size+(j-1)] + a0[(i+1)*size+(j+1)]) * fac;
      }
    }
    upc_forall (i = 1; i < n+1; i++; &a0[i*size+i]) {
      for (j = 1; j < n+1; j++) {
	a0[i*size+j] = a1[i*size+j];
      }
    }

    upc_barrier;

  } /* end iteration loop */

  clock_gettime(CLOCK, &end);
  
  if(MYTHREAD == 0){
    elapsed_time_hr(start, end, "Stencil - 9 point");

    /* Free malloc'd memory to prevent leaks */
    upc_free(a0);
    upc_free(a1);
  }

}


void stencil5(unsigned int size){

  int i, j, iter;
  int n = size-2;
  int array_size = size*size;
  double fac = 1.0/8;
  
  /* Work buffers */
  shared double *a0 = (shared double*)upc_all_alloc(THREADS, sizeof(double)*array_size);
  shared double *a1 = (shared double*)upc_all_alloc(THREADS, sizeof(double)*array_size);
  
  if(a0==NULL||a1==NULL){
    /* Something went wrong in the memory allocation here, fail gracefully */
    if(MYTHREAD==0) printf("5-point Stencil Error: Unable to allocate memory\n");
  }
  
  if (MYTHREAD == 0){
    printf("Running with %d UPC thread(s):\n\n",THREADS);
  }
  
  struct timespec start,end;
  
  /* zero all of array (including halos) */
  upc_forall (i = 0; i < array_size; i++; &a0[i]) {
    a0[i] = 0.0;
  }
  
  /* use random numbers to fill interior */
  upc_forall (i = 1; i < n+1; i++; &a0[i]) {
    for (j = 1; j < n+1; j++){
      a0[i*size+j] = (double) rand()/ (double)(1.0 + RAND_MAX);
    }
  }
  
  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(n+1);
  upc_barrier_timer();

  upc_barrier;

  /* run main computation */
  clock_gettime(CLOCK, &start);
  for (iter = 0; iter < REPS; iter++) {
    
    upc_forall (i = 1; i < n+1; i++; &a1[i*size+i]) {
      for (j = 1; j < n+1; j++){
	a1[i*size+j] = (a0[i*size+(j-1)] + a0[i*size+(j+1)] 
			+ a0[(i-1)*size+j] + a0[(i+1)*size+j]) * fac;
      }
    }
 
    upc_forall(i = 1; i < n+1; i++; &a0[i*size+i]) {
      for (j = 1; j < n+1; j++) {
	a0[i*size+j] = a1[i*size+j];
      }
    }
    
    upc_barrier;
  }
  
  clock_gettime(CLOCK, &end);

  if(MYTHREAD == 0){
    elapsed_time_hr(start, end, "Stencil - 5 point");
    
    /* Free malloc'd memory to prevent leaks */
    upc_free(a0);
    upc_free(a1);
  }
  
}
