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

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include <upc.h>

#include "level1.h"
#include "utils.h"
#include "matrix_utils.h"


/*
 * Vector dot product, integers
 *
 * result = result + v1_i * v2_i
 *
 * Input: size of the vectors (in number of elements)
 * Output: dot product
 *
 */
int int_dot_product(unsigned int size){

  int i;
  int local_size = size / THREADS;

  /* create two vectors */
  shared int *v1 = (shared int *)upc_all_alloc(THREADS, local_size * sizeof(int));
  shared int *v2 = (shared int *)upc_all_alloc(THREADS, local_size * sizeof(int));

  /* local result array and final result variable */
  shared unsigned int *tmp_result = (shared unsigned int *)upc_all_alloc(THREADS, sizeof(int));
  static shared unsigned int result = 0;

  if(v1 == NULL || v2 == NULL){
    if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the two arrays.\n");
    return 0;
  }

  struct timespec start, end;

  tmp_result[MYTHREAD] = 0;

  /* fill vectors with random integer values */
  if(MYTHREAD == 0){
    for(i=0; i<local_size*THREADS; i++){
      v1[i] = KISS;
      v2[i] = KISS;
    }
  }
 
  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size*THREADS);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  /* perform dot product */

  upc_forall(i=0; i<local_size*THREADS; i++; &v1[i]){
    tmp_result[MYTHREAD] = tmp_result[MYTHREAD] + v1[i] * v2[i];
  }

  upc_barrier;

  if (MYTHREAD==0){
    for(i=0;i<THREADS;i++){
      result = result + tmp_result[i];
    }
  }

  upc_barrier;
  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){

    /* print result so compiler does not throw it away */
    printf("Dot product result: %d\n", result);

    elapsed_time_hr(start, end, "Integer dot product.");

    upc_free(v1);
    upc_free(v2);
  }

  return 0;

}

/*
 * Vector dot product, floats
 *
 * result = result + v1_i * v2_i
 *
 * Input: size of the vectors (in number of elements)
 * Output: dot product
 *
 */
int float_dot_product(unsigned int size){

  int i;
  int local_size = size / THREADS;

  /* create three vectors */
  shared float *v1 = (shared float *)upc_all_alloc(THREADS, local_size * sizeof(float));
  shared float *v2 = (shared float *)upc_all_alloc(THREADS, local_size * sizeof(float));

  /* local result array and final result variable */
  shared float *tmp_result = (shared float *)upc_all_alloc(THREADS, sizeof(float));
  static shared float result = 0.0;

  if(v1 == NULL || v2 == NULL){
	  if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the two arrays.\n");
    return 0;
  }

  struct timespec start, end;

  tmp_result[MYTHREAD] = 0.0;

  /* fill vectors with random floats */
  if(MYTHREAD == 0){
    for(i=0; i<local_size*THREADS; i++){
      v1[i] = UNI;
      v2[i] = UNI;
    }
  }

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size*THREADS);
  upc_barrier_timer();
  upc_barrier;

  clock_gettime(CLOCK, &start);

  /* perform dot product */

  upc_forall(i=0; i<local_size*THREADS; i++; &v1[i]){
    tmp_result[MYTHREAD] = tmp_result[MYTHREAD] + v1[i] * v2[i];
  }

  upc_barrier;
  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
	/* print result so compiler does not throw it away */
	printf("Dot product result: %f\n", result);

	elapsed_time_hr(start, end, "Float dot product.");

	upc_free(v1);
	upc_free(v2);
  }
  return 0;

}


/*
 * Vector dot product, doubles
 *
 * result = result + v1_i * v2_i
 *
 * Input: size of the vectors (in number of elements)
 * Output: dot product
 *
 */
int double_dot_product(unsigned int size){

  int i;
  int local_size = size / THREADS;

  /* create three vectors */
  shared double *v1 = (shared double *)upc_all_alloc(THREADS, local_size * sizeof(double));
  shared double *v2 = (shared double *)upc_all_alloc(THREADS, local_size * sizeof(double));

  /* local result array and final result variable */
  shared double *tmp_result = (shared double *)upc_all_alloc(THREADS, sizeof(double));
  static shared double result = 0.0;

  if(v1 == NULL || v2 == NULL){
	if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the two arrays.\n");
    return 0;
  }

  struct timespec start, end;

  tmp_result[MYTHREAD] = 0.0;

  /* fill vectors with random doubles */
  if(MYTHREAD == 0){
    for(i=0; i<local_size*THREADS; i++){
      v1[i] = VNI;
      v2[i] = VNI;
    }
  }
  
  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size*THREADS);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  /* perform dot product */

  upc_forall(i=0; i<local_size*THREADS; i++; &v1[i]){
    tmp_result[MYTHREAD] = tmp_result[MYTHREAD] + v1[i] * v2[i];
  }

  upc_barrier;
  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
    /* print result so compiler does not throw it away */
	printf("Dot product result: %f\n", result);

    elapsed_time_hr(start, end, "Double dot product.");

    upc_free(v1);
    upc_free(v2);
  }
  return 0;

}


/* Vector scalar multiplication, integers    */
/* v_i = a * v1_i                     */
int int_scalar_mult(unsigned int size){

  int i;
  int local_size = size/THREADS;

  /* create vector and scalar */
  shared int *v = (shared int *)upc_all_alloc(THREADS, local_size * sizeof(int));
  static shared unsigned int a = 0;

  if(v == NULL){
    if (MYTHREAD==0) printf("Out Of Memory: could not allocate space for the array.\n");
    return 0;
  }

  struct timespec start, end;

  /* fill vector with random ints */
  if(MYTHREAD == 0){
    for(i=0; i<local_size*THREADS; i++){
      v[i] = KISS;
    }

    /* assign random int value */
    a = KISS;
  }
  
  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  /* perform scalar product */
  upc_forall(i=0; i<local_size*THREADS; i++; &v[i]){
    v[i] = a * v[i];
  }
  
  upc_barrier;

  clock_gettime(CLOCK, &end);

  if(MYTHREAD == 0){
    /* print result so compiler does not throw it away */
    printf("Scalar product result: %d\n", v[0]);
    
    elapsed_time_hr(start, end, "Int scalar multiplication.");
    
    upc_free(v);
  }

  return 0;

}

/* Vector scalar multiplication, floats    */
/* v_i = a * v1_i                     */
int float_scalar_mult(unsigned int size){

  int i;
  int local_size = size/THREADS;

  /* create vector and scalar */
  shared float *v = (shared float *)upc_all_alloc(THREADS, local_size * sizeof(float));
  static shared float a = 0.0;

  if(v == NULL){
    if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the array.\n");
    return 0;
  }

  struct timespec start, end;

  /* fill vector with random floats */
  if (MYTHREAD == 0){
    for(i=0; i<size; i++){
      v[i] = UNI;
    }  

    /* assign random float value */
    a = UNI;

  }

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  /* perfom scalar product */
  upc_forall(i=0; i<local_size*THREADS; i++; &v[i]){
    v[i] = a * v[i];
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD == 0){
    /* print result so compiler does not throw it away */
    printf("Scalar product result: %f\n", v[0]);
    
    elapsed_time_hr(start, end, "Float scalar multiplication.");
    
    upc_free(v);
  }

  return 0;

}

/* Vector scalar multiplication, doubles    */
/* v_i = a * v1_i                     */
int double_scalar_mult(unsigned int size){

  int i;
  int local_size = size / THREADS;

  /* create vector and scalar */
  shared double *v = (shared double *)upc_all_alloc(THREADS, local_size * sizeof(double));
  static shared double a = 0.0;

  if(v == NULL){
    if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the array.\n");
    return 0;
  }

  struct timespec start, end;

  if (MYTHREAD==0){
    /* fill vector with random doubles */
    for(i=0; i<local_size*THREADS; i++){
      v[i] = VNI;
    }
    
    /* assign random double value */
    a = VNI;
  }
  
  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  /* perform scalar product */
  upc_forall(i=0; i<local_size*THREADS; i++; &v[i]){
    v[i] = a * v[i];
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
    /* print result so compiler does not throw it away */
    printf("Scalar product result: %f\n", v[0]);
    
    elapsed_time_hr(start, end, "Double scalar multiplication.");
    
    upc_free(v);
  }
  return 0;

}

/*
 * compute the Euclidean norm of an int vector
 */
int int_norm(unsigned int size){

  int i;
  int local_size = size/THREADS;

  shared int *v = (shared int *)upc_all_alloc(THREADS, local_size * sizeof(int));
  shared int *part_sum = (shared int *)upc_all_alloc(THREADS, sizeof(int));
  static shared int sum = 0, norm = 0;

  if(v == NULL){
    if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the array.\n");
    return 0;
  }

  struct timespec start, end;

  if (MYTHREAD == 0){
    /* fill vector with random floats */
    for(i=0; i<size; i++){
      v[i] = UNI;
    }
  }

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  part_sum[MYTHREAD] = 0;

  upc_barrier;

  clock_gettime(CLOCK, &start);

  upc_forall (i=0; i<size; i++; &v[i]){
    part_sum[MYTHREAD] = part_sum[MYTHREAD] + (v[i]*v[i]);
  }

  upc_barrier;

  if(MYTHREAD==0){
    for(i=0;i<THREADS;i++){
      sum = sum + part_sum[i];
    }
    norm = sqrt(sum);
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
    elapsed_time_hr(start, end, "Integer vector norm.");

    /* print result so compiler does not throw it away */
    printf("Norm = %d\n", norm);

    upc_free(v);
  }

  return 0;
}

/*
 * compute the Euclidean norm of a float vector
 */
int float_norm(unsigned int size){

  int i;
  int local_size = size/THREADS;

  shared float *v = (shared float *)upc_all_alloc(THREADS, local_size * sizeof(float));
  shared float *part_sum = (shared float *)upc_all_alloc(THREADS, sizeof(float));
  static shared float sum = 0.0, norm = 0.0;

  if(v == NULL){
    if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the array.\n");
    return 0;
  }

  struct timespec start, end;

  if (MYTHREAD == 0){
    /* fill vector with random floats */
    for(i=0; i<size; i++){
      v[i] = UNI;
    }
  }

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  part_sum[MYTHREAD] = 0.0;

  upc_barrier;

  clock_gettime(CLOCK, &start);

  upc_forall (i=0; i<size; i++; &v[i]){
    part_sum[MYTHREAD] = part_sum[MYTHREAD] + (v[i]*v[i]);
  }

  upc_barrier;

  if(MYTHREAD==0){
    for(i=0;i<THREADS;i++){
      sum = sum + part_sum[i];
    }
    norm = sqrt(sum);
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
    elapsed_time_hr(start, end, "Float vector norm.");
    
    /* print result so compiler does not throw it away */
    printf("Norm = %f\n", norm);
    
    upc_free(v);
  }

  return 0;
}

/*
 * compute the Euclidean norm of a float vector
 */
int double_norm(unsigned int size){

  int i;
  int local_size = size/THREADS;

  shared double *v = (shared double *)upc_all_alloc(THREADS, local_size * sizeof(double));
  shared double *part_sum = (shared double *)upc_all_alloc(THREADS, sizeof(double));
  static shared double sum = 0.0, norm = 0.0;

  if(v == NULL){
    if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the array.\n");
    return 0;
  }

  struct timespec start, end;

  if (MYTHREAD == 0){
    /* fill vector with random floats */
    for(i=0; i<size; i++){
      v[i] = UNI;
    }
  }

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  part_sum[MYTHREAD] = 0.0;

  upc_barrier;

  clock_gettime(CLOCK, &start);

  upc_forall (i=0; i<size; i++; &v[i]){
    part_sum[MYTHREAD] = part_sum[MYTHREAD] + (v[i]*v[i]);
  }

  upc_barrier;

  if(MYTHREAD==0){
    for(i=0;i<THREADS;i++){
      sum = sum + part_sum[i];
    }
    norm = sqrt(sum);
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
    elapsed_time_hr(start, end, "Double vector norm.");

    /* print result so compiler does not throw it away */
    printf("Norm = %f\n", norm);

    upc_free(v);
  }

  return 0;
}

/*
 *
 * Compute vector-scalar product
 * AXPY, integers
 *
 * y = a * x + y
 *
 * Naive implementation
 *
 */
int int_axpy(unsigned int size){

  int i;
  int local_size = size / THREADS;

  static shared int a;

  shared int *x = (shared int *)upc_all_alloc(THREADS, local_size * sizeof(int));
  shared int *y = (shared int *)upc_all_alloc(THREADS, local_size * sizeof(int));

  if(x == NULL || y == NULL){
    if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the two arrays.\n");
    return 0;
  }

  if(MYTHREAD==0){

    a = KISS;
    
    /* fill x and y vectors with random ints */
    for(i=0; i<size; i++){
      x[i] = KISS;
      y[i] = KISS;
    }
  }

  struct timespec start, end;

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  upc_forall(i=0; i<size; i++; &y[i]){
    y[i] = a * x[i] + y[i];
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD == 0){

    elapsed_time_hr(start, end, "Int AXPY.");

    /* print some of the result so compiler does not throw it away */
    printf("APXY result = %d\n", y[0]);

    upc_free(x);
    upc_free(y);

  }

  return 0;
}

/*
 *
 * Compute vector-scalar product
 * AXPY, floats
 *
 * y = a * x + y
 *
 * Naive implementation
 *
 */
int float_axpy(unsigned int size){

  int i;
  int local_size = size / THREADS;

  static shared float a;

  shared float *x = (shared float *)upc_all_alloc(THREADS, local_size * sizeof(float));
  shared float *y = (shared float *)upc_all_alloc(THREADS, local_size * sizeof(float));

  if(x == NULL || y == NULL){
	if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the two arrays.\n");
    return 0;
  }

  if(MYTHREAD==0){

	a = UNI;

    /* fill x and y vectors with random ints */
    for(i=0; i<size; i++){
      x[i] = UNI;
      y[i] = UNI;
    }

  }

  struct timespec start, end;

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  upc_forall(i=0; i<local_size*THREADS; i++; &y[i]){
    y[i] = a * x[i] + y[i];
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
	elapsed_time_hr(start, end, "Float AXPY.");

    /* print some of the result so compiler does not throw it away */
    printf("APXY result = %f\n", y[0]);

    upc_free(x);
    upc_free(y);
  }

  return 0;
}

/*
 *
 * Compute vector-scalar product
 * AXPY, doubles
 *
 * y = a * x + y
 *
 * Naive implementation
 *
 */
int double_axpy(unsigned int size){

  int i;
  int local_size = size / THREADS;

  static shared double a;

  shared double *x = (shared double *)upc_all_alloc(THREADS, local_size * sizeof(double));
  shared double *y = (shared double *)upc_all_alloc(THREADS, local_size * sizeof(double));

  if(x == NULL || y == NULL){
	if (MYTHREAD == 0) printf("Out Of Memory: could not allocate space for the two arrays.\n");
    return 0;
  }

  if(MYTHREAD==0){

    a = VNI;

    /* fill x and y vectors with random ints */
    for(i=0; i<size; i++){
      x[i] = VNI;
      y[i] = VNI;
    }

  }

  struct timespec start, end;

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  upc_forall(i=0; i<local_size*THREADS; i++; &y[i]){
    y[i] = a * x[i] + y[i];
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
	elapsed_time_hr(start, end, "Double AXPY.");

    /* print some of the result so compiler does not throw it away */
    printf("APXY result = %f\n", y[0]);

    upc_free(x);
    upc_free(y);
  }
  return 0;
}

/*
 * Dense Matrix-Vector product, integers
 *
 * y = A * x
 * where A is a square matrix
 *
 * Input:  number of elements in vectors and of rows/cols
 *         in matrix specified as number of ints
 *
 */
int int_dmatvec_product(unsigned int size){

  int i,j;
  int r1,r2;
  int local_size = size/THREADS;
  int mat_size = pow(local_size*THREADS, 2);
  printf("local size: %ld, matrix size: %ld\n", local_size, mat_size);

  /* create two vectors */
  shared int *x = (shared int *)upc_all_alloc(THREADS, local_size * sizeof(int));
  shared int *y = (shared int *)upc_all_alloc(THREADS, local_size * sizeof(int));

  /* create matrix - represented as a 1D array */
  shared int *A = (shared int *)upc_all_alloc(THREADS, (local_size*local_size*THREADS) * sizeof(int));

  //  printf("size of A = %d\n", sizeof(A));

  if(x == NULL || y == NULL || A == NULL){
    if (MYTHREAD==0) printf("Out Of Memory: could not allocate space for the vectors and matrix.\n");
    return 0;
  }

  struct timespec start, end;

  r1 = KISS;
  r2 = KISS;
  
  if(MYTHREAD==0){
    /* fill vector x and matrix A with random integer values */
    for(i=0; i<local_size*THREADS; i++){
      x[i] = r1;
    }
    for(i=0; i<mat_size; i++){
      A[i] = r2;
      /* if(i>=536872378){ */
      /* 	printf("size of A[%d] = %d\n", i, sizeof(A[i])); */
      /* 	printf("thread affinity A[%d]= %d\n", i, upc_threadof(&A[i])); */
      /* } */
    }
  }

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  /* perform matrix-vector product */
  upc_forall(i=0; i<local_size*THREADS; i++; &y[i]){
    for(j=0; j<local_size*THREADS; j++){
      y[i] = y[i] + A[(i*local_size*THREADS)+j] * x[j];
    }
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
    elapsed_time_hr(start, end, "Integer Dense Matrix-Vector product.");

    /* print result so compiler does not throw it away */
    printf("Result vector y[0] = %d\n", y[0]);
    
    upc_free(x);
    upc_free(y);
    upc_free(A);
  }
  return 0;

}

/*
 * Dense Matrix-Vector product, floats
 *
 * y = A * x
 * where A is a square matrix
 *
 * Input:  number of elements in vectors and of rows/cols
 *         in matrix specified as number of floats
 *
 */
int float_dmatvec_product(unsigned int size){

  int i,j;
  float r1,r2;
  int local_size = size/THREADS;
  int mat_size = pow(local_size*THREADS, 2);

  /* create two vectors */
  shared float *x = (shared float *)upc_all_alloc(THREADS, local_size * sizeof(float));
  shared float *y = (shared float *)upc_all_alloc(THREADS, local_size * sizeof(float));

  /* create matrix - represented as a 1D array */
  shared float *A = (shared float *) upc_all_alloc(THREADS*local_size, THREADS*local_size * sizeof(float));

  if(x == NULL || y == NULL || A == NULL){
	if (MYTHREAD==0) printf("Out Of Memory: could not allocate space for the vectors and matrix.\n");
    return 0;
  }

  struct timespec start, end;

  r1 = UNI;
  r2 = UNI;

  if(MYTHREAD==0){
    /* fill vector x and matrix A with random integer values */
    for(i=0; i<local_size*THREADS; i++){
      x[i] = r1;
    }
    for(i=0; i<mat_size; i++){
      A[i] = r2;
    }
  }

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  /* perform matrix-vector product */
  upc_forall(i=0; i<local_size*THREADS; i++; &y[i]){
    for(j=0; j<local_size*THREADS; j++){
      y[i] = y[i] + A[(i*local_size*THREADS)+j] * x[j];
    }
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
	elapsed_time_hr(start, end, "Float Dense Matrix-Vector product.");

    /* print result so compiler does not throw it away */
    printf("Result vector y[0] = %f\n", y[0]);

    upc_free(x);
    upc_free(y);
    upc_free(A);

  }

  return 0;

}

/*
 * Dense Matrix-Vector product, doubles
 *
 * y = A * x
 * where A is a square matrix
 *
 * Input:  number of elements in vectors and of rows/cols
 *         in matrix specified as number of floats
 *
 */
int double_dmatvec_product(unsigned int size){

  int i,j;
  double r1,r2;
  int local_size = size/THREADS;
  int mat_size = pow(local_size*THREADS, 2);

  /* create two vectors */
  shared double *x = (shared double *)upc_all_alloc(THREADS, local_size * sizeof(double));
  shared double *y = (shared double *)upc_all_alloc(THREADS, local_size * sizeof(double));

  /* create matrix - represented as a 1D array */
  shared double *A = (shared double *) upc_all_alloc(THREADS*local_size, THREADS*local_size * sizeof(double));


  if(x == NULL || y == NULL || A == NULL){
    printf("Out Of Memory: could not allocate space for the vectors or the matrix.\n");
    return 0;
  }


  struct timespec start, end;

  r1 = VNI;
  r2 = VNI;

  if(MYTHREAD==0){
    /* fill vector x and matrix A with random integer values */
    for(i=0; i<local_size*THREADS; i++){
      x[i] = r1;
    }
    for(i=0; i<mat_size; i++){
      A[i] = r2;
    }
  }

  /* measuring upc_forall and upc_barrier overheads */
  upc_loop_timer_nop(local_size);
  upc_barrier_timer();

  upc_barrier;

  clock_gettime(CLOCK, &start);

  /* perform matrix-vector product */
  upc_forall(i=0; i<local_size*THREADS; i++; &y[i]){
    for(j=0; j<local_size*THREADS; j++){
      y[i] = y[i] + A[(i*local_size*THREADS)+j] * x[j];
    }
  }

  upc_barrier;

  clock_gettime(CLOCK, &end);

  if (MYTHREAD==0){
	elapsed_time_hr(start, end, "Dense Matrix-Vector product.");

    /* print result so compiler does not throw it away */
    printf("Result vector y[0] = %f\n", y[0]);

    upc_free(x);
    upc_free(y);
    upc_free(A);
  }

  return 0;

}

int float_spmatvec_product(unsigned long r){

  int m, n, nz;
  int *row_idx, *col_idx;
  float *values;
  float *x,*b;

  struct timespec start,end;

  char *filename = "matrix_in.csr";

  int i,j,rep;

  FILE *f;
  char line[64];

  if(r==ULONG_MAX) r=1000;

  if ((f = fopen(filename, "r")) == NULL) {
    printf ("can't open file <%s> \n", filename);
    exit(1);
  }

  fgets(line, sizeof(line), f);
  sscanf(line, "%d %d %d", &nz, &n, &m);

  printf("Number of elements of values and col_idx: %d; number of values in row_idx: %d\n", nz, m);

  row_idx = malloc(m * sizeof(int));
  col_idx = malloc(n * sizeof(int));
  values = malloc(nz * sizeof(float));
  x = malloc((m-1) * sizeof(float));
  b = malloc((m-1) * sizeof(float));

  if (!row_idx || !col_idx || !values || !x || !b){
    printf ("cannot allocate memory for sparse matrix and vectors\n");
    exit(1);
  }
  else{
    printf("memory allocated\n");
  }

  // fill values
  for(i=0;i<nz;i++){
    fgets(line, sizeof(line), f);
    sscanf(line, "%f", &values[i]);
  }

  // fill col_idx
  for(i=0;i<nz;i++){
    fgets(line, sizeof(line), f);
    sscanf(line, "%d", &col_idx[i]);
  }

  // fill row_idx
  for(i=0;i<m;i++){
    fgets(line, sizeof(line), f);
    sscanf(line, "%d", &row_idx[i]);
  }

  for (i=0; i<m-1; i++)  {
    x[i]=i+1.5; // give basic values to vector x
  }

  clock_gettime(CLOCK, &start);

  for(rep=0;rep<r;rep++){
    /* Ax=b */
    for(i=0; i<m-1; i++){
      for(j=row_idx[i];j<row_idx[i+1];j++){
        b[i] = b[i] + values[j]* x[col_idx[j]];
      }
    }
  }

  clock_gettime(CLOCK, &end);

  elapsed_time_hr(start, end, "Sparse DMVs.");

  free(x);
  free(b);
  free(row_idx);
  free(col_idx);
  free(values);

  fclose(f);

  return 0;
}


int double_spmatvec_product(unsigned long r){


  /* static shared int m; */
  /* static shared int n; */
  /* static shared int nz; */

  shared [] int *row_idx;
  shared [] int *col_idx;
  shared [] double *values;
  shared [] double *x;
  shared [] double *b;

  shared [] int *m_s;
  shared [] int *n_s;
  shared [] int *nz_s;

  m_s = (shared [] int *)upc_all_alloc(1, sizeof(int));
  n_s = (shared [] int *)upc_all_alloc(1, sizeof(int));
  nz_s = (shared [] int *)upc_all_alloc(1, sizeof(int));

  int m, n, nz;

  int i,j,rep;

  /* int *row_idx, *col_idx; */
  /* double *values; */
  /* double *x,*b; */

  struct timespec start,end;
  FILE *f;
  char line[64];
  if(r==ULONG_MAX) r=1000;

  if(MYTHREAD==0){


    char *filename = "matrix_in.csr";

    if ((f = fopen(filename, "r")) == NULL) {
      printf ("can't open file <%s> \n", filename);
      exit(1);
    }


    fgets(line, sizeof(line), f);
    int nz_t, m_t, n_t;
    sscanf(line, "%d %d %d", &nz_t, &n_t, &m_t);
    nz_s[0] = nz_t;
    m_s[0] = m_t;
    n_s[0] = n_t;

  }
  upc_barrier;
  nz = nz_s[0];
  m = m_s[0];
  n = n_s[0];

  printf("[%d] Number of elements of values and col_idx: %d; number of values in row_idx: %d\n", MYTHREAD, nz, m);



  upc_barrier;

  // Make all shared space
  /* row_idx = malloc(m * sizeof(int)); */
  /* col_idx = malloc(n * sizeof(int)); */
  /* values = malloc(nz * sizeof(double)); */
  /* x = malloc((m-1) * sizeof(double)); */
  /* b = malloc((m-1) * sizeof(double)); */
  row_idx = (shared [] int *)upc_all_alloc(1, m*sizeof(int));
  col_idx = (shared [] int *)upc_all_alloc(1, n*sizeof(int));
  values = (shared [] double *)upc_all_alloc(1, nz*sizeof(double));
  x = (shared [] double *)upc_all_alloc(1, (m-1)*sizeof(double));
  b = (shared [] double *)upc_all_alloc(1, (m-1)*sizeof(double));

  if (!row_idx || !col_idx || !values || !x || !b){
    printf ("cannot allocate memory for sparse matrix and vectors\n");
    exit(1);
  }
  else{
    printf("memory allocated\n");
  }

  upc_barrier;

  if(MYTHREAD==0){
    double t_d = 0;
    // fill values
    for(i=0;i<nz;i++){
      fgets(line, sizeof(line), f);
      sscanf(line, "%lf", &t_d);
      values[i] = t_d;
    }

    int t_i = 0;
    // fill col_idx
    for(i=0;i<nz;i++){
      fgets(line, sizeof(line), f);
      sscanf(line, "%d", &t_i);
      /* sscanf(line, "%d", &col_idx[i]); */
      col_idx[i] = t_i;
    }

    // fill row_idx
    for(i=0;i<m;i++){
      fgets(line, sizeof(line), f);
      sscanf(line, "%d", &t_i);
      row_idx[i] = t_i;
    }
    fclose(f);

  }
  // Use affinity of x
  upc_forall(i=0; i<m-1; i++;i)  {
    x[i]=i+1.5; // give basic values to vector x
    b[i] = 0;
  }

  upc_barrier;
  if(MYTHREAD==0){
    clock_gettime(CLOCK, &start);
  }

  /* Main algorithm loop */
  for(rep=0;rep<r;rep++){
    /* Ax=b */
    upc_forall(i=0; i<m-1; i++;i){
      for(j=row_idx[i];j<row_idx[i+1];j++){
        b[i] = b[i] + values[j]* x[col_idx[j]];
      }
    }
  }

  if(MYTHREAD==0){
    clock_gettime(CLOCK, &end);
    elapsed_time_hr(start, end, "Sparse DMVs.");
  }
  upc_all_free(x);
  upc_all_free(b);
  upc_all_free(row_idx);
  upc_all_free(col_idx);
  upc_all_free(values);



  return 0;
}

int double_spgemm(unsigned long r){

  int *m, *n, *nz;
  int *row_csr_idx, *col_csr_idx;
  int *row_csc_idx, *col_csc_idx;
  double *A_csr, *B_csc, *C; // matrices

  char *filename = "matrix_in.txt";
  char line[64];

  int i,j,k;
  unsigned long rep = 0;
  int count=0, nz_count = 0;

  struct timespec start,end;

  if(r==ULONG_MAX) r=100;

  m = malloc(sizeof(int));
  n = malloc(sizeof(int));
  nz = malloc(sizeof(int));

  get_matrix_size(filename, m, n, nz);
  printf("NZ = %d, M = %d, N = %d\n", *nz, *m, *n);


  row_csr_idx = malloc(*nz * sizeof(int));
  col_csr_idx = malloc(*nz * sizeof(int));
  A_csr = malloc(*nz * sizeof(double));
  row_csc_idx = malloc(*nz * sizeof(int));
  col_csc_idx = malloc(*nz * sizeof(int));
  B_csc = malloc(*nz * sizeof(double));
  C = calloc((*m * *n), sizeof(double));

  if (!row_csr_idx || !col_csr_idx || !A_csr ||  !row_csc_idx || !col_csc_idx || !B_csc || !C){
    printf ("cannot allocate memory for %d, %d, %d sparse matrices and vector\n", *m, *n, *nz);
    exit(1);
  }
  else{
    printf("memory allocated\n");
  }

  if(file_exists("matrix_in.csr")){

    printf("File exists.\n");

    clock_gettime(CLOCK, &start);

    FILE *f;
    f = fopen("matrix_in.csr", "r");
    /* printf("File open.\n"); */
    fgets(line, sizeof(line), f);
    sscanf(line, "%d %d %d", nz, n, m);

    // fill A_csr
    for(i=0;i<*nz;i++){
      fgets(line, sizeof(line), f);
      sscanf(line, "%lf", &A_csr[i]);
    }

    // fill col_csr_idx
    for(i=0;i<*nz;i++){
      fgets(line, sizeof(line), f);
      sscanf(line, "%d", &col_csr_idx[i]);
    }

    // fill row_csr_idx
    for(i=0;i<*m;i++){
      fgets(line, sizeof(line), f);
      sscanf(line, "%d", &row_csr_idx[i]);
    }

    fclose(f);
    clock_gettime(CLOCK, &end);
    elapsed_time_hr(start, end, "Read in CSR file");

  }
  else{
    clock_gettime(CLOCK, &start);
    mm_to_csr(filename, *m, *n, *nz, row_csr_idx, col_csr_idx, A_csr);
    clock_gettime(CLOCK, &end);
    elapsed_time_hr(start, end, "MM to CSR conversion");
  }

  *m = *m - 1;
  *n = *m;

  col_csc_idx[0]=0;

  /* create B_csc from A_csr */
  for(i=0; i<*m; i++){

    nz_count=0;

    for(j=0; j<*m; j++){

      for(k=row_csr_idx[j]; k<row_csr_idx[j+1]; k++){

        if(col_csr_idx[k]==i){
          B_csc[count]=A_csr[nz_count];
          row_csc_idx[count]=j;
          //      printf("j = %d and k = %d\n", j, k);
          count++;
        }
        nz_count++;
        //      printf("nz = %d\n", nz_count);
      }
      col_csc_idx[i+1]=count;
    }
  }

  clock_gettime(CLOCK, &start);

  for(rep=0; rep<r; rep++){

    /* AB=C */

    for(j=0; j<*n; j++){ // cols

      /* temporary padded column vector */
      double *temp_vec = malloc(*nz * sizeof(double));

      for(k = col_csc_idx[j]; k < col_csc_idx[j+1]; k++){
        temp_vec[row_csc_idx[k]]=B_csc[k];
      }

      /* spgemm */
      for(i=0; i<*m; i++){ // rows

        for (k = row_csr_idx[i]; k < row_csr_idx[i+1]; k++){
          C[j + i * *m] = C[j + i * *m] + A_csr[k] * temp_vec[col_csr_idx[k]];
        }

      }

      /* free temporary vector */
      free(temp_vec);

    }

  }

  clock_gettime(CLOCK, &end);

  elapsed_time_hr(start, end, "Sparse DGEMMs");

  /* free memory*/
  free(m);
  free(n);
  free(nz);

  free(row_csr_idx);
  free(col_csr_idx);
  free(row_csc_idx);
  free(col_csc_idx);
  free(A_csr);
  free(B_csc);
  free(C);

  return 0;
}
