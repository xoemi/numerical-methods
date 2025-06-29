#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int task3 (void);
double middle (double a, double b);
double phi_v(int i, int Tr, double x, double y, int *triangles, double *vertices);
double phi_im(int i, int Tr, double x, double y, int *triangles, double *i_middle, int *i_edges, double *vertices);
double phi_bm(int i, int Tr, double x, double y, int *triangles, double *b_middle, int *b_edges, double *vertices);
int point_in_triangle (double x, double y, int j, double *vertices, int *triangles);
double phi(int i, double x, double y, double x_i, double y_i, double *vertices, int *triangles, int A, int B, int f);
void get_transposed_matrix (int m, int n, double *M, double *MT);
void print_matrix (int m, int n, double *A);
void print_vec (int n, double *A);
void matrix_product (int m, int n, double *a, double *b, double *res);
void matrix_vec_product (int m, int n, double *matrix, double *vector, double *res);
int solve_system (int n, double *a, double *b, int *order, double *res);
void swap_col (double *a, int n, int i, int j);
void swap_row_in_vec (double *a, int i, int j);
void swap_row (double *a, int n, int i, int j);
