#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int task1 (int Nx, int Ny);
void triangulation(FILE *output, double Lx, double Ly, int Nx, int Ny);
void vertex_coordinates(FILE *output, double Lx, double Ly, int Nx, int Ny);
void triangle_vertices(FILE *output, int Nx, int Ny);
void interior_edge_endpoints(FILE *output, int Nx, int Ny);
void boundary_edge_endpoints(FILE *output, int Nx, int Ny);