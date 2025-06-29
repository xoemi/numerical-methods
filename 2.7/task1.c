#include "task1.h"
#include "f.h"

int task1 (int Nx, int Ny) {
FILE *output;
double Lx = 1.0, Ly = 1.0; // given rectangle 

    output = fopen("triangulation.txt", "w");
    if(!output) {
        printf("cannot open\n");
        return 1;
    }
    
    triangulation(output, Lx, Ly, Nx, Ny);

fclose(output);
return 0;
}


void triangulation(FILE *output, double Lx, double Ly, int Nx, int Ny) {

    fprintf(output, "%d\n", Nx); // N for task 2
    fprintf(output, "%d\n", (Nx + 1) * (Ny + 1)); // number of vertices
    fprintf(output, "%d\n", 2 * Nx * Ny); // number of triangles
    fprintf(output, "%d\n", 3 * Nx  * Ny - (Nx + Ny)); // number of interior edges
    fprintf(output, "%d\n", 2 * (Nx + Ny)); // number of boundary edges

    vertex_coordinates(output, Lx, Ly, Nx, Ny);
    triangle_vertices(output, Nx, Ny);
    interior_edge_endpoints(output, Nx, Ny);
    boundary_edge_endpoints(output, Nx, Ny);

}


void vertex_coordinates(FILE *output, double Lx, double Ly, int Nx, int Ny) {
double x = Lx / Nx, y = Ly / Ny;

    for (int i = 0; i <= Nx; i ++)
        for (int j = 0; j <= Ny; j ++) 
            fprintf(output, "%d %lf %lf\n", j + i * (Ny + 1) + 1, i * x, j * y);

}


void triangle_vertices(FILE *output, int Nx, int Ny) {

    for (int i = 0; i < Nx; i ++)
        for (int j = 0; j < Ny; j ++) {
            fprintf(output, "%d %d %d %d\n", 2 * (j + i * Ny) + 1, j + i * (Ny + 1) + 1, j + i * (Ny + 1) + 2, j + (i + 1) * (Ny + 1) + 1 ); 
            fprintf(output, "%d %d %d %d\n", 2 * (j + i * Ny) + 2, j + (i + 1) * (Ny + 1) + 1, j + (i + 1) * (Ny + 1) + 2, j + i * (Ny + 1) + 2);
        }

}


void interior_edge_endpoints(FILE *output, int Nx, int Ny) {
int edge_num = 1;

    for(int i = 0; i < Nx; i ++)
        for(int j = 0; j < Ny; j++)
            fprintf(output, "%d %d %d\n", edge_num ++, j  + (i + 1) * (Ny + 1) + 1,  j + i * (Ny + 1) + 2);

    for(int i = 1; i < Nx; i ++)
        for (int j = 0; j < Ny; j++)
            fprintf(output, "%d %d %d\n", edge_num ++, j + i * (Ny + 1) + 1,  j + i * (Ny + 1) + 2);

    for(int i = 0; i < Nx; i ++)
        for (int j = 1; j < Ny; j++)
            fprintf(output, "%d %d %d\n", edge_num ++,  j + i * (Ny + 1) + 1,  j + (i + 1) * (Ny + 1) + 1);
    
}


void boundary_edge_endpoints(FILE *output, int Nx, int Ny) {
int edge_num = 1;

    for (int j = 0; j < Ny; j ++) 
        fprintf(output, "%d %d %d\n", edge_num ++,  j + 1,  j + 2);

    for (int j = 0; j < Ny; j ++) 
        fprintf(output, "%d %d %d\n", edge_num ++,  j + Nx * (Ny + 1) + 1,  j + Nx * (Ny + 1) + 2);

    for (int i = 0; i < Nx; i ++) 
        fprintf(output, "%d %d %d\n", edge_num ++,  i * (Ny + 1) + 1,  (i + 1) * (Ny + 1) + 1);

    for (int i = 0; i < Nx; i ++) 
        fprintf(output, "%d %d %d\n", edge_num ++,  i * (Ny + 1) + (Nx + 1),  (i + 1) * (Ny + 1) + (Nx + 1));

}
