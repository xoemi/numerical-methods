#include "task2.h"
#include "f.h"

int task2 (void) {
FILE *input;
FILE *output;
double sum = 0.0;
int N;
int V; // number of vertices
int T; // number of triangles
int I_E; // number of interior edges
int B_E; // number of boundary edges
int A, B, C;
double A_x, A_y, B_x, B_y, C_x, C_y, M1_x, M1_y, M2_x, M2_y, M3_x, M3_y;
double *vertices;
int *triangles;

    input = fopen("triangulation.txt", "r");
    if(!input) {
        printf("cannot open\n");
        return 1;
    }

    output = fopen("int_2D.txt", "a");
    if(!output) {
        fclose(input);
        printf("cannot open\n");
        return 1;
    }

    if (!(fscanf(input, "%d", &N)) || !(fscanf(input, "%d", &V)) || !(fscanf(input, "%d", &T)) || !(fscanf(input, "%d", &I_E)) || !(fscanf(input, "%d", &B_E))) {
        fclose(input);
        fclose(output);
        printf("cannot read n\n");
        return 2;
    }

    vertices = (double *)malloc(3 * V * sizeof(double));
    if (!vertices) {
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }
    
    triangles = (int *)malloc(4 * T * sizeof(int));
    if (!triangles) {
        free(vertices);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    for (int i = 0; i < 3 * V; i ++)
        if (!(fscanf(input, "%lf", vertices + i) == 1)) {
            free(vertices);
            free(triangles);            
            fclose(input);
            fclose(output);
            printf("cannot read\n");
            return 3;
        }

    for (int i = 0; i < 4 * T; i ++)
        if (!(fscanf(input, "%d", triangles + i) == 1)) {
            free(vertices);
            free(triangles);            
            fclose(input);
            fclose(output);
            printf("cannot read\n");
            return 3;
        }

    for (int i = 0; i < T; i ++) {
        A = triangles[1 + i * 4];
        A_x = vertices[1 + (A - 1) * 3];
        A_y = vertices[2 + (A - 1) * 3];

        B = triangles[2 + i * 4];
        B_x = vertices[1 + (B - 1) * 3];
        B_y = vertices[2 + (B - 1) * 3];

        C = triangles[3 + i * 4];
        C_x = vertices[1 + (C - 1) * 3];
        C_y = vertices[2 + (C - 1) * 3];

        M1_x = middle(A_x, B_x);
        M1_y = middle(A_y, B_y);

        M2_x = middle(A_x, C_x);
        M2_y = middle(A_y, C_y);

        M3_x = middle(B_x, C_x);
        M3_y = middle(B_y, C_y);

        sum += triangle_area(A_x, A_y, B_x, B_y, C_x, C_y) * (f(M1_x, M1_y) + f(M2_x, M2_y) + f(M3_x, M3_y)) / 3.0;
    }

    fprintf(output, "%15.15lf %15.15lf\n", log(N), log(fabs(sum - 23.0 / 45.0)));

free(vertices);
free(triangles);
fclose(input);
fclose(output);
return 0;
}

double middle (double a, double b) {
    
    return 0.5 * (a + b);

}

double triangle_area (double x1, double y1, double x2, double y2, double x3, double y3) {

    return 0.5 * fabs((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3));

}

