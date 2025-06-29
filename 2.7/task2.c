#include "task2.h"
#include "f.h"

int task2 (void) {
FILE *input;
FILE *output;
double P1[2] = {0};
double P2[2] = {0};
int A, B;
int N;
int V; // number of vertices
int T; // number of triangles
int I_E; // number of interior edges
int B_E; // number of boundary edges
double max_dif = 0.0;
double *vertices;
int *triangles;
int *i_edges;
int *b_edges;
double *i_middle;
double *b_middle;

    input = fopen("triangulation.txt", "r");
    if(!input) {
        printf("cannot open\n");
        return 1;
    }

    output = fopen("res2.txt", "w");
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

    i_edges = (int *)malloc(3 * I_E * sizeof(int));
    if (!i_edges) {
        free(vertices);
        free(triangles);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    b_edges = (int *)malloc(3 * B_E * sizeof(int));
    if (!b_edges) {
        free(vertices);
        free(triangles);
        free(i_edges);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    i_middle = (double *)malloc(2 * I_E * sizeof(double));
    if (!i_middle) {
        free(vertices);
        free(triangles);
        free(i_edges);
        free(b_edges);  
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    b_middle = (double *)malloc(2 * B_E * sizeof(double));
    if (!b_middle) {
        free(vertices);
        free(triangles);
        free(i_edges);
        free(b_edges);
        free(i_edges);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }


    for (int i = 0; i < 3 * V; i ++)
        if (!(fscanf(input, "%lf", vertices + i) == 1)) {
            free(vertices);
            free(triangles);
            free(i_edges);
            free(b_edges);    
            free(i_middle);
            free(b_middle);        
            fclose(input);
            fclose(output);
            printf("cannot read\n");
            return 3;
        }

    for (int i = 0; i < 4 * T; i ++)
        if (!(fscanf(input, "%d", triangles + i) == 1)) {
            free(vertices);
            free(i_edges);
            free(b_edges);
            free(i_middle);
            free(b_middle);    
            free(triangles);            
            fclose(input);
            fclose(output);
            printf("cannot read\n");
            return 3;
        }

    for (int i = 0; i < 3 * I_E; i ++)
        if (!(fscanf(input, "%d", i_edges + i) == 1)) {
            free(vertices);
            free(triangles);            
            fclose(input);
            free(i_edges);
            free(b_edges);
            free(i_middle);
            free(b_middle);    
            fclose(output);
            printf("cannot read\n");
            return 3;
        }

    for (int i = 0; i < 3 * B_E; i ++)
        if (!(fscanf(input, "%d", b_edges + i) == 1)) {
            free(vertices);
            free(i_edges);
            free(b_edges);
            free(i_middle);
            free(b_middle);    
            free(triangles);            
            fclose(input);
            fclose(output);
            printf("cannot read\n");
            return 3;
        }

    for (int i = 0; i < I_E; i ++) {
        A = i_edges[1 + i * 3];
        B = i_edges[2 + i * 3];

        P1[0] = vertices[1 + (A - 1) * 3];
        P1[1] = vertices[2 + (A - 1) * 3];

        P2[0] = vertices[1 + (B - 1) * 3];
        P2[1] = vertices[2 + (B - 1) * 3];

        i_middle[i * 2] = middle(P1[0], P2[0]);
        i_middle[1 + i * 2] = middle(P1[1], P2[1]);
    }

    for (int i = 0; i < B_E; i ++) {
        A = b_edges[1 + i * 3];
        B = b_edges[2 + i * 3];

        P1[0] = vertices[1 + (A - 1) * 3];
        P1[1] = vertices[2 + (A - 1) * 3];

        P2[0] = vertices[1 + (B - 1) * 3];
        P2[1] = vertices[2 + (B - 1) * 3];

        b_middle[i * 2] = middle(P1[0], P2[0]);
        b_middle[1 + i * 2] = middle(P1[1], P2[1]);
    }

    for (int i = 0; i <= N; i ++) {
        for (int j = 0; j <= N; j ++) {
            double x = 1.0 * i / N;
            double y = 1.0 * j / N;
            double sum = 0.0;
            for (int k = 0; k < V; k ++) {
                sum += f(vertices[1 + 3 * k], vertices[2 + 3 * k]) * phi_v(k, T, x, y, triangles, vertices);
            }
            for (int k = 0; k < I_E; k ++) {
                sum += f(i_middle[2 * k], i_middle[1 + 2 * k]) * phi_im(k, T, x, y, triangles, i_middle, i_edges, vertices);
            }
            for (int k = 0; k < B_E; k ++) { 
                sum += f(b_middle[2 * k], b_middle[1 + 2 * k]) * phi_bm(k, T, x, y, triangles, b_middle, b_edges, vertices);
            }

            if (fabs(sum - f(x, y)) >= max_dif)
                max_dif = fabs(sum - f(x, y));
            fprintf(output, "%lf %lf %lf %lf\n", x, y, f(x, y), sum);
        }
    }
    printf("max_dif_on_vertices = %e\n", max_dif);
    max_dif = 0.0;

    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            double x = 1.0 * (i + 0.5) / N;
            double y = 1.0 * (j + 0.5) / N;
            double sum = 0.0;
            for (int k = 0; k < V; k ++) {
                sum += f(vertices[1 + 3 * k], vertices[2 + 3 * k]) * phi_v(k, T, x, y, triangles, vertices);
            }
            for (int k = 0; k < I_E; k ++) {
                sum += f(i_middle[2 * k], i_middle[1 + 2 * k]) * phi_im(k, T, x, y, triangles, i_middle, i_edges, vertices);
            }
            for (int k = 0; k < B_E; k ++) { 
                sum += f(b_middle[2 * k], b_middle[1 + 2 * k]) * phi_bm(k, T, x, y, triangles, b_middle, b_edges, vertices);
            }

            if (fabs(sum - f(x, y)) >= max_dif)
                max_dif = fabs(sum - f(x, y));
            fprintf(output, "%lf %lf %lf %lf\n", x, y, f(x, y), sum);
        }
    }
    printf("max_dif_not_on_vertices = %e\n", max_dif);


free(vertices);
free(i_edges);
free(b_edges);
free(i_middle);
free(b_middle);    
free(triangles);
fclose(input);
fclose(output);
return 0;
}


double middle (double a, double b) {
    
    return 0.5 * (a + b);

}


double phi_v(int i, int Tr, double x, double y, int *triangles, double *vertices) {
int f = 0;
int T = 0;

    for (int j = 0; j < Tr; j ++) { 
        if (triangles[1 + j * 4] == (i + 1) || triangles[2 + j * 4] == (i + 1) || triangles[3 + j * 4] == (i + 1)) {
            if (point_in_triangle(x, y, j, vertices, triangles)) {
                T = triangles[j * 4];
                break;
        }
        }
    }
    
if (T == 0)
    return 0;

return phi(T, x, y, vertices[1 + i * 3], vertices[2 + i * 3], vertices, triangles, i + 1, i + 1, f);
}


double phi_im(int i, int Tr, double x, double y, int *triangles, double *i_middle, int *i_edges, double *vertices) {
int A, B, f = 1;
int T = 0;

    A = i_edges[1 + i * 3];
    B = i_edges[2 + i * 3];

    for (int j = 0; j < Tr; j ++) { 
         if ((triangles[1 + j * 4] == A || triangles[2 + j * 4] == A || triangles[3 + j * 4] == A) && (triangles[1 + j * 4] == B || triangles[2 + j * 4] == B || triangles[3 + j * 4] == B)) {
            if (point_in_triangle(x, y, j, vertices, triangles)) {
                T = triangles[j * 4];
                break;
        }
        }
    }

if (T == 0)
    return 0;

return phi(T, x, y, i_middle[i * 2], i_middle[1 + i * 2], vertices, triangles, A, B, f);
}


double phi_bm(int i, int Tr, double x, double y, int *triangles, double *b_middle, int *b_edges, double *vertices) {
int A, B, f = 1;
int T = 0;

    A = b_edges[1 + i * 3];
    B = b_edges[2 + i * 3];

    for (int j = 0; j < Tr; j ++) { 
         if ((triangles[1 + j * 4] == A || triangles[2 + j * 4] == A || triangles[3 + j * 4] == A) && (triangles[1 + j * 4] == B || triangles[2 + j * 4] == B || triangles[3 + j * 4] == B)) {
            if (point_in_triangle(x, y, j, vertices, triangles)) {
                T = triangles[j * 4];
                break;
        }
        }
    }

if (T == 0)
    return 0;

return phi(T, x, y, b_middle[i * 2], b_middle[1 + i * 2], vertices, triangles, A, B, f);
}


int point_in_triangle (double x, double y, int j, double *vertices, int *triangles) {
double P1[2] = {vertices[1 + 3 * (triangles[1 + 4 * j] - 1)], vertices[2 + 3 * (triangles[1 + 4 * j] - 1)]};
double P2[2] = {vertices[1 + 3 * (triangles[2 + 4 * j] - 1)], vertices[2 + 3 * (triangles[2 + 4 * j] - 1)]};
double P3[2] = {vertices[1 + 3 * (triangles[3 + 4 * j] - 1)], vertices[2 + 3 * (triangles[3 + 4 * j] - 1)]};
double d1 = 0.0, d2 = 0.0, d3 = 0.0;

    d1 = (x - P1[0]) * (P2[1] - P1[1]) - (y - P1[1]) * (P2[0] - P1[0]);
    d2 = (x - P2[0]) * (P3[1] - P2[1]) - (y - P2[1]) * (P3[0] - P2[0]);
    d3 = (x - P3[0]) * (P1[1] - P3[1]) - (y - P3[1]) * (P1[0] - P3[0]);

    if (((d1 <= 0) && (d2 <= 0) && (d3 <= 0)) || ((d1 >= 0) && (d2 >= 0) && (d3 >= 0)))
        return 1;

return 0;
}


double phi(int i, double x, double y, double x_i, double y_i, double *vertices, int *triangles, int A, int B, int f) {
double P1[2] = {vertices[1 + 3 * (triangles[1 + 4 * (i - 1)] - 1)], vertices[2 + 3 * (triangles[1 + 4 * (i - 1)] - 1)]};
double P2[2] = {vertices[1 + 3 * (triangles[2 + 4 * (i - 1)] - 1)], vertices[2 + 3 * (triangles[2 + 4 * (i - 1)] - 1)]};
double P3[2] = {vertices[1 + 3 * (triangles[3 + 4 * (i - 1)] - 1)], vertices[2 + 3 * (triangles[3 + 4 * (i - 1)] - 1)]};
double P4[2] = {middle(P1[0], P2[0]), middle(P1[1], P2[1])};
double P5[2] = {middle(P2[0], P3[0]), middle(P2[1], P3[1])};
double P6[2] = {middle(P1[0], P3[0]), middle(P1[1], P3[1])};
double psi1[2], psi2[2], psi3[2], psi4[2], psi5[2], psi6[2];
double eta1[2], eta2[2], eta3[2], eta4[2], eta5[2], eta6[2];  //0 for x, y; 1 for x_i, y_i
int k = 0;

    psi1[0] = (x - P2[0]) * (P3[1] - P2[1]) - (y - P2[1]) * (P3[0] - P2[0]);
    psi2[0] = (x - P1[0]) * (P3[1] - P1[1]) - (y - P1[1]) * (P3[0] - P1[0]);
    psi3[0] = (x - P1[0]) * (P2[1] - P1[1]) - (y - P1[1]) * (P2[0] - P1[0]);
    psi4[0] = (x - P4[0]) * (P6[1] - P4[1]) - (y - P4[1]) * (P6[0] - P4[0]);
    psi5[0] = (x - P4[0]) * (P5[1] - P4[1]) - (y - P4[1]) * (P5[0] - P4[0]);
    psi6[0] = (x - P5[0]) * (P6[1] - P5[1]) - (y - P5[1]) * (P6[0] - P5[0]);

    psi1[1] = (x_i - P2[0]) * (P3[1] - P2[1]) - (y_i - P2[1]) * (P3[0] - P2[0]);
    psi2[1] = (x_i - P1[0]) * (P3[1] - P1[1]) - (y_i - P1[1]) * (P3[0] - P1[0]);
    psi3[1] = (x_i - P1[0]) * (P2[1] - P1[1]) - (y_i - P1[1]) * (P2[0] - P1[0]);
    psi4[1] = (x_i - P4[0]) * (P6[1] - P4[1]) - (y_i - P4[1]) * (P6[0] - P4[0]);
    psi5[1] = (x_i - P4[0]) * (P5[1] - P4[1]) - (y_i - P4[1]) * (P5[0] - P4[0]);
    psi6[1] = (x_i - P5[0]) * (P6[1] - P5[1]) - (y_i - P5[1]) * (P6[0] - P5[0]);

    eta1[0] = psi1[0] * psi4[0];
    eta2[0] = psi2[0] * psi5[0];
    eta3[0] = psi3[0] * psi6[0];
    eta4[0] = psi1[0] * psi2[0];
    eta5[0] = psi2[0] * psi3[0];
    eta6[0] = psi1[0] * psi3[0];

    eta1[1] = psi1[1] * psi4[1];
    eta2[1] = psi2[1] * psi5[1];
    eta3[1] = psi3[1] * psi6[1];
    eta4[1] = psi1[1] * psi2[1];
    eta5[1] = psi2[1] * psi3[1];
    eta6[1] = psi1[1] * psi3[1];
    
    if (f == 0) {
        if (A == (triangles[1 + 4 * (i - 1)])) 
            k = 1;
        if (A == (triangles[2 + 4 * (i - 1)])) 
            k = 2;
        if (A == (triangles[3 + 4 * (i - 1)])) 
            k = 3;
    }
    else {
        if ((A == (triangles[1 + 4 * (i - 1)]) && B == (triangles[2 + 4 * (i - 1)])) || (A == (triangles[2 + 4 * (i - 1)]) && B == (triangles[1 + 4 * (i - 1)])))
            k = 4;
        if ((A == (triangles[1 + 4 * (i - 1)]) && B == (triangles[3 + 4 * (i - 1)])) || (A == (triangles[3 + 4 * (i - 1)]) && B == (triangles[1 + 4 * (i - 1)])))
            k = 6;
        if ((A == (triangles[3 + 4 * (i - 1)]) && B == (triangles[2 + 4 * (i - 1)])) || (A == (triangles[2 + 4 * (i - 1)]) && B == (triangles[3 + 4 * (i - 1)])))
            k = 5;
    }
    
    switch(k) {
        case 1: return (eta1[0] / eta1[1]);
        case 2: return (eta2[0] / eta2[1]);
        case 3: return (eta3[0] / eta3[1]);
        case 4: return (eta4[0] / eta4[1]);
        case 5: return (eta5[0] / eta5[1]);
        case 6: return (eta6[0] / eta6[1]);
    }

return 0;
}
