#include "task3.h"
#include "f.h"

int task3 (void) {
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
double *M;
double *MT;
double *MTb;
double *b;
double *P;
double *c;
int *order;
int n = 0, m = 0;

    input = fopen("triangulation.txt", "r");
    if(!input) {
        printf("cannot open\n");
        return 1;
    }

    output = fopen("res3.txt", "w");
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
        free(i_middle);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    M = (double *)malloc((2 * N + 1) * (2 * N + 1) * (V + I_E + B_E) * sizeof(double));
    if (!M) {
        free(vertices);
        free(triangles);
        free(i_edges);
        free(b_edges);
        free(i_middle);
        free(b_middle);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    b = (double *)malloc((2 * N + 1) * (2 * N + 1) * sizeof(double));
    if (!b) {
        free(vertices);
        free(triangles);
        free(i_edges);
        free(b_edges);
        free(i_middle);
        free(b_middle);
        free(M);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    MT = (double *)malloc((2 * N + 1) * (2 * N + 1) * (V + I_E + B_E) * sizeof(double));
    if (!M) {
        free(vertices);
        free(triangles);
        free(i_edges);
        free(b_edges);
        free(i_middle);
        free(b_middle);
        free(b);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    P = (double *)malloc((2 * N + 1) * (2 * N + 1) * (V + I_E + B_E) * sizeof(double));
    if (!M) {
        free(vertices);
        free(triangles);
        free(i_edges);
        free(b_edges);
        free(i_middle);
        free(b_middle);
        free(b);
        free(MT);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    c = (double *)malloc((2 * N + 1) * (2 * N + 1) * sizeof(double));
    if (!c) {
        free(vertices);
        free(triangles);
        free(i_edges);
        free(b_edges);
        free(i_middle);
        free(b_middle);
        free(M);
        free(MT);
        free(P);
        free(b);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    order = (int *)malloc((2 * N + 1) * (2 * N + 1) * sizeof(int));
    if (!c) {
        free(vertices);
        free(triangles);
        free(i_edges);
        free(b_edges);
        free(i_middle);
        free(b_middle);
        free(M);
        free(MT);
        free(P);
        free(b);
        free(c);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }

    MTb = (double *)malloc((V + I_E + B_E) * sizeof(double));
    if (!c) {
        free(vertices);
        free(triangles);
        free(i_edges);
        free(b_edges);
        free(i_middle);
        free(b_middle);
        free(M);
        free(MT);
        free(P);
        free(b);
        free(c);
        free(order);
        fclose(input);
        fclose(output);
        printf("cannot allocate\n");
        return 3;
    }


    for (int i = 0; i < (2 * N + 1) * (2 * N + 1); i ++)
        order[i] = i;

    for (int i = 0; i < 3 * V; i ++)
        if (!(fscanf(input, "%lf", vertices + i) == 1)) {
            free(vertices);
            free(triangles);
            free(i_edges);
            free(b_edges);    
            free(i_middle);
            free(b_middle);  
            free(M);    
            free(b);  
            free(MT);
            free(P);
            free(c);
            free(MTb);
            free(order);
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
            free(M);   
            free(b);  
            free(MT);
            free(P);
            free(c);
            free(MTb);
            free(order);
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
            free(M);      
            free(b); 
            free(MT);
            free(c);
            free(order);
            free(P);
            free(MTb);
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
            free(M);       
            free(b);   
            free(MT);
            free(P); 
            free(c);
            free(MTb);
            free(order);
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

    for (double i = 0; i <= N; i += 0.5) {
        n = 0;
        for (double j = 0; j <= N; j += 0.5) {
            double x = 1.0 * i / N;
            double y = 1.0 * j / N;
            int t = 0;

            b[(m * (2 * N + 1) + n)] = f(x, y);

            for (int k = 0; k < V; k ++) {
                M[(V + I_E + B_E) * (m * (2 * N + 1) + n) + t] =  phi_v(k, T, x, y, triangles, vertices);
                t ++;
            }
            for (int k = 0; k < I_E; k ++) {
                M[(V + I_E + B_E) * (m * (2 * N + 1) + n) + t] =  phi_im(k, T, x, y, triangles, i_middle, i_edges, vertices);
                t ++;
            }
            for (int k = 0; k < B_E; k ++) { 
                M[(V + I_E + B_E) * (m * (2 * N + 1) + n) + t] =  phi_bm(k, T, x, y, triangles, b_middle, b_edges, vertices);
                t ++;
            }
            n ++;
        }
        m ++;
    }

    get_transposed_matrix((2 * N + 1) * (2 * N + 1), (V + I_E + B_E), M, MT);

    matrix_product((V + I_E + B_E), (2 * N + 1) * (2 * N + 1), MT, M, P);
    //printf("\n");
    //print_matrix((2 * N + 1) * (2 * N + 1), (2 * N + 1) * (2 * N + 1), P);

    matrix_vec_product((V + I_E + B_E), (2 * N + 1) * (2 * N + 1), MT, b, MTb);
    solve_system((V + I_E + B_E), P, MTb, order, c);

    //printf("\n");
    //print_vec((2 * N + 1) * (2 * N + 1), c);
    
    for (double i = 0; i <= N; i += 0.5) {
        for (double j = 0; j <= N; j += 0.5) {
            double x = 1.0 * i / N;
            double y = 1.0 * j / N;
            double sum = 0.0;
            int t = 0;

            for (int k = 0; k < V; k ++) {
                sum += c[t] * phi_v(k, T, x, y, triangles, vertices);
                t ++;
            }
            for (int k = 0; k < I_E; k ++) {
                sum += c[t] * phi_im(k, T, x, y, triangles, i_middle, i_edges, vertices);
                t ++;
            }
            for (int k = 0; k < B_E; k ++) { 
                sum += c[t] * phi_bm(k, T, x, y, triangles, b_middle, b_edges, vertices);
                t ++;
            }
            if (fabs(sum - f(x, y)) >= max_dif)
                max_dif = fabs(sum - f(x, y));
            fprintf(output, "%lf %lf %lf %lf\n", x, y, f(x, y), sum);
        }
        
    }
    printf("max_dif = %e\n", max_dif);


free(vertices);
free(i_edges);
free(b_edges);
free(i_middle);
free(b_middle);    
free(triangles);
free(M);     
free(b); 
free(MT);
free(MTb);
free(P);
free(c);
free(order);
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

void get_transposed_matrix (int m, int n, double *M, double *MT) {

    for (int i = 0; i < m; i ++)
        for (int j = 0; j < n; j ++) {
            MT[j * m + i] = M[i * n + j];
        }


}

void print_matrix (int m, int n, double *A) {

    for (int i = 0; i < m; i ++) {
        for (int j = 0; j < n; j ++)
            printf("%lf ", A[i * n + j]);
        printf("\n");
    }
}


void print_vec (int n, double *A) {

    for (int j = 0; j < n; j ++)
        printf("%e\n", A[j]);

}


void matrix_product (int m, int n, double *a, double *b, double *res) {
   
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            res[i * n + j] = 0.0;
            for (int k = 0; k < n; k++) {
                res[i* n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }

}


void matrix_vec_product (int m, int n, double *matrix, double *vector, double *res) {

    for (int i = 0; i < m; i++) {
        res[i] = 0.0;
        for (int j = 0; j < n; j++) {
            res[i] += matrix[i * n + j] * vector[j];
        }
    }

}


int solve_system (int n, double *a, double *b, int *order, double *res) {
double tmp;
double tmp_;
double max = 0.0;
int i_max = 0, j_max = 0; 

    for (int i = 0; i < n; i ++) {
        max = fabs(a[i + i * n]);
        i_max = i;
        j_max = i;
        for (int f = i; f < n; f ++) 
            for (int g = i + 1; g < n; g ++) 
                if (fabs(a[g + f * n]) > max) {
                    max = fabs(a[g + f * n]);
                    i_max = f;
                    j_max = g;
                }
        swap_row(a, n, i_max, i); 
        swap_row_in_vec(b, i_max, i); 
        swap_col(a, n, j_max, i);
        tmp_ = order[i];                           
        order[i] = order[j_max];
        order[j_max] = tmp_;

        if (fabs(a[i*n + i]) < 1e-15) {
            printf("cant help you");
            return 1;   
        }
        tmp = a[i + i * n];
        for (int t = i + 1; t < n; t ++) 
            a[t + i * n] /= tmp;
        b[i] /= tmp;

        for (int k = i + 1; k < n; k ++) {
            tmp = a[i + k * n];
            for (int t = i + 1; t < n; t ++)
                a[t + k * n] -= tmp * a[t + i * n];
            b[k] -= b[i] * tmp;
        }
    }

    for (int i = n - 1; i >= 0; i --) {
        tmp = 0.0;
        for (int t = i + 1; t < n; t ++) 
            tmp += a[t + i * n] * res[order[t]];
        res[order[i]] = b[i] - tmp;
        
    }

return 0;
}


void swap_col (double *a, int n, int i, int j) {
double tmp;
double *a1 = a + i;
double *a2 = a + j;
    
    if (i == j) return;
    
    for (int k = 0; k < n; k ++) {
        tmp = a1[k * n];
        a1[k * n] = a2[k * n];
        a2[k * n] = tmp;
    }
    
}



void swap_row_in_vec (double *a, int i, int j) {
double tmp;
    
    if (i == j) return;
    
    tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;

}


void swap_row (double *a, int n, int i, int j) {
double tmp;
double *a1 = a + i * n;
double *a2 = a + j * n;
    
    if (i == j) return;
    
    for (int k = 0; k < n; k ++) {
        tmp = a1[k];
        a1[k] = a2[k];
        a2[k] = tmp;
    }

}
