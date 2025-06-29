#include "task1.h"

int main (int argc, char *argv[]) {
int Nx, Ny;

    if (!(argc == 3 && sscanf(argv[1], "%d", &Nx) == 1 && sscanf(argv[2], "%d", &Ny) == 1 && Nx >= 1 && Ny >= 1)) {
        printf("\nusage:\n\n %s Nx, Ny >= 1\n", argv[0]);
        return 1;
    }

    if (task1(Nx, Ny))
        return 2;

return 0; 
}
