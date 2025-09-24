#include <stdio.h>
#include <math.h>
#define TOL 1e-6
#define MAX_ITERATION 10000000

int main() {
    int n,h,k,T_AMBIENT,T_LeftWall,heatFlux;
    T_AMBIENT=300;
    heatFlux=200;
    h=25;
    k=0.1;
    T_LeftWall=(T_AMBIENT)-(heatFlux/h); // whatever heat flux outside the wall it gets transferred by convection to the wall
    printf("Enter the number of grid-points: ");
    scanf("%d", &n);
    
    float L = 0.1;
    float dx = L / (n + 1);  // Length of each partition.
    float T[n + 2]; // Array for temperature values
    float error;

    // Setting The Boundary Conditions
    T[0] =T_LeftWall; // Temperature at the left face
    T[n + 1] = 50; // Temperature at the right face

    // Initializing the temperature array at nodes
    for (int i = 1; i <= n; i++) {
        T[i] = 0;
    }

    // Using Gauss-Seidel method:
    int iter = 0;
    while (iter < MAX_ITERATION) {
        error = 0;
        for (int i = 1; i <= n; i++) {
            float TNEW = 0.5 * (T[i - 1] + T[i + 1]);
            float residue = fabs(TNEW - T[i]);
            T[i] = TNEW; // Update the temperature value
            
            if (residue > error) {
                error = residue;
            }
        }

        if (error < TOL) {
            break;
        }
        iter++;
    }

    float x = 0;
    printf("X        TEMPERATURE\n");
    for (int i = 0; i <= n + 1; i++) {
        printf("%.3fm        %.3f\n", x, T[i]);
        x += dx;
    }

    // Saving the results to a text file
    FILE *fptr = fopen("q2.txt", "w");
    if (fptr == NULL) {
        printf("Error opening file!\n");
        return 1; // Return an error code
    }

    x = 0; // Reset x for file writing
    fprintf(fptr, "X(in m)       TEMPERATURE\n");
    for (int i = 0; i <= n + 1; i++) {
        fprintf(fptr, "%.3fm        %.3f\n", x, T[i]);
        x += dx;
    }

    fclose(fptr);
    printf("Data successfully written to q2.txt\n");

    return 0;
}