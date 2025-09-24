#include <stdio.h>

// Function prototype
void solveTemperature(double *T, int N, double Q, double k, double L, double T0, double TN);

int main() {
    // Parameters
    int N = 6; // Number of points (including boundaries)
    double L = 0.02; // Thickness of the plate in meters
    double k = 0.5; // Thermal conductivity in W/m-K
    double Q = 1e6; // Heat generation in W/m^3
    double T0 = 100; // Temperature at x = 0 in °C
    double TN = 200; // Temperature at x = L in °C

    double T[N]; // Temperature array

    // Solve the system
    solveTemperature(T, N, Q, k, L, T0, TN);

    // Open file for writing
    FILE *file = fopen("numerical_solution.txt", "w");
    if (file == NULL) {
        printf("Error opening file for writing.\n");
        return 1;
    }

    // Write the results to the file
    fprintf(file, "Numerical Temperature Distribution:\n");
    for (int i = 0; i < N; i++) {
        fprintf(file, "T[%d] = %.2f °C\n", i, T[i]);
    }

    // Close the file
    fclose(file);

    return 0;
}

void solveTemperature(double *T, int N, double Q, double k, double L, double T0, double TN) {
    double dx = L / (N - 1); // Distance between points
    double alpha = Q * dx * dx / k;

    // Coefficients of the matrix
    double A[N][N];
    double B[N];

    // Initialize the matrix and the right-hand side vector
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = 0;
        }
        B[i] = 0;
    }

    // Boundary conditions
    T[0] = T0;
    T[N - 1] = TN;

    // Set up the matrix for interior points
    for (int i = 1; i < N - 1; i++) {
        A[i][i - 1] = 1;
        A[i][i] = -2;
        A[i][i + 1] = 1;
        B[i] = -alpha;
    }

    // Set boundary conditions in the matrix
    B[0] = T0;
    B[N - 1] = TN;

    // Solve the system using Gaussian elimination
    for (int i = 1; i < N - 1; i++) {
        for (int j = i + 1; j < N - 1; j++) {
            double factor = A[j][i - 1] / A[i][i];
            for (int k = i; k < N; k++) {
                A[j][k] -= factor * A[i][k];
            }
            B[j] -= factor * B[i];
        }
    }

    // Back-substitution
    for (int i = N - 2; i > 0; i--) {
        for (int j = i + 1; j < N - 1; j++) {
            B[i] -= A[i][j] * T[j];
        }
        T[i] = B[i] / A[i][i];
    }
}
