#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_N 1000

// Function prototypes
void read_input(int *N, double *dx, double *T_left, double *T_right, double *tolerance);
void solve_fd(int N, double dx, double T_left, double T_right, double tolerance, double *T);
void write_output(int N, double dx, double *T);
void plot_results(int N, double dx, double *T, double T_left, double T_right);

int main() {
    int N;
    double dx;
    double T_left, T_right;
    double tolerance;
    double T[MAX_N];

    // Read input
    read_input(&N, &dx, &T_left, &T_right, &tolerance);

    // Solve the finite difference problem
    solve_fd(N, dx, T_left, T_right, tolerance, T);

    // Write output
    write_output(N, dx, T);

    // Plot results
    plot_results(N, dx, T, T_left, T_right);

    return 0;
}

void read_input(int *N, double *dx, double *T_left, double *T_right, double *tolerance) {
    FILE *file = fopen("input.dat", "r"); // No path, will look in the current directory
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    fscanf(file, "%d", N);
    fscanf(file, "%lf", dx);
    fscanf(file, "%lf", T_left);
    fscanf(file, "%lf", T_right);
    fscanf(file, "%lf", tolerance);
    fclose(file);
}

void solve_fd(int N, double dx, double T_left, double T_right, double tolerance, double *T) {
    double *T_new = (double *)malloc(N * sizeof(double));
    int i;
    double max_error;

    // Initialize boundary conditions and initial guess
    T[0] = T_left;
    T[N-1] = T_right;
    for (i = 1; i < N-1; i++) {
        T[i] = 0.0;  // Initial guess
    }

    // Main iterative loop
    do {
        max_error = 0.0;
        for (i = 1; i < N-1; i++) {
            T_new[i] = 0.5 * (T[i-1] + T[i+1]);
            max_error = fmax(max_error, fabs(T_new[i] - T[i]));
        }
        // Update temperature array
        for (i = 1; i < N-1; i++) {
            T[i] = T_new[i];
        }
    } while (max_error > tolerance);

    free(T_new);
}

void write_output(int N, double dx, double *T) {
    FILE *file = fopen("output.txt", "w"); // No path, will save in the current directory
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "Position (m)    Temperature (C)\n");
    for (int i = 0; i < N; i++) {
        fprintf(file, "%.2f          %.2f\n", i * dx, T[i]);
    }
    fclose(file);
}

void plot_results(int N, double dx, double *T, double T_left, double T_right) {
    FILE *file = fopen("analytical_solution.txt", "w"); // No path, will save in the current directory
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "Position (m)    Numerical Temp (C)    Analytical Temp (C)\n");
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        double L = (N - 1) * dx;  // Correct domain length
        double T_analytical = T_left + (T_right - T_left) * x / L;
        fprintf(file, "%.2f          %.2f          %.2f\n", x, T[i], T_analytical);
    }
    fclose(file);
}    
