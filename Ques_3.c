#include <stdio.h>
#include <math.h>

int max_iterations = 1000;
double tolerance = 1e-6;

// Function to calculate flux term
double Flux(double V)
{
    return (V - pow(V, 2)) / 2.0;
}

// Predictor-Corrector Scheme
double PredictorCorrector(double V[], double time_step, double space_step)
{
    double alpha = 0.001 * time_step / pow(space_step, 2), error_array[1000], V_intermediate[41], V_updated[41];
    int iteration;

    // Set boundary conditions
    V_intermediate[0] = V[0];
    V_intermediate[40] = V[40];
    V_updated[0] = V[0];
    V_updated[40] = V[40];

    for (iteration = 0; iteration < max_iterations; iteration++)
    {
        // Predictor step
        for (int i = 1; i < 40; i++)
            V_intermediate[i] = V[i] - (time_step / space_step) * (Flux(V[i + 1]) - Flux(V[i])) 
                                + alpha * (Flux(V[i + 1]) - 2 * Flux(V[i]) + Flux(V[i - 1]));
        
        // Corrector step
        for (int i = 1; i < 40; i++)
            V_updated[i] = 0.5 * (V[i] + V_intermediate[i] 
                                  - (time_step / space_step) * (Flux(V_intermediate[i]) - Flux(V_intermediate[i - 1])) 
                                  + alpha * (V_intermediate[i + 1] - 2 * V_intermediate[i] + V_intermediate[i - 1]));

        // Calculate error
        error_array[iteration] = 0.0;
        for (int j = 1; j < 40; j++)
            error_array[iteration] += pow(V_updated[j] - V[j], 2);

        error_array[iteration] = sqrt(error_array[iteration]);

        // Update V values
        for (int j = 1; j < 40; j++)
            V[j] = V_updated[j];

        if (error_array[iteration] < tolerance)
            break;
    }

    return iteration;
}

// Initialize boundary conditions
double InitializeBoundary(double V[], double space_step)
{
    for (int i = 0; i < 41; i++)
        V[i] = 0.5 * (1 + tanh(250.0 * (space_step * i - 20.0)));
    
    return 0;
}

int main()
{
    double V_half_dt[41], V_full_dt[41], V_exact[41], space_step, domain_length = 40.0;
    space_step = domain_length / 40.0;

    // Initialize boundary conditions
    InitializeBoundary(V_half_dt, space_step);
    InitializeBoundary(V_full_dt, space_step);
    InitializeBoundary(V_exact, space_step);

    // Perform calculations for time step = 0.5
    int iterations = PredictorCorrector(V_half_dt, 0.5, space_step);
    printf("\nIterations required for time step 0.5: %d", iterations);

    // Perform calculations for time step = 1.0
    iterations = PredictorCorrector(V_full_dt, 1.0, space_step);
    printf("\nIterations required for time step 1.0: %d", iterations);

    // Output the results to a file
    FILE *output_file = fopen("Output_Q3.txt", "w");
    for (int i = 0; i < 41; i++)
        fprintf(output_file, "%lf\t%lf\t%lf\t%lf\n", i * space_step, V_half_dt[i], V_full_dt[i], V_exact[i]);
    fclose(output_file);

    return 0;
}