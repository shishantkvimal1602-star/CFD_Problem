#include <stdio.h>
#include <math.h>

int max_iterations = 1000;
double tolerance = 1e-6;

double FTCS(double V[], double V_new[21])  // Changed from UpwindScheme to FTCS
{
    double coef = 0.4, diffusion = 0.1, error_values[1000];
    int n_iter;

    for(n_iter = 0; n_iter < max_iterations; n_iter++)
    {
        for(int i = 1; i < 20; i++)
            V_new[i] = V[i] - 0.5 * coef * (V[i+1] - V[i-1]) + diffusion * (V[i+1] - 2 * V[i] + V[i-1]);

        error_values[n_iter] = 0.0;
        for(int j = 1; j < 20; j++)
            error_values[n_iter] += pow(V_new[j] - V[j], 2);

        error_values[n_iter] = sqrt(error_values[n_iter]);

        for(int j = 1; j < 20; j++)
            V[j] = V_new[j];

        if(error_values[n_iter] < tolerance)
            break;
    }

    return n_iter;
}

double MacCormack(double V[], double V_new[])  // Changed from PredictorCorrector to MacCormack
{
    double coef = 0.4, diffusion = 0.1, error_values[1000], intermediate[21];
    int n_iter;

    intermediate[0] = V[0];
    intermediate[20] = V[20];

    for(n_iter = 0; n_iter < max_iterations; n_iter++)
    {
        for(int i = 1; i < 20; i++)
            intermediate[i] = V[i] - coef * (V[i+1] - V[i]) + diffusion * (V[i+1] - 2 * V[i] + V[i-1]);

        for(int i = 1; i < 20; i++)
            V_new[i] = 0.5 * (V[i] + intermediate[i] - coef * (intermediate[i] - intermediate[i-1]) + diffusion * (intermediate[i+1] - 2 * intermediate[i] + intermediate[i-1]));

        error_values[n_iter] = 0.0;
        for(int j = 1; j < 20; j++)
            error_values[n_iter] += pow(V_new[j] - V[j], 2);

        error_values[n_iter] = sqrt(error_values[n_iter]);

        for(int j = 1; j < 20; j++)
            V[j] = V_new[j];

        if(error_values[n_iter] < tolerance)
            break;
    }

    return n_iter;
}

double AnalyticalSolution(double V[], double dx)
{
    double coef = 0.4, diffusion = 0.1;

    for(int i = 0; i < 21; i++)
        V[i] = 100 * (1 - exp((coef / (diffusion * dx)) * (dx * i - 1)) / (1 - exp(-(coef / (diffusion * dx)))));
}

double ApplyBoundaryConditions(double V[], double V_new[])
{
    for(int i = 0; i < 21; i++)
    {
        if(i == 0)
        {
            V[i] = 100.0;
            V_new[i] = 100.0;
        }
        else
        {
            V[i] = 0.0;
            V_new[i] = 0.0;
        }
    }
    return 0;
}

int main()
{
    double V_ftcs[21], V_maccormack[21], V_analytical[21], V_new[21];
    double dx = 1.0 / 20.0;

    // Apply boundary conditions and solve using FTCS method
    ApplyBoundaryConditions(V_ftcs, V_new);
    int iterations = FTCS(V_ftcs, V_new);
    printf("\nNumber of iterations taken to converge for FTCS: %d", iterations);

    // Apply boundary conditions and solve using MacCormack method
    ApplyBoundaryConditions(V_maccormack, V_new);
    iterations = MacCormack(V_maccormack, V_new);
    printf("\nNumber of iterations taken to converge for MacCormack: %d", iterations);

    // Analytical solution
    AnalyticalSolution(V_analytical, dx);

    // Write the results to file
    FILE *output_file = fopen("Output_Q2.txt", "w");
    for(int i = 0; i < 21; i++)
        fprintf(output_file, "%lf\t%lf\t%lf\t%lf\n", i * dx, V_ftcs[i], V_maccormack[i], V_analytical[i]);
    fclose(output_file);

    return 0;
}