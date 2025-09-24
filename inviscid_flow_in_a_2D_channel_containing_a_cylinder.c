#include <stdio.h>
#include <math.h>

int maxIterations = 1000;
double tolerance = 1e-6;

int main()
{
    double potential[100][100], potential_new[100][100], residuals[1000];

    // Initialize potential arrays to zero
    for (int i = 0; i < 21; i++)
        for (int j = 0; j < 21; j++)
        {
            potential[i][j] = 0.0;
            potential_new[i][j] = 0.0;
        }

    // Boundary conditions for x direction
    for (int i = 1; i < 21; i++)
    {
        potential[i][0] = 0.0;
        potential[i][20] = 0.1;
        potential_new[i][0] = 0.0;
        potential_new[i][20] = 0.1;
    }

    // Boundary conditions for y direction
    for (int j = 0; j < 21; j++)
    {
        potential[0][j] = 0.1 * j / 20.0;
        potential_new[0][j] = 0.1 * j / 20.0;
    }

    int iteration;
    for (iteration = 0; iteration < maxIterations; iteration++)
    {
        residuals[iteration] = 0.0;

        // Update potential inside the domain
        for (int i = 1; i < 20; i++)
            for (int j = 1; j < 20; j++)
            {
                if (i > 15 && j < 5)
                {
                    if ((i == 16 && j == 3) || (i == 16 && j == 4) || (i == 17 && j == 4))
                        continue;
                }
                potential_new[i][j] = (potential_new[i - 1][j] + potential[i + 1][j] + potential_new[i][j - 1] + potential[i][j + 1]) / 4.0;
            }

        // Calculate residuals and update potential array
        for (int i = 1; i < 21; i++)
            for (int j = 1; j < 21; j++)
                residuals[iteration] += pow(potential_new[i][j] - potential[i][j], 2);

        residuals[iteration] = sqrt(residuals[iteration]);

        for (int j = 5; j < 21; j++)
            potential_new[20][j] = potential_new[19][j];

        for (int i = 1; i < 21; i++)
            for (int j = 1; j < 21; j++)
                potential[i][j] = potential_new[i][j];

        if (residuals[iteration] < tolerance)
            break;
    }

    // Calculate velocities at the cross-section
    double velocities[16];
    for (int i = 4; i < 20; i++)
        velocities[i] = (potential[20][i + 1] - potential[20][i]) / 0.00533;

    // Calculate coefficient of pressure
    double Cp = 1 - pow(velocities[4] / 1.0, 2);
    printf("\nPressure Coefficient: %lf", Cp);

    // Calculate flow rate across the section
    double flowRate = 0.0;
    for (int i = 4; i < 20; i++)
        flowRate += 0.00533 * velocities[i];

    printf("\nFlow Rate: %lf", flowRate);

    // Write potential values to file
    FILE *file_output = fopen("Output_PSI_Q1.txt", "w");
    for (int j = 20; j >= 0; j--)
    {
        for (int i = 0; i < 21; i++)
            fprintf(file_output, "%lf\t", potential[i][j]);
        fprintf(file_output, "\n");
    }
    printf("\nTotal iterations: %d", iteration);
    fclose(file_output);

    // Write residuals to file
    file_output = fopen("Residual.txt", "w");
    for (int i = 0; i < iteration; i++)
        fprintf(file_output, "%d\t%lf\n", i + 1, log(residuals[i]));
    fclose(file_output);

    // Write velocities at the cross-section to file
    file_output = fopen("Velocity_CD.txt", "w");
    for (int i = 5; i < 20; i++)
        fprintf(file_output, "%d\t%lf\n", (i - 5), velocities[i]);
    fclose(file_output);

    return 0;
}
