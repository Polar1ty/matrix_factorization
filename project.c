// Program implements different matrix decompositions
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void header(void)
{
    printf("###-------------------------------------------###\n");
    printf("                  Course project\n");
    printf("              Programming in C language\n");
    printf("           Made by: Alexey Baida KM - 01\n");
    printf("");
    printf("###-------------------------------------------###\n");
    printf("Task: Implement the following types of matrix decomposition\n");
    printf("    - LU decomposition\n");
    printf("    - Rank factorization\n");
    printf("    - Kholetsky's decomposition\n");
    printf("    - QR decomposition\n\n");
}

void print_notation(void)
{
    printf("\nBefore entering the matrix note that:\n");
    printf("    - LU decomposition requires square matrix\n");
    printf("    - Kholetsky's decomposition requires symmetric non-negative matrix\n\n");
    printf("Enter the matrix manually or generate randomly?\n0 - randomly\n1 - manually\n");
}

// Log user moves
void log_print(char arr[])
{
    FILE *f;
    f = fopen("user_actions.log", "a+"); // a+ (create + append)
    if (f == NULL) 
    {
        printf("Log file was not created\n");
    }
    fprintf(f, "%s", arr);
    fclose(f);
}

// Check correctness of entered integer, returns valid integer
int check_int(void)
{
    int num;
    char term;
    while (1)
    {
        if (scanf("%d%c", &num, &term) != 2 || term != '\n')
        {
            printf("Invalud input. Please enter an integer value\n");
            int c;
            while (('\n' != (c=fgetc(stdin))) && (c != EOF));   // clear up to end of line
        }
        else
        {
            break;
        }
    }
    return num;
}

// User enter a nrow*ncol matrix
void enter_manually(int nrows, int ncols, float matrix[][ncols])
{
    int i, j;
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            printf("Enter [%d][%d] element of matrix: ", i + 1, j + 1);
            scanf("%f", &matrix[i][j]);
        }
    }
}

// Generates nrow*ncols matrix
void generate_randomly(int nrows, int ncols, float matrix[][ncols])
{
    int i, j;
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            //Generates pseudorandom float between 0.1 and 9.9
            matrix[i][j] = (float) rand() / RAND_MAX * (9.9 - 0.1) + 0.1;
        }
    }
}

// Prints given matrix
void show_matrix(int nrows, int ncols, float matrix[][ncols])
{
    int i, j;

    for (i = 0; i < nrows; i++)
    {
        printf("  |");
        for (j = 0; j < ncols; j++)
        {
            if (j == ncols - 1)
            {
                printf("%f", matrix[i][j]);
            }
            else
            {
                printf("%f\t", matrix[i][j]);
            }
        }
        printf("|\n");
    }
}

// Impelements LU decomposition, outputs results
void lu_decomposition(int nrows, int ncols, float matrix[][ncols])
{
    // L - Lower triangle matrix
    float L[nrows][ncols];
    // U - Upper triangle matrix
    float U[nrows][ncols];
    // Iterators for cycles
    int i, j, k;

    printf("\nLU decomposition:  A = L * U (where A - entered matrix, L - lower triangle matrix, U - upper triangle matrix)\n");

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            U[i][j] = 0;
            L[i][j] = 0;
            L[i][i] = 1;
        }
    }

    // Variables to store summations
    float s, s1;
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            if (i <= j)
            {
                s = 0;
                for (k = 0; k <= i - 1; k++)
                {
                    s += L[i][k] * U[k][j];
                }
                U[i][j] = matrix[i][j] - s;
            }
            else if (i > j)
            {
                s1 = 0;
                for (k = 0; k <= j - 1; k++)
                {
                    s1 += L[i][k] * U[k][j];
                }
                L[i][j] = (matrix[i][j] - s1) / U[j][j];
            }
        }
    }

    // Show Lower matrix
    printf("L matrix:\n");
    show_matrix(nrows, ncols, L);

    // Show Upper matrix
    printf("U matrix:\n");
    show_matrix(nrows, ncols, U);
}

// Impelements Kholetsky's decomposition, outputs results
void Kholetskys_decomposition(int nrows, int ncols, float matrix[][ncols])
{
    // Iterators for cycles
    int i, j, p;
    // L - Lower triangle matrix
    float L[nrows][ncols];
    // L_transpose - transposed matrix L
    float L_transpose[nrows][ncols];
    // Variables to store summations
    float s, s1;

    printf("\nKholetskys_decomposition:  A = L * L_transpose (where A - entered matrix, L - lower triangle matrix)\n");

    L[0][0] = sqrt(matrix[0][0]);
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            s = 0;
            for (p = 0; p <= i - 1; p++)
            {
                s += L[i][p] * L[i][p];
            }
            L[i][i] = sqrt(matrix[i][i] - s);

            if (j < i)
            {
                s1 = 0;
                for (p = 0; p <= j - 1; p++)
                {
                    s1 += L[i][p] * L[j][p];
                }
                L[i][j] = (matrix[i][j] - s1) / L[j][j];
            }
        }
    }

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            L_transpose[j][i] = L[i][j];
        }
    }

    printf("L matrix:\n");
    show_matrix(nrows, ncols, L);
    printf("L transposed:\n");
    show_matrix(nrows, ncols, L_transpose);
}

void calculate_U(int nrows, int ncols, float matrix[][ncols], float U[][ncols])
{
    int i, j, k, a;
    /* Calculates U matrix */

    // Clone 1-st vector of A matrix into U
    for (i = 0; i < nrows; i++)
    {
        U[i][0] = matrix[i][0];
    }

    for (k = 1; k < ncols; k++)
    {
        // Initialize sum vector
        float summ[nrows];
        for (a = 0; a < nrows; a++)
        {
            summ[a] = 0;
        }
        // Calculates sum of projections
        for (j = 0; j <= k - 1; j++)
        {
            float sk, sk1, proj[nrows];
            //proj uj xk = (xk,uj)/(uj,uj)*uj
            sk = 0;
            for (a = 0; a < nrows; a++)
            {
                sk += matrix[a][k] * U[a][j];
            }
            sk1 = 0;
            for (a = 0; a < nrows; a++)
            {
                sk1 += U[a][j] * U[a][j];
            }
            for (a = 0; a < nrows; a++)
            {
                proj[a] = (sk / sk1) * U[a][j];
            }
            for (a = 0; a < nrows; a++)
            {
                summ[a] += proj[a];
            }
        }
        // Performs Ak vector minus sum of projections
        for (a = 0; a < nrows; a++)
        {
            U[a][k] = matrix[a][k] - summ[a];
        }
    }
    /* End of calculation U matrix */
}

// Impelements QR decomposition, outputs results
void qr_decomposition(int nrows, int ncols, float matrix[][ncols], float U[][ncols])
{
    int i, j, a, k;
    float Q[nrows][ncols];
    float Q_transpose[nrows][ncols];
    float R[nrows][ncols];

    printf("\nQR - decomposition:  A = Q * R (where A - entered matrix, Q - orthogonal matrix, R - upper trianle matrix)\n");


    /* Calculates Q matrix (based on U matrix) */
    for (k = 0; k < ncols; k++)
    {
        float sum_of_squares, norm;
        sum_of_squares = 0;
        for (a = 0; a < nrows; a++)
        {
            sum_of_squares += U[a][k] * U[a][k];
        }
        norm = sqrt(sum_of_squares);

        for (a = 0; a < nrows; a++)
        {
            Q[a][k] = U[a][k] / norm;
        }
    }

    /* Calculates Q_transposed*/
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            Q_transpose[j][i] = Q[i][j];
        }
    }

    /* Calculates R matrix as product of Q_transposed and A(entered matrix) */
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            R[i][j] = 0;
            for (k = 0; k < nrows; k++)
            {
                R[i][j] += Q_transpose[i][k] * matrix[k][j];
            }
        }
    }

    printf("Q matrix:\n");
    show_matrix(nrows, ncols, Q);
    printf("R matrix:\n");
    show_matrix(nrows, ncols, R);
}

// decomposition menu
int decomposition_menu(int nrows, int ncols, float matrix[][ncols])
{
    int answer;
    answer = 5;
    printf("\nChoose the decomposition method you want to test:\n1 - LU decomposition\n");
    printf("2 - Kholetsky's decomposition\n3 - QR decomposition\n4 - exit to main menu\n");
    answer = check_int();
    if (answer == 1)
    {
        // Log and perform LU factorization
        lu_decomposition(nrows, ncols, matrix);
        log_print("User performed LU decomposition\n");
    }
    else if (answer == 2)
    {
        // Log and perform Kholetsky's factorization
        Kholetskys_decomposition(nrows, ncols, matrix);
        log_print("User performed Kholetsky's decomposition\n");
    }
    else if (answer == 3)
    {
        // Log and perform QR factorization
        float U[nrows][ncols];
        calculate_U(nrows, ncols, matrix, U);
        qr_decomposition(nrows, ncols, matrix, U);
        log_print("User performed QR decomposition\n");
    }
    return answer;
}

int main(void)
{
    int answer;
    answer = 5;  // any non zero value :)

    system("cls");
    header();  // print header

    log_print("\n---START OF LOGGING---\n");
    while (1)
    {
        printf("Enter 1 - for start testing the program, Any number - for exit:\n");
        answer = check_int();
        if (answer == 1)
        {
            log_print("User entered program\n");
            while (1)
            {
                print_notation();
                answer = check_int();

                int nrows, ncols, i, j;
                printf("Enter number of rows:\n");
                nrows = check_int();
                printf("Enter number of columns:\n");
                ncols = check_int();
                float matrix[nrows][ncols];

                // Generate randomly
                if (answer == 0)
                {
                    generate_randomly(nrows, ncols, matrix);
                }
                // Enter manually
                else if (answer == 1)
                {
                    enter_manually(nrows, ncols, matrix);
                }

                log_print("User entered matrix correct\n");
                // Show entered matrix
                printf("Entered matrix:\n");
                show_matrix(nrows, ncols, matrix);

                answer = decomposition_menu(nrows, ncols, matrix);
                if (answer == 4)
                {
                    // Log and perform exit to main menu
                    log_print("User exited to main menu\n");
                    break;
                }
                printf("Test another matrix?\n0 - Yes\n1 - No (return to main menu)\n");
                answer = check_int();
                if (answer == 1)
                {
                    log_print("User exited to main menu\n");
                    break;
                }
            }
        }
        else
        {
            log_print("User exited program\n");
            break;
        }
    }
    printf("---Thanks for using!---\n");
    log_print("---END OF LOGGING---\n");
}