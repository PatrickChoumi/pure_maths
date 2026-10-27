#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#define CLEAR_SCREEN() system("cls || clear")
#define EPSILON 1e-12

// Prototypes de fonctions
long double** create_matrix(int rows, int cols);
void free_matrix(long double** matrix, int rows);
void read_matrix(long double** matrix, int rows, int cols, const char* name);
void print_matrix(long double** matrix, int rows, int cols, const char* title);
int read_int(const char* prompt);
long double read_long_double(const char* prompt);

// Opérations matricielles
void matrix_addition();
void matrix_multiplication();
void matrix_transpose();
void matrix_determinant();
void matrix_inverse();
void matrix_trace();
void row_echelon_form();
void random_matrix_mode();

// Fonctions mathématiques
long double det(long double** matrix, int n);
long double power(long double x, int n);
long double** inverse(long double** matrix, int n);
long double fabsl(long double x);
long double** matrix_echelon(long double** matrix, int rows, int cols);
long double** swap_rows(long double** matrix, int row1, int row2);
long double cofactor(long double** matrix, int n, int row, int col);

// Gestion de l'interface
void show_main_menu();
void matrix_operations();

/*------------------------------------------------
----------------- IMPLÉMENTATION -----------------
--------------------------------------------------*/

long double** create_matrix(int rows, int cols) {
    long double** matrix = (long double**)malloc(rows * sizeof(long double*));
    if (!matrix) {
        fprintf(stderr, "Erreur d'allocation mémoire\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        matrix[i] = (long double*)malloc(cols * sizeof(long double));
        if (!matrix[i]) {
            fprintf(stderr, "Erreur d'allocation mémoire\n");
            exit(EXIT_FAILURE);
        }
    }
    return matrix;
}

void free_matrix(long double** matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

int read_int(const char* prompt) {
    int value;
    printf("%s", prompt);
    while (scanf("%d", &value) != 1 || value <= 0) {
        printf("Entrée invalide! Réessayez: ");
        while (getchar() != '\n');
    }
    return value;
}

long double read_long_double(const char* prompt) {
    long double value;
    printf("%s", prompt);
    while (scanf("%Lf", &value) != 1) {
        printf("Entrée invalide! Réessayez: ");
        while (getchar() != '\n');
    }
    return value;
}

void read_matrix(long double** matrix, int rows, int cols, const char* name) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%s[%d][%d] = ", name, i+1, j+1);
            matrix[i][j] = read_long_double("");
        }
    }
}

void print_matrix(long double** matrix, int rows, int cols, const char* title) {
    printf("\n%s (%dx%d):\n", title, rows, cols);
    for (int i = 0; i < rows; i++) {
        printf("| ");
        for (int j = 0; j < cols; j++) {
            printf("%12.6Lf ", matrix[i][j]);
        }
        printf("|\n");
    }
    printf("\n");
}

/*----------------- OPÉRATIONS MATHÉMATIQUES -----------------*/

long double** swap_rows(long double** matrix, int row1, int row2) {
    long double* temp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = temp;
    return matrix;
}

long double power(long double x, int n){
    if(n == 0) return 1;
    if(n < 0) return 1/power(x,-n);
    if(n%2 == 0){
        long double temp = power(x,n/2);
        return temp * temp;
    }else{
        return x * power(x,n-1);
    }
}


long double fabsl(long double x){
    if(x >= 0){
        return x;
    }else{
        return -x;
    }
}

long double** matrix_echelon(long double** matrix, int rows, int cols) {
    int pivot = 0;
    
    for (int col = 0; col < cols && pivot < rows; col++) {
        int max_row = pivot;
        for (int i = pivot + 1; i < rows; i++) {
            if (fabsl(matrix[i][col]) > fabsl(matrix[max_row][col])) {
                max_row = i;
            }
        }

        if (fabsl(matrix[max_row][col]) < EPSILON) continue;

        if (max_row != pivot) {
            swap_rows(matrix, pivot, max_row);
        }

        for (int i = pivot + 1; i < rows; i++) {
            long double factor = matrix[i][col] / matrix[pivot][col];
            for (int j = col; j < cols; j++) {
                matrix[i][j] -= factor * matrix[pivot][j];
            }
        }
        pivot++;
    }
    return matrix;
}

long double det(long double** matrix, int n) {
    long double determinant = 1.0;
    int sign = 1;
    long double** copy = create_matrix(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            copy[i][j] = matrix[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            if (fabsl(copy[k][i]) > fabsl(copy[max_row][i])) {
                max_row = k;
            }
        }

        if (max_row != i) {
            swap_rows(copy, i, max_row);
            sign *= -1;
        }

        if (fabsl(copy[i][i]) < EPSILON) {
            determinant = 0.0;
            break;
        }

        determinant *= copy[i][i];

        for (int k = i + 1; k < n; k++) {
            long double factor = copy[k][i] / copy[i][i];
            for (int j = i; j < n; j++) {
                copy[k][j] -= factor * copy[i][j];
            }
        }
    }

    free_matrix(copy, n);
    return sign * determinant;
}

long double cofactor(long double** matrix, int n, int row, int col) {
    long double** minor = create_matrix(n-1, n-1);
    int mi = 0, mj = 0;
    
    for (int i = 0; i < n; i++) {
        if (i == row) continue;
        mj = 0;
        for (int j = 0; j < n; j++) {
            if (j == col) continue;
            minor[mi][mj++] = matrix[i][j];
        }
        mi++;
    }
    
    long double result = det(minor, n-1) * power(-1, row+col);
    free_matrix(minor, n-1);
    return result;
}

long double** inverse(long double** matrix, int n) {
    long double determinant = det(matrix, n);
    if (fabsl(determinant) < EPSILON) return NULL;
    
    long double** adjugate = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            adjugate[j][i] = cofactor(matrix, n, i, j) / determinant;
        }
    }
    return adjugate;
}

/*----------------- OPÉRATIONS MATRICIELLES -----------------*/

void matrix_addition() {
    CLEAR_SCREEN();
    printf("\t[1] ADDITION DE MATRICES\n");
    
    int rows = read_int("Nombre de lignes : ");
    int cols = read_int("Nombre de colonnes : ");
    
    long double** A = create_matrix(rows, cols);
    read_matrix(A, rows, cols, "A");
    
    long double** B = create_matrix(rows, cols);
    read_matrix(B, rows, cols, "B");
    
    long double** result = create_matrix(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    
    print_matrix(result, rows, cols, "Résultat de l'addition");
    
    free_matrix(A, rows);
    free_matrix(B, rows);
    free_matrix(result, rows);
}

void matrix_multiplication() {
    CLEAR_SCREEN();
    printf("\t[2] MULTIPLICATION DE MATRICES\n");
    
    int n = read_int("Lignes de la première matrice (n) : ");
    int m = read_int("Colonnes de la première matrice (m) : ");
    int p = read_int("Colonnes de la deuxième matrice (p) : ");
    
    long double** A = create_matrix(n, m);
    read_matrix(A, n, m, "A");
    
    long double** B = create_matrix(m, p);
    read_matrix(B, m, p, "B");
    
    long double** C = create_matrix(n, p);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < m; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    
    print_matrix(C, n, p, "Résultat de la multiplication");
    
    free_matrix(A, n);
    free_matrix(B, m);
    free_matrix(C, n);
}

void matrix_transpose() {
    CLEAR_SCREEN();
    printf("\t[5] TRANSPOSÉE DE MATRICE\n");
    
    int rows = read_int("Nombre de lignes : ");
    int cols = read_int("Nombre de colonnes : ");
    
    long double** matrix = create_matrix(rows, cols);
    read_matrix(matrix, rows, cols, "A");
    
    long double** transpose = create_matrix(cols, rows);
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            transpose[i][j] = matrix[j][i];
        }
    }
    
    print_matrix(transpose, cols, rows, "Matrice transposée");
    
    free_matrix(matrix, rows);
    free_matrix(transpose, cols);
}

void matrix_determinant() {
    CLEAR_SCREEN();
    printf("\t[3] DÉTERMINANT\n");
    
    int n = read_int("Taille de la matrice carrée (n) : ");
    long double** matrix = create_matrix(n, n);
    read_matrix(matrix, n, n, "A");
    
    long double determinant = det(matrix, n);
    printf("\nDéterminant : %.4Lf\n", determinant);
    
    free_matrix(matrix, n);
}

void matrix_inverse() {
    CLEAR_SCREEN();
    printf("\t[4] INVERSE DE MATRICE\n");
    
    int n = read_int("Taille de la matrice carrée (n) : ");
    long double** matrix = create_matrix(n, n);
    read_matrix(matrix, n, n, "A");
    
    long double** inv = inverse(matrix, n);
    if (!inv) {
        printf("Matrice non inversible (déterminant nul)!\n");
    } else {
        print_matrix(inv, n, n, "Matrice inverse");
        free_matrix(inv, n);
    }
    
    free_matrix(matrix, n);
}

void matrix_trace() {
    CLEAR_SCREEN();
    printf("\t[6] TRACE DE MATRICE\n");
    
    int rows = read_int("Nombre de lignes : ");
    int cols = read_int("Nombre de colonnes : ");
    int dim = rows < cols ? rows : cols;
    
    long double** matrix = create_matrix(rows, cols);
    read_matrix(matrix, rows, cols, "A");
    
    long double trace = 0.0;
    for (int i = 0; i < dim; i++) {
        trace += matrix[i][i];
    }
    
    printf("\nTrace : %.4Lf\n", trace);
    free_matrix(matrix, rows);
}

void row_echelon_form() {
    CLEAR_SCREEN();
    printf("\t[7] FORME ÉCHELONNÉE\n");
    
    int rows = read_int("Nombre de lignes : ");
    int cols = read_int("Nombre de colonnes : ");
    
    long double** matrix = create_matrix(rows, cols);
    read_matrix(matrix, rows, cols, "A");
    
    long double** echelon_matrix = create_matrix(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            echelon_matrix[i][j] = matrix[i][j];
        }
    }
    
    matrix_echelon(echelon_matrix, rows, cols);
    print_matrix(echelon_matrix, rows, cols, "Forme échelonnée");
    
    free_matrix(matrix, rows);
    free_matrix(echelon_matrix, rows);
}

void random_matrix_mode() {
    CLEAR_SCREEN();
    printf("\t[8] MODE ALÉATOIRE\n");
    
    int n = read_int("Taille de la matrice carrée (n) : ");
    long double** matrix = create_matrix(n, n);
    
    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = (rand() % 2000 - 1000) / 10.0L;
        }
    }
    
    print_matrix(matrix, n, n, "Matrice aléatoire");
    printf("Déterminant : %.4Lf\n", det(matrix, n));
    long double** inv = inverse(matrix, n);
    if (!inv) {
        printf("Matrice non inversible (déterminant nul)!\n");
    } else {
        print_matrix(inv, n, n, "Matrice inverse");
        free_matrix(inv, n);
    }
    free_matrix(matrix, n);
}