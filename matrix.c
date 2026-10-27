#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void clear_screen(){
    #ifdef _WIN32
        system("cls");
    #else
        system("clear");
    #endif
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

//fonction pour creer une matrice carree
long double** create_matrix(int n, int m){
    long double** T = (long double**)malloc(n*sizeof(long double*));
    
    for(int i = 0; i<n; i++){
        T[i] = (long double*)malloc(m*sizeof(long double));
    }
    return T;
    
}

void print_matrix(long double** A, int n, int m){
    for(int i=0; i<n; i++){
        printf("\t");
        for(int j=0; j<m; j++){
            printf("A[%d][%d] = %.10Lf\t\t",i+1,j+1,A[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");    
}

//fonction pour liberer la memoire occupee par la matrice
void free_matrix(long double** T, int n){
    for(int i = 0; i<n; i++){
        free(T[i]);
    }
    free(T);
}

//fonction additionner deux matrices
void addition_matrix(){
    clear_screen();
    int n,m;
    printf("\t[1]. Addition de deux matrices\n");
    printf("Entrez la  taille n x m des matrices\n");
    printf("n = ");
    scanf("%d",&n);
    printf("m = ");
    scanf("%d",&m);
    printf("\tDonnees de la premiere matrice\n\n");

    long double** A = create_matrix(n,m);
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            printf("\tA[%d][%d] = ",i+1,j+1);
            scanf("%Lf",&A[i][j]);
        }
    }
    printf("\n\n");

    print_matrix(A,n,m);
    printf("\tDonnees de la deuxieme matrice\n\n");

    long double** B = create_matrix(n,m);
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            printf("\tB[%d][%d] = ",i+1,j+1);
            scanf("%Lf",&B[i][j]);
        }
    }
    printf("\n\n");
    for(int i=0; i<n; i++){
        printf("\t");
        for(int j=0; j<m; j++){
            printf("B[%d][%d] = %.10Lf\t\t",i+1,j+1,B[i][j]);
        }
        printf("\n");
    }
    printf("\nResultat : \n\n");
    long double** C = create_matrix(n,m);
    for(int i=0; i<n; i++){
        printf("\t");
        for(int j=0; j<m; j++){
            C[i][j] = A[i][j] + B[i][j];
            printf("C[%d][%d] = %.10Lf\t\t",i+1,j+1,C[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    free_matrix(A,n);
    free_matrix(B,n);
    free_matrix(C,n);
}

void multi_matrix(){
    clear_screen();
    int n,m,p;
    printf("\t[2]. Multiplication de deux matrices\n");
    printf("Entrez la  taille n x m de la premiere matrice \n");
    printf("n = ");
    scanf("%d",&n);
    printf("m = ");
    scanf("%d",&m);
    printf("\tDonnees de la premiere matrice\n\n");

    long double** A = create_matrix(n,m);
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            printf("\tA[%d][%d] = ",i+1,j+1);
            scanf("%Lf",&A[i][j]);
        }
    }
    printf("\n\n");

    print_matrix(A,n,m);
    printf("Entrez la  taille m x p de la deuxieme matrice \n");
    printf("m = %d\n",m);
    printf("p = ");
    scanf("%d",&p);
    printf("\tDonnees de la deuxieme matrice\n\n");

    long double** B = create_matrix(m,p);
    for(int i=0; i<m; i++){
        for(int j=0; j<p; j++){
            printf("\tB[%d][%d] = ",i+1,j+1);
            scanf("%Lf",&B[i][j]);
        }
    }
    printf("\n\n");
    for(int i=0; i<m; i++){
        printf("\t");
        for(int j=0; j<p; j++){
            printf("B[%d][%d] = %.10Lf\t\t",i+1,j+1,B[i][j]);
        }
        printf("\n");
    }
    printf("\nResultat : \n\n");
    long double** C = create_matrix(n,p);
    for(int i=0; i<n; i++){
        printf("\t");
        for(int j=0; j<p; j++){
            C[i][j] = 0;
            for(int k=0; k<m; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
            printf("C[%d][%d] = %.10Lf\t\t",i+1,j+1,C[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    free_matrix(A,n);
    free_matrix(B,n);
    free_matrix(C,n);
}

//fonction pour echanger des lignes lors de l'echelonnement
long double** echanger_lignes(long double** T, int a, int b){
    long double* temp = T[a];
    T[a] = T[b];
    T[b] = temp;
    return T;
}

//fonction pour echelonner la matrice
long double** echelon(long double** T, int n, int m){
    for(int i = 0; i<n ; i++){
        if(T[i][i] != 0){
            for(int k=i+1; k<n; k++){
                long double temp = T[k][i];
                for(int j=0; j<m; j++){
                    T[k][j] -= (temp/T[i][i])*T[i][j];
                }
            }
        }else{
            for(int k = i+1; k<n; k++){
                if(T[k][i] != 0){
                    T = echanger_lignes(T,i,k);
                    i--;
                }
                break;
            }
        }
    }
    return T;
}

// Calcule effectivement le determinant
long double det(long double** T, int n){
    switch (n)
    {
    case 1:
        return T[0][0];
        break;
    case 2:
        return T[0][0]*T[1][1] - T[0][1]*T[1][0];
        break;
    default:
        T = echelon(T,n,n);
        int temp = 1;
        for(int i = 0; i<n; i++){
            temp *= T[i][i];
        }
        return temp;
        break;
    }
}

void determinant_matrix(){
    clear_screen();
    int n;
    printf("\t[3]. Determinant d'une matrice\n");
    printf("Entrez la  taille n de la matrice carre \n");
    printf("n = ");
    scanf("%d",&n);
    printf("\tDonnees de la matrice\n\n");

    long double** A = create_matrix(n,n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("\tA[%d][%d] = ",i+1,j+1);
            scanf("%Lf",&A[i][j]);
        }
    }
    printf("\n\n");
    printf("\tMatrice actuelle\n\n");

    print_matrix(A,n,n);

    long double** B = echelon(A,n,n);
    printf("\tMatrice echelonnee\n\n");
    for(int i=0; i<n; i++){
        printf("\t");
        for(int j=0; j<n; j++){
            printf("A[%d][%d] = %.10Lf\t\t",i+1,j+1,B[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    printf("\nResultat : \n\n");

    printf("\tdet(A) = %.10Lf",det(A,n));
    //petit test de la fonction echelon avant de l'implementer
    /*A=echelon(A,n,n);
    for(int i=0; i<n; i++){
    printf("\t");
        for(int j=0; j<n; j++){
            printf("A[%d][%d] = %.10Lf\t\t",i+1,j+1,A[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");*/

    printf("\n\n");
    free_matrix(A,n);
}

long double cofacteur(long double** T,int n, int i, int j){
    long double** temp = create_matrix(n-1, n-1);
    int x = 0, y = 0;
    for(int a = 0; a < n; a++){
        for(int b = 0; b < n; b++){
            if(a != i && b != j){
                temp[x][y++] = T[a][b];
                if(y == n - 1){
                    y = 0;
                    x++;
                }
            } 
        }
    }
    return det(temp,n-1) * power(-1,i+j);
}

long double** inverse(long double** T,int n){
   
    long double** C = create_matrix(n,n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            C[i][j] =  cofacteur(T,n,j,i);

        }
    }
    long double t = det(T,n);
    if(t == 0) return NULL;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            C[i][j] *= 1/t;

        }
    }
    return C;
}

void inverse_matrix(){
    clear_screen();

    int n;
    printf("\t[4]. Inverse d'une matrice\n");
    printf("Entrez la  taille n de la matrice carre \n");
    printf("n = ");
    scanf("%d",&n);
    printf("\tDonnees de la matrice\n\n");

    long double** A = create_matrix(n,n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("\tA[%d][%d] = ",i+1,j+1);
            scanf("%Lf",&A[i][j]);
        }
    }
    printf("\n\n");
    printf("\tMatrice actuelle\n\n");

    for(int i=0; i<n; i++){
        printf("\t");
        for(int j=0; j<n; j++){
            printf("A[%d][%d] = %.10Lf\t\t",i+1,j+1,A[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");

    A = inverse(A,n);

    printf("\nResultat : \n\n");

    if(A == NULL) {
        printf("Matrice non inversible car son determinant est nul !!!\n\n");
        return;
    }

    for(int i=0; i<n; i++){
        printf("\t");
        for(int j=0; j<n; j++){
            printf("A[%d][%d] = %.10Lf\t\t",i+1,j+1,A[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    free_matrix(A,n);
}



void trace_matrix(){
    clear_screen();
    int n,m;
    printf("\t[6]. Trace d'une matrice\n");
    printf("Entrez la  taille n x m d'une matrice \n");
    printf("n = ");
    scanf("%d",&n);
    printf("m = ");
    scanf("%d",&m);
    printf("\tDonnees de la matrice\n\n");

    long double** A = create_matrix(n,m);
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            printf("\tA[%d][%d] = ",i+1,j+1);
            scanf("%Lf",&A[i][j]);
        }
    }
    printf("\n\n");
    printf("\tMatrice actuelle\n\n");

    print_matrix(A,n,m);

    int a = n <= m ? n : m;
    long double tr = 0;
    for(int i = 0; i< a; i++){
        tr += A[i][i];
    }
    
    printf("\nResultat : \n\n");

    printf("\tTr(A) = %.10Lf",tr);
    

    printf("\n\n");
    free_matrix(A,n);
}

void echelon_matrix(){
    clear_screen();
    int n,m;
    printf("\t[7]. Echelonnee d'une matrice\n");
    printf("Entrez la  taille n x m des matrices\n");
    printf("n = ");
    scanf("%d",&n);
    printf("m = ");
    scanf("%d",&m);
    printf("\tDonnees de la premiere matrice\n\n");

    long double** A = create_matrix(n,m);
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            printf("\tA[%d][%d] = ",i+1,j+1);
            scanf("%Lf",&A[i][j]);
        }
    }
    printf("\n\n");

    printf("\tMatrice actuelle\n\n");

    print_matrix(A,n,m);

    A = echelon(A,n,m);

    printf("\nResultat : \n\n");

    
    print_matrix(A,n,m);
    free_matrix(A,n);
}

void matrice_aleatoire(long double** A, int n, int m){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            int t = 1 + rand()%100;
            A[i][j] = (-500 + rand()%1000)/t;
        }
    }
}

void mode_aleatoire(){
    clear_screen();

    int n;
    srand(time(NULL));
    printf("\t[18037]. Mode Aleatoire\n");
    printf("Entrez la  taille n de la matrice carre \n");
    printf("n = ");
    scanf("%d",&n);
    printf("\tDonnees de la matrice\n\n");

    long double** A = create_matrix(n,n);
    matrice_aleatoire(A,n,n);
    printf("\n\n");
    printf("\tMatrice actuelle\n\n");

    for(int i=0; i<n; i++){
        printf("|");
        for(int j=0; j<n; j++){
            printf("  %.3Lf",A[i][j]);
        }
        printf("|\n");
    }
    printf("\n\n");

    printf("\tdet(A) = %.25Lf",det(A,n));

    A = inverse(A,n);

    printf("\nResultat : \n\n");

    if(A == NULL) {
        printf("Matrice non inversible car son determinant est nul !!!\n\n");
        return;
    }

    /*for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("A[%d][%d] = %.25Lf  \n",i+1,j+1,A[i][j]);
        }
    }*/
    print_matrix(A,n,n);
    printf("\n\n");
    free_matrix(A,n);
}

//menu
void menu_matrix(){
    int option_matrix;
    clear_screen();
    printf("\t\tMenu principal\n\n");
    printf("\t1. Addition de deux matrices\n");
    printf("\t2. Multiplication de deux matrices\n");
    printf("\t3. Determinant d'une matrice\n");
    printf("\t4. Inverse d'une matrice\n");
    printf("\t5. Transposee d'une matrice\n");
    printf("\t6. Trace d'une matrice\n");
    printf("\t7. Echelonnee d'une matrice\n");
    printf("\t8. Quitter le programme\n");
    printf("\n\n\t18037.Mode Aleatoire \n");

    printf("\tChoix = ");
    scanf("%d",&option_matrix);

    switch (option_matrix)
    {
    case 1:
        addition_matrix();
        break;
    case 2:
        multi_matrix();
        break;

    case 3:
        determinant_matrix();
        break;
    
    case 4:
        inverse_matrix();
        break;

    case 5:
        //inverse_matrix();
        break;

    case 6:
        trace_matrix();
        break;
    case 7:
        echelon_matrix();
        break;
    case 18037:
        mode_aleatoire();
        break;
    default:
        break;
    }
}

void matrix(){
    int choix_matrix = 1;
    
    while (choix_matrix == 1)
    {
        menu_matrix();
        printf("Continuer?\n");
        printf("[1]. Retour au menu principal\n");
        printf("[Autre valeur]. Quitter le programme\n");
        printf("Choix = ");
        scanf("%d",&choix_matrix);
    }
    printf("Exit...:) See you soon\n\n");   
}

int main(){
    matrix();
    return 0;
}