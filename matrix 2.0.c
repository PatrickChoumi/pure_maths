#include "pure_maths.h"


void show_main_menu() {
    CLEAR_SCREEN();
    printf("\tCALCULATRICE MATRICIELLE\n\n");
    printf("\t1. Addition de matrices\n");
    printf("\t2. Multiplication de matrices\n");
    printf("\t3. Déterminant\n");
    printf("\t4. Matrice inverse\n");
    printf("\t5. Transposée\n");
    printf("\t6. Trace\n");
    printf("\t7. Forme échelonnée\n");
    printf("\t8. Mode aléatoire\n");
    printf("\t9. Quitter\n");
}

void matrix_operations() {
    int choice;
    do {
        show_main_menu();
        choice = read_int("\nVotre choix : ");
        
        switch(choice) {
            case 1: matrix_addition(); break;
            case 2: matrix_multiplication(); break;
            case 3: matrix_determinant(); break;
            case 4: matrix_inverse(); break;
            case 5: matrix_transpose(); break;
            case 6: matrix_trace(); break;
            case 7: row_echelon_form(); break;
            case 8: random_matrix_mode(); break;
            case 9: 
                printf("Exit...:) See you soon\n\n"); 
                exit(0);
            default: printf("Option invalide!\n");
        }
        
        if(choice != 9) {
            printf("\nAppuyez sur Entrée pour continuer...");
            while(getchar() != '\n');
            getchar();
        }
    } while(choice != 9);
}

int main() {
    srand(time(NULL));
    matrix_operations();
    return 0;
}