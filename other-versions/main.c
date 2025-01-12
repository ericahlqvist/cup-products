
#include <pari/pari.h>
#include <stdio.h>
#include <time.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#include "misc_functions.h"
#include "tests.h"
#include "ext_and_aut.h"
#include "artin_symbol.h"
#include "find_cup_matrix.h"
#include "test_artin.h"


int
main (int argc, char *argv[])	  
{
    //--------
    clock_t start = clock();
    //--------
    

    int p_int, my_int, p_rk, i;

    int min;
    int sec;
    int msec;
    
    pari_init(4000000000,500000);
    // printf("Initial adress: %ld\n", avma);
    // pari_sp limit = stack_lim(avma, 1);
    
    GEN p = cgeti(DEFAULTPREC);
    GEN s = pol_x(fetch_user_var("s"));
    GEN K, f, Kcyc, p_ClFld_pol, D_prime_vect, D, J_vect, Ja_vect;

    p = gp_read_str(argv[1]);
    p_int = atoi(argv[1]);
    my_int = atoi(argv[2]);

    D = stoi(my_int);
    D_prime_vect = gel(factor(D), 1);

    // Define K.pol
    f = gsubgs(gsqr(s), my_int);
    printf("\n");

    // Define base field K
    K = Buchall(f, nf_FORCE, DEFAULTPREC);
    //tu = bnf_get_tuU(K);
    Kcyc = bnf_get_cyc(K);
    p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);

    J_vect = my_find_p_gens(K, p);
    pari_printf("J_vect: %Ps\n\n", J_vect);
    Ja_vect = my_find_Ja_vect(K, J_vect);
    pari_printf("Ja_vect: %Ps\n\n", Ja_vect);
    // if (my_is_p_torsion(K,J_vect, p))
    // {
    //     printf(ANSI_COLOR_GREEN "p-torsion test passed\n\n" ANSI_COLOR_RESET);
    // }
    // else {
    //     printf(ANSI_COLOR_RED "p-torsion test  failed\n\n" ANSI_COLOR_RESET);
    // }
    

    p_rk = glength(J_vect)-1;
    printf("p-rank: %d\n\n", p_rk);


    pari_printf("p_int: %d\n\nmy_int: %d\n\nK_cyc: %Ps\n\nK_basis: %Ps\n\n", p_int, my_int, Kcyc, nf_get_zk(bnf_get_nf(K)));

    GEN K_ext = my_ext(K, p_ClFld_pol, my_int, s, p, D_prime_vect, p_rk);

    
    //Defines a matrix over F_2 with index (i*k, j) corresponding to 
    //< x_i\cup x_k, (a_j, J_j)>
    GEN cup_matrix = my_cup_matrix(K_ext, K, p, p_int, p_rk, J_vect, Ja_vect);
    int k;
    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (i=1; i<p_rk+1; ++i) {
        for (k=1; k<p_rk+1; ++k) {
            pari_printf(ANSI_COLOR_CYAN "(%d,%d)  :  %Ps\n\n" ANSI_COLOR_RESET, i,k, gel(cup_matrix, (i-1)*p_rk+k));
        } 
    }
    printf("\n\n");

    // int cup_det_12 = smodis(gsub(gmul(gmael2(cup_matrix, 1,1), gmael2(cup_matrix, 2,2)), gmul(gmael2(cup_matrix, 1,2), gmael2(cup_matrix, 2,1))), p_int);
    // int cup_det_13 = smodis(gsub(gmul(gmael2(cup_matrix, 1,1), gmael2(cup_matrix, 3,2)), gmul(gmael2(cup_matrix, 3,1), gmael2(cup_matrix, 1,2))), p_int);
    // int cup_det_23 = smodis(gsub(gmul(gmael2(cup_matrix, 2,1), gmael2(cup_matrix, 3,2)), gmul(gmael2(cup_matrix, 3,1), gmael2(cup_matrix, 2,2))), p_int);

    // FILE *fptr;
    // fptr = fopen(file_name, "w");

    // if (my_QV_equal0(gel(cup_matrix, 1)) && my_QV_equal0(gel(cup_matrix, 2)) && my_QV_equal0(gel(cup_matrix, 3)))
    // {
    //     printf(ANSI_COLOR_GREEN "Rank 0 \n\n" ANSI_COLOR_RESET);
    //     pari_fprintf(fptr, "{\"p\": \"%d\", \"D\": \"%Ps\", \"K-cyc\": \"%Ps\", \"Lx-cyc\": \"%Ps\", \"Ly-cyc\": \"%Ps\", \"M-rk\": \"0\", \"CM\": \"%Ps\"},\n", cup_matrix);
    // }
    // else if (cup_det_12==1 || cup_det_13==1 || cup_det_23==1)
    // {
    //     printf(ANSI_COLOR_YELLOW "Rank 2 \n\n" ANSI_COLOR_RESET);
    //     pari_fprintf(fptr, "{\"p\": \"%d\", \"D\": \"%Ps\", \"K-cyc\": \"%Ps\", \"Lx-cyc\": \"%Ps\", \"Ly-cyc\": \"%Ps\", \"M-rk\": \"2\", \"CM\": \"%Ps\"},\n", cup_matrix);
    // }
    // else {
    //     printf(ANSI_COLOR_YELLOW "Rank 1 \n\n" ANSI_COLOR_RESET);
    //     pari_fprintf(fptr, "{\"p\": \"%d\", \"D\": \"%Ps\", \"K-cyc\": \"%Ps\", \"Lx-cyc\": \"%Ps\", \"Ly-cyc\": \"%Ps\", \"M-rk\": \"1\", \"CM\": \"%Ps\"},\n", cup_matrix);
    // }
    
    // fclose(fptr);

    printf(ANSI_COLOR_GREEN "Done! \n \n" ANSI_COLOR_RESET);

    // Close pari
    pari_close();

    //--------
    // Compute the time the whole program took
    clock_t duration = (clock()-start) / 1000;
    msec = duration%1000000;
    sec = (duration/1000)%60;
    min = duration/60000;

    printf (ANSI_COLOR_YELLOW "Runtime: %d min, %d,%d sec\n\n" ANSI_COLOR_RESET, min, sec, msec);
    
    //-----------
    return 0;
}