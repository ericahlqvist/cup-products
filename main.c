
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
#include "find_orbit.h"
#include "ext_and_aut.h"
#include "find_variables.h"
#include "find_basis.h"
#include "artin_symbol.h"
#include "find_cup_matrix.h"

int
main (int argc, char *argv[])	  
{
    //--------
    clock_t start = clock();
    //--------
    // char prime_str = argv[1];
    // char my_det = argv[2];
    //char swap_str[100];
    

    int p_int, my_int, i;

    int min;
    int sec;
    int msec;
    
    pari_init(4000000000,500000);
    // printf("Initial adress: %ld\n", avma);
    // pari_sp limit = stack_lim(avma, 1);
    
    GEN p = cgeti(DEFAULTPREC);
    GEN s = pol_x(fetch_user_var("s"));
    GEN K, f, Kcyc, p_ClFld_pol, D_prime_vect, D;

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
    Kcyc = bnf_get_cyc(K);
    p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);

    pari_printf("p_int: %d\n\nmy_int: %d\n\nK_cyc: %Ps\n\nK_basis: %Ps\n\n", p_int, my_int, Kcyc, nf_get_zk(bnf_get_nf(K)));

    GEN K_ext_aut   =   my_ext(K, p_ClFld_pol, my_int, s, p, D_prime_vect, 0);
    GEN LxAbs       =   gel(K_ext_aut, 1);
    GEN LxRel       =   gel(K_ext_aut, 2);
    GEN LyAbs       =   gel(K_ext_aut, 3);
    GEN LyRel       =   gel(K_ext_aut, 4);
    GEN sigma_x     =   gel(K_ext_aut, 5);
    GEN sigma_y     =   gel(K_ext_aut, 6);

    GEN Lx_cyc = bnf_get_cyc(LxAbs);
    GEN Ly_cyc = bnf_get_cyc(LyAbs);

    GEN J_vect = my_find_p_gens(K, p);
    
    GEN T_x = rnfisnorminit(K, rnf_get_pol(LxRel), 1);
    
    GEN x_basis = my_find_basis_2(LxAbs, LxRel, K, sigma_x, p, J_vect, T_x);
    GEN I_vect = my_find_I_from_basis(x_basis);

    GEN cup_matrix = my_cup_matrix(LxAbs, LxRel, LyAbs, LyRel, K, sigma_x, sigma_y, p, J_vect, I_vect, T_x, p_int);

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