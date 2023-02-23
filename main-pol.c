
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
#include "artin_symbol.h"
#include "test_artin.h"
#include "ext_and_aut2.h"
#include "find_cup_matrix.h"
#include "find_cup_matrix2.h"



int
main (int argc, char *argv[])	  
{
    //--------
    clock_t start = clock();
    //--------
    

    int p_int, p_rk, r_rk;

    int min;
    int sec;
    int msec;
    
    pari_init(10000000000,500000);
    // printf("Initial adress: %ld\n", avma);
    // pari_sp limit = stack_lim(avma, 1);
    
    GEN p = cgeti(DEFAULTPREC);
    GEN s = pol_x(fetch_user_var("s"));
    // GEN x = pol_x(fetch_user_var("x"));
    GEN K, f, Kcyc, p_ClFld_pol, J_vect, Ja_vect, D, D_prime_vect;

    p = gp_read_str(argv[1]);
    p_int = atoi(argv[1]);
    
    f = gp_read_str(argv[2]);
    D = gp_read_str(argv[3]);
    D_prime_vect = gel(factor(D), 1);

    pari_printf("POL: %Ps\n", f);

    // Define K.pol
    printf("\n");
    
    // Define base field K
    K = Buchall(f, nf_FORCE, DEFAULTPREC);
    //tu = bnf_get_tuU(K);
    Kcyc = bnf_get_cyc(K);
    pari_printf("K cyc: %Ps\n\n", Kcyc);
    pari_printf("Discriminant: %Ps\n\n", D);
    p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);
    // pari_printf("p Cl Fld: %Ps\n\n", p_ClFld_pol);
    // pari_printf("Fund units: %Ps\n", bnf_get_fu(K));
    // pari_printf("Tors unit: %Ps\n\n", algtobasis(K, bnf_get_tuU(K)));
    // pari_printf("Tors unit 2: %Ps\n\n", nfpow(K, bnf_get_tuU(K), gen_2));
    // pari_printf("Tors unit 3: %Ps\n\n", nfpow(K, bnf_get_tuU(K), stoi(3)));

    // GEN clf_pol = gsubstpol(polredabs0(mkvec2(bnrclassfield(K, p, 2, DEFAULTPREC), D_prime_vect),0), x, s);
    // // GEN clf_pol = bnrclassfield(K, p, 2, DEFAULTPREC);
    // pari_printf("abs pol: %Ps\n\n", clf_pol);
    // my_unramified_p_extensions(K, p, D_prime_vect);
    // my_unramified_p_extensions_with_trivial_action(K, p, D_prime_vect);
    
    J_vect = my_find_p_gens(K, p);
    // J_vect = bnf_get_gen(K);
    //pari_printf("J_vect: %Ps\n\n", J_vect);
    
    GEN units_mod_p = my_find_units_mod_p(K, p);
    printf("Nr of units mod p: %ld\n", glength(units_mod_p));

    r_rk = glength(J_vect)+glength(units_mod_p);
    printf("r-rank: %d\n\n", r_rk);

    // Ja_vect = my_find_Ja_vect(K, J_vect, p);
    
    //pari_printf("Ja_vect: %Ps\n\n", Ja_vect);
    p_rk = glength(J_vect);
    printf("p-rank: %d\n\n", p_rk);



    pari_printf("p_int: %d\n\nmy_pol: %Ps\n\nK_cyc: %Ps\n\n", p_int, f, Kcyc);

    GEN K_ext = my_ext(K, p_ClFld_pol, s, p, p_rk, D_prime_vect);
    printf("Extensions found\n\n");
    
    Ja_vect = my_find_Ja_vect_modified(gel(gel(K_ext, 1), 1), gel(gel(K_ext, 1), 2), K, gel(gel(K_ext, 1), 3), J_vect, units_mod_p, p);

    GEN p_power_units = my_find_p_power_units(K, units_mod_p, p);
    
    // if (my_is_p_torsion(K,J_vect, p))
    // {
    //     printf(ANSI_COLOR_GREEN "p-torsion test passed\n\n" ANSI_COLOR_RESET);
    // }
    // else {
    //     printf(ANSI_COLOR_RED "p-torsion test  failed\n\n" ANSI_COLOR_RESET);
    //     pari_close();
    //     exit(0);
    // }
    
    
    
    //Defines a matrix over F_2 with index (i*k, j) corresponding to 
    //< x_i\cup x_k, (a_j, J_j)>

    // my_cup_matrix(K_ext, K, p, p_int, p_rk, J_vect, Ja_vect, units_mod_p, r_rk);
    // my_cup_matrix_transpose(K_ext, K, p, p_int, p_rk, J_vect, Ja_vect, units_mod_p, r_rk);

    //-----Faster version--------
    // my_cup_matrix_2(K_ext, K, p, p_int, p_rk, J_vect, units_mod_p, r_rk);
    my_cup_matrix_2_transpose(K_ext, K, p, p_int, p_rk, J_vect, Ja_vect, units_mod_p, r_rk, p_power_units);
    
    printf("\n\n");
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