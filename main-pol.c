
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
    

    pari_printf("POL: %Ps\n", f);
    
    // Define K.pol
    printf("\n");
    // GEN f_red = polredabs0(f, 0);
    // pari_printf("pol_red: %Ps\n", f_red);
    // Define base field K
    K = Buchall(f, nf_FORCE, DEFAULTPREC);
    // MEDDEFAULTPREC, BIGDEFAULTPREC
    //K = Buchall_param(f, 1.5,1.5,4, nf_FORCE, DEFAULTPREC);
    D = nf_get_disc(bnf_get_nf(K));
    //D = gp_read_str(argv[3]);
    D_prime_vect = gel(factor(D), 1);
    //tu = bnf_get_tuU(K);
    Kcyc = bnf_get_cyc(K);
    int Knr = itos(bnf_get_no(K));

    if (Knr%p_int != 0)
    {
        pari_printf("%Ps does not divide the class number %d\n", p, Knr);
        pari_printf("Disc: %Ps\n", D);
        pari_close();
        exit(0);
    }
    

    pari_printf("K cyc: %Ps\n\n", Kcyc);
    pari_printf("Discriminant: %Ps\n\n", D);
    p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);
    pari_printf("p Cl Fld: %Ps\n\n", p_ClFld_pol);
    // pari_printf("Fund units: %Ps\n", bnf_get_fu(K));
    // pari_printf("Tors unit: %Ps\n\n", algtobasis(K, bnf_get_tuU(K)));
    // pari_printf("Tors unit 2: %Ps\n\n", nfpow(K, bnf_get_tuU(K), gen_2));
    // pari_printf("Tors unit 3: %Ps\n\n", nfpow(K, bnf_get_tuU(K), stoi(3)));

    // GEN clf_pol = polredabs0(mkvec2(bnrclassfield(K, p, 2, DEFAULTPREC), D_prime_vect),0);
    // pari_printf("H fld pol: %Ps\n\n", clf_pol);
    // // GEN clf_pol = bnrclassfield(K, p, 2, DEFAULTPREC);
    // GEN LAB = Buchall(clf_pol, nf_FORCE, DEFAULTPREC);
    // pari_printf("L cyc: %Ps\n\n", bnf_get_cyc(LAB));
    // my_unramified_p_extensions(K, p, D_prime_vect);
    // my_unramified_p_extensions_with_trivial_action(K, p, D_prime_vect);
    
    J_vect = my_find_p_gens(K, p);
    p_rk = glength(J_vect);
    printf("p-rank: %d\n\n", p_rk);

    if (p_rk<2)
    {
        printf("p-rank less than 2 --> finite tower\n\n");
        
        pari_close();
        exit(0);
    }

    GEN units_mod_p = my_find_units_mod_p(K, p);
    printf("Nr of units mod p: %ld\n", glength(units_mod_p));

    r_rk = glength(J_vect)+glength(units_mod_p);
    printf("r-rank: %d\n\n", r_rk);
 
    pari_printf("p_int: %d\n\nmy_pol: %Ps\n\nK_cyc: %Ps\n\n", p_int, f, Kcyc);

    GEN K_ext = my_ext(K, p_ClFld_pol, s, p, p_rk, D_prime_vect);
    printf("Extensions found\n\n");
    //my_norms(K, K_ext, p);
    //my_diffs(K_ext, p);
    // my_matrices(K_ext, p);
    // pari_close();
    // exit(0);
    // For Ja_vect, choose an extension gel(K_ext, 3) with big class group. Might need to be change once known.
    //int ext_nr = 1;
    //Ja_vect = my_find_Ja_vect_modified(gel(gel(K_ext, ext_nr), 1), gel(gel(K_ext, ext_nr), 2), K, gel(gel(K_ext, ext_nr), 3), J_vect, units_mod_p, p);
    Ja_vect = my_find_Ja_vect(K, J_vect, p, units_mod_p);
    pari_printf("Ja_vect: %Ps\n\n", Ja_vect);
    // GEN p_power_units = my_find_p_power_units(K, units_mod_p, p);
    // // GEN p_power_units_2 = my_get_unit_group(K, units_mod_p, p);
    // // GEN p_power_units = gconcat(p_power_units_1, p_power_units_2);
    // printf("Unit group, size: %ld\n\n", glength(p_power_units));
    
    
    
    //Defines a matrix over F_2 with index (i*k, j) corresponding to 
    //< x_i\cup x_k, (a_j, J_j)>

    // my_cup_matrix(K_ext, K, p, p_int, p_rk, J_vect, Ja_vect, units_mod_p, r_rk);
    // my_cup_matrix_transpose(K_ext, K, p, p_int, p_rk, J_vect, Ja_vect, units_mod_p, r_rk);

    //-----Faster version--------
    // my_cup_matrix_2(K_ext, K, p, p_int, p_rk, J_vect, units_mod_p, r_rk);
    //my_cup_matrix_2_transpose(K_ext, K, p, p_int, p_rk, J_vect, Ja_vect, units_mod_p, r_rk, p_power_units);
    my_cup_matrix_3(K_ext, K, p, p_int, p_rk, Ja_vect, r_rk);


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