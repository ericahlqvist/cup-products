
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

#include "misc_functions2.h"
#include "tests.h"
#include "artin_symbol.h"
//#include "test_artin.h"
#include "ext_and_aut2.h"
#include "find_cup_matrix3.h"



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
    GEN K, f, Kcyc, p_ClFld_pol, J_vect, D, D_prime_vect;

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
    
    
    J_vect = my_find_p_gens(K, p);
    p_rk = glength(J_vect);
    printf("p-rank: %d\n\n", p_rk);

    GEN units_mod_p = my_find_units_mod_p(K, p);
    printf("Nr of units mod p: %ld\n", glength(units_mod_p));

    r_rk = glength(J_vect)+glength(units_mod_p);
    printf("r-rank: %d\n\n", r_rk);
 
    pari_printf("p_int: %d\n\nmy_pol: %Ps\n\nK_cyc: %Ps\n\n", p_int, f, Kcyc);

    GEN K_ext = my_ext(K, p_ClFld_pol, s, p, p_rk, D_prime_vect);
    printf("Extensions found\n\n");
    
    //Defines a vector over F_2 with index i*j*k corresponding to 
    //< x_i \cup x_j \cup x_k, -1>

    int mat_rk = my_cup_matrix_12(K_ext, K, p, p_int, p_rk, r_rk, p_ClFld_pol);


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
    return mat_rk;
}