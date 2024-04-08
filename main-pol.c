
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
#include "test_artin.h"
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
    GEN K, f, Kcyc, p_ClFld_pol, J_vect, Ja_vect, D, D_prime_vect;

    // Read the prime number p from arguments
    p = gp_read_str(argv[1]);
    p_int = atoi(argv[1]);
    
    // Read the defining polynomial for K
    f = gp_read_str(argv[2]);
    pari_printf("POL: %Ps\n\n", f);
    
    //--------------------------------------------------
    // Define base field K
    K = Buchall(f, nf_FORCE, DEFAULTPREC);
    //--------------------------------------------------

    // Other possible precisions
    // MEDDEFAULTPREC, BIGDEFAULTPREC

    //K = Buchall_param(f, 1.5,1.5,4, nf_FORCE, DEFAULTPREC);

    //--------------------------------------------------
    // Discriminant
    D = nf_get_disc(bnf_get_nf(K));
    pari_printf("Discriminant: %Ps\n\n", D);

    // Factor discriminant
    D_prime_vect = gel(factor(D), 1);
    //--------------------------------------------------

    //tu = bnf_get_tuU(K);

    //--------------------------------------------------
    // Class group of K (cycle type)
    Kcyc = bnf_get_cyc(K);
    pari_printf("K cyc: %Ps\n\n", Kcyc);

    // Class numer of K
    int Knr = itos(bnf_get_no(K));

    // Test if p divides the class number. If not, then H^1(X, Z/pZ) = 0 and there is nothing to compute. 
    if (Knr%p_int != 0)
    {
        pari_printf("%Ps does not divide the class number %d\n", p, Knr);
        pari_printf("Disc: %Ps\n", D);
        pari_close();
        exit(0);
    }
    //--------------------------------------------------

    //--------------------------------------------------
    // Define polynomials for the generating fields for the part of the Hilbert class field corresp to Cl(K)/p. 
    p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);
    pari_printf("p Cl Fld: %Ps\n\n", p_ClFld_pol);
    //--------------------------------------------------

    // pari_printf("Fund units: %Ps\n", bnf_get_fu(K));
    // pari_printf("Tors unit: %Ps\n\n", algtobasis(K, bnf_get_tuU(K)));
    // pari_printf("Tors unit 2: %Ps\n\n", nfpow(K, bnf_get_tuU(K), gen_2));
    // pari_printf("Tors unit 3: %Ps\n\n", nfpow(K, bnf_get_tuU(K), stoi(3)));

    // GEN clf_pol = polredabs0(mkvec2(bnrclassfield(K, p, 2, DEFAULTPREC), D_prime_vect),0);
    // pari_printf("H fld pol: %Ps\n\n", clf_pol);
    // // // GEN clf_pol = bnrclassfield(K, p, 2, DEFAULTPREC);
    // GEN LAB = Buchall(clf_pol, nf_FORCE, DEFAULTPREC);
    // pari_printf("L cyc: %Ps\n\n", bnf_get_cyc(LAB));

    //-----------------------------
    // Uncomment this if you just want the polynomials of the unramified deg p extensions

    // my_unramified_p_extensions(K, p, D_prime_vect);
    
    // pari_close();
    // exit(0);

    //------------------------------

    // my_unramified_p_extensions_with_trivial_action(K, p, D_prime_vect);
    
    //--------------------------------------------------
    // Find generators for the p-torsion of the class group
    J_vect = my_find_p_gens(K, p);
    p_rk = glength(J_vect);
    printf("p-rank: %d\n\n", p_rk);
    //--------------------------------------------------

    // // 6,7
    // if (p_rk<2)
    // {
    //     printf("p-rank less than 2 --> finite tower\n\n");
        
    //     pari_close();
    //     exit(0);
    // }

    //--------------------------------------------------
    // find generators for the group of units modulo p
    GEN units_mod_p = my_find_units_mod_p(K, p);
    printf("Nr of units mod p: %ld\n", glength(units_mod_p));

    // Define r_rk -- the rank of H^2(X, Z/pZ)
    r_rk = glength(J_vect)+glength(units_mod_p);
    printf("r-rank: %d\n\n", r_rk);
    //--------------------------------------------------

    //--------------------------------------------------
    // Define the extensions generating the p-part of the Hilbert class field corresponding to CL(K)/p
    GEN K_ext = my_ext(K, p_ClFld_pol, s, p, p_rk, D_prime_vect);
    printf("Extensions found\n\n");
    //--------------------------------------------------

    //--------------------------------------------------
    // Find generators for H^1(X, mu_p), which is dual to H^2(X, Z/pZ)
    Ja_vect = my_find_Ja_vect(K, J_vect, p, units_mod_p);
    //pari_printf("Ja_vect: %Ps\n\n", Ja_vect);
    //--------------------------------------------------

    // GEN p_power_units = my_find_p_power_units(K, units_mod_p, p);
    // // GEN p_power_units_2 = my_get_unit_group(K, units_mod_p, p);
    // // GEN p_power_units = gconcat(p_power_units_1, p_power_units_2);
    // printf("Unit group, size: %ld\n\n", glength(p_power_units));
    
    
    //--------------------------------------------------
    // Defines a matrix over F_p with index (i*k, j) corresponding to 
    // < x_i\cup x_k, (a_j, J_j) > if i is not equal to j and
    // < B(x_i), (a_j, J_j) > if i=j. 
    // Here < - , - > denotes the Artin--Verdier pairing, which may be computed using our cup product formula and the Artin symbol. 
    int mat_rk = my_relations(K_ext, K, p, p_int, p_rk, Ja_vect, r_rk);
    //printf("\n");

    //--------------------------------------------------

    //printf(ANSI_COLOR_GREEN "Done! \n \n" ANSI_COLOR_RESET);

    // Close pari
    pari_close();

    // //--------
    // // Compute the time the whole program took to run
    // clock_t duration = (clock()-start) / 1000;
    // msec = duration%1000000;
    // sec = (duration/1000)%60;
    // min = duration/60000;

    // printf (ANSI_COLOR_YELLOW "Runtime: %d min, %d,%d sec\n\n" ANSI_COLOR_RESET, min, sec, msec);
    // //-----------
    return mat_rk;
}