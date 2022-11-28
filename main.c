
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