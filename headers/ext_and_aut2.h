/*
Generates all extensions of base together with all necessary automorphisms 
*/

GEN my_ext_old(GEN base, GEN base_clf, int disc, GEN s, GEN p, GEN D_prime_vect, int swap) 
{   
    
    printf("Finding extensions... \n\n");
    GEN x, y;

    x = pol_x(fetch_user_var("x"));
    y = pol_x(fetch_user_var("y"));
    

    GEN base_ext = cgetg(6, t_VEC);
    pari_printf("Base_clf: %Ps\n\n", base_clf);

    // pari_close();
    // exit(0);
    GEN p1, p2;
    if (!swap)
    {
        p1 = gel(base_clf, 1);
        p2 = gel(base_clf, 2);
    }
    else
    {
        p1 = gel(base_clf, 2);
        p2 = gel(base_clf, 1);
    }
    
    
    GEN q1 = gsubstpol(p1, x, y);
    GEN q2 = gsubstpol(p2, x, y);
    
    /* Define LxRel */
    GEN p1red = rnfpolredbest(base, mkvec2(q1, D_prime_vect), 0);
    GEN p2red = rnfpolredbest(base, mkvec2(q2, D_prime_vect), 0);
    
    // GEN p1red = q1;
    // GEN p2red = q2;
    printf("p1red: ");
    output(p1red);
    printf("\n\n");
    printf("D_prime_vect: ");
    output(D_prime_vect);
    
    GEN LxRel = rnfinit(base, p1red);
    printf("---> LxRel <--- \n");
    

    GEN LxAbs = Buchall(rnf_get_polabs(LxRel), nf_FORCE, DEFAULTPREC);
    printf("---> LxAbs <--- \n\n");
    output(rnf_get_polabs(LxRel));
    output(bnf_get_cyc(LxAbs));
    pari_printf("Disc LxAbs: %Ps\n\n", nf_get_disc(bnf_get_nf(LxAbs)));

    output(p2red);
    // p2red = rnfpolredabs(base, p2red, 0);
    // output(p2red);
    GEN LyRel = rnfinit(base, p2red);
    printf("\n---> LyRel <--- \n");
    
    
    GEN LyAbs = Buchall(rnf_get_polabs(LyRel), nf_FORCE, DEFAULTPREC);
    printf("---> LyAbs <--- \n\n");
    output(rnf_get_polabs(LyRel));
    printf("\n");
    pari_printf("Lx_cyc: %Ps\n\n", bnf_get_cyc(LxAbs));
    pari_printf("Ly_cyc: %Ps\n\n", bnf_get_cyc(LyAbs));


    GEN s_lift_x = rnfeltup0(LxRel, s, 1);
    GEN s_lift_y = rnfeltup0(LyRel, s, 1);

    GEN cx = galoisconj(LxAbs, NULL);
    pari_printf("GALOIS: %Ps\n\n", cx);
    GEN cy = galoisconj(LyAbs, NULL);

    
    int i;
    GEN sigma_x;
    for (i = 1; i < glength(cx)+1; ++i)
    {
        if ( (!my_QV_equal(algtobasis(LxAbs,gel(cx, i)), algtobasis(LxAbs,y))) && my_QV_equal(galoisapply(LxAbs, gel(cx,i), s_lift_x), s_lift_x)) 
        {
            sigma_x = algtobasis(LxAbs, gel(cx, i));
            // sigma_x = galoisapply(LxAbs, sigma_x, sigma_x);
            // printf(ANSI_COLOR_GREEN "sigma_x found %d\n\n" ANSI_COLOR_RESET, i);
            pari_printf("sigma_x: %Ps, sigma_x^2: %Ps, sigma_x^3: %Ps\n\n", sigma_x, galoisapply(LxAbs, sigma_x, sigma_x), basistoalg(LxAbs, galoisapply(LxAbs, sigma_x, galoisapply(LxAbs, sigma_x, sigma_x))));
            break;
        }
    }

    printf(ANSI_COLOR_CYAN "---> sigma_x <--- \n \n" ANSI_COLOR_RESET);

    GEN sigma_y = pol_x(fetch_user_var("sigma_y"));
    for (i = 1; i < glength(cy)+1; ++i)
    {
        if ( (!my_QV_equal(algtobasis(LyAbs,gel(cy, i)), algtobasis(LyAbs,y))) && my_QV_equal(galoisapply(LyAbs, gel(cy,i), s_lift_y), s_lift_y)) 
        {
            sigma_y = algtobasis(LyAbs, gel(cy, i));
            // sigma_y = galoisapply(LyAbs, sigma_y, sigma_y);
            // printf(ANSI_COLOR_GREEN "sigma_y found %d\n\n" ANSI_COLOR_RESET, i);
            // pari_printf("sigma_y: %Ps, sigma_y^2: %Ps, sigma_y^3: %Ps\n\n", sigma_y, galoisapply(LyAbs, sigma_y, sigma_y), basistoalg(LyAbs, galoisapply(LyAbs, sigma_y, galoisapply(LyAbs, sigma_y, sigma_y))));
            break;
        }
    }

    printf(ANSI_COLOR_CYAN "---> sigma_y <--- \n \n" ANSI_COLOR_RESET);

    gel(base_ext, 1) = LxAbs;
    gel(base_ext, 2) = LxRel;
    gel(base_ext, 3) = LyAbs;
    gel(base_ext, 4) = LyRel;
    gel(base_ext, 5) = sigma_x;
    gel(base_ext, 6) = sigma_y;





    return base_ext;
}

GEN my_ext(GEN base, GEN base_clf, GEN s, GEN p, int p_rk, GEN D_prime_vect) 
{   
    pari_sp av = avma;
    printf("\n--------------------------\nStart: my_ext\n--------------------------\n\n");
    GEN x, y, p1, q1, p1red, Lrel, Labs, s_lift_x, cx, sigma;

    x = pol_x(fetch_user_var("x"));
    y = pol_x(fetch_user_var("y"));

    GEN base_ext = cgetg(p_rk+1,t_VEC);

    // printf("base l: %ld\n", glength(base_ext));
    // pari_printf("Base_clf: %Ps\n\n", base_clf);
    char filename[100];
    
    int i, j;
    for (i=1; i<p_rk+1; ++i) {
        p1 = gel(base_clf, i);
        q1 = gsubstpol(p1, x, y);
        
        /* Define Lrel/Labs */
        p1red = rnfpolredbest(base, mkvec2(q1, D_prime_vect), 0);
        // p1red = q1;
        printf("Reduced polynomial for relative extension found\n");
        Lrel = rnfinit(base, p1red);
        printf("Lrel found\n");
        pari_printf("Abs pol: %Ps\n", rnf_get_polabs(Lrel));
        
        Labs = Buchall(rnf_get_polabs(Lrel), nf_FORCE, DEFAULTPREC);
        
        // May change flag to nf_FORCE
        // Other precisions: MEDDEFAULTPREC, BIGDEFAULTPREC
        //Labs = Buchall_param(rnf_get_polabs(Lrel), 1.5,1.5,4, nf_FORCE, DEFAULTPREC);
        printf("Labs found\n");
        pari_printf("L_cyc[%d]: %Ps\n", i, bnf_get_cyc(Labs));
        pari_printf("rel_pol[%d]: %Ps\n", i, p1red);
        pari_printf("\nabs pol: %Ps\n\n", gsubstpol(rnf_get_polabs(Lrel),y, s));

        s_lift_x = rnfeltup0(Lrel, s, 1);
        //pari_printf("\nLift: %Ps\n\n", s_lift_x);
        cx = galoisconj(Labs, NULL);
        //pari_printf("\nGalois conj: %Ps\n\n", cx);
        if (glength(cx)==nf_get_degree(bnf_get_nf(Labs)))
        {
            printf(ANSI_COLOR_GREEN "\n------------------------\nLabs is Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
        }
        else {
            printf(ANSI_COLOR_RED "\n------------------------\nLabs is not Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
        }
        

        for (j = 1; j < glength(cx)+1; ++j)
        {
            if ( (!my_QV_equal(algtobasis(Labs,gel(cx, j)), algtobasis(Labs,y))) && my_QV_equal(galoisapply(Labs, gel(cx,j), s_lift_x), s_lift_x)) 
            {
                sigma = algtobasis(Labs, gel(cx, j));
                //pari_printf("sigma: %Ps\n\n", sigma);
                break;
            }
        }
        //printf(ANSI_COLOR_CYAN "---> sigma <--- \n \n" ANSI_COLOR_RESET);
        
    
        gel(base_ext, i) = mkvec3(Labs, Lrel, sigma);

        sprintf(filename, "/Users/eric/Documents/Matematik/cup-products/large-fields/test/ext_%d", i);
        writebin(filename, gel(base_ext, i));
        
        
        // my_test_artin_symbol (Labs, Lrel, base, itos(p), sigma);
    }
    printf(ANSI_COLOR_GREEN "Extensions and generators found\n----------------------------\n\n" ANSI_COLOR_RESET);
    base_ext = gerepilecopy(av, base_ext);
    printf("\n--------------------------\nEnd: my_ext\n--------------------------\n\n");
    return base_ext;
}

GEN my_ext_from_file(GEN base, char *fields[], GEN base_clf, GEN s, GEN p, int p_rk, GEN D_prime_vect) 
{   
    pari_sp av = avma;
    printf("\n--------------------------\nStart: my_ext_from_file\n--------------------------\n\n");
    GEN x, y, p1, q1, p1red, Lrel, Labs, s_lift_x, cx, sigma;

    x = pol_x(fetch_user_var("x"));
    y = pol_x(fetch_user_var("y"));

    GEN base_ext = zerovec(p_rk);
    // printf("base l: %ld\n", glength(base_ext));
    // pari_printf("Base_clf: %Ps\n\n", base_clf);

    int i, j, check_sigma=0;
    for (i=1; i<p_rk+1; ++i) {

        Labs = gp_read_file(fields[i-1]);
        // gel(gel(Labs, 7), 1) = gsubstpol(nf_get_pol(bnf_get_nf(Labs)), gpolvar(nf_get_pol(bnf_get_nf(Labs))), y);
        printf("\n\n------------------------------------\n\nLabs found\n");
        pari_printf("L_cyc[%d]: %Ps\n", i, bnf_get_cyc(Labs));

        // Uncomment to test if the GRH can be removed
        // setalldebug(1);
        // if (bnfcertify0(Labs, 0))
        // {
        //     printf(ANSI_COLOR_GREEN "\n------------------------\nGRH removed\n------------------------\n\n" ANSI_COLOR_RESET);
        // }
        
        // pari_printf("rel_pol[%d]: %Ps\n", i, p1red);
        pari_printf("\nabs pol: %Ps\n\n", nf_get_pol(bnf_get_nf(Labs)));
        pari_printf("Variable: %Ps\n\n", gpolvar(nf_get_pol(bnf_get_nf(Labs))));

        p1 = gel(base_clf, i);
        q1 = gsubstpol(p1, x, y);
        
        /* Define Lrel/Labs */
        
        p1red = rnfpolredbest(base, mkvec2(q1, D_prime_vect), 0);
        
        //p1red = q1;
        
        
        printf("\n\n------------------------------------\n\nReduced polynomial for relative extension found\n");
        Lrel = rnfinit(base, p1red);
        printf("Lrel found\n");
        pari_printf("Abs pol: %Ps\n------------------------------------\n", rnf_get_polabs(Lrel));
        
        s_lift_x = rnfeltup0(Lrel, s, 1);
        pari_printf("\nLift: %Ps\n\n%Ps\n\n", s_lift_x, basistoalg(Labs, s_lift_x));
        cx = galoisconj0(Labs, 0, NULL, DEFAULTPREC);
        if (glength(cx)==nf_get_degree(bnf_get_nf(Labs)))
        {
            printf(ANSI_COLOR_GREEN "\n------------------------\nLabs is Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
        }
        else {
            printf(ANSI_COLOR_RED "\n------------------------\nLabs is not Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
        }
        //printf(ANSI_COLOR_YELLOW "\nGalois group size: %ld\n\n" ANSI_COLOR_RESET, glength(cx));
        
        for (j = 1; j < lg(cx); ++j)
        {
            pari_printf("sigma s[%d]: %Ps\n", j, galoisapply(Labs, gel(cx,j), s_lift_x));
            pari_printf("s[%d]: %Ps\n\n", j, s_lift_x);

            if ( (!my_QV_equal(algtobasis(Labs,gel(cx, j)), algtobasis(Labs,y))) && my_QV_equal(galoisapply(Labs, gel(cx,j), s_lift_x), s_lift_x)) 
            {
                sigma = algtobasis(Labs, gel(cx, j));
                pari_printf("sigma: %Ps\n\n", sigma);
                check_sigma=1;
                break;
            }
        }
        if (!check_sigma)
        {
            printf(ANSI_COLOR_RED "\nERROR: No sigma found\n\n" ANSI_COLOR_RESET);
            pari_close();
            exit(111);
        }
        
        //printf(ANSI_COLOR_CYAN "---> sigma <--- \n \n" ANSI_COLOR_RESET);
        
        

        gel(base_ext, i) = mkvec3(Labs, Lrel, sigma);
        // my_test_artin_symbol (Labs, Lrel, base, itos(p), sigma);
    }
    base_ext = gerepilecopy(av, base_ext);
    printf("\n--------------------------\nEnd: my_ext_from_file\n--------------------------\n\n");
    return base_ext;
}

