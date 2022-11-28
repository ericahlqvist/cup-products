/*
Generates all extensions of base together with all necessary automorphisms 
*/

GEN my_ext(GEN base, GEN base_clf, int disc, GEN s, GEN p, GEN D_prime_vect, int swap) 
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
    output(p1red);
    
    GEN LxRel = rnfinit(base, mkvec2(p1red, D_prime_vect));
    printf("---> LxRel <--- \n");
    

    GEN LxAbs = Buchall(rnf_get_polabs(LxRel), nf_FORCE, DEFAULTPREC);
    printf("---> LxAbs <--- \n\n");
    output(rnf_get_polabs(LxRel));
    output(bnf_get_cyc(LxAbs));
    pari_printf("Disc LxAbs: %Ps\n\n", nf_get_disc(bnf_get_nf(LxAbs)));

    output(p2red);
    // p2red = rnfpolredabs(base, p2red, 0);
    // output(p2red);
    GEN LyRel = rnfinit(base, mkvec2(p2red, D_prime_vect));
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
            pari_printf("sigma_y: %Ps, sigma_y^2: %Ps, sigma_y^3: %Ps\n\n", sigma_y, galoisapply(LyAbs, sigma_y, sigma_y), basistoalg(LyAbs, galoisapply(LyAbs, sigma_y, galoisapply(LyAbs, sigma_y, sigma_y))));
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

