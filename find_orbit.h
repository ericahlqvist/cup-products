GEN my_find_norm_group (GEN LxAbs, GEN LxRel, GEN K, int p) {

    GEN L_clgp = bnf_get_gen(LxAbs);// my_get_clgp(LxAbs);
    int clnr = glength(L_clgp);
    GEN norm_basis_non_red = zerovec(clnr);
    GEN norm_group = mkvec(idealhnf0(K, gen_1, NULL));
    GEN lifts = mkvec(idealhnf0(LxAbs, gen_1, NULL));
    GEN test_ideal;

    int i, j, check;
    for (i = 1; i < clnr+1; i++)
    {
        gel(norm_basis_non_red, i) = idealred(K, rnfidealnormrel(LxRel, gel(L_clgp, i)));
        check = 1;
        
        for (j = 1; j < glength(norm_group)+1; j++)
        {
            test_ideal = idealdiv(K, gel(norm_basis_non_red, i), gel(norm_group, j));
            if (!my_QV_equal0(my_is_principal_mod_p (K, test_ideal, stoi(p))))
            {
                check = 0;
            }
        }
            
        if (check)
        {
            norm_group = gconcat(norm_group, mkvec(gel(norm_basis_non_red, i)));
            lifts = gconcat(lifts, mkvec(gel(L_clgp, i)));
            output(norm_group);
            
        }
        
    }

    return mkvec2(norm_group, lifts);
}

GEN my_find_p_lifts (GEN LxAbs, GEN LxRel, GEN K, GEN p) {
    GEN gens = my_find_p_gens(K, p);
    GEN i_x_Cl_K_p = zerovec(9);

    gel(i_x_Cl_K_p, 1) = idealhnf0(LxAbs, gen_1, NULL);
    gel(i_x_Cl_K_p, 2) = idealred(LxAbs, rnfidealup0(LxRel, gel(gens, 1), 1));
    gel(i_x_Cl_K_p, 3) = idealred(LxAbs, idealpow(LxAbs, gel(i_x_Cl_K_p, 2), gen_2));
    gel(i_x_Cl_K_p, 4) = idealred(LxAbs, rnfidealup0(LxRel, gel(gens, 2), 1));
    gel(i_x_Cl_K_p, 5) = idealred(LxAbs, idealpow(LxAbs, gel(i_x_Cl_K_p, 4), gen_2));

    gel(i_x_Cl_K_p, 6) = idealred(LxAbs, idealmul(LxAbs, gel(i_x_Cl_K_p, 2), gel(i_x_Cl_K_p, 4)));
    gel(i_x_Cl_K_p, 7) = idealred(LxAbs, idealmul(LxAbs, gel(i_x_Cl_K_p, 2), gel(i_x_Cl_K_p, 5)));
    gel(i_x_Cl_K_p, 8) = idealred(LxAbs, idealmul(LxAbs, gel(i_x_Cl_K_p, 3), gel(i_x_Cl_K_p, 4)));
    gel(i_x_Cl_K_p, 9) = idealred(LxAbs, idealmul(LxAbs, gel(i_x_Cl_K_p, 3), gel(i_x_Cl_K_p, 5)));


    return i_x_Cl_K_p;
}

int my_same_orbit (GEN LxAbs, GEN I1, GEN I2, GEN H) {
    GEN diff = idealred(LxAbs, idealdiv(LxAbs, I1, I2));
    GEN test_ideal;
    int i, same = 0;
    for (i = 1; i < glength(H)+1; i++)
    {
        test_ideal = idealred(LxAbs, idealmul(LxAbs, diff, gel(H, i)));
        if (my_SQ_MAT_equal(idealhnf0(LxAbs, gen_1, NULL), test_ideal))
        {
            same = 1; 
        }
        
    }
    return same;
}

GEN my_lift_Cl_K (GEN LxAbs, GEN LxRel, GEN K) {
    GEN i_x_Cl_K = my_get_clgp(K);
    GEN i_x_Cl_K_red = mkvec(idealhnf0(LxAbs, gen_1, NULL));
    GEN test_ideal;
    int i, j, check;
    for (i = 1; i < glength(i_x_Cl_K)+1; i++)
    {
        gel(i_x_Cl_K, i) = idealred(LxAbs, rnfidealup0(LxRel, gel(i_x_Cl_K, i), 1));

        check = 1;
        
        for (j = 1; j < glength(i_x_Cl_K_red)+1; j++)
        {
            test_ideal = idealred(LxAbs, idealdiv(LxAbs, gel(i_x_Cl_K, i), gel(i_x_Cl_K_red, j)));
            if (my_QV_equal0(gel(bnfisprincipal0(LxAbs, test_ideal, 1), 1)))
            {
                check = 0;
            }
        }
            
        if (check)
        {
            i_x_Cl_K_red = gconcat(i_x_Cl_K_red, mkvec(gel(i_x_Cl_K, i)));
        }

    }
    printf("i_x_Cl_K_red: ");
    printf("%ld", glength(i_x_Cl_K_red));
    printf("\n");
    
    return i_x_Cl_K_red;
}

GEN my_find_orbit (GEN LxAbs, GEN LxRel, GEN K, GEN J2, GEN sigma_x, GEN p) {
    GEN test_ideal, J2_prime;
    GEN Nx_Cl_Lx_mod_p = my_find_norm_group(LxAbs, LxRel, K, itos(p));
    gel(Nx_Cl_Lx_mod_p, 1) = gconcat(gel(Nx_Cl_Lx_mod_p, 1), mkvec(idealred(K, idealpow(K, gmael2(Nx_Cl_Lx_mod_p, 1, 2), gen_2))));
    gel(Nx_Cl_Lx_mod_p, 2) = gconcat(gel(Nx_Cl_Lx_mod_p, 2), mkvec(idealpow(LxAbs, gmael2(Nx_Cl_Lx_mod_p, 2, 2), gen_2)));
    int i, j, k;

    GEN orbits = gel(Nx_Cl_Lx_mod_p, 2);
    GEN Nx_J2 = rnfidealnormrel(LxRel, idealred(LxAbs, J2));
    for (i = 1; i < glength(gel(Nx_Cl_Lx_mod_p, 1))+1; i++)
    {
        test_ideal = idealdiv(K, Nx_J2, gmael2(Nx_Cl_Lx_mod_p, 1, i));
        if (my_is_principal_mod_p (K, test_ideal, p))
            {
                J2_prime = gel(orbits, i);
                break;
            }
    }
    

    GEN i_x_Cl_K = my_lift_Cl_K(LxAbs, LxRel, K);
    // GEN i_x_Cl_K_p = my_find_p_lifts(LxAbs, LxRel, K, p);
    GEN L_clgp = my_get_clgp(LxAbs);
    
    for (j = 1; j < glength(L_clgp)+1; j++) {
        for (k = 1; k < glength(i_x_Cl_K)+1; k++) {
            pari_sp ltop = avma;
            printf("%d, %d\n", j, k);
            test_ideal = idealred(LxAbs, idealmul(LxAbs, idealdiv(LxAbs, J2, J2_prime), idealmul(LxAbs, my_1MS_ideal(LxAbs, sigma_x, gel(L_clgp, j)), gel(i_x_Cl_K, k))));
            
            if (my_QV_equal0(bnfisprincipal0(LxAbs, test_ideal, 0)))
            {
                printf("Orbit found! \n\n");
                printf("Current adress: %ld\n", avma);
                return mkvec2(gel(orbits, i), gel(L_clgp, j));
                
            }
            pari_sp lbot = avma;
            // if (my_SQ_MAT_equal(test_ideal, idealhnf(LxAbs, gen_1)))
            // {
            //     printf("Orbit found! \n\n");
            //     printf("Current adress: %ld\n", avma);
            //     return mkvec2(gel(orbits, i), gel(L_clgp, j));
                
            // }
            gerepile(ltop, lbot, NULL);
        }
    }
    
    printf(ANSI_COLOR_RED "No orbit found\n\n" ANSI_COLOR_RESET);
    pari_close();
    exit(0);
}