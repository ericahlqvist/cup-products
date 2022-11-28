
GEN my_find_I_from_basis(GEN my_basis) {

    GEN I_vect = zerovec(glength(my_basis));

    int i;
    for (i = 1; i < glength(I_vect)+1; ++i) {
        gel(I_vect, i) = gel(gel(my_basis, i), 4);
    }

    return I_vect;
}

GEN my_find_a_from_basis(GEN my_basis) {

    GEN a_vect = zerovec(glength(my_basis));

    int i;
    for (i = 1; i < glength(a_vect)+1; ++i) {
        gel(a_vect, i) = gel(gel(my_basis, i), 2);
    }

    return a_vect;
}


GEN my_find_J (GEN my_basis) {

    GEN J_vect = zerovec(glength(my_basis));

    int i;
    for (i = 1; i < glength(J_vect)+1; ++i) {
        gel(J_vect, i) = gel(gel(my_basis, i), 3);
    }

    return J_vect;
}


GEN my_find_I_prime (GEN LxAbs, GEN LxRel, GEN K, GEN I_vect, GEN sigma_x, GEN p) {
    GEN I_prime_vect = zerovec(glength(I_vect)), modifiers;
    int i;
    for (i = 1; i < glength(I_vect)+1; i++)
    {
        modifiers = my_find_orbit(LxAbs, LxRel, K, gel(I_vect, i), sigma_x, p);
        printf("Orbit: ");
        output(gel(modifiers, 1));
        printf("\n\n");
        gel(I_prime_vect, i) = gel(modifiers, 2);
    }
    
    return I_prime_vect;
}




GEN my_find_I2 (GEN LyAbs, GEN LyRel, GEN K, GEN sigma_y, GEN a2_vect, GEN J_vect, int p) {
    // printf(ANSI_COLOR_CYAN "my_find_I2\n\n" ANSI_COLOR_RESET);
    GEN I2_vect = zerovec(glength(J_vect));
    GEN div_a2 = pol_x(fetch_user_var("div_a2"));
    // GEN rel_ideal = pol_x(fetch_user_var("relative"));
    GEN iJ = pol_x(fetch_user_var("iJ"));
    GEN iJ_div_a2 = pol_x(fetch_user_var("iJ_div_a2"));
    
    int i;
    for (i = 1; i < glength(a2_vect)+1; i++)
    {
        div_a2 = idealhnf0(LyAbs, gel(a2_vect, i), NULL);
        printf("Computing iI\n\n");
        iJ = rnfidealup0(LyRel, gel(J_vect, i),1);
        printf("Computing div(b1) + iI\n\n");
        iJ_div_a2 = idealmul(LyAbs, iJ, div_a2);
        
        // printf("Computing relative div(b1) + iI\n\n");
        // rel_ideal = rnfidealabstorel(LyRel, iJ_div_a2);
        // printf("Computing N (div(b1) + iI)\n\n");
        // pari_printf(ANSI_COLOR_CYAN "\nN (div(b1) + iI): %Ps\n\n" ANSI_COLOR_RESET, rnfidealnormrel(LyRel, rel_ideal));
        // printf(ANSI_COLOR_YELLOW "\n----------\n" ANSI_COLOR_RESET);
        
        if (my_SQ_MAT_equal(iJ_div_a2, idealhnf0(LyAbs, gen_1, NULL)))
        {
            gel(I2_vect, i) = idealhnf0(LyAbs, gen_1, NULL);
        }
        else {
            gel(I2_vect, i) = my_find_H90_ideal(LyAbs, LyRel, K, iJ_div_a2, sigma_y, p);
        }
        
        if (my_SQ_MAT_equal(iJ_div_a2, my_SM1_ideal(LyAbs, sigma_y, gel(I2_vect, i))))
        {
            printf(ANSI_COLOR_GREEN "I2[%d] is an H90 companion\n\n" ANSI_COLOR_RESET, i);
        }
        else {
            outmat(idealfactor(LyAbs, iJ_div_a2));
            outmat(iJ_div_a2);
            outmat(my_SM1_ideal(LyAbs, sigma_y, gel(I2_vect, i)));
            printf(ANSI_COLOR_RED "ERROR - I2[%d] incorrect\n\n" ANSI_COLOR_RESET, i);
        }
        
    }
    
    return I2_vect;
}


GEN my_find_I_prime_2 (GEN LxAbs, GEN LxRel, GEN K, GEN I_vect, GEN sigma_x, GEN p, GEN T) {
    pari_sp av = avma;
    GEN I_prime_vect = zerovec(glength(I_vect)), norms, I_prime, u, div_u, i_J_prime;
    int i;
    for (i = 1; i < glength(I_vect)+1; i++)
    {
        norms = my_is_principal_mod_p(K, rnfidealnormrel(LxRel, rnfidealabstorel(LxRel, gel(I_vect, i))), p);
        u = rnfeltreltoabs(LxRel, gel(rnfisnorm(T, gel(norms,1), 0),1));
        div_u = idealhnf0(LxAbs, u, NULL);
        i_J_prime = rnfidealup0(LxRel, gel(norms,2), 1); 
        I_prime = my_find_H90_ideal(LxAbs, LxRel, K, idealdiv(LxAbs, idealmul(LxAbs, div_u, i_J_prime), gel(I_vect, i)), sigma_x, itos(p));
        gel(I_prime_vect, i) = I_prime;
    }
    I_prime_vect = gerepilecopy(av, I_prime_vect);
    return I_prime_vect;
}