GEN my_find_uB (GEN Lrel, GEN K, GEN J, GEN I, int p_int) 
{
    
    GEN class_group = my_get_clgp(K);
    int clnr = itos(bnf_get_no(K));
    GEN p = stoi(p_int);
    GEN exp = stoi(p_int*(p_int-1)/2);
    //pari_printf("%Ps, %Ps\n", J, idealhnf(K, gen_1));
    GEN ideal_2 = idealpow(K, J, exp);
    
    GEN norm = rnfidealnormrel(Lrel, I);
    
    GEN prod = idealmul(K, norm, ideal_2);
    GEN test, pB, test_ideal, B;
    
    int condition = 1;

    int n;
    for (n = 1; n < clnr-1; ++n) {
        B = gel(class_group, n);
        pB = idealpow(K, B, p);
        test_ideal = idealdiv(K, prod, pB);
        test = bnfisprincipal0(K, test_ideal, 3);
        
        if (my_QV_equal0(gel(test, 1))) {
            printf(ANSI_COLOR_GREEN "Found u, B\n\n" ANSI_COLOR_RESET);
            pari_printf("u: %Ps\n\nB: %Ps\n\n", gel(test, 2), B);
            condition = 0;
            break;
        }
    }
    if (condition) {
        printf(ANSI_COLOR_RED "ERROR: No u, B found" ANSI_COLOR_RESET);
        pari_close();
        exit(0);
    }
    GEN output = mkvec2(gel(test, 2), B);
    return output;
}

//---------------------
// A vector of elements (a,J,I) with (N_x(a), J) forming a basis for Z^1/B^1.
//--------------------

GEN my_find_basis (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN J_vect, GEN p) 
{
    pari_sp av = avma;
    int i;
    GEN ext_gens = mkvec2(rnfidealup0(Lrel, gel(J_vect, 1), 1), rnfidealup0(Lrel, gel(J_vect, 2), 1));
    GEN my_basis = zerovec(2), aJI_vect = zerovec(2);
    int nr_comp = glength(J_vect);
    int equal;
    for (i = 1; i < nr_comp+1; ++i)
    {
        GEN aI = my_find_I(Labs, K, sigma, gel(ext_gens,i));
        
        GEN aJI = mkvec3(gel(aI, 1),gel(J_vect, i),gel(aI, 2));
        
        // Test relation

        equal = my_SQ_MAT_equal(idealmul(Labs, idealhnf0(Labs, gel(aI, 1), NULL), idealmul(Labs, gel(ext_gens,i), my_1MS_ideal(Labs, sigma, gel(aI, 2)))), idealhnf0(Labs, gen_1, NULL));

        if (equal)
        {
            printf (ANSI_COLOR_GREEN "Test passed\n" ANSI_COLOR_RESET);
        }
        else {
            printf(ANSI_COLOR_RED "ERROR - test failed for (a,J,I)\n\n" ANSI_COLOR_RED);
            pari_close();
            exit(0);
        }


        //-------------
        gel(aJI_vect, i) = aJI;
    } 

    // Modify the basis to satify relation. See the document "Split H^2" for details
    GEN uB, elem1, elem2, ideal1, ideal2, Gamma_a, i_x_u, i_x_B;
    
    for (i = 1; i < 3; i++) {
        
        uB = my_find_uB(Lrel, K, gel(gel(aJI_vect, i), 2), gel(gel(aJI_vect, i), 3), itos(p));
        
        Gamma_a = my_Gamma(Labs, sigma, gel(gel(aJI_vect, i), 1), itos(p));
        i_x_u = rnfeltup0(Lrel, gel(uB,1), 1);
        i_x_B = rnfidealup0(Lrel, gel(uB,2), 1);
        
        elem1 = nfdiv(Labs, Gamma_a, i_x_u);
        elem2 = gel(gel(aJI_vect, i), 1);
        ideal1 = gel(gel(aJI_vect, i), 2);
        ideal2 = idealdiv(Labs, gel(gel(aJI_vect, i), 3), i_x_B);
        
        gel(my_basis, i) = mkvec4(elem1, elem2, ideal1, ideal2);
        
    }

    my_basis = gerepilecopy(av, my_basis);
    return my_basis;
}


// ----- Basis version 2 ---------
// Faster but requires rnfisnorm

GEN my_find_basis_2 (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN p, GEN J_vect, GEN T) {
    pari_sp av = avma;
    int p_int = itos(p);
    // GEN T = rnfisnorminit(K, rnf_get_pol(LxRel), 1);
    GEN my_basis = zerovec(2), aJI_vect = zerovec(2);
    GEN Na_vect = zerovec(2);
    GEN a_vect = zerovec(2);
  
    GEN test;
    printf("Finding basis\n\n");
    int i;
    for (i = 1; i < 3; i++)
    {
        printf("%d\n\n", i);
        test = bnfisprincipal0(K, idealpow(K, gel(J_vect, i), p), 3);
        if (!my_QV_equal0(gel(test,1)))
        {
            printf(ANSI_COLOR_RED "ERROR in my_find_basis_2 \n\n" ANSI_COLOR_RESET);
            pari_close();
            exit(0);
        }
        
        gel(Na_vect, i) = nfinv(K, gel(test, 2));
        
        gel(a_vect, i) = algtobasis(Labs, rnfeltreltoabs(Lrel, gel(rnfisnorm(T, gel(Na_vect, i), 0),1)));
    }
    
    pari_printf(ANSI_COLOR_YELLOW "BASIS: \n\nNa_vect: %Ps\n\nJ_vect: %Ps\n\n" ANSI_COLOR_RESET, Na_vect, J_vect);

    // Given a_2 and J, finds I_2 such that div(a_2)+i(J)+(1-sigma)I_2 = 0
    GEN I_vect = my_find_I2 (Labs, Lrel, K, sigma, a_vect, J_vect, p_int);
    
    for (i = 1; i < 3; i++) {
        gel(I_vect, i) = idealred(Labs, gel(I_vect, i));
        gel(a_vect, i) = nfinv(Labs, gel(bnfisprincipal0(Labs, idealmul(Labs, rnfidealup0(Lrel, gel(J_vect, i), 1), my_1MS_ideal(Labs, sigma, gel(I_vect, i))), 3), 2));
        gel(aJI_vect, i) = mkvec3(gel(a_vect, i), gel(J_vect, i), gel(I_vect, i));

    }
    // Modify the basis to satify relation. See the document "Split H^2" for details
    GEN uB, elem1, elem2, ideal1, ideal2, Gamma_a, i_x_u, i_x_B;
    
    for (i = 1; i < 3; i++) {
        
        uB = my_find_uB(Lrel, K, gel(gel(aJI_vect, i), 2), gel(gel(aJI_vect, i), 3), p_int);
        
        Gamma_a = my_Gamma(Labs, sigma, gel(gel(aJI_vect, i), 1), p_int);
        i_x_u = rnfeltup0(Lrel, gel(uB,1), 1);
        i_x_B = rnfidealup0(Lrel, gel(uB,2), 1);
        
        elem1 = nfdiv(Labs, Gamma_a, i_x_u);
        elem2 = gel(gel(aJI_vect, i), 1);
        ideal1 = gel(gel(aJI_vect, i), 2);
        ideal2 = idealdiv(Labs, gel(gel(aJI_vect, i), 3), i_x_B);
        
        gel(my_basis, i) = mkvec4(elem1, elem2, ideal1, ideal2);
        
    }

    my_basis = gerepilecopy(av, my_basis);
    return my_basis;
}