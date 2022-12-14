//Defines a matrix over F_2 with index (i*k, j) corresponding to 
//< x_i\cup x_k, (a_j, J_j)>

GEN my_cup_matrix (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN J_vect, GEN Ja_vect, GEN units_mod_2)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = (p_rk*p_rk);
    GEN cup_matrix = zerovec(nr_col);
    
    int i, j, k;
    for (i=1; i<(p_rk*p_rk+1); ++i) {
        gel(cup_matrix, i) = zerovec(p_rk+glength(units_mod_2));
    }
    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        if (p_int ==2) {
            I_vect = my_find_I_vect2(Labs_cup, Lrel_cup, K, sigma_cup, J_vect, Ja_vect, units_mod_2, p_int);
        }
        else {
            I_vect = my_find_I_vect(Labs_cup, Lrel_cup, K, sigma_cup, J_vect);
        }
        
        
        //pari_printf("I_vect: %Ps\n\n", I_vect);

        // //Test that we get a correct I
        // if (my_test_H90_ideal(Labs_cup, Lrel_cup, K, sigma_cup, I_vect, J_vect))
        // {
        //     printf(ANSI_COLOR_GREEN "H90 test passed\n\n" ANSI_COLOR_RESET);
        // }
        // else {
        //     printf(ANSI_COLOR_RED "H90 test  failed\n\n" ANSI_COLOR_RESET);
        // }
        // //----- test done
        int added;
        if (p_int == 2) {
            added = 1+glength(units_mod_2);
        }
        else {
            added = 1;
        }
        for (j=1; j<p_rk+added; ++j) {
            // evaluate cup on j:th basis class (a,J) 
           
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            if (j <= glength(J_vect)) {
                NIpJ = idealmul(K, rnfidealnormrel(Lrel_cup, I_rel), gel(J_vect, j));
            }
            else {
                NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
            }
            
            
            for (k=1; k<p_rk+1; ++k) {
                
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                sigma = gel(gel(K_ext, k), 3);
                
                gmael2(cup_matrix, (i-1)*p_rk+k,j) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int));
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }
    return cup_matrix;
}    