//Defines a matrix over F_2 with index (i*k, j) corresponding to 
//< x_i\cup x_k, (a_j, J_j)>

void my_cup_matrix (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN J_vect, GEN Ja_vect, GEN units_mod_p, int r_rk)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = (p_rk*p_rk);
    GEN cup_matrix = zerovec(nr_col);
    
    int i, j, k;
    for (i=1; i<(p_rk*p_rk+1); ++i) {
        gel(cup_matrix, i) = zerovec(p_rk+glength(units_mod_p));
    }
    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        I_vect = my_find_I_vect2(Labs_cup, Lrel_cup, K, sigma_cup, J_vect, Ja_vect, units_mod_p, p_int);
        
        
        printf("I_vect nr: %d\n\n", i);

        //Test that we get a correct I
        if (my_test_H90_ideal(Labs_cup, Lrel_cup, K, sigma_cup, I_vect, J_vect))
        {
            printf(ANSI_COLOR_GREEN "H90 test passed\n\n" ANSI_COLOR_RESET);
        }
        else {
            printf(ANSI_COLOR_RED "H90 test  failed\n\n" ANSI_COLOR_RESET);
        }
        //----- test done
        int added = 1+glength(units_mod_p);
        
        for (j=1; j<p_rk+added; ++j) {
            // evaluate cup on j:th basis class (a,J) 
           
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            printf("I, %d to rel\n\n", j);
            if (j <= glength(J_vect)) {
                if (p_int == 2) {
                    NIpJ = idealmul(K, gel(J_vect, j), rnfidealnormrel(Lrel_cup, I_rel));
                }
                else {
                    NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
                }
                
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
  
    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (i=1; i<p_rk+1; ++i) {
        for (k=i; k<p_rk+1; ++k) {
            pari_printf(ANSI_COLOR_CYAN "(%d,%d)  :  %Ps\n\n" ANSI_COLOR_RESET, i,k, gel(cup_matrix, (i-1)*p_rk+k));
        } 
    }
    printf(ANSI_COLOR_YELLOW "rank of image: ");
    pari_printf(ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, FpM_rank((ZM_copy(cup_matrix)), p));

    printf("\n\n");

    char letters[] = "abcdefghijklmnopqr";
    printf(ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
    printf("{");
    for (j=1; j<r_rk+1; j++) {
        for (i=1; i<p_rk+1; ++i) {
            for (k=i; k<p_rk+1; k++) {
                if (gequal1(gel(gel(cup_matrix, (i-1)*p_rk+k), j))) {
                    if (i==k) {
                        printf("%c^2", letters[i-1]);
                    }
                    else {
                        printf("(%c,%c)", letters[i-1],letters[k-1]);
                    }
                }
            }
        }
        if (j==r_rk) {
            printf("}\n");
        }
        else {
            printf(", ");
        }
    }
    
}    

void my_cup_matrix_transpose (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN J_vect, GEN Ja_vect, GEN units_mod_p, int r_rk)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = (p_rk*(p_rk+1)/2);
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<nr_row+1; ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
    }
    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        I_vect = my_find_I_vect2(Labs_cup, Lrel_cup, K, sigma_cup, J_vect, Ja_vect, units_mod_p, p_int);
        
        
        printf("I_vect nr: %d\n\n", i);

        //Test that we get a correct I
        if (my_test_H90_ideal(Labs_cup, Lrel_cup, K, sigma_cup, I_vect, J_vect))
        {
            printf(ANSI_COLOR_GREEN "H90 test passed\n\n" ANSI_COLOR_RESET);
        }
        else {
            printf(ANSI_COLOR_RED "H90 test  failed\n\n" ANSI_COLOR_RESET);
        }
        //----- test done
        
        for (j=1; j<r_rk+1; ++j) {
            // evaluate cup on j:th basis class (a,J) 
           
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            printf("I, %d to rel\n\n", j);
            if (j <= glength(J_vect)) {
                if (p_int == 2) {
                    NIpJ = idealmul(K, gel(J_vect, j), rnfidealnormrel(Lrel_cup, I_rel));
                }
                else {
                    NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
                }
            }
            else {
                NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
            }
            
            
            for (k=i; k<p_rk+1; ++k) {
                
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                sigma = gel(gel(K_ext, k), 3);
                
                gmael2(cup_matrix,j,(2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int));
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }
  
    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<r_rk+1; ++j) {

        pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix, j));

    }
  
    printf(ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(cup_matrix)), p);
    pari_printf(ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    if (mat_rk > 0) {
        GEN cup_hnf = FpM_red(hnf((ZM_copy(cup_matrix))),p);
        printf(ANSI_COLOR_YELLOW "Hermite normal form:  \n\n" ANSI_COLOR_RESET);
        for (i=1;i<glength(cup_hnf)+1;++i) {
            pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_hnf, i));
        }
        
        printf("\n\n");
        
        char letters[] = "abcdefghijklmnopqr";
        printf(ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
        int hnf_r_rk = glength(cup_hnf);
        printf("{");
        for (j=1; j<hnf_r_rk+1; j++) {
            for (i=1; i<p_rk+1; ++i) {
                for (k=i; k<p_rk+1; k++) {
                    if (gequal1(gel(gel(cup_hnf, j), (2*p_rk-(i-2))*(i-1)/2+k-(i-1)))) {
                        if (i==k) {
                            printf("%c^2", letters[i-1]);
                        }
                        else {
                            printf("(%c,%c)", letters[i-1],letters[k-1]);
                        }
                    }
                }
            }
            if (j==hnf_r_rk) {
                printf("}\n");
            }
            else {
                printf(", ");
            }
        }
    }
}    