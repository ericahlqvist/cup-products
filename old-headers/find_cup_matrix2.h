//Defines a matrix over F_2 with index (i*k, j) corresponding to 
//< x_i\cup x_k, (a_j, J_j)>

void my_cup_matrix_2 (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN J_vect, GEN units_mod_p, int r_rk, GEN p_power_units)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = (p_rk*(p_rk+1)/2);
    GEN cup_matrix = zerovec(nr_col);
    
    int i, j, k;
    for (i=1; i<(nr_col+1); ++i) {
        gel(cup_matrix, i) = zerovec(p_rk+glength(units_mod_p));
    }

    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        I_vect = my_find_I_vect3(Labs_cup, Lrel_cup, K, sigma_cup, J_vect, units_mod_p, p_int, p_power_units);
        
        printf("I_vect nr: %d\n\n", i);

        //Test that we get a correct I
        if (my_test_H90_ideal(Labs_cup, Lrel_cup, K, sigma_cup, I_vect, J_vect))
        {
            printf(ANSI_COLOR_GREEN "H90 test passed\n\n" ANSI_COLOR_RESET);
        }
        else {
            printf(ANSI_COLOR_RED "H90 test  failed\n\n" ANSI_COLOR_RESET);
            pari_close();
            exit(0);
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
                gmael2(cup_matrix, (2*p_rk-(i-2))*(i-1)/2+k-(i-1),j) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int));
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }

    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (i=1; i<p_rk+1; ++i) {
        for (k=i; k<p_rk+1; ++k) {
            pari_printf(ANSI_COLOR_CYAN "(%d,%d)  :  %Ps\n\n" ANSI_COLOR_RESET, i,k, gel(cup_matrix, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)));
        } 
    }
    
    pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, ZM_hnf((ZM_copy(cup_matrix))));

    printf("\n\n");
    char letters[] = "abcdefghijklmnopqr";
    printf(ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
    printf("{");
    for (j=1; j<r_rk+1; j++) {
        for (i=1; i<p_rk+1; ++i) {
            for (k=i; k<p_rk+1; k++) {
                if (gequal1(gel(gel(cup_matrix, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)), j))) {
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

void my_cup_matrix_2_transpose (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN J_vect, GEN Ja_vect, GEN units_mod_p, int r_rk, GEN p_power_units)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = (p_rk*(p_rk+1)/2);
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<(nr_row+1); ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
    }

    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        I_vect = my_find_I_vect3(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, units_mod_p, p_int, p_power_units);
        
        printf("I_vect nr: %d\n\n", i);

        // //Test that we get a correct I
        // if (my_test_H90_ideal(Labs_cup, Lrel_cup, K, sigma_cup, I_vect, J_vect))
        // {
        //     printf(ANSI_COLOR_GREEN "H90 test passed\n\n" ANSI_COLOR_RESET);
        // }
        // else {
        //     printf(ANSI_COLOR_RED "H90 test  failed\n\n" ANSI_COLOR_RESET);
        //     pari_close();
        //     exit(0);
        // }
        // //----- test done
        
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
                gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred(K,NIpJ), p_int, sigma), p_int));
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }

    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<nr_row+1; ++j) {
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

void my_cup_matrix_2_transpose_script (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN J_vect, GEN Ja_vect, GEN units_mod_p, int r_rk, GEN p_power_units)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = (p_rk*(p_rk+1)/2);
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<(nr_row+1); ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
    }

    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        I_vect = my_find_I_vect3(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, units_mod_p, p_int, p_power_units);
        
        //printf("I_vect nr: %d\n\n", i);

        // //Test that we get a correct I
        // if (my_test_H90_ideal(Labs_cup, Lrel_cup, K, sigma_cup, I_vect, J_vect))
        // {
        //     printf(ANSI_COLOR_GREEN "H90 test passed\n\n" ANSI_COLOR_RESET);
        // }
        // else {
        //     printf(ANSI_COLOR_RED "H90 test  failed\n\n" ANSI_COLOR_RESET);
        //     pari_close();
        //     exit(0);
        // }
        // //----- test done
        
        for (j=1; j<r_rk+1; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            //printf("I, %d to rel\n\n", j);
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
                gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred(K,NIpJ), p_int, sigma), p_int));
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }

    // printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    // for (j=1; j<nr_row+1; ++j) {
    //     pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix, j));
    // }
    
    
    printf(ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(cup_matrix)), p);
    pari_printf(ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    FILE *fptr;
    fptr = fopen("data/polynomials/ranks-rd-40-80-ov7.json", "a");

    pari_fprintf(fptr, "{\"pol\": \"%Ps\", \"cl_gp\": \"%Ps\", \"rel_rk\": \"%d\", \"mat_rank\": \"%ld\"},\n", nf_get_pol(bnf_get_nf(K)), bnf_get_cyc(K), r_rk, mat_rk);

    fclose(fptr);

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

void my_cup_matrix_3 (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = (p_rk*(p_rk+1)/2);
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<(nr_row+1); ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
    }

    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        printf(ANSI_COLOR_MAGENTA "-----------\n\n\nStarting i = %d\n\n\n-------------\n" ANSI_COLOR_RESET, i);
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        //I_vect = my_find_I_vect_full(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, p_int);
        I_vect = my_H90_vect(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, stoi(p_int));
        pari_printf(ANSI_COLOR_GREEN "-----------\n\n\nI_vect nr: %d\n\n\n-------------\n%Ps\n" ANSI_COLOR_RESET, i, I_vect);

        // //Test that we get a correct I
        // if (my_test_H90_ideal(Labs_cup, Lrel_cup, K, sigma_cup, I_vect, my_find_p_gens(K,  stoi(p_int))))
        // {
        //     printf(ANSI_COLOR_GREEN "H90 test passed\n\n" ANSI_COLOR_RESET);
        // }
        // else {
        //     printf(ANSI_COLOR_RED "H90 test  failed\n\n" ANSI_COLOR_RESET);
        //     pari_close();
        //     exit(0);
        // }
        // //----- test done
        
        for (j=1; j<r_rk+1; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            printf(ANSI_COLOR_MAGENTA "-----------\n\n\nStarting [i,j] = [%d, %d]\n\n\n-------------\n" ANSI_COLOR_RESET,i, j);
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            printf("I, %d to rel\n\n", j);
            
            if (p_int == 2) {
                NIpJ = idealmul(K, gel(gel(Ja_vect, j), 2), rnfidealnormrel(Lrel_cup, I_rel));
            }
            else {
                NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
            }
            
            for (k=i; k<p_rk+1; ++k) {
                printf(ANSI_COLOR_YELLOW "----------------\n\n\n\nStarting [i,j,k] = [%d, %d, %d]\n\n\n\n----------------\n" ANSI_COLOR_RESET, i,j,k);
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                sigma = gel(gel(K_ext, k), 3);
                gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred(K,NIpJ), p_int, sigma), p_int));
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }

    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<nr_row+1; ++j) {
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

void my_cup_matrix_3_script (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk, GEN ext_rks, GEN ext_cyc, GEN ext_pol)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = (p_rk*(p_rk+1)/2);
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<(nr_row+1); ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
    }

    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        I_vect = my_find_I_vect_full(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, p_int);
        
        printf("I_vect nr: %d\n\n", i);

        // //Test that we get a correct I
        // if (my_test_H90_ideal(Labs_cup, Lrel_cup, K, sigma_cup, I_vect, J_vect))
        // {
        //     printf(ANSI_COLOR_GREEN "H90 test passed\n\n" ANSI_COLOR_RESET);
        // }
        // else {
        //     printf(ANSI_COLOR_RED "H90 test  failed\n\n" ANSI_COLOR_RESET);
        //     pari_close();
        //     exit(0);
        // }
        // //----- test done
        
        for (j=1; j<r_rk+1; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            printf("I, %d to rel\n\n", j);
            
            if (p_int == 2) {
                NIpJ = idealmul(K, gel(gel(Ja_vect, j), 2), rnfidealnormrel(Lrel_cup, I_rel));
            }
            else {
                NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
            }
            
            for (k=i; k<p_rk+1; ++k) {
                
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                sigma = gel(gel(K_ext, k), 3);
                gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred(K,NIpJ), p_int, sigma), p_int));
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }

    

    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<nr_row+1; ++j) {
        pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix, j));
    }
    
    
    printf(ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(cup_matrix)), p);
    pari_printf(ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    FILE *fptr;
    fptr = fopen("data/discriminants/ranks-[2,2,2,2]_4_mod_16.json", "a");

    pari_fprintf(fptr, "{\"pol\": \"%Ps\", \"K_cyc\": \"%Ps\", \"ext_cyc\": \"%Ps\", \"ext_pol\": \"%Ps\"},\n", nf_get_pol(bnf_get_nf(K)), bnf_get_cyc(K), ext_cyc, ext_pol);

    fclose(fptr);

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

void my_quaternion_matrix (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup, J;
    int nr_col = (p_rk*(p_rk+1)/2);
    int nr_row = r_rk;
    int cup, B_x, B_y, sum;
    GEN cup_matrix = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<(nr_row+1); ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
    }

    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        I_vect = my_H90_vect(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, stoi(p_int));
        
        printf("I_vect nr: %d\n\n", i);

        // //Test that we get a correct I
        // if (my_test_H90_ideal(Labs_cup, Lrel_cup, K, sigma_cup, I_vect, J_vect))
        // {
        //     printf(ANSI_COLOR_GREEN "H90 test passed\n\n" ANSI_COLOR_RESET);
        // }
        // else {
        //     printf(ANSI_COLOR_RED "H90 test  failed\n\n" ANSI_COLOR_RESET);
        //     pari_close();
        //     exit(0);
        // }
        // //----- test done
        
        for (j=1; j<r_rk+1; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            J = gel(gel(Ja_vect, j), 2);
            printf("I, %d to rel\n\n", j);
            
            if (p_int == 2) {
                NIpJ = idealmul(K, J, rnfidealnormrel(Lrel_cup, I_rel));
            }
            else {
                NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
            }
            
            for (k=i; k<p_rk+1; ++k) {
                
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                sigma = gel(gel(K_ext, k), 3);

                cup = smodis(my_Artin_symbol(Labs, Lrel, K, idealred(K,NIpJ), p_int, sigma), p_int);
                B_y = smodis(my_Artin_symbol(Labs, Lrel, K, J, p_int, sigma), p_int);
                B_x = smodis(my_Artin_symbol(Labs_cup, Lrel_cup, K, J, p_int, sigma_cup), p_int);

                sum = cup+B_x+B_y;

                gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(sum % p_int);
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }

    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<nr_row+1; ++j) {
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
        
        // char letters[] = "abcdefghijklmnopqr";
        // printf(ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
        // int hnf_r_rk = glength(cup_hnf);
        // printf("{");
        // for (j=1; j<hnf_r_rk+1; j++) {
        //     for (i=1; i<p_rk+1; ++i) {
        //         for (k=i; k<p_rk+1; k++) {
        //             if (gequal1(gel(gel(cup_hnf, j), (2*p_rk-(i-2))*(i-1)/2+k-(i-1)))) {
        //                 if (i==k) {
        //                     printf("%c^2", letters[i-1]);
        //                 }
        //                 else {
        //                     printf("(%c,%c)", letters[i-1],letters[k-1]);
        //                 }
        //             }
        //         }
        //     }
        //     if (j==hnf_r_rk) {
        //         printf("}\n");
        //     }
        //     else {
        //         printf(", ");
        //     }
        // }
    }
    
} 