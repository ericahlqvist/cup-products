//Defines a matrix over F_2 with index (i*k, j) corresponding to 
//< x_i\cup x_k, (a_j, J_j)>


int my_cup_matrix_3 (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk)
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
    for (i=1; i<=p_rk; ++i) {
        printf(ANSI_COLOR_MAGENTA "-----------\n\n\nStarting i = %d\n\n\n-------------\n" ANSI_COLOR_RESET, i);
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        //I_vect = my_find_I_vect_full(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, p_int);
        I_vect = my_H90_vect(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, stoi(p_int), 1);
        if (gequal(I_vect,stoi(-1)))
        {
            I_vect = my_find_I_vect_full(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, p_int);
            // return 111;
        }
        
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
        
        for (j=1; j<=r_rk; ++j) {
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
            
            for (k=i; k<=p_rk; ++k) {
                printf(ANSI_COLOR_YELLOW "----------------\n\n\n\nStarting [i,j,k] = [%d, %d, %d]\n\n\n\n----------------\n" ANSI_COLOR_RESET, i,j,k);
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                sigma = gel(gel(K_ext, k), 3);
                gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred0(K,NIpJ, NULL), p_int, sigma), p_int));
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }

    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<=nr_row; ++j) {
        pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix, j));
    }
    
    
    printf(ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(cup_matrix)), p);
    pari_printf(ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    if (mat_rk > 0) {
        GEN cup_hnf = FpM_red(hnf((ZM_copy(cup_matrix))),p);
        printf(ANSI_COLOR_YELLOW "Hermite normal form:  \n\n" ANSI_COLOR_RESET);
        for (i=1;i<lg(cup_hnf);++i) {
            pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_hnf, i));
        }
        
        printf("\n\n");
        
        char letters[] = "abcdefghijklmnopqr";
        printf(ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
        int hnf_r_rk = glength(cup_hnf);
        printf("{");
        for (j=1; j<=hnf_r_rk; j++) {
            for (i=1; i<=p_rk; ++i) {
                for (k=i; k<=p_rk; k++) {
                    if (!gequal0(gel(gel(cup_hnf, j), (2*p_rk-(i-2))*(i-1)/2+k-(i-1)))) {
                        if (i==k && p_int==2) {
                            pari_printf("%c^(%d*%Ps)", letters[i-1], p_int, gel(gel(cup_hnf, j), (2*p_rk-(i-2))*(i-1)/2+k-(i-1)));
                        }
                        if (i<k) {
                            pari_printf("(%c,%c)^%Ps", letters[i-1],letters[k-1], gel(gel(cup_hnf, j), (2*p_rk-(i-2))*(i-1)/2+k-(i-1)));
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
    return mat_rk;
} 

int my_relations (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk)
{
    printf("\n--------------------------\nStart: my_relations\n--------------------------\n\n");
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;

    //-------------------------------------------------
    // Define a matrix with r_rk rows and columns indexed by (i, k) for 1 <= i <= k <= p_rk. 

    int nr_col = (p_rk*(p_rk+1)/2), nr_col_full = p_rk*p_rk;
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row), cup_matrix_full = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<=nr_row; ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
        gel(cup_matrix_full, j) = zerovec(nr_col_full);
    }
    //-------------------------------------------------

    // <x_i cup x_k, (a_j, J_j)> = < x_k , ... > ( = Artin symbol) ----->  i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<=p_rk; ++i) {
        //printf(ANSI_COLOR_MAGENTA "-----------\n\n\nStarting i = %d\n\n\n-------------\n" ANSI_COLOR_RESET, i);
        Labs_cup = gmael(K_ext, i, 1);
        Lrel_cup = gmael(K_ext, i, 2);
        sigma_cup = gmael(K_ext, i, 3);

        // Artin symbol test
        //my_test_artin_symbol(Labs_cup, Lrel_cup, K, p_int, sigma_cup);
        // my_test_artin_on_norms(Labs_cup, Lrel_cup, K, p_int, sigma_cup);
        // I_vect corresp. to i:th extension
        // I_vect = my_find_I_vect_full(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, p_int);

        I_vect = my_H90_vect(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, stoi(p_int), 1);
        if (gequal(I_vect,gen_m1))
        {
            I_vect = my_find_I_vect_full(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, p_int);
            // return 111;
        }
        
        //pari_printf(ANSI_COLOR_GREEN "-----------\n\n\nI_vect nr: %d\n\n\n-------------\n%Ps\n" ANSI_COLOR_RESET, i, I_vect);

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
        
        for (j=1; j<=r_rk; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            // printf(ANSI_COLOR_MAGENTA "-----------\n\n\nStarting [i,j] = [%d, %d]\n\n\n-------------\n" ANSI_COLOR_RESET,i, j);
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            // printf("I, %d to rel\n\n", j);
            
            if (p_int == 2) {
                NIpJ = idealmul(K, gel(gel(Ja_vect, j), 2), rnfidealnormrel(Lrel_cup, I_rel));
            }
            else {
                NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
            }
            
            for (k=1; k<=p_rk; ++k) {
                // printf(ANSI_COLOR_YELLOW "----------------\n\n\n\nStarting [i,j,k] = [%d, %d, %d]\n\n\n\n----------------\n" ANSI_COLOR_RESET, i,j,k);
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                sigma = gel(gel(K_ext, k), 3);
                if (i<k) {
                    gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred0(K,NIpJ, NULL), p_int, sigma), p_int));
                    gmael2(cup_matrix_full, j, p_rk*(k-1)+i) = gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1));
                }
                if (i==k) {
                    gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, gel(gel(Ja_vect, j), 2), p_int, sigma), p_int));
                    gmael2(cup_matrix_full, j, p_rk*(k-1)+i) = gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1));
                }
                if (i>k) {
                    gmael2(cup_matrix_full, j, p_rk*(k-1)+i) = stoi(smodis(gneg(my_Artin_symbol(Labs, Lrel, K, idealred0(K,NIpJ, NULL), p_int, sigma)), p_int));
                }
                
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }
    // ------------------------
    //Check symmetry:
    printf(ANSI_COLOR_YELLOW "Full Cup Matrix:  (for (i,k) with i>k we have - (x_i cup x_k) instead of x_i cup x_k in order to see the symmetry).\n\n" ANSI_COLOR_RESET);
    for (j=1; j<nr_row+1; ++j) {
        pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix_full, j));
    }
    int sym_check = 1;
    GEN fail;
    for (i = 1; i <= p_rk; i++)
    {
        for (j = 1; j < r_rk; j++)
        {
            for (k = 1; k <= p_rk; k++)
            {
                if (!gequal(gmael2(cup_matrix_full, j, p_rk*(k-1)+i), gmael2(cup_matrix_full, j, p_rk*(i-1)+k)))
                {
                    // gel(fail, 1) = stoi(i);
                    // gel(fail, 2) = stoi(k);
                    fail = mkvec2s(i, k);
                    
                    sym_check = 0;
                    break;
                }
                
            }
        }
    }
    switch (sym_check)
    {
    case 1:
        printf(ANSI_COLOR_GREEN "\nSymmetry test passed\n\n" ANSI_COLOR_RESET);
        break;
    case 0:
        pari_printf(ANSI_COLOR_RED "\nSymmetry test failed at: %Ps\n\n" ANSI_COLOR_RESET, fail);
        pari_close();
        exit(111);
        break;
    default:
        break;
    }
    
    printf("-----------------------------------------------------------------------------------------------\n\n");
    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<=nr_row; ++j) {
        pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix, j));
    }
    printf(ANSI_COLOR_RED "For indices (i,i) we have B(x_i) instead of x_i cup x_i to get the correct presentation of Q_2.\n\n" ANSI_COLOR_RESET);
    
    printf(ANSI_COLOR_CYAN "This determines the second quotient Q_2 for the lower p-central series and not only the Zassenhaus quotient ZQ_2.\n\n" ANSI_COLOR_RESET);
    printf("-----------------------------------------------------------------------------------------------\n\n");
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
        printf("\"\"\"[");
        for (j=1; j<=hnf_r_rk; j++) {
            for (i=1; i<=p_rk; ++i) {
                for (k=i; k<=p_rk; k++) {
                    if (!gequal0(gel(gel(cup_hnf, j), (2*p_rk-(i-2))*(i-1)/2+k-(i-1)))) {
                        if (i==k) {
                            printf("%c^%ld", letters[i-1], p_int*smodis(gneg(gmael2(cup_hnf, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1))), p_int));
                        }
                        if (i<k) {
                            printf("%c_%c^%ld", letters[i-1],letters[k-1], smodis(gneg(gmael2(cup_hnf, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1))), p_int));
                        }
                    }
                }
            }
            if (j==hnf_r_rk) {
                printf("]\"\"\"\n");
            }
            else {
                printf(", ");
            }
        }
    }
    printf("\n--------------------------\nEnd: my_relations\n--------------------------\n\n");
    printf("-----------------------------------------------------------------------------------------------\n\n");
    
    return mat_rk;
}
void my_cup_matrix_and_relators (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = (p_rk*(p_rk+1)/2);
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<=nr_row; ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
    }

    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<=p_rk; ++i) {
        //printf(ANSI_COLOR_MAGENTA "-----------\n\n\nStarting i = %d\n\n\n-------------\n" ANSI_COLOR_RESET, i);
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        //I_vect = my_find_I_vect_full(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, p_int);
        I_vect = my_H90_vect(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, stoi(p_int), 1);
        //pari_printf(ANSI_COLOR_GREEN "-----------\n\n\nI_vect nr: %d\n\n\n-------------\n%Ps\n" ANSI_COLOR_RESET, i, I_vect);

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
        
        for (j=1; j<=r_rk; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            //printf(ANSI_COLOR_MAGENTA "-----------\n\n\nStarting [i,j] = [%d, %d]\n\n\n-------------\n" ANSI_COLOR_RESET,i, j);
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            //printf("I, %d to rel\n\n", j);
            
            if (p_int == 2) {
                NIpJ = idealmul(K, gel(gel(Ja_vect, j), 2), rnfidealnormrel(Lrel_cup, I_rel));
            }
            else {
                NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
            }
            
            for (k=i; k<=p_rk; ++k) {
                //printf(ANSI_COLOR_YELLOW "----------------\n\n\n\nStarting [i,j,k] = [%d, %d, %d]\n\n\n\n----------------\n" ANSI_COLOR_RESET, i,j,k);
                // take Artin symbol with resp. to k:th extension
                Labs = gmael(K_ext, k, 1);
                Lrel = gmael(K_ext, k, 2);
                sigma = gmael(K_ext, k, 3);
                gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred0(K,NIpJ, NULL), p_int, sigma), p_int));
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }

    

    // printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    // for (j=1; j<nr_row+1; ++j) {
    //     pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix, j));
    // }
    
    
    // printf(ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(cup_matrix)), p);
    //pari_printf(ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    FILE *fptr;
    fptr = fopen("data/discriminants/relators-[2,2,2,2]_4_mod_16.json", "a");


    if (mat_rk > 0) {
        GEN cup_hnf = FpM_red(hnf((ZM_copy(cup_matrix))),p);
        // printf(ANSI_COLOR_YELLOW "Hermite normal form:  \n\n" ANSI_COLOR_RESET);
        // for (i=1;i<glength(cup_hnf)+1;++i) {
        //     pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_hnf, i));
        // }
        
        //printf("\n\n");
        
        char letters[] = "abcdefghijklmnopqr";
        //printf(ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
        int hnf_r_rk = glength(cup_hnf);
        pari_fprintf(fptr, "{\"pol\": \"%Ps\",\"2-rk\": \"%d\", \"relations\": \"[", nf_get_pol(bnf_get_nf(K)), hnf_r_rk);
        
        for (j=1; j<=hnf_r_rk; j++) {
            for (i=1; i<=p_rk; ++i) {
                for (k=i; k<=p_rk; k++) {
                    if (gequal1(gel(gel(cup_hnf, j), (2*p_rk-(i-2))*(i-1)/2+k-(i-1)))) {
                        if (i==k) {
                            pari_fprintf(fptr, "%c^2*", letters[i-1]);
                            
                        }
                        else {
                            pari_fprintf(fptr, "%c_%c*", letters[i-1],letters[k-1]);
                        }
                    }
                }
            }
            if (j==hnf_r_rk) {
                pari_fprintf(fptr, "]\"},\n");
            }
            else {
                pari_fprintf(fptr, ", ");
            }
        }
    }
    fclose(fptr);
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
        I_vect = my_H90_vect(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, stoi(p_int), 1);
        
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
            J = gel(gel(Ja_vect, j), 2);
            //printf("I, %d to rel\n\n", j);
            
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

                cup = smodis(my_Artin_symbol(Labs, Lrel, K, idealred0(K,NIpJ, NULL), p_int, sigma), p_int);
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

int my_cup_matrix_12 (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, int r_rk, GEN p_ClFld_pol)
{
    // This versions assumes that p=2 and that the torsion units mod 2 is generated by -1
    GEN NIpJ, Labs_x, Lrel_x, sigma_x, Labs_y, Lrel_y, sigma_y, Labs_z, Lrel_z, sigma_z, b, a, J, I, c, class_group;
    GEN x = pol_x(fetch_user_var("x"));
    // int vec_len = (p_rk-1)*p_rk*(p_rk+1)/6;
    GEN cup_matrix = zerovec(0);
    
    
    int i, j, k;

    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        printf(ANSI_COLOR_MAGENTA "-----------\n\n\nStarting i = %d\n\n\n-------------\n" ANSI_COLOR_RESET, i);
        Labs_x = gel(gel(K_ext, i), 1);
        Lrel_x = gel(gel(K_ext, i), 2);
        sigma_x = gel(gel(K_ext, i), 3);
        
        for (j=i; j<p_rk+1; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            printf(ANSI_COLOR_MAGENTA "-----------\n\n\nStarting [i,j] = [%d, %d]\n\n\n-------------\n" ANSI_COLOR_RESET,i, j);
            Labs_y = gel(gel(K_ext, j), 1);
            Lrel_y = gel(gel(K_ext, j), 2);
            sigma_y = gel(gel(K_ext, j), 3);
            
            for (k=i; k<p_rk+1; ++k) {
                printf(ANSI_COLOR_YELLOW "----------------\n\n\n\nStarting [i,j,k] = [%d, %d, %d]\n\n\n\n----------------\n" ANSI_COLOR_RESET, i,j,k);
                // take Artin symbol with resp. to k:th extension
                Labs_z = gel(gel(K_ext, k), 1);
                Lrel_z = gel(gel(K_ext, k), 2);
                sigma_z = gel(gel(K_ext, k), 3);

                
                pari_printf("rel_pol: %Ps\n", gel(p_ClFld_pol, i));
                c = nfmul(Labs_x, stoi(-1), nfdiv(Labs_x, nfadd(Labs_x, galoisapply(Labs_x, sigma_x, rnfeltreltoabs(Lrel_x, x)), rnfeltreltoabs(Lrel_x, x)), gen_2));
                b = nfadd(Labs_x, c, algtobasis(Labs_x, rnfeltreltoabs(Lrel_x, x)));
                pari_printf("b: %Ps\n", basistoalg(Labs_x, b));
                pari_printf("sigma x: %Ps\n", galoisapply(Labs_x, sigma_x, rnfeltreltoabs(Lrel_x, x)));
                pari_printf("x+sigma x: %Ps\n", nfadd(Labs_x, galoisapply(Labs_x, sigma_x, rnfeltreltoabs(Lrel_x, x)), rnfeltreltoabs(Lrel_x, x)));
                pari_printf("(1-s)b: %Ps\n", my_1MS_elt(Labs_x, sigma_x, b));
                if (!gequal(basistoalg(Labs_x, my_1MS_elt(Labs_x, sigma_x, b)), stoi(-1)))
                {
                    printf("(1-s)b isn't -1\n");
                    pari_close();
                    exit(111);
                }
                
                a = rnfeltdown0(Lrel_x, rnfeltabstorel(Lrel_x, nfpow(Labs_x, b, stoi(-2))), 1);
                J = rnfidealdown(Lrel_x, rnfidealabstorel(Lrel_x, idealhnf0(Labs_x, b, NULL)));
                I = my_H90_12(Labs_y, Lrel_y, K, sigma_y, a, J, p, 1);
                if (gequal(I, stoi(-1)))
                {
                    class_group = my_get_clgp(Labs_y);
                    I = gel(my_find_I_full(Labs_y, Lrel_y, K, sigma_y, idealhnf0(Labs_x, b, NULL), a, class_group, p_int),2);
                }
                

                NIpJ = idealmul(K, J, rnfidealnormrel(Lrel_y, rnfidealabstorel(Lrel_y, I)));
                
                cup_matrix = shallowconcat(cup_matrix, mkvec(stoi(smodis(my_Artin_symbol(Labs_z, Lrel_z, K, idealred0(K,NIpJ, NULL), p_int, sigma_z), p_int))));
                
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }

    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, cup_matrix);
    
    
    long mat_rk;
    printf(ANSI_COLOR_YELLOW "rank: ");
    if (my_QV_equal0(cup_matrix))
    {
        mat_rk = 0;
    }
    else {
        mat_rk = 1;
    }
    
    pari_printf(ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    
    return mat_rk;
} 

int my_massey_matrix (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk, int n)
{
    GEN NIpJ, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup, I_prime_vect;
    int nr_col = p_rk*p_rk;
    int nr_row = r_rk;
    GEN massey_matrix = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<(nr_row+1); ++j) {
        gel(massey_matrix, j) = zerovec(nr_col);
    }

    //pari_printf("cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        printf("Start round %d/%d\n\n", i, p_rk);
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        I_prime_vect = my_H90_vect(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, stoi(p_int), n);
        printf("I'_vect found\n\n");

        for (j=1; j<r_rk+1; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            // pari_printf("I'_vect[%d]: %Ps\n\n", j, gel(I_prime_vect, j));
            // printf("j=[%d]\n\n", j);
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_prime_vect, j));
            //printf("I, %d to rel\n\n", j);
            
            if (p_int == 3) {
                NIpJ = idealmul(K, gel(gel(Ja_vect, j),2), rnfidealnormrel(Lrel_cup, I_rel));
            }
            else {
                NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
            }
            
            
            for (k=1; k<p_rk+1; ++k) {
                printf(ANSI_COLOR_GREEN "Start: [%d,%d,%d]\n" ANSI_COLOR_RESET, i,j,k);
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                sigma = gel(gel(K_ext, k), 3);
                gmael2(massey_matrix, j, p_rk*(i-1)+k) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred(K,NIpJ), p_int, sigma), p_int));
                printf(ANSI_COLOR_GREEN "End: [%d,%d,%d]\n\n" ANSI_COLOR_RESET, i,j,k);
                //pari_printf("ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }
    printf(ANSI_COLOR_MAGENTA "\n-------------------------------------------------------\n %d-fold Massey products of the form < x, x, ..., x, y >\n-------------------------------------------------------\n\n" ANSI_COLOR_RESET, n+1);
    pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, massey_matrix);
    printf(ANSI_COLOR_YELLOW "Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<nr_row+1; ++j) {
        pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(massey_matrix, j));
    }
    
    
    printf(ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(massey_matrix)), p);
    pari_printf(ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    // FILE *fptr;
    // fptr = fopen("data/polynomials/ranks_S3_[2,2]", "a");

    // pari_fprintf(fptr, "pol: %Ps,    rel_rk: %d,     mat_rank: %ld\n", nf_get_pol(bnf_get_nf(K)), r_rk, mat_rk);

    // fclose(fptr);

    if (mat_rk > 0) {
        GEN massey_hnf = FpM_red(hnf((ZM_copy(massey_matrix))),p);
        printf(ANSI_COLOR_YELLOW "Hermite normal form:  \n\n" ANSI_COLOR_RESET);
        for (i=1;i<glength(massey_hnf)+1;++i) {
            pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(massey_hnf, i));
        }
        
        printf("\n\n");
        
        char letters[] = "abcdefghijklmnopqr";
        printf(ANSI_COLOR_YELLOW "Massey relations:  \n\n" ANSI_COLOR_RESET);
        int hnf_r_rk = glength(massey_hnf);
        printf("{");
        for (j=1; j<hnf_r_rk+1; j++) {
            for (i=1; i<p_rk+1; ++i) {
                for (k=1; k<p_rk+1; k++) {
                    if (!gequal0(gel(gel(massey_hnf, j), p_rk*(i-1)+k))) {
                        if (p_int==3) {
                            if (i==k) {
                                pari_printf("%c^3*%Ps", letters[i-1], gel(gel(massey_hnf, j), p_rk*(i-1)+k));
                            }
                            else {
                                if (i<k)
                                {
                                    pari_printf("((%c,%c), %c)^%Ps", letters[i-1],letters[k-1], letters[i-1], gneg(gel(gel(massey_hnf, j), p_rk*(i-1)+k)));
                                }
                                else {
                                    pari_printf("((%c,%c), %c)^%Ps", letters[i-1],letters[k-1], letters[i-1], gel(gel(massey_hnf, j), p_rk*(i-1)+k));
                                }
                                
                            }
                        }   
                        else {
                            if (i<k)
                            {
                                pari_printf("((%c,%c), %c)^%Ps", letters[i-1],letters[k-1], letters[i-1], gneg(gel(gel(massey_hnf, j), p_rk*(i-1)+k)));
                            }
                            if (i>k) {
                                pari_printf("((%c,%c), %c)^%Ps", letters[i-1],letters[k-1], letters[i-1], gel(gel(massey_hnf, j), p_rk*(i-1)+k));
                            }
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
        printf("\n\n");
    }

    int is_zero;
    GEN next_col = zerocol(r_rk);
    if (mat_rk > 0) {
        for (i=1; i<p_rk+1; ++i) {
            for (k=1; k<p_rk+1; ++k) {
                is_zero = 1;
                for (j=1; j<r_rk+1; ++j) {
                    if (!gequal0(gmael2(massey_matrix, j, p_rk*(i-1)+k)))
                    {
                        is_zero = 0;
                        break;
                    }
                }
                //printf("< %d, %d > is zero: %d\n", i, k, is_zero);
                if (i!=k && is_zero)
                {
                    printf("Start round %d/%d\n\n", i, p_rk);
                    Labs_cup = gel(gel(K_ext, i), 1);
                    Lrel_cup = gel(gel(K_ext, i), 2);
                    sigma_cup = gel(gel(K_ext, i), 3);

                    I_prime_vect = my_H90_vect(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, stoi(p_int), n+2);
                    printf("I'_vect found\n\n");
   
                    // take Artin symbol with resp. to k:th extension
                    Labs = gel(gel(K_ext, k), 1);
                    Lrel = gel(gel(K_ext, k), 2);
                    sigma = gel(gel(K_ext, k), 3);
                    for (j=1; j<r_rk+1; ++j) {
                        I_rel = rnfidealabstorel(Lrel_cup, gel(I_prime_vect, j));
                        //printf("I, %d to rel\n\n", j);
                        
                        if (p_int == 3) {
                            NIpJ = idealmul(K, gel(gel(Ja_vect, j),2), rnfidealnormrel(Lrel_cup, I_rel));
                        }
                        else {
                            NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
                        }
                        gel(next_col, j) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred(K,NIpJ), p_int, sigma), p_int));
                    }
                    printf(ANSI_COLOR_YELLOW "\n-----------------------------------------\n%d-fold Massey < x, x, ..., x, y >: \n-----------------------------------------\n\n" ANSI_COLOR_RESET, n+3); 
                    pari_printf(ANSI_COLOR_CYAN "%Ps\n\n-------------------------\n\n" ANSI_COLOR_RESET, next_col);
                }
                
            }
        }
        
    }
    return mat_rk;
} 


//------------------------------------------------
// Parallel version of my_relations
//------------------------------------------------

/*
GP;install("compute_my_relations","Gp","./stdin.so");
*/
// Define the function prototype
int my_relations_par(GEN K_ext, GEN K, GEN p, int p_rk, GEN Ja_vect, int r_rk);

// Wrapper function for parallel computation
GEN compute_my_relations(long i, GEN args) {
    printf("\n--------------------------\nStart: compute_my_relations\n--------------------------\n\n");

    GEN K_ext   = gel(args, 1);
    GEN K       = gel(args, 2);
    GEN p       = gel(args, 3); 
    GEN gp_rk   = gel(args, 4);
    GEN Ja_vect = gel(args, 5);
    GEN gr_rk   = gel(args, 6);

    pari_sp av = avma;
    GEN result = zerovec(2); // A vector to store results
    
    int p_rk = itos(gp_rk), r_rk = itos(gr_rk);
    
    
    GEN s_cup_matrix = zerovec(r_rk);
    GEN s_cup_matrix_full = zerovec(r_rk);

    GEN Labs_cup, Lrel_cup, sigma_cup, I_vect, I_rel, NIpJ, Labs, Lrel, sigma;
    long j, k;

    for (j = 1; j <= r_rk; ++j) {
        gel(s_cup_matrix, j) = zerovec((p_rk * (p_rk + 1)) / 2);
        gel(s_cup_matrix_full, j) = zerovec(p_rk * p_rk);
    }

    Labs_cup = gmael(K_ext, i, 1);
    Lrel_cup = gmael(K_ext, i, 2);
    sigma_cup = gmael(K_ext, i, 3);
    
    // Artin symbol test
    //my_test_artin_symbol(Labs_cup, Lrel_cup, K, itos(p), sigma_cup);
    I_vect = my_H90_vect(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, p, 1);
    if (gequal(I_vect, gen_m1)) {
        I_vect = my_find_I_vect_full(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, itos(p));
    }

    for (j = 1; j < r_rk + 1; ++j) {
        
        I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
        if (itos(p) == 2) {
            NIpJ = idealmul(K, gel(gel(Ja_vect, j), 2), rnfidealnormrel(Lrel_cup, I_rel));
        } else {
            NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
        }
        for (k = 1; k < p_rk + 1; ++k) {
            
            Labs = gel(gel(K_ext, k), 1);
            Lrel = gel(gel(K_ext, k), 2);
            sigma = gel(gel(K_ext, k), 3);
            if (i < k) {
                
                gmael2(s_cup_matrix, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, idealred0(K, NIpJ, NULL), itos(p), sigma), itos(p)));
                gmael2(s_cup_matrix_full, j, p_rk * (k - 1) + i) = gmael2(s_cup_matrix, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1));
                
            }
            if (i == k) {
                
                gmael2(s_cup_matrix, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1)) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, gel(gel(Ja_vect, j), 2), itos(p), sigma), itos(p)));
                gmael2(s_cup_matrix_full, j, p_rk * (k - 1) + i) = gmael2(s_cup_matrix, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1));
                
            }
            if (i > k) {
                
                gmael2(s_cup_matrix_full, j, p_rk * (k - 1) + i) = stoi(smodis(gneg(my_Artin_symbol(Labs, Lrel, K, idealred0(K, NIpJ, NULL), itos(p), sigma)), itos(p)));
                
            }
        }
    }

    gel(result, 1) = s_cup_matrix;
    
    gel(result, 2) = s_cup_matrix_full;
    
    printf("\n--------------------------\nEnd: compute_my_relations\n--------------------------\n\n");
    result = gerepilecopy(av, result);
    return result;
}

int my_relations_par(GEN K_ext, GEN K, GEN p, int p_rk, GEN Ja_vect, int r_rk) {
    printf("\n--------------------------\nStart: my_relations_par\n--------------------------\n\n");
    //GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;

    int nr_col = (p_rk * (p_rk + 1) / 2), nr_col_full = p_rk * p_rk;
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row), cup_matrix_full = zerovec(nr_row);
    
    for (int j = 1; j < nr_row + 1; ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
        gel(cup_matrix_full, j) = zerovec(nr_col_full);
    }

    long i, taskid, pending;
    GEN done, res, res_full;

    struct pari_mt pt;
    
    GEN args = mkvecn(6, K_ext, K, p, stoi(p_rk), Ja_vect, stoi(r_rk));
    //pari_printf("strtoclosure: %Ps\n", strtoclosure("_worker", 1, args));
    mt_queue_start(&pt, strtoclosure("_worker", 1, args));
    
    for (i = 1; i <= p_rk || pending; i++)
    { 
        printf("for i = %ld <= %d\n", i, p_rk);
        mt_queue_submit(&pt, i, i<=p_rk? mkvecs(i): NULL);
        done = mt_queue_get(&pt, &taskid, &pending);
        if (done) {
            printf("Done i = %ld \n", i);
            res = gel(done, 1);
            res_full = gel(done, 2);
            
            for (int j = 1; j <= r_rk; ++j) {
                for (int k = 1; k <= nr_col; ++k) {
                    
                    gmael2(cup_matrix, j, k) = gadd(gmael2(cup_matrix, j, k), gmael2(res, j, k));
                }
                for (int k = 1; k <= nr_col_full; ++k) {
                    gmael2(cup_matrix_full, j, k) = gadd(gmael2(cup_matrix_full, j, k), gmael2(res_full, j, k));
                }
            }
        }
        
    }
    mt_queue_end(&pt); /* end parallelism */


    printf(ANSI_COLOR_YELLOW "\n\nFull Cup Matrix:  (for (i,k) with i>k we have - (x_i cup x_k) instead of x_i cup x_k in order to see the symmetry).\n\n" ANSI_COLOR_RESET);
    for (int j = 1; j < nr_row + 1; ++j) {
        pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix_full, j));
    }

    int sym_check = 1;
    GEN fail;
    for (i = 1; i <= p_rk; i++) {
        for (int j = 1; j < r_rk; j++) {
            for (int k = 1; k <= p_rk; k++) {
                if (!gequal(gmael2(cup_matrix_full, j, p_rk * (k - 1) + i), gmael2(cup_matrix_full, j, p_rk * (i - 1) + k))) {
                    fail = mkvec2s(i, k);
                    sym_check = 0;
                    break;
                }
            }
        }
    }
    switch (sym_check) {
    case 1:
        printf(ANSI_COLOR_GREEN "\nSymmetry test passed\n\n" ANSI_COLOR_RESET);
        break;
    case 0:
        pari_printf(ANSI_COLOR_RED "\nSymmetry test failed at: %Ps\n\n" ANSI_COLOR_RESET, fail);
        pari_close();
        exit(111);
        break;
    default:
        break;
    }
    
    printf("-----------------------------------------------------------------------------------------------\n\n");
    printf(ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (int j = 1; j <= nr_row; ++j) {
        pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix, j));
    }
    printf(ANSI_COLOR_RED "For indices (i,i) we have B(x_i) instead of x_i cup x_i to get the correct presentation of Q_2.\n\n" ANSI_COLOR_RESET);
    printf(ANSI_COLOR_CYAN "This determines the second quotient Q_2 for the lower p-central series and not only the Zassenhaus quotient ZQ_2.\n\n" ANSI_COLOR_RESET);
    printf("———————————————————————————————–\n\n");
    printf(ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(cup_matrix)), p);
    pari_printf(ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    if (mat_rk > 0) {
        GEN cup_hnf = FpM_red(hnf((ZM_copy(cup_matrix))), p);
        printf(ANSI_COLOR_YELLOW "Hermite normal form:  \n\n" ANSI_COLOR_RESET);
        for (int i = 1; i < glength(cup_hnf) + 1; ++i) {
            pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_hnf, i));
        }
        
        printf("\n\n");
        
        char letters[] = "abcdefghijklmnopqr";
        printf(ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
        int hnf_r_rk = glength(cup_hnf);
        printf("\"\"\"[");
        for (int j = 1; j <= hnf_r_rk; j++) {
            for (int i = 1; i <= p_rk; ++i) {
                for (int k = i; k <= p_rk; k++) {
                    if (!gequal0(gel(gel(cup_hnf, j), (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1)))) {
                        if (i == k) {
                            printf("%c^%ld", letters[i - 1], itos(p) * smodis(gneg(gmael2(cup_hnf, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1))), itos(p)));
                        }
                        if (i < k) {
                            printf("%c_%c^%ld", letters[i - 1], letters[k - 1], smodis(gneg(gmael2(cup_hnf, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1))), itos(p)));
                        }
                    }
                }
            }
            if (j == hnf_r_rk) {
                printf("]\"\"\"\n");
            } else {
                printf(", ");
            }
        }
    }
    printf("\n--------------------------\nEnd: my_relations_par\n--------------------------\n\n");
    printf("-----------------------------------------------------------------------------------------------\n\n");

    return mat_rk;
}