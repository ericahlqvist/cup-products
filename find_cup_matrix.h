GEN my_cup_matrix (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN J_vect)
{
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;
    int nr_col = p_rk^2;
    GEN cup_matrix = cgetg(nr_col, t_VEC);
    
    int i, j, k;
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; i++) {
        
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);

        // I_vect corresp. to i:th extension
        I_vect = my_find_I_vect(Labs_cup, Lrel_cup, K, sigma_cup, J_vect);
        
        for (j=1; j<p_rk+1; ++j) {
            // evaluate cup on j:th ideal J 
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            NIpJ = idealmul(K, rnfidealnormrel(Lrel_cup, I_rel), gel(J_vect, j));

            for (k=1; k<p_rk+1; ++k) {
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                sigma = gel(gel(K_ext, k), 3);
                gmael2(cup_matrix, i*j,k) = stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int));
            }
        }
        
    }
    return cup_matrix;
}    