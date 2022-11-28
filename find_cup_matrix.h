GEN my_cup_matrix (GEN LxAbs, GEN LxRel, GEN LyAbs, GEN LyRel, GEN K, GEN sigma_x, GEN sigma_y, GEN p, GEN J_vect, GEN I_vect, GEN T_x, int p_int)
{
    GEN cup_matrix = mkvec2(zerovec(2),zerovec(2));
    GEN NIpJ, I_rel;
    int i;
    for (i=1; i<glength(J_vect)+1; i++) {
        I_rel = rnfidealabstorel(LxRel, gel(I_vect, i));
        NIpJ = idealmul(K, rnfidealnormrel(LxRel, I_rel), gel(J_vect, i));

        gmael2(cup_matrix, 1,i) = stoi(smodis(my_Artin_symbol(LyAbs, LyRel, K, NIpJ, p_int, sigma_y), p_int));

        gmael2(cup_matrix, 2,i) = stoi(smodis(my_Artin_symbol(LxAbs, LxRel, K, NIpJ, p_int, sigma_x), p_int));
    }
    return cup_matrix;
}    