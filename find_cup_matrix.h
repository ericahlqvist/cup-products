GEN my_cup_matrix (GEN LxAbs, GEN LxRel, GEN LyAbs, GEN LyRel, GEN K, GEN sigma_x, GEN sigma_y, GEN p, GEN J_vect, GEN Ix_vect, GEN Iy_vect, int p_int)
{
    GEN cup_matrix = mkvec3(zerovec(2),zerovec(2),zerovec(2));
    GEN NIpJ_x, NIpJ_y, Ix_rel, Iy_rel;
    int i;
    for (i=1; i<glength(J_vect)+1; i++) {
        Ix_rel = rnfidealabstorel(LxRel, gel(Ix_vect, i));
        Iy_rel = rnfidealabstorel(LyRel, gel(Iy_vect, i));
        NIpJ_x = idealmul(K, rnfidealnormrel(LxRel, Ix_rel), gel(J_vect, i));
        NIpJ_y = idealmul(K, rnfidealnormrel(LyRel, Iy_rel), gel(J_vect, i));

        gmael2(cup_matrix, 1,i) = stoi(smodis(my_Artin_symbol(LyAbs, LyRel, K, NIpJ_x, p_int, sigma_y), p_int));

        gmael2(cup_matrix, 2,i) = stoi(smodis(my_Artin_symbol(LxAbs, LxRel, K, NIpJ_x, p_int, sigma_x), p_int));

        gmael2(cup_matrix, 3,i) = stoi(smodis(my_Artin_symbol(LyAbs, LyRel, K, NIpJ_y, p_int, sigma_y), p_int));
    }
    return cup_matrix;
}    