
GEN my_xvector (int l, int x) {
    GEN vec = zerovec(l);
    int i; 
    GEN a = stoi(x);
    for (i = 1; i < l+1; i++)
    {
        gel(vec, i) = a;
    }
    return vec;
    
}

GEN my_xcol (int l, int x) {
    GEN vec = zerocol(l);
    int i; 
    GEN a = stoi(x);
    for (i = 1; i < l+1; i++)
    {
        gel(vec, i) = a;
    }
    return vec;
    
}

GEN my_int_to_frac_vec (GEN v) {
    GEN vec = zerovec(glength(v));
    GEN frac = cgetg(3, t_FRAC);
    int i;
    for (i = 1; i < glength(v)+1; ++i) {
        gel(frac, 1) = gel(v,i);
        gel(frac, 2) = gen_1;
        gel(vec,i) = frac;
    }
    return vec;
}

GEN my_QC_add (GEN v1, GEN v2) {
    if (!(glength(v1)==glength(v2))) {
        printf(ANSI_COLOR_RED "ERROR in my_ZC_add: vectors has different lengths\n\n" ANSI_COLOR_RESET);
        pari_printf("v1: %Ps\n\nv2: %Ps\n\n", v1, v2);
        exit(111);
    }
    GEN sum = zerocol(glength(v1));
    int i;
    for (i = 1; i < glength(v1)+1; ++i) {
        gel(sum,i) = gadd(gel(v1,i), gel(v2,i));
    }
    return sum;
}

int my_QV_equal1 (GEN v) {
    int output = 1;
    int i;
    if (!gequal1(gel(v,1))) {
            output = 0;
        }
    for (i=2; i<glength(v)+1; ++i) {
        if (!gequal0(gel(v,i))) {
            output = 0;
        }
    }
    return output;
}

int my_QV_equal0 (GEN v) {
    int output = 1;
    int i;
    for (i=1; i<glength(v)+1; ++i) {
        if (!gequal0(gel(v,i))) {
            output = 0;
        }
    }
    return output;
}

int my_QV_equal0_mod_p (GEN v, int p) {
    int output = 1;
    int i;
    for (i=1; i<glength(v)+1; ++i) {
        if (!(itos(gel(v,i))%p == 0)) {
            output = 0;
        }
    }
    return output;
}

int my_QV_equal (GEN v1, GEN v2) {
    int output = 1;
    int i;
    for (i=1; i<glength(v1)+1; ++i) {
        if (!gequal(gel(v1,i), gel(v2,i))) {
            output = 0;
        }
    }
    return output;
}

int my_SQ_MAT_equal (GEN M1, GEN M2) {
    // printf(ANSI_COLOR_CYAN "my_SQ_MAT_equal\n\n" ANSI_COLOR_RESET);
    int output = 1;
    int i;
    int j;
    
    // outmat(M1);
    // outmat(M2);
    for (i = 1; i < glength(gel(M1, 1))+1; ++i) {

        for (j = 1; j < glength(gel(M1, 1))+1; j++)
        {

            if (!gequal(gel(gel(M1, i), j), gel(gel(M2, i), j))) {
                
                return 0;
            }
        }
        
    }
    return output;
}

int my_SQ_MAT_equal0 (GEN M) {
    // printf(ANSI_COLOR_CYAN "my_SQ_MAT_equal\n\n" ANSI_COLOR_RESET);
    int output = 1;
    int i;
    int j;
    
    for (i = 1; i < glength(gel(M, 1))+1; ++i) {
        for (j = 1; j < glength(gel(M, 1))+1; j++)
        {
            if (!gequal0(gel(gel(M, i), j))) {
                output = 0;
                break;
            }
        }
        
    }
    return output;
}

GEN my_norm (GEN Labs, GEN Lrel, GEN L, GEN elem, GEN sigma, int p) {
    int i;
    GEN new_elem = algtobasis(Labs, gen_1);
    for (i=0; i<p; ++i) {
        new_elem = basistoalg(Labs, nfmul(Labs, new_elem, elem));
        
        elem = galoisapply(Labs, sigma, elem);
        // pari_printf("new_elem: %Ps\n\nelem: %Ps\n\n", new_elem, elem);
    }
    
    GEN norm = algtobasis(L, rnfeltdown(Lrel, rnfeltabstorel(Lrel, new_elem)));
    
    return norm;
}

GEN my_norm_ideal (GEN Labs, GEN Lrel, GEN L, GEN I, GEN sigma, int p) {
    pari_sp av = avma;
    int i;
    GEN new_ideal = algtobasis(Labs, gen_1);
    for (i=0; i<p; ++i) {
        new_ideal = basistoalg(Labs, nfmul(Labs, new_ideal, I));
        
        I = galoisapply(Labs, sigma, I);
        // pari_printf("new_elem: %Ps\n\nelem: %Ps\n\n", new_elem, elem);
    }
    
    GEN norm = rnfidealdown(Lrel, rnfidealabstorel(Lrel, new_ideal));
    // printf("Hoooj\n\n");
    norm = gerepilecopy(av, norm);
    return norm;
}

GEN my_i (GEN Labs, GEN Lrel, GEN sigma, GEN elem) {
    GEN lift = rnfeltup0(Lrel, elem, 1);
    output(lift);
    GEN mod_lift = galoisapply(Labs, sigma, lift);
    output(mod_lift);
    printf("\n");
    return mod_lift;
}

GEN my_i_ideal (GEN Labs, GEN Lrel, GEN sigma, GEN ideal) {
    GEN lift = rnfidealup0(Lrel, ideal, 1);
    GEN mod_lift = galoisapply(Labs, sigma, lift);
    return mod_lift;
}

/*------------------------------------
 The function sigma_x-1 on ideals
------------------------------------
* L - a number field
* sigma - An auto nfalgtobasis(nf, c[2]) where c = nfgaloisconj(nf);
* I - an ideal of nf given as a matrix in hnf

* output - (\sigma_x-1)(I)

DEFINITION:
-----------------------*/

GEN my_SM1_ideal (GEN L, GEN sigma, GEN I)
{
    return idealmul(L, galoisapply(L, sigma, I), idealinv(L, I));
    // return idealred(L, idealmul(L, galoisapply(L, sigma, I), idealinv(L, I)));
} 

GEN my_1MS_ideal (GEN L, GEN sigma, GEN I) 
{
    return idealmul(L, I, idealinv(L, galoisapply(L, sigma, I)));
    // return idealred(L, idealmul(L, I, idealinv(L, galoisapply(L, sigma, I))));
} 

/*
SM1 for elements
*/
GEN my_SM1_elt (GEN L, GEN sigma, GEN elem) 
{
    return nfdiv(L, galoisapply(L, sigma, elem), elem);
}

GEN my_1MS_elt (GEN L, GEN sigma, GEN elem) 
{
    return nfdiv(L, elem, galoisapply(L, sigma, elem));
}

GEN my_Gamma (GEN Labs, GEN sigma, GEN l, int p_int) 
{
    pari_sp av = avma;
    int n;
    int i;

    GEN Gamma_l = algtobasis(Labs, gen_1);
    for (n = 1; n < p_int; ++n)
    {
        GEN elem_n = algtobasis(Labs, l);
        GEN n_gen = stoi(n);
        elem_n = nfpow(Labs, elem_n, n_gen);
        for (i = 1; i < n+1; ++i)
        {
            elem_n = galoisapply(Labs, sigma, elem_n);
        }     
        Gamma_l = nfmul(Labs, Gamma_l, elem_n);
    }
    Gamma_l = nfinv(Labs, Gamma_l);
    Gamma_l = gerepilecopy(av , Gamma_l);
    return Gamma_l;
}

GEN my_Gamma_ideal (GEN Labs, GEN sigma, GEN I, int p_int) 
{
    pari_sp av = avma;
    int n;
    int i;

    GEN Gamma_l = idealhnf0(Labs, gen_1, NULL);
    for (n = 1; n < p_int; ++n)
    {
        GEN elem_n = gcopy(I);
        GEN n_gen = stoi(n);
        elem_n = idealpow(Labs, elem_n, n_gen);
        for (i = 1; i < n+1; ++i)
        {
            elem_n = galoisapply(Labs, sigma, elem_n);
        }     
        Gamma_l = idealmul(Labs, Gamma_l, elem_n);
    }
    Gamma_l = idealinv(Labs, Gamma_l);
    Gamma_l = gerepilecopy(av, Gamma_l);
    return Gamma_l;
}

/* The function p-N*/
GEN my_pmN (GEN Labs, GEN Lrel, GEN z, GEN p)  
{
    GEN elem1 = nfpow(Labs, z, p);
    GEN elem2 = rnfeltup(Lrel, rnfeltnorm(Lrel, z));
    return nfdiv(Labs, elem1, elem2);
}

GEN my_tau (GEN Labs, GEN a, GEN b, GEN sigma)
{
    return nfmul(Labs, a, galoisapply(Labs, sigma, b));
}

GEN my_trace_tau (GEN Labs, GEN a, GEN b, GEN sigma, int p)
{
    pari_sp av = avma;
    GEN new_b = algtobasis(Labs, b);
    GEN trace = my_int_to_frac_vec(algtobasis(Labs, gen_0));
    // pari_printf("Trace: %Ps\n\n", trace);
    trace = my_QC_add(trace, new_b);
    // pari_printf("Trace: %Ps\n\n", trace);
    
    int i;
    for (i = 1; i<p; ++i) {
        new_b = my_tau(Labs, a, new_b, sigma);
        // pari_printf("new_b: %Ps\n\n", new_b);
        trace = my_QC_add(trace, new_b);
        // pari_printf("Trace: %Ps\n\n", trace);
    }
    trace = gerepilecopy(av, trace);
    return trace;
}


GEN my_find_primes_in_factorization(GEN LyAbs, GEN factorization) {
    pari_sp av = avma;
    int l = glength(gel(factorization, 1));
    GEN primes = zerovec(l);
    GEN es = zerovec(l);
    int i;
    for (i = 1; i < l+1; i++)
    {
        gel(primes,i) = idealhnf0(LyAbs, gel(gel(gel(factorization, 1), i), 1), gel(gel(gel(factorization, 1), i), 2));
        gel(es, i) = gel(gel(factorization, 2), i);
    }
    GEN primes_and_es = mkvec2(primes, es);
    primes_and_es = gerepilecopy(av, primes_and_es);
    return primes_and_es;
}

GEN my_1MS_operator (GEN L, GEN sigma) {
    printf("my_1MS_operator\n");
    pari_sp av = avma;
    GEN cyc = bnf_get_cyc(L), col; 
    GEN gens = bnf_get_gen(L);
    int l = glength(cyc);
    GEN M = zeromatcopy(l,l);
    int i;
    int j;
    for (i = 1; i < l+1; i++)
    {
        col = bnfisprincipal0(L, my_1MS_ideal(L, sigma, gel(gens, i)), 0);
        for (j = 1; j < l+1; j++)
        {
            gel(col, j) = modii(gel(col, j), gel(cyc, j));
        }
        gel(M, i) = col;
    }
    //pari_printf("my_1MS_operator: %Ps\n", M);
    M = gerepilecopy(av, M);
    return M;
}

GEN my_1MS_elt_operator (GEN L, GEN sigma) {
    printf("my_1MS_elt_operator\n");
    pari_sp av = avma;
    GEN int_basis = nf_get_zk(bnf_get_nf(L)), col, basis_vec;
    //pari_printf("int_basis: %Ps\n", int_basis);
    int l = glength(int_basis);
    GEN M = zeromatcopy(l,l);
    int i;
    
    for (i = 1; i < l+1; i++)
    {
        basis_vec = algtobasis(L, gel(int_basis, i));
        col = algtobasis(L, my_1MS_elt(L, sigma, basis_vec));
        //pari_printf("col: %Ps\n", col);
        gel(M, i) = col;
    }
    M = gerepilecopy(av, M);
    return M;
}


GEN my_find_ext_ranks(GEN K_ext) {
    int l = glength(K_ext);
    GEN ext_rks = zerovec(l);

    int i;

    for (i = 1; i < l+1; i++)
    {
        gel(ext_rks, i) = stoi(glength(bnf_get_cyc(gel(gel(K_ext, i), 1))));
    }
    
    return vecsort0(ext_rks, NULL, 0);
}

GEN my_find_ext_cyc(GEN K_ext) {
    int l = glength(K_ext);
    GEN ext_cyc = zerovec(l);

    int i;

    for (i = 1; i < l+1; i++)
    {
        gel(ext_cyc, i) = bnf_get_cyc(gel(gel(K_ext, i), 1));
    }
    
    return ext_cyc;
}

int my_is_p_power (GEN vect, int p_int, int order) {
    int l = glength(vect), check_p = 1, i;

    if (l==0) {
        return 0;
    }
    
    if (itos(gel(vect, l))!=(p_int%order) && itos(gel(vect, l))%p_int!=0) {
            return 0;
        }
    for (i=1;i<l;++i) {
        if (itos(gel(vect, i))%p_int!=0) {
            check_p = 0;
            break;
        }
    } 

    return check_p;
}

GEN my_find_ext_pol(GEN K_ext) {
    int l = glength(K_ext);
    GEN ext_pol = zerovec(l);

    int i;

    for (i = 1; i < l+1; i++)
    {
        gel(ext_pol, i) = nf_get_pol(bnf_get_nf(gel(gel(K_ext, i), 1)));
    }
    
    return ext_pol;
}

GEN my_Gal_rel (GEN Labs, GEN Lrel, GEN K, GEN sigma, int p) {
    
    GEN Gal_rel = zerovec(p);
    GEN current_sigma = sigma;
    int i;
    
    for (i = 1; i < p+1; i++)
    {
        gel(Gal_rel, i) = current_sigma;
        current_sigma = galoisapply(Labs, sigma, current_sigma); 
        
    }

    return Gal_rel;
}

// Creates a set (vector) of of exponents in bijection with Cl(L) 
GEN my_get_vect (int n, GEN cyc)
{
    int b = itos(gel(cyc, n+1));
    GEN next_vect;
    
    if (n > 0) {
        GEN prev_vect = my_get_vect(n-1, cyc);
        
        int l = glength(prev_vect);
        next_vect = zerovec(b*l);
        //pari_printf("%Ps\n", next_vect);
        
        //printf("%ld\n", glength(next_vect));
        int i;
        for (i = 0; i < b*l; ++i) {
            double num = i+1;
            double den = b;
            int first_index = ceil(num/den);
            //printf("%d\n",first_index);
            gel(next_vect, i+1) = gconcat(gel(prev_vect, first_index), mkvec(stoi((i+1)%b)));
            
        }

    }
    else {
        
        next_vect = zerovec(b);

        int i;
        for (i = 0; i<b; ++i) {
            gel(next_vect, i+1) = stoi(i);
        }
    }
    return next_vect;
}

GEN my_get_clgp (GEN K)
{
    pari_sp av = avma;
    printf("-------\n\nComputing the class group\n\n---------\n");
    GEN Kcyc = bnf_get_cyc(K);
    GEN Kgen = bnf_get_gen(K);
    GEN class_number = bnf_get_no(K);
    int clnr = itos(class_number);
    int nr_comp = glength(Kcyc);
    GEN class_group_exp;
    int n;
    if (nr_comp > 1) {
        class_group_exp = my_get_vect( nr_comp - 1, Kcyc );
    }
    else {
        class_group_exp = zerovec(clnr);
        for (n=0; n<clnr; ++n) {
            gel(class_group_exp, n+1) = mkvec(stoi(n));
        }
    }
    // pari_printf("cyc: %Ps\n\n", Kcyc);
    // pari_printf("cl_gp_exp: %Ps\n\n", class_group_exp);
    GEN class_group = zerovec(clnr);
    GEN current_I, exponents, pow;
    
    
    for (n = 1; n < clnr + 1; ++n) {
        exponents = gel(class_group_exp, n);
        
        int i;
        current_I = idealhnf0(K, gen_1, NULL);
        for ( i = 1; i < nr_comp + 1; ++i ) {
            
            pow = idealpow(K, gel(Kgen, i), gel(exponents, i));
            
            current_I = idealmul(K, current_I, pow);
            // printf("ideal norm: ");
            // output(idealnorm(K, current_I));
            // printf("\n");
        }
        if (n%1000 == 0) {
            printf("%d/%d\n", n, clnr);
        }
        
        gel(class_group, n) = idealred0(K, current_I, NULL); 
    }
    
    class_group = gerepilecopy(av, class_group);
    return class_group;
}

GEN my_ideal_from_exp (GEN K, GEN exp) {
    //printf("\nmy_ideal_from_exp\n");
    pari_printf("exp: %Ps\n", exp);
    pari_sp av = avma;
    GEN ideal = idealhnf0(K, gen_1, NULL);
    GEN gens = bnf_get_gen(K);
    int l = glength(exp);
    int i;
    for (i = 1; i < l+1; i++)
    {
        printf("%d/%d\n", i, l);
        ideal = idealred0(K, idealmul(K, ideal, idealred0(K, idealpow(K, gel(gens, i), gel(exp, i)), NULL)), NULL);
        //printf("%d/%d\n", i, l);
    }
    

    ideal = gerepilecopy(av, ideal);
    return ideal;

}

GEN my_vect_from_exp (GEN basis, GEN exp) {
    //printf("\nmy_vect_from_exp\n");
    pari_sp av = avma;
    int n = glength(gel(basis, 1));
    GEN vect = zerocol(n);
    //pari_printf("exp: %Ps\n", exp);
    int l = glength(exp);

    int i;
    for (i = 1; i < l+1; i++)
    {
        
        vect = ZC_add(vect, ZC_Z_mul(gel(basis, i), gel(exp, i)));
    }
    
    vect = gerepilecopy(av, vect);
    return vect;

}

GEN my_get_unit_group (GEN K, GEN unit_gens, GEN p)
{
    pari_sp av = avma;
    printf("-------\n\nComputing the unit group\n\n---------\n");
    int nr_comp = glength(unit_gens);
    GEN cyc = zerovec(nr_comp);

    int i;
    
    for (i = 1; i < nr_comp+1; i++)
    {
        gel(cyc, i) = p;
    }
    
    int nr = pow(itos(p), nr_comp);
    
    GEN group_exp;
    int n;
    if (nr_comp > 1) {
        group_exp = my_get_vect( nr_comp - 1, cyc );
    }
    else {
        group_exp = zerovec(nr);
        for (n=0; n<nr; ++n) {
            gel(group_exp, n+1) = mkvec(stoi(n));
        }
    }
    
    // pari_printf("cyc: %Ps\n\n", Kcyc);
    // pari_printf("cl_gp_exp: %Ps\n\n", class_group_exp);
    GEN group = zerovec(nr);
    GEN current_u, exponents, pow;
    
    
    for (n = 1; n < nr + 1; ++n) {
        exponents = gel(group_exp, n);
        
        int i;
        current_u = gen_1;
        
        for ( i = 1; i < nr_comp + 1; ++i ) {
            
            pow = nfpow(K, gel(unit_gens, i), gel(exponents, i));
            
            current_u = nfmul(K, current_u, pow);
            // printf("ideal norm: ");
            // output(idealnorm(K, current_I));
            // printf("\n");
        }
        
        if (n%1000 == 0) {
            printf("%d/%d\n", n, nr);
        }
        
        gel(group, n) = algtobasis(K, current_u); 
    }
    
    group = gerepilecopy(av, group);
    return group;
}



//------------------------------
// Returns an ideal I in Div(L) such that iJ = (1-sigma)I in Cl(L)
//------------------------------ 
GEN my_H90 (GEN L, GEN iJ, GEN sigma) {
    printf("\nmy_H90\n");
    pari_sp av = avma;
    GEN H90_ideal, M, B, E, D;


    // Cycle type of Cl(L) as a column vector
    D = gtocol(bnf_get_cyc(L));
    
    // l x l zero matrix
    // M = zeromatcopy(l,l);
    // finding the matrix M consisting of exponents for the (1-sigma)g_i, Cl(L) = < g_1, ..., g_n > 
    // that is, finding the matrix corresponding to the linear operator (1-sigma)
    M = my_1MS_operator(L, sigma);
    //D = zerocol(l);

    // finding exponents for iJ when written as a products of ideals in gens
    B = bnfisprincipal0(L, iJ, 0);

    
    // int i;
    // for (i = 1; i < l+1; i++)
    // {
    //     gel(M, i) = bnfisprincipal0(L, idealred0(L, my_1MS_ideal(L, sigma, gel(gens, i)), NULL), 0);
    // }

    // Gauss: Solving the system M*X = B (i.e., finding I s.t. (1-sigma)I = iJ)
    E = matsolvemod(M,D,B,0);
    
    // Defining the ideal I from the exponents E
    H90_ideal = my_ideal_from_exp(L, E);
    printf(ANSI_COLOR_GREEN "H90_ideal found\n" ANSI_COLOR_RESET);

    // Test that nothing went wrong
    GEN Itest = idealred(L,idealdiv(L, iJ, my_1MS_ideal(L, sigma, H90_ideal)));
    GEN test_vec = bnfisprincipal0(L, Itest, 0);
    
    if (!my_QV_equal0(test_vec)) {
        printf(ANSI_COLOR_RED "Problem in my_H90\n" ANSI_COLOR_RESET);
        pari_close();
        exit(0);
    }

    H90_ideal = gerepilecopy(av, H90_ideal);
    return H90_ideal;
}

GEN my_H90_exp (GEN L, GEN iJ, GEN sigma) {
    printf("\nmy_H90_exp\n");
    pari_sp av = avma;
    GEN H90_exp, gens, M, B, D;
    gens = bnf_get_gen(L);
    D = gtocol(bnf_get_cyc(L));
    int l = glength(gens);
    
    M = zeromatcopy(l,l);
    //D = zerocol(l);
    // finding exponents for iJ 
    B = bnfisprincipal0(L, iJ, 0);

    // finding the matrix M consisting of exponents for the (1-s)g_i, Cl(L) = < g_1, ..., g_n > 
    int i;
    for (i = 1; i < l+1; i++)
    {
        gel(M, i) = bnfisprincipal0(L, idealred0(L, my_1MS_ideal(L, sigma, gel(gens, i)), NULL), 0);
    }
    // Gauss
    // pari_printf("M: %Ps\n", M);
    // pari_printf("B: %Ps\n", B);
    //printf("Gauss start\n");
    H90_exp = matsolvemod(M,D,B,0);
    //printf("Gauss done\n");

    H90_exp = gerepilecopy(av, H90_exp);
    return H90_exp;
}


GEN my_find_units_mod_p (GEN K, GEN p) {
    
    GEN fund_units = bnf_get_fu(K);
    GEN tors_unit = bnf_get_tuU(K);
    int tors_gp_ord = bnf_get_tuN(K);

    GEN units_mod_p = fund_units;
    // pari_printf("A: %Ps\n", gpow(p, stoi(4), DEFAULTPREC));
    // pari_printf("B: %Ps\n", nfpow(K, tors_unit, gpow(p, stoi(4), DEFAULTPREC)));
    // pari_printf("tors_unit: %Ps\n", tors_unit);
    // pari_printf("fund_units: %Ps\n", fund_units);
    // printf("torsion order: %d\n\n", tors_gp_ord);

    if (tors_gp_ord%itos(p)==0) {
        units_mod_p = shallowconcat(units_mod_p, mkvec(tors_unit));
    }
    
    // printf("Units:\n\n");
    // for (i=1;i<glength(units_mod_p)+1;++i) {
    //     pari_printf("u %d: %Ps\n", i, algtobasis(K, gel(units_mod_p, i)));
    //     pari_printf("u^-1 %d: %Ps\n", i, nfinv(K, algtobasis(K, gel(units_mod_p, i))));
    // }
    // printf("\n");

    // ad-hoc modification
    // gel(units_mod_p, 2) = nfmul(K, gel(units_mod_p, 1), gel(units_mod_p, 4));
    // gel(units_mod_p, 3) = nfdiv(K, gel(units_mod_p, 1), gel(units_mod_p, 4));

    return units_mod_p;
}

GEN my_norm_operator (GEN Labs, GEN Lrel, GEN K, GEN p) {
    pari_sp av = avma;
    GEN M, g, g_rel, Ng, Lgens;
    int Klength, Llength, i;
    Lgens = shallowconcat(bnf_get_fu(Labs), bnf_get_tuU(Labs));
    Llength = glength(bnf_get_fu(Labs))+1;
    Klength = glength(bnf_get_fu(K));
    M = zeromatcopy(Klength, Llength);

    for (i = 1; i < Llength+1; i++)
    {
        g = gel(Lgens, i);
        g_rel = rnfeltabstorel(Lrel, g);
        Ng = rnfeltnorm(Lrel, g_rel);
        gel(M, i) = bnfisunit0(K, Ng, NULL);
    }
    

    M = gerepilecopy(av, M);
    return M;
}

GEN my_find_p_gens (GEN K, GEN p)
{
    pari_sp av = avma;
    GEN my_n = bnf_get_cyc(K);
    int l = glength(my_n);
    
    GEN p_gens;
    p_gens = zerovec(l);
    GEN current_gen;
    
    GEN my_gens = bnf_get_gen(K);
    //pari_printf("my_gens: %Ps\n\n", my_gens);
    int i;
    for (i = 1; i < l+1; ++i)
    {
        //pari_printf("my_n, my_n/p, gen^my_n: %Ps, %Ps, %Ps\n\n", gel(my_n, i), gdiv(gel(my_n, i),p), gel(bnfisprincipal0(K, idealpow(K, gel(my_gens, i), gel(my_n, i)), 1), 1));
        current_gen = idealpow(K,gel(my_gens, i), gdiv(gel(my_n, i),p));
        //pari_printf("current_gen: %Ps\n\n", current_gen);
        gel(p_gens, i) = idealred0(K, current_gen, NULL);
    }
    p_gens = gerepilecopy(av, p_gens);
    return p_gens;
}


void my_unramified_p_extensions(GEN K, GEN p, GEN D_prime_vect) {
    int i;
    GEN s = pol_x(fetch_user_var("s"));
    GEN x = pol_x(fetch_user_var("x"));
    GEN index = mkvec(p);
    // GEN R = nfsubfields0(clf_pol,4,1);
    // pari_printf("subgrouplist: %Ps\n", subgrouplist0(bnf_get_cyc(K), mkvec(gpow(p,gen_2, DEFAULTPREC)), 0));

    GEN R = bnrclassfield(K, subgrouplist0(bnf_get_cyc(K), index, 2), 2, DEFAULTPREC);
    //GEN Rsq = bnrclassfield(K, subgrouplist0(bnf_get_cyc(K), mkvec(gpow(p,gen_2, DEFAULTPREC)), 2), 2, DEFAULTPREC);
    //GEN Rcb = bnrclassfield(K, subgrouplist0(bnf_get_cyc(K), mkvec(gpow(p,stoi(3), DEFAULTPREC)), 2), 2, DEFAULTPREC);
    GEN abs_pol;
    printf("[p]-extensions:\n\n");
    for (i=1;i<glength(R)+1;i++) {
        //printf("i: %d\n", i);
        abs_pol = polredabs0(mkvec2(gel(R, i), D_prime_vect), 0);
        //abs_pol = polredabs(gel(R, i));
        pari_printf("%Ps, cyc: %Ps\n", gsubstpol(abs_pol, x, s), bnf_get_cyc(Buchall(abs_pol, nf_FORCE, DEFAULTPREC)));
        
    }
    printf("\n");
    // printf("[p,p]-extensions:\n\n");
    // for (i=1;i<glength(Rsq)+1;i++) {
    //     abs_pol = polredabs0(mkvec2(gel(Rsq, i), D_prime_vect), 0);
    //     //abs_pol = polredabs(gel(R, i));
    //     pari_printf("%Ps, cyc: %Ps\n", gsubstpol(abs_pol, x, s), bnf_get_cyc(Buchall(abs_pol, nf_FORCE, DEFAULTPREC)));
    // }
    // printf("\n");
    // printf("[p,p,p]-extensions:\n\n");
    // for (i=1;i<glength(Rcb)+1;i++) {
    //     abs_pol = polredabs0(mkvec2(gel(Rcb, i), D_prime_vect), 0);
    //     //abs_pol = polredabs(gel(R, i));
    //     pari_printf("%Ps, cyc: %Ps\n", gsubstpol(abs_pol, x, s), bnf_get_cyc(Buchall(abs_pol, nf_FORCE, DEFAULTPREC)));
    // }
    // printf("\n");
}


// Returns vector of tuples (a,J) with div(a)+pJ = 0.
GEN my_find_Ja_vect(GEN K, GEN J_vect, GEN p, GEN units_mod_p) {
    pari_sp av = avma;
    int l = glength(J_vect);
    int r_rk = l + glength(units_mod_p);
    GEN Ja_vect = zerovec(r_rk);
    GEN a;
    int i;
    
    for (i=1; i<l+1; ++i) {
        a = nfinv(K, gel(bnfisprincipal0(K, idealpow(K, gel(J_vect, i), p), 1), 2));
        gel(Ja_vect, i) = mkvec2(a, gel(J_vect, i));
    }
    
    for (i = l+1; i < r_rk+1; i++)
    {
        gel(Ja_vect, i) = mkvec2(gel(units_mod_p, i-l), idealhnf0(K, gen_1, NULL));
    }
    
    Ja_vect = gerepilecopy(av, Ja_vect);
    return Ja_vect;
}



GEN my_get_sums (GEN basis) {
    printf("my_get_sums\n");
    pari_sp av = avma;
    int l = glength(basis), i, ord = 2;
    int n = pow(ord, l);
    GEN gp = zerovec(n);
    GEN cyc = zerovec(l);
    
    if (l<2)
    {
        gp = shallowconcat(mkvec(zerocol(glength(gel(basis, 1)))), basis);
        //pari_printf(ANSI_COLOR_MAGENTA "gp: %Ps\n" ANSI_COLOR_RESET, gp);
        gp = gerepilecopy(av, gp);
        return gp;
    }
    

    else {
        for (i = 1; i < l+1; i++)
        {
            gel(cyc, i) = stoi(ord);
        }
        
        GEN exp = my_get_vect( l - 1, cyc );
        //pari_printf(ANSI_COLOR_MAGENTA "exp_list: %Ps\n" ANSI_COLOR_RESET, exp);
        for (i = 1; i < n+1; i++)
        {
            gel(gp, i) = my_vect_from_exp(basis, gel(exp, i));
        }
    }
    // pari_printf(ANSI_COLOR_RED "Basis: %Ps\n" ANSI_COLOR_RESET, basis);
    // pari_printf(ANSI_COLOR_RED "Sums: %Ps\n" ANSI_COLOR_RESET, gp);
    gp = gerepilecopy(av, gp);
    return gp;
}

//-------------------------------------------
// The function my_H90_vect finds for each (a, J) in Ja_vect, a fractional ideal I in Div(L) such that 
// i(J) = (1-sigma)I + div(t), where t in L^x satisfies N(t)*a = 1.  
// Returns: a vector of these I's 
//-------------------------------------------
GEN my_H90_vect (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN Ja_vect, GEN p) {
    printf(ANSI_COLOR_CYAN "\nmy_H90_vect\n" ANSI_COLOR_RESET);
    pari_sp av = avma;
    int r_rk = glength(Ja_vect), f, j, done = 0;
    GEN I_vect = zerovec(r_rk), a, iJ, F, ker_T, ker_T_basis, ker_sol, F_ker_T, t, t_rel, Nt, diff, exp, is_princ, is_norm, cyc;
    cyc = gtocol(bnf_get_cyc(Labs));
    int i;
    
    for (i = 1; i < r_rk+1; ++i)
    {
        done = 0;
        a = gel(gel(Ja_vect, i), 1);
        iJ = rnfidealup0(Lrel, gel(gel(Ja_vect, i),2), 1);

        // Find I s.t. (1-sigma)I = iJ in Cl(L)
        F = my_H90(Labs, iJ, sigma);
        // Then (1-sigma)I+div(t) = iJ for some t in L^x. However, it might be the case that N(t)*a is not 1. 
        // We know from theory that there should be an I and a t s.t (1-sigma)I+div(t) = iJ and N(t)*a=1. 
        // If I' and t' are such, then I-I' is in ker (1-sigma) so in order to find the correct I' we need to go through the kernel of (1-sigma) and for each I'' in there, take I+I'' and check f the corresponding t' satisfies N(t')*a = 1.   

        // Find generators for the kernel of (1-sigma)
        ker_sol = matsolvemod(my_1MS_operator(Labs, sigma), cyc, zerocol(glength(cyc)), 1);
        ker_T_basis = gtovec(gel(ker_sol, 2));
        // printf(ANSI_COLOR_GREEN "------------------------\n\n\nker_T_basis size: %ld\n\n\n------------------------\n" ANSI_COLOR_RESET, glength(ker_T_basis));
        
        // If the kernel is trivial, then we are done. 
        if (glength(ker_T_basis)==0)
        {
            gel(I_vect, i) = F;
            done = 1;
        }
        else {
            // Generate part of the kernel (or more precisely, the corresponding exponents)
            ker_T = my_get_sums(ker_T_basis);

            // WARNING: The next line makes it faster but might not find the answer. Switch back to the above ker_T!
            //ker_T = shallowconcat(mkvec(zerocol(glength(gel(ker_T_basis, 1)))), ker_T_basis);
            
            //pari_printf(ANSI_COLOR_CYAN "ker_T: %Ps\n" ANSI_COLOR_RESET, ker_T);
            // pari_printf("ker_T[2]: %Ps\n", gel(ker_T, 2));
            f = glength(ker_T);
            printf("Searching a chunk of ker (1-sigma) of size: %d\n", f);
            /*
                Find F_ker_T --------------------
            */
            
            for (j = 1; j < f+1; j++)
            {
                printf("\nSearching: %d/%d\n", j, f);
                // pari_printf("Adding the exp: %Ps\n", gel(ker_T, j));
                // pari_printf("And the ideal: %Ps\n", my_ideal_from_exp(Labs, gel(ker_T, j)));

                // This is our I+I'' as explained above
                F_ker_T = idealred0(Labs, idealmul(Labs, F, idealred0(Labs, my_ideal_from_exp(Labs, gel(ker_T, j)), NULL)), NULL);
                //pari_printf("F_ker_T: %Ps\n", F_ker_T);

                // Now find the corresponding t
                is_princ = bnfisprincipal0(Labs, idealdiv(Labs, iJ, my_1MS_ideal(Labs, sigma, F_ker_T)), nf_GEN);

                // Sanity check
                if (!my_QV_equal0(gel(is_princ, 1)))
                {
                    printf(ANSI_COLOR_RED "Problem in my_H90_vect\n" ANSI_COLOR_RESET);
                    pari_close();
                    exit(111);
                }
                
                // The corresponding t
                t = gel(is_princ, 2);
                //pari_printf("t: %Ps\n", t);
                //printf("length t: %ld\n", glength(t));

                // In case precision is too low, t may not be given, but instead a "place holder" []^~.
                if (glength(t)==0)
                {
                    return stoi(-1);
                }
                
                // In relative coordinates
                t_rel = rnfeltabstorel(Lrel, t);

                // N(t)
                Nt = rnfeltnorm(Lrel, t_rel);

                // N(t)*a
                diff = nfmul(K, Nt, a);

                // Since div(a)+pJ = 0 and div(N(t))-pJ = 0, we get that div(N(t)*a) = 0 and therefore N(t)*a must be a unit.
                // Next we find its exponents in terms of the fixed generators of the unit group 
                exp = bnfisunit0(K, diff, NULL);
                pari_printf(ANSI_COLOR_MAGENTA "exp: %Ps\n" ANSI_COLOR_RESET, exp);

                // Check if N(t)*a is the norm of a unit. If it is, we may modify t by this unit without effecting the equality 
                // (1-sigma)I=iJ and hence we are done. Returns zero if not a norm (which should never happen).
                is_norm = matsolvemod(my_norm_operator(Labs, Lrel, K, p), zerocol(glength(exp)), gtocol(exp), 0);
                pari_printf(ANSI_COLOR_CYAN "is_norm: %Ps\n" ANSI_COLOR_RESET, is_norm);

                // if N(t)*a = 1 or if N(t)*a is a norm
                if (my_QV_equal0(exp) || !gequal0(is_norm))
                {
                    //pari_printf(ANSI_COLOR_CYAN "Norm of elt with exp: %Ps\n" ANSI_COLOR_RESET, is_norm);
                    // w = gel(bnf_get_fu(Labs), 1);
                    // w_rel = rnfeltabstorel(Lrel, w);
                    // Nw = rnfeltnorm(Lrel, w_rel);
                    // pari_printf("Norm w: %Ps\n", Nw);
                    gel(I_vect, i) = F_ker_T;
                    printf(ANSI_COLOR_GREEN "\nI found\n\n" ANSI_COLOR_RESET);
                    done = 1;
                    break;
                }
            }
        }
        // If some of the I's were never found, then we return -1 and another (slower) function will take over.  
        if (!done)
        {
            printf(ANSI_COLOR_RED "my_H90_vect ended with a problem \n" ANSI_COLOR_RESET);
            return stoi(-1);
            // pari_close();
            // exit(111);
        }
    } 
    
    I_vect = gerepilecopy(av, I_vect);
    return I_vect;
}

GEN my_H90_12 (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN a, GEN J, GEN p) {
    printf(ANSI_COLOR_CYAN "\nmy_H90_vect\n" ANSI_COLOR_RESET);
    pari_sp av = avma;
    int f, j, done = 0;
    GEN I, iJ, F, ker_T, ker_T_basis, ker_sol, F_ker_T, t, t_rel, Nt, diff, exp, is_princ, is_norm, cyc;
    cyc = gtocol(bnf_get_cyc(Labs));
    
    

    done = 0;
    iJ = rnfidealup0(Lrel, J, 1);
    F = my_H90(Labs, iJ, sigma);
    // pari_printf("F: %Ps\n", F);
    ker_sol = matsolvemod(my_1MS_operator(Labs, sigma), cyc, zerocol(glength(cyc)), 1);
    ker_T_basis = gtovec(gel(ker_sol, 2));
    printf(ANSI_COLOR_GREEN "------------------------\n\n\nker_T_basis size: %ld\n\n\n------------------------\n" ANSI_COLOR_RESET, glength(ker_T_basis));
    
    if (glength(ker_T_basis)==0)
    {
        I = F;
        done = 1;
    }
    else {
        ker_T = my_get_sums(ker_T_basis);

        // WARNING: This makes it faster but might not find the answer. Switch back to the above ker_T
        //ker_T = shallowconcat(mkvec(zerocol(glength(gel(ker_T_basis, 1)))), ker_T_basis);
        
        //pari_printf(ANSI_COLOR_CYAN "ker_T: %Ps\n" ANSI_COLOR_RESET, ker_T);
        // pari_printf("ker_T[2]: %Ps\n", gel(ker_T, 2));
        f = glength(ker_T);
        printf("Searching a chunk of ker (1-sigma) of size: %d\n", f);
        /*
            Find F_ker_T --------------------
        */
        
        for (j = 1; j < f+1; j++)
        {
            printf("\nSearching: %d/%d\n", j, f);
            // pari_printf("Adding the exp: %Ps\n", gel(ker_T, j));
            // pari_printf("And the ideal: %Ps\n", my_ideal_from_exp(Labs, gel(ker_T, j)));
            F_ker_T = idealred0(Labs, idealmul(Labs, F, idealred0(Labs, my_ideal_from_exp(Labs, gel(ker_T, j)), NULL)), NULL);
            //pari_printf("F_ker_T: %Ps\n", F_ker_T);
            is_princ = bnfisprincipal0(Labs, idealdiv(Labs, iJ, my_1MS_ideal(Labs, sigma, F_ker_T)), 1);
            if (!my_QV_equal0(gel(is_princ, 1)))
            {
                printf(ANSI_COLOR_RED "Problem in my_H90_vect\n" ANSI_COLOR_RESET);
                pari_close();
                exit(111);
            }
            
            t = gel(is_princ, 2);
            printf("t: %ld\n", glength(t));
            if (glength(t)==0)
            {
                return stoi(-1);
            }
            
            t_rel = rnfeltabstorel(Lrel, t);
            Nt = rnfeltnorm(Lrel, t_rel);
            diff = nfmul(K, Nt, a);
            exp = bnfisunit0(K, diff, NULL);
            // pari_printf(ANSI_COLOR_YELLOW "exp: %Ps\n" ANSI_COLOR_RESET, exp);
            // pari_printf("M: %Ps\n", my_norm_operator(Labs, Lrel, K, p));
            // pari_printf("D: %Ps\n", my_xcol(glength(exp), itos(p)));
            // pari_printf("B: %Ps\n", gtocol(exp));
            /* check if exp lies in the image of the operator associated to N: (O_L)^x -> (O_K)^x */
            // if (my_QV_equal0(my_norm_operator(Labs, Lrel, K, p)))
            // {
            //     gel(I_vect, i) = F_ker_T;
            //     done = 1;
            //     break;
            // }
            is_norm = matsolvemod(my_norm_operator(Labs, Lrel, K, p), zerocol(glength(exp)), gtocol(exp), 0);
            //pari_printf(ANSI_COLOR_CYAN "is_norm: %Ps\n" ANSI_COLOR_RESET, is_norm);
            if (my_QV_equal0(exp) || !gequal0(is_norm))
            {
                //pari_printf(ANSI_COLOR_CYAN "Norm of elt with exp: %Ps\n" ANSI_COLOR_RESET, is_norm);
                // w = gel(bnf_get_fu(Labs), 1);
                // w_rel = rnfeltabstorel(Lrel, w);
                // Nw = rnfeltnorm(Lrel, w_rel);
                // pari_printf("Norm w: %Ps\n", Nw);
                I = F_ker_T;
                done = 1;
                break;
            }
        }
    }
    if (!done)
    {
        printf(ANSI_COLOR_RED "my_H90_vect ended with a problem \n" ANSI_COLOR_RESET);
        return stoi(-1);
        // pari_close();
        // exit(111);
    }
    
    
    I = gerepilecopy(av, I);
    return I;
}

GEN my_H90_exp_vect (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN Ja_vect, GEN p) {
    // This function is not ready to be used
    printf(ANSI_COLOR_CYAN "\nmy_H90_exp_vect\n" ANSI_COLOR_RESET);
    pari_sp av = avma;
    int r_rk = glength(Ja_vect), f, j, done = 0;
    GEN I_vect = zerovec(r_rk), a, iJ, F, ker_T, ker_T_basis, ker_sol, F_ker_T, t, t_rel, Nt, diff, exp, is_princ, is_norm, cyc;
    cyc = gtocol(bnf_get_cyc(Labs));
    int i;
    
    for (i = 1; i < r_rk+1; ++i)
    {
        done = 0;
        a = gel(gel(Ja_vect, i), 1);
        iJ = rnfidealup0(Lrel, gel(gel(Ja_vect, i),2), 1);
        F = my_H90_exp(Labs, iJ, sigma);
        // pari_printf("F: %Ps\n", F);
        ker_sol = matsolvemod(my_1MS_operator(Labs, sigma), cyc, zerocol(glength(cyc)), 1);
        ker_T_basis = gtovec(gel(ker_sol, 2));
        printf(ANSI_COLOR_GREEN "------------------------\n\n\nker_T_basis size: %ld\n\n\n------------------------\n" ANSI_COLOR_RESET, glength(ker_T_basis));
        
        if (glength(ker_T_basis)==0)
        {
            gel(I_vect, i) = F;
            done = 1;
        }
        else {
            ker_T = my_get_sums(ker_T_basis);

            // WARNING: This makes it faster but might not find the answer. Switch back to the above ker_T
            //ker_T = shallowconcat(mkvec(zerocol(glength(gel(ker_T_basis, 1)))), ker_T_basis);
            
            //pari_printf(ANSI_COLOR_CYAN "ker_T: %Ps\n" ANSI_COLOR_RESET, ker_T);
            // pari_printf("ker_T[2]: %Ps\n", gel(ker_T, 2));
            f = glength(ker_T);
            printf("Searching a chunk of ker (1-sigma) of size: %d\n", f);
            
            /*
                Find F_ker_T --------------------
            */
            
            for (j = 1; j < f+1; j++)
            {
                printf("%d/%d\n", j, f);
                // pari_printf("Adding the exp: %Ps\n", gel(ker_T, j));
                // pari_printf("And the ideal: %Ps\n", my_ideal_from_exp(Labs, gel(ker_T, j)));
                F_ker_T = ZC_add(F, gel(ker_T, j));
                //pari_printf("F_ker_T: %Ps\n", F_ker_T);

                // Here we should check, using linear algebra--------- 
                is_princ = bnfisprincipal0(Labs, idealdiv(Labs, iJ, my_1MS_ideal(Labs, sigma, F_ker_T)), 1);
                //-----------------------------------
                if (!my_QV_equal0(gel(is_princ, 1)))
                {
                    printf(ANSI_COLOR_RED "Problem in my_H90_exp_vect\n" ANSI_COLOR_RESET);
                    pari_close();
                    exit(111);
                }
                
                t = gel(is_princ, 2);
                
                t_rel = rnfeltabstorel(Lrel, t);
                Nt = rnfeltnorm(Lrel, t_rel);
                diff = nfmul(K, Nt, a);
                exp = bnfisunit0(K, diff, NULL);
                // pari_printf(ANSI_COLOR_YELLOW "exp: %Ps\n" ANSI_COLOR_RESET, exp);
                // pari_printf("M: %Ps\n", my_norm_operator(Labs, Lrel, K, p));
                // pari_printf("D: %Ps\n", my_xcol(glength(exp), itos(p)));
                // pari_printf("B: %Ps\n", gtocol(exp));
                /* check if exp lies in the image of the operator associated to N: (O_L)^x -> (O_K)^x */
                // if (my_QV_equal0(my_norm_operator(Labs, Lrel, K, p)))
                // {
                //     gel(I_vect, i) = F_ker_T;
                //     done = 1;
                //     break;
                // }
                is_norm = matsolvemod(my_norm_operator(Labs, Lrel, K, p), zerocol(glength(exp)), gtocol(exp), 0);
                //pari_printf(ANSI_COLOR_CYAN "is_norm: %Ps\n" ANSI_COLOR_RESET, is_norm);
                if (my_QV_equal0(exp) || !gequal0(is_norm))
                {
                    //pari_printf(ANSI_COLOR_CYAN "Norm of elt with exp: %Ps\n" ANSI_COLOR_RESET, is_norm);
                    // w = gel(bnf_get_fu(Labs), 1);
                    // w_rel = rnfeltabstorel(Lrel, w);
                    // Nw = rnfeltnorm(Lrel, w_rel);
                    // pari_printf("Norm w: %Ps\n", Nw);
                    gel(I_vect, i) = F_ker_T;
                    done = 1;
                    break;
                }
            }
        }
        if (!done)
        {
            printf(ANSI_COLOR_RED "my_H90_exp_vect ended with a problem \n" ANSI_COLOR_RESET);
            pari_close();
            exit(111);
        }
        
         
        
    } 
    printf("Finished search\n");
    I_vect = gerepilecopy(av, I_vect);
    return I_vect;
}

void my_compute_unit_norms(GEN Lrel, GEN K, GEN units_group) {
    
    GEN rel_unit, norm;
    int l = glength(units_group);
    
    int i;
    for (i = 1; i < l+1; i++)
    {
        rel_unit = rnfeltabstorel(Lrel, gel(units_group, i));
        norm = rnfeltnorm(Lrel, rel_unit);

        pari_printf("norm[%d]: %Ps\n", i, bnfisunit0(K, norm, NULL));

    }
} 

void my_norms(GEN K, GEN K_ext, GEN p) {
    GEN Labs, Lrel, unit_group;
    int l = glength(K_ext);
    int i;
    for (i = 1; i < l+1; i++)
    {
        Labs = gel(gel(K_ext, i), 1);
        unit_group = my_find_units_mod_p(Labs, p);
        Lrel = gel(gel(K_ext, i), 2);
        my_compute_unit_norms(Lrel, K, unit_group);
    }
    
}

void my_compute_unit_diffs(GEN Labs, GEN sigma, GEN units_group) {
    
    GEN diff;
    int l = glength(units_group);
    
    int i;
    for (i = 1; i < l+1; i++)
    {
        diff = my_1MS_elt(Labs, sigma, gel(units_group, i));

        pari_printf("diff[%d]: %Ps\n", i, bnfisunit0(Labs, diff, NULL));

    }
} 

void my_diffs(GEN K_ext, GEN p) {
    GEN Labs, unit_group, sigma;
    int l = glength(K_ext);
    int i;
    for (i = 1; i < l+1; i++)
    {
        Labs = gel(gel(K_ext, i), 1);
        sigma = gel(gel(K_ext, i), 3);
        unit_group = my_find_units_mod_p(Labs, p);
       
        my_compute_unit_diffs(Labs, sigma, unit_group);
    }
    
}

void my_compute_matrices(GEN Labs, GEN sigma, GEN units_group) {
    
    GEN image;
    int l = glength(units_group);
    
    int i;
    printf("1:\n\n");
    for (i = 1; i < l+1; i++)
    {
        image = bnfisunit0(Labs, gel(units_group, i), NULL);

        pari_printf("%Ps\n", image);

    }
    printf("\nsigma:\n\n");
    for (i = 1; i < l+1; i++)
    {
        image = bnfisunit0(Labs, galoisapply(Labs, sigma, gel(units_group, i)), NULL);

        pari_printf("%Ps\n", image);

    }
    printf("\nsigma^2:\n\n");
    for (i = 1; i < l+1; i++)
    {
        image = bnfisunit0(Labs, galoisapply(Labs, sigma, galoisapply(Labs, sigma, gel(units_group, i))), NULL);

        pari_printf("%Ps\n", image);

    }
} 

void my_matrices(GEN K_ext, GEN p) {
    GEN Labs, unit_group, sigma;
    int l = glength(K_ext);
    int i;
    for (i = 1; i < l+1; i++)
    {
        Labs = gel(gel(K_ext, i), 1);
        sigma = gel(gel(K_ext, i), 3);
        unit_group = my_find_units_mod_p(Labs, p);
        printf("Field ext [%d]:\n\n", i);
        my_compute_matrices(Labs, sigma, unit_group);
    }
    
}

GEN my_find_I_full (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN i_xJ, GEN a, GEN class_group, int p_int)
{
    pari_sp av = avma;
    GEN diff, current_I, test_vec, class_number = bnf_get_no(Labs), Itest, t, t_rel, norm_t, exponents, norms, L_units, L_unit_group, ret;
    int clnr = itos(class_number), n, roots_of_unity_nr = bnf_get_tuN(K), l;
    
    
    // for (n = 1; n < clnr + 1; ++n) {
    //     printf("%d/%d\n", n, clnr);
    //     current_I = gel(class_group, n);
    //     Itest = idealdiv(Labs, i_xJ, my_1MS_ideal(Labs, sigma, current_I));
    //     test_vec = bnfisprincipal0(Labs, Itest, 1);

    //     if (my_QV_equal0(gel(test_vec, 1)))
    //     {
    //         t = gel(test_vec, 2);
    //         t_rel = rnfeltabstorel(Lrel, t);
    //         norm_t = rnfeltnorm(Lrel, t_rel);
    //         exponents = bnfisunit(K, nfmul(K, norm_t, a));
    //         pari_printf("exponents: %Ps\n", exponents);
    //         if (my_is_p_power(exponents, p_int, roots_of_unity_nr))
    //         {
    //             printf("I found\n\n");
    //             ret = mkvec2(gel(test_vec, 2), current_I);
    //             ret = gerepilecopy(av, ret);
    //             return ret;
    //         }
            
    //     } 

    // }

    //printf("Entering the slow round\n\n");
    printf("roots of unity L: %d\n\n", roots_of_unity_nr);
    L_units = my_find_units_mod_p(Labs, stoi(p_int));
    printf("nr of L_units mod p: %ld\n\n", glength(L_units));
    L_unit_group = my_get_unit_group(Labs, L_units, stoi(p_int));
    printf("size of L_unit_group: %ld\n\n", glength(L_unit_group));
    
    l = glength(L_unit_group);
    norms = zerovec(l);
    
    for (n = 1; n < l+1; n++)
    {
        gel(norms, n) = algtobasis(K, rnfeltnorm(Lrel, rnfeltabstorel(Lrel, gel(L_unit_group, n))));
    }
    // pari_close();
    // exit(0);
    int i;
    for (n = 1; n < clnr + 1; ++n) {
        printf("%d/%d\n", n, clnr);
        current_I = gel(class_group, n);
        Itest = idealdiv(Labs, i_xJ, my_1MS_ideal(Labs, sigma, current_I));
        test_vec = bnfisprincipal0(Labs, Itest, nf_GEN);

        if (my_QV_equal0(gel(test_vec, 1)))
        {
            t = gel(test_vec, 2);
            pari_printf("t: %Ps\n", t);
            printf("length t: %ld\n", glength(t));
            if (glength(t)!=0)
            {

                t_rel = rnfeltabstorel(Lrel, t);
                norm_t = rnfeltnorm(Lrel, t_rel);
                diff = nfmul(K, norm_t, a);
                for ( i = 1; i < l+1; i++)
                {
                    exponents = bnfisunit0(K, nfdiv(K, diff, gel(norms, i)), NULL);
                    if (l%1000==0)
                    {
                        pari_printf("exponents[%d]: %Ps\n", i, exponents);
                    }
                    
                    if (my_is_p_power(exponents, p_int, roots_of_unity_nr))
                    {
                        printf(ANSI_COLOR_GREEN "\nI found\n\n" ANSI_COLOR_RESET);
                        ret = mkvec2(gel(test_vec, 2), current_I);
                        ret = gerepilecopy(av, ret);
                        return ret;
                    }
                    
                }
            }
                
            
        } 

    }
    printf(ANSI_COLOR_RED "my_find_I_full ended with a problem\n" ANSI_COLOR_RESET);
   
    printf(ANSI_COLOR_RED "No I found in my_find_I\n\n" ANSI_COLOR_RESET);
    pari_close();
    exit(111);

    
    return 0;
}

GEN my_find_I_vect_full (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN Ja_vect, int p_int) {
    pari_sp av = avma;
    int r_rk = glength(Ja_vect);
    GEN I_vect = zerovec(r_rk), a, i_xJ;
   
    GEN class_group = my_get_clgp (Labs);
   
    int i;
    
    for (i = 1; i < r_rk+1; ++i)
    {
        a = gel(gel(Ja_vect, i), 1);
        i_xJ = rnfidealup0(Lrel, gel(gel(Ja_vect, i),2), 1);
        gel(I_vect, i) = gel(my_find_I_full(Labs, Lrel, K, sigma, i_xJ, a, class_group, p_int), 2);
    } 
    
    I_vect = gerepilecopy(av, I_vect);
    return I_vect;
}
