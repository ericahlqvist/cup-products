

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
        exit(0);
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

GEN my_find_prime_vect(GEN LyAbs, GEN sigma_y, GEN p_1, int p) {
    // printf(ANSI_COLOR_CYAN "my_find_prime_vect\n\n" ANSI_COLOR_RESET);
    pari_sp av = avma;
    GEN prime_vect = zerovec(p);
    gel(prime_vect, 1) = idealhnf0(LyAbs, p_1, NULL);

    int i;
    
    GEN new_p = idealmul(LyAbs, idealhnf0(LyAbs, gen_1, NULL), p_1);
    for (i = 1; i < p; i++)
    {
        
        new_p = galoisapply(LyAbs, sigma_y, new_p);
        gel(prime_vect, i+1) = new_p;
    }
    prime_vect = gerepilecopy(av, prime_vect);
    return prime_vect;
}

GEN my_find_e_vect(GEN LyAbs, GEN sigma_y, GEN prime_vect, GEN primes, GEN es, int p) {
    // printf(ANSI_COLOR_CYAN "my_find_e_vect\n\n" ANSI_COLOR_RESET);
    pari_sp av = avma;
    GEN e_vect = zerovec(p);
    GEN current_e = pol_x(fetch_user_var("current_e"));
    GEN current_prime = pol_x(fetch_user_var("current_prime"));

    int i;
    int j;
    
    for (i = 1; i < glength(es)+1; i++)
    {
        current_e = gel(es,i);
        current_prime = gel(primes, i);
        for (j = 1; j < p+1; j++)
        {
            if (my_SQ_MAT_equal(current_prime, gel(prime_vect, j)))
            {
                // output(current_e);
                gel(e_vect, j) = current_e;
            }
        } 
    }
    e_vect = gerepilecopy(av, e_vect);
    return e_vect;
}

GEN my_find_primes_under(GEN LyRel, GEN K, GEN prime_vect) {
    pari_sp av = avma;
    GEN primes_under = mkvec(idealhnf0(K, rnfidealdown(LyRel, rnfidealabstorel(LyRel, gel(prime_vect, 1))), NULL));
    int l = glength(prime_vect);
    // GEN primes_under_unsorted = zerovec(l);
    GEN current_prime = pol_x(fetch_user_var("current_prime"));
    
    int i;
    int j;
    int test;
    
    for (i = 2; i < l+1; i++)
    {
        current_prime = idealhnf0(K, rnfidealdown(LyRel, rnfidealabstorel(LyRel, gel(prime_vect, i))), NULL);
        
        test = 1;
        for (j = 1; j < glength(primes_under)+1; j++) {

            if (my_SQ_MAT_equal(gel(primes_under, j), current_prime)) {
                test = 0;
            }
        }

        if (test) {
            primes_under = gconcat(primes_under, mkvec(current_prime));
        }
        
    }
    // printf("Primes under: ");
    // output(primes_under);
    // printf("\n\n");
    primes_under = gerepilecopy(av, primes_under);
    return primes_under;
}

GEN my_sort_primes_and_es(GEN LyRel, GEN K, GEN primes_and_es_in_factorization, GEN primes_under) {
    pari_sp av = avma;
    int l = glength(primes_under);
    GEN prime_vect = gel(primes_and_es_in_factorization, 1);
    GEN e_vect = gel(primes_and_es_in_factorization, 2);
    GEN primes_and_es_sorted = zerovec(l);
    GEN current_prime = pol_x(fetch_user_var("current_prime"));

    int i;
    int j;
    for (i = 1; i < l+1; i++)
    {
        gel(primes_and_es_sorted, i) = mkvec(gen_0);
        for (j = 1; j < glength(prime_vect)+1; j++) {
            current_prime = idealhnf0(K, rnfidealdown(LyRel, rnfidealabstorel(LyRel, gel(prime_vect, j))), NULL);
            if (my_SQ_MAT_equal(current_prime, gel(primes_under, i)))
            {
                if (!my_QV_equal0(gel(primes_and_es_sorted, i))) {
                    gel(gel(primes_and_es_sorted, i), 1) = gconcat(gel(gel(primes_and_es_sorted, i), 1), mkvec(gel(prime_vect, j))); 
                    gel(gel(primes_and_es_sorted, i), 2) = gconcat(gel(gel(primes_and_es_sorted, i), 2), mkvec(gel(e_vect, j))); 
                }
                else {
                    gel(primes_and_es_sorted, i) = mkvec2(mkvec(gel(prime_vect, j)), mkvec(gel(e_vect, j)));
                }
            }
            
        }
    }
    // printf("Primes and Es sorted: ");
    // output(primes_and_es_sorted);
    // printf("\n\n");
    primes_and_es_sorted = gerepilecopy(av, primes_and_es_sorted);
    return primes_and_es_sorted;
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

GEN my_find_H90_ideal_single_prime (GEN LyAbs, GEN LyRel, GEN K, GEN primes, GEN es, GEN sigma_y, int p) {
    pari_sp av = avma;
    
    GEN p_1 = gel(primes,1);
    GEN prime_vect = my_find_prime_vect(LyAbs, sigma_y, p_1, p);
    GEN e_vect = my_find_e_vect(LyAbs, sigma_y, prime_vect, primes, es, p);
    
    
    GEN H90_ideal = idealhnf0(LyAbs, gen_1, NULL);
    GEN new_ideal = idealhnf0(LyAbs, gen_1, NULL);
    int i;
    int j;
    
    for (i = 2; i < glength(prime_vect)+1; i++)
    {
        new_ideal = idealhnf0(LyAbs, gen_1, NULL);
        for (j = 1; j < i; j++)
        {
            // printf("p[%d]", j);
            new_ideal = idealmul(LyAbs, new_ideal, gel(prime_vect, j));
        }
        // printf("^e[%d]", i);
        // printf("\n\n");
        
        new_ideal = idealpow(LyAbs, new_ideal, gel(e_vect, i));
        H90_ideal = idealmul(LyAbs, H90_ideal, new_ideal);
    }
    H90_ideal = gerepilecopy(av, H90_ideal);
    return H90_ideal;
}



GEN my_find_H90_ideal (GEN LyAbs, GEN LyRel, GEN K, GEN iJ_div_a2, GEN sigma_y, int p) {
    pari_sp av = avma;
    GEN factorization = idealfactor(LyAbs, iJ_div_a2);
    
    GEN primes_and_es_in_factorization = my_find_primes_in_factorization(LyAbs, factorization);
    //pari_printf("primes_and_es: %Ps\n\n", primes_and_es_in_factorization);
    // output(primes_and_es_in_factorization);
    GEN prime_vect = gel(primes_and_es_in_factorization,1);
    //pari_printf("prime_vect: %Ps\n\n", prime_vect);
    GEN primes_under = my_find_primes_under(LyRel, K, prime_vect);
    GEN primes_and_es_by_primes_under = my_sort_primes_and_es(LyRel, K, primes_and_es_in_factorization, primes_under);
    //printf("Factorization for H90 done\n\n");
    

    GEN H90_ideal = idealhnf0(LyAbs, gen_1, NULL);
    GEN new_ideal = idealhnf0(LyAbs, gen_1, NULL);
    int i;
    
    for (i = 1; i < glength(primes_under)+1; i++)
    {
        new_ideal = my_find_H90_ideal_single_prime(LyAbs, LyRel, K, gel(gel(primes_and_es_by_primes_under, i), 1), gel(gel(primes_and_es_by_primes_under, i), 2), sigma_y, p);
        // printf("my_find_H90_ideal_single_prime[%d] done", i);
        H90_ideal = idealmul(LyAbs, H90_ideal, new_ideal);
        //H90_ideal = idealred(LyAbs, idealmul(LyAbs, H90_ideal, new_ideal));
    }
    
    // if (!my_SQ_MAT_equal(my_SM1_ideal(LyAbs, sigma_y, H90_ideal), idealred(LyAbs, iJ_div_a2))) {
    //     printf("Wrong H90:\n\n");
    //     pari_printf("%Ps\n\n", my_SM1_ideal(LyAbs, sigma_y, H90_ideal));
    //     pari_printf("%Ps\n\n", idealred(LyAbs, iJ_div_a2));
    //     pari_close();
    //     exit(0);
    // }
    if (!my_SQ_MAT_equal(my_SM1_ideal(LyAbs, sigma_y, H90_ideal), iJ_div_a2)) {
        printf("Wrong H90:\n\n");
        pari_printf("%Ps\n\n", my_SM1_ideal(LyAbs, sigma_y, H90_ideal));
        pari_printf("%Ps\n\n", iJ_div_a2);
        pari_close();
        exit(0);
    }
    else {
        printf(ANSI_COLOR_GREEN "------------\nH90 TEST PASSED\n-----------\n" ANSI_COLOR_RESET);
    }
    H90_ideal = gerepilecopy(av, H90_ideal);
    return H90_ideal;
}

GEN my_1MS_operator (GEN L, GEN sigma) {
    pari_sp av = avma;
    GEN cyc = bnf_get_cyc(L); 
    GEN gens = bnf_get_gen(L);
    int l = glength(cyc);
    GEN M = zeromatcopy(l,l);
    int i;
    for (i = 1; i < l+1; i++)
    {
        gel(M, i) = bnfisprincipal0(L, my_1MS_ideal(L, sigma, gel(gens, i)), 0);
    }
    
    M = gerepilecopy(av, M);
    return M;
}


GEN my_norm_operator (GEN Labs, GEN Lrel, GEN K, GEN p) {
    pari_sp av = avma;
    GEN M, Lgens, Kgens, g, g_rel, Ng;
    int Klength, Llength, i;
    Lgens = my_find_units_mod_p(Labs, p);
    Kgens = my_find_units_mod_p(K, p);
    Llength = glength(Lgens);
    Klength = glength(Kgens);
    M = zeromatcopy(Klength, Llength);

    for (i = 1; i < Llength+1; i++)
    {
        g = gel(Lgens, i);
        g_rel = rnfeltabstorel(Lrel, g);
        Ng = rnfeltnorm(Lrel, g_rel);
        gel(M, i) = bnfisunit(K, Ng);
    }
    

    M = gerepilecopy(av, M);
    return M;
}

GEN my_H90 (GEN L, GEN iJ, GEN sigma) {
    pari_sp av = avma;
    GEN H90_ideal, gens, M, D, B, E, Itest, test_vec;
    gens = bnr_get_gen(L);
    int l = glength(gens);
    
    M = zeromatcopy(l,l);
    D = zerovec(l);

    // finding exponents for iJ 
    B = bnfisprincipal0(L, iJ, 0);

    // finding the matrix M consisting of exponents for the (1-s)g_i, Cl(L) = < g_1, ..., g_n > 
    int i;
    for (i = 1; i < l+1; i++)
    {
        gel(M, i) = bnfisprincipal0(L, my_1MS_ideal(L, sigma, gel(gens, i)), 0);
    }
    // Gauss
    E = matsolvemod(M,D,B,0);

    H90_ideal = idealhnf0(L, gen_1, NULL);
    for (i = 1; i < l+1; i++)
    {
        H90_ideal = idealmul(L, H90_ideal, idealpow(L, gel(gens, i), gel(E, i)));
       
    }

    H90_ideal = gerepilecopy(av, H90_ideal);
    return H90_ideal;
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

GEN my_factor_I_by_primes_below (GEN Labs, GEN Lrel, GEN K, GEN I) {
    pari_sp av = avma;
    GEN factorization = idealfactor(Labs, I);
    GEN primes_and_es_in_factorization = my_find_primes_in_factorization(Labs, factorization);
    GEN prime_vect = gel(primes_and_es_in_factorization,1);
    // GEN e_vect = gel(primes_and_es_in_factorization,2);
    GEN primes_under = my_find_primes_under(Lrel, K, prime_vect);
    GEN primes_and_es_by_primes_under = my_sort_primes_and_es(Lrel, K, primes_and_es_in_factorization, primes_under);
    // pari_printf("e_vect: %Ps\n\n", e_vect);
    // pari_printf("prime_vect: %Ps\n\n", prime_vect);

    primes_and_es_by_primes_under = gerepilecopy(av, primes_and_es_by_primes_under);
    return primes_and_es_by_primes_under;
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

GEN my_ideal_class_representative (GEN K, GEN class_group, GEN ideal) {
    GEN rep, test;

    int i;
    for (i = 1; i < glength(class_group)+1; i++)
    {
        rep = gel(class_group, i);
        test = gel(bnfisprincipal0(K, idealdiv(K, ideal, rep), 1), 1);
        if (my_QV_equal0(test))
        {
            break;
        }
        
    }
    return rep;
}

GEN my_reduce_ideal_by_ideals_below (GEN Labs, GEN Lrel, GEN K, GEN ideal) {

    GEN red_ideal = ideal;
    GEN to_be_removed, current_prime;
    GEN factorization = idealfactor(Labs, ideal);

    int i;
    for (i = 1; i < glength(gel(factorization, 1))+1; i++)
    {
        // Remove inert ideals
        if (pr_is_inert(gmael2(factorization, 1, i)))
        {
            current_prime = idealhnf0(Labs, gmael3(factorization, 1, i, 1), gmael3(factorization, 1, i, 2));
            to_be_removed = idealpow(Labs, current_prime, gmael2(factorization, 2, i));
            red_ideal = idealdiv(Labs, red_ideal, to_be_removed);
            printf(ANSI_COLOR_MAGENTA "Removed inert ideal\n\n" ANSI_COLOR_RESET);
        }
        
    }
    

    return red_ideal;
}

// GEN my_represent_ideal_class_by_prime (GEN K, GEN ideal) {
//     GEN list_of_primes = primes_interval(gen_1, stoi(1000000));
//     printf("Hej\n\n");
    
//     GEN Q = nfinit(mkpoln(0,gen_0), DEFAULTPREC);
    
//     GEN Krel = rnfinit(Q, nf_get_pol(bnf_get_nf(K)));
    
//     GEN prime_representative, factorization, primes_and_es, primes, test, test_ideal;
//     int verif = 1;

//     int i;
//     int j;
//     for (i = 1; i < glength(list_of_primes)+1; i++)
//     {
//         factorization = idealfactor(K, rnfidealup0(Krel, gel(list_of_primes, i), 1));
//         primes_and_es = my_find_primes_in_factorization(K, factorization);
//         primes = gel(primes_and_es, 1);
//         for (j = 1; j < glength(primes)+1; j++)
//         {
//             test_ideal = idealdiv(K, ideal, gel(primes, j));
//             test = gel(bnfisprincipal0(K, test_ideal, 1), 1);
//             if (my_QV_equal0(test))
//             {
//                 prime_representative = gel(primes, j);
//                 verif = 0;
//                 break;
//             }
            
//         }
        
//     }
    
//     if (verif)
//     {
//         printf(ANSI_COLOR_RED "ERROR: No prime representative for ideal class found\n\n" ANSI_COLOR_RESET);
//         printf("Proceeds with I as it was\n\n");
//         prime_representative = ideal;
//     }
    

//     return prime_representative;
// }
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
        
        gel(class_group, n) = idealred(K, current_I); 
    }
    
    class_group = gerepilecopy(av, class_group);
    return class_group;
}

GEN my_get_clgp_p (GEN K, GEN p) {
    GEN clgp = my_get_clgp(K);
    int i; 
    for (i=1;i<glength(clgp)+1;i++) {
        gel(clgp, i) = idealpow(K, gel(clgp, i), p);
    }
    return clgp;
}

GEN my_is_principal_mod_p (GEN K, GEN I, GEN p) {
    
    GEN clgp = my_get_clgp(K);
    int clnr = glength(clgp);
    GEN p_J, is_pr;
    GEN test_ideal, test_vec;
    int i;
    for (i = 1; i < clnr+1; i++)
    {
        p_J = idealred(K, idealpow(K, gel(clgp, i), p));
        test_ideal = idealdiv(K, I, p_J);
        is_pr = bnfisprincipal0(K, test_ideal, 1);
        test_vec = gel(is_pr, 1);
        //output(test_vec);
        if (my_QV_equal0(test_vec))
        {
            printf("Norms found\n\n");
            output(gel(is_pr, 2));
            return mkvec2(gel(is_pr, 2), gel(clgp, i));
        }
        

    }
    

    return mkvec2(gen_0, gen_0);
}

int my_is_principal_mod_N (GEN Labs, GEN Lrel, GEN K, GEN I) {
    
    GEN clgp = my_get_clgp(Labs);
    int clnr = glength(clgp);
    GEN N;
    GEN test_ideal, test_vec;
    int i;
    for (i = 1; i < clnr+1; i++)
    {
        N = idealred(K, rnfidealnormrel(Lrel, rnfidealabstorel(Lrel, gel(clgp, i))));
        test_ideal = idealmul(K, I, N);
        test_vec = bnfisprincipal0(K, test_ideal, 0);
        // output(bnfisprincipal0(K, test_ideal, 0));
        if (my_QV_equal0(test_vec))
        {
            return 1;
            break;
        }
        

    }
    
    return 0;
}

GEN my_reduce_H1_mod (GEN LxAbs, GEN sigma_x, GEN v, GEN elem, GEN p) {

    GEN new_v = zerovec(4);
    gel(new_v, 1) = nfdiv(LxAbs, gel(v, 1), nfpow(LxAbs, elem, p));
    gel(new_v, 2) = nfmul(LxAbs, gel(v, 2), my_SM1_elt(LxAbs, sigma_x, elem));
    gel(new_v, 3) = gel(v, 3);
    gel(new_v, 4) = idealmul(LxAbs, idealhnf0(LxAbs, elem, NULL), gel(v, 4));

    return new_v;
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
    printf("torsion order: %d\n\n", tors_gp_ord);

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

GEN my_find_p_gens (GEN K, GEN p)
{
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
        gel(p_gens, i) = idealred(K, current_gen);
    }
    
    return p_gens;
}

int my_action_on_clgp (GEN L, GEN sigma, GEN p) {
    GEN gens = bnf_get_gen(L);
    GEN clgpxp = my_get_clgp_p(L, p);
    GEN test_ideal, test_vec;

    
    int i;
    int j;
    int check = 0;
    GEN ideal;
    for (i=1;i<glength(gens)+1;i++) {
        ideal = my_SM1_ideal(L, sigma, gel(gens, i));
        for (j=1; j<glength(clgpxp)+1;j++) {
            test_ideal = idealred(L, idealmul(L, ideal, gel(clgpxp, j)));
            test_vec = bnfisprincipal0(L, test_ideal, 1);
            if (my_QV_equal0(gel(test_vec, 1))) {
                check = 1;
                break;
            }
            // if (my_SQ_MAT_equal(test_ideal, unit)) {
            //     check = 1;
            //     break;
            // }
        }
        if (check) {
            printf(ANSI_COLOR_GREEN "ok\n" ANSI_COLOR_RESET);
        }
        else {
            printf(ANSI_COLOR_RED "-\n" ANSI_COLOR_RESET);
        }
        check = 0;
    }
    printf("\n\n");
    return 0;
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

void my_unramified_p_extensions_with_trivial_action(GEN K, GEN p, GEN D_prime_vect) {
    GEN x, y, p1, q1, p1red, Lrel, Labs, s_lift_x, cx, sigma;

    x = pol_x(fetch_user_var("x"));
    y = pol_x(fetch_user_var("y"));

    GEN s = pol_x(fetch_user_var("s"));
    GEN index = mkvec(p);
    // GEN R = nfsubfields0(clf_pol,4,1);
    // pari_printf("subgrouplist: %Ps\n", subgrouplist0(bnf_get_cyc(K), mkvec(gpow(p,gen_2, DEFAULTPREC)), 0));

    GEN R = bnrclassfield(K, subgrouplist0(bnf_get_cyc(K), index, 2), 1, DEFAULTPREC);
    pari_printf("R: %Ps\n", R);
    GEN Rsq = bnrclassfield(K, subgrouplist0(bnf_get_cyc(K), mkvec(gpow(p,gen_2, DEFAULTPREC)), 2), 1, DEFAULTPREC);
    pari_printf("Rsq: %Ps\n", Rsq);
    GEN Rcb = bnrclassfield(K, subgrouplist0(bnf_get_cyc(K), mkvec(gpow(p,stoi(3), DEFAULTPREC)), 2), 1, DEFAULTPREC);
    pari_printf("Rcb: %Ps\n", Rcb);
    
    printf("[p]-extensions:\n\n");
    int i, j;
    for (i=1; i<glength(R)+1; ++i) {
        p1 = gel(R, i);
        q1 = gsubstpol(p1, x, y);

        /* Define Lrel/Labs */
        p1red = rnfpolredbest(K, mkvec2(q1, D_prime_vect), 0);
        Lrel = rnfinit(K, p1red);
        //printf("Lrel found\n");
        Labs = Buchall(polredabs0(mkvec2(rnf_get_polabs(Lrel), D_prime_vect), 0), nf_FORCE, DEFAULTPREC);
        //printf("Labs found\n");
        
        pari_printf("abs pol: %Ps, ", polredabs0(mkvec2(gsubstpol(rnf_get_polabs(Lrel),y, s), D_prime_vect), 0));
        pari_printf("L_cyc[%d]: %Ps\n", i, bnf_get_cyc(Labs));

        s_lift_x = rnfeltup0(Lrel, s, 1);
        cx = galoisconj(Labs, NULL);

        for (j = 1; j < glength(cx)+1; ++j)
        {
            if ( (!my_QV_equal(algtobasis(Labs,gel(cx, j)), algtobasis(Labs,y))) && my_QV_equal(galoisapply(Labs, gel(cx,j), s_lift_x), s_lift_x)) 
            {
                sigma = algtobasis(Labs, gel(cx, j));
                //pari_printf("sigma: %Ps\n\n", sigma);
                break;
            }
        }
        //printf(ANSI_COLOR_CYAN "---> sigma <--- \n \n" ANSI_COLOR_RESET);
        
        my_action_on_clgp(Labs, sigma, p);

        //Should throw away some trash here to free memory
        //my_test_artin_symbol (Labs, Lrel, base, itos(p), sigma);
    }
    printf("\n");
    printf("[p,p]-extensions:\n\n");
    for (i=1; i<glength(Rsq)+1; ++i) {
        p1 = gel(Rsq, i);
        q1 = gsubstpol(p1, x, y);

        /* Define Lrel/Labs */
        p1red = rnfpolredbest(K, mkvec2(q1, D_prime_vect), 0);
        Lrel = rnfinit(K, p1red);
        //printf("Lrel found\n");
        Labs = Buchall(rnf_get_polabs(Lrel), nf_FORCE, DEFAULTPREC);
        //printf("Labs found\n");
        
        pari_printf("abs pol: %Ps, ", gsubstpol(rnf_get_polabs(Lrel),y, s));
        pari_printf("L_cyc[%d]: %Ps\n", i, bnf_get_cyc(Labs));

        s_lift_x = rnfeltup0(Lrel, s, 1);
        cx = galoisconj(Labs, NULL);

        for (j = 1; j < glength(cx)+1; ++j)
        {
            if ( (!my_QV_equal(algtobasis(Labs,gel(cx, j)), algtobasis(Labs,y))) && my_QV_equal(galoisapply(Labs, gel(cx,j), s_lift_x), s_lift_x)) 
            {
                sigma = algtobasis(Labs, gel(cx, j));
                //pari_printf("sigma: %Ps\n\n", sigma);
                break;
            }
        }
        //printf(ANSI_COLOR_CYAN "---> sigma <--- \n \n" ANSI_COLOR_RESET);
        
        my_action_on_clgp(Labs, sigma, p);

        //Should throw away some trash here to free memory
        //my_test_artin_symbol (Labs, Lrel, base, itos(p), sigma);
    }
    printf("\n");
    printf("[p,p,p]-extensions:\n\n");
    for (i=1; i<glength(Rcb)+1; ++i) {
        p1 = gel(Rcb, i);
        q1 = gsubstpol(p1, x, y);

        /* Define Lrel/Labs */
        p1red = rnfpolredbest(K, mkvec2(q1, D_prime_vect), 0);
        Lrel = rnfinit(K, p1red);
        //printf("Lrel found\n");
        Labs = Buchall(rnf_get_polabs(Lrel), nf_FORCE, DEFAULTPREC);
        //printf("Labs found\n");
        
        pari_printf("abs pol: %Ps, ", gsubstpol(rnf_get_polabs(Lrel),y, s));
        pari_printf("L_cyc[%d]: %Ps\n", i, bnf_get_cyc(Labs));

        s_lift_x = rnfeltup0(Lrel, s, 1);
        cx = galoisconj(Labs, NULL);

        for (j = 1; j < glength(cx)+1; ++j)
        {
            if ( (!my_QV_equal(algtobasis(Labs,gel(cx, j)), algtobasis(Labs,y))) && my_QV_equal(galoisapply(Labs, gel(cx,j), s_lift_x), s_lift_x)) 
            {
                sigma = algtobasis(Labs, gel(cx, j));
                //pari_printf("sigma: %Ps\n\n", sigma);
                break;
            }
        }
        //printf(ANSI_COLOR_CYAN "---> sigma <--- \n \n" ANSI_COLOR_RESET);
        
        my_action_on_clgp(Labs, sigma, p);

        //Should throw away some trash here to free memory
        //my_test_artin_symbol (Labs, Lrel, base, itos(p), sigma);
    }
    printf("\n");
}

GEN my_find_p_power_tors_units(GEN K, GEN p) {
    
    if (itos(p)!=2) {
        return mkvec(algtobasis(K, gen_1));
    }

    else {
        return mkvec4(algtobasis(K, gen_1), algtobasis(K, nfpow(K, bnf_get_tuU(K), p)), algtobasis(K, nfpow(K, bnf_get_tuU(K), gmul(p,p))), algtobasis(K, nfpow(K, bnf_get_tuU(K), gmul(p,gmul(p,p)))));
    }
}

GEN my_find_p_power_units(GEN K, GEN unit_group, GEN p) {
    // This function must be improved to return more ("small") powers
    
    int cyc_order = glength(unit_group);
    if (cyc_order==0) {
        return mkvec(algtobasis(K, gen_1));
    }

    if (cyc_order == 1) {
        return mkvec2(algtobasis(K, gen_1), algtobasis(K, nfpow(K, gel(unit_group, 1), p)));
    }

    int order = pow(pow(itos(p),3),cyc_order);
    //printf("order: %d\ncyc_order: %d\n", order, cyc_order);
    GEN cyc = zerovec(cyc_order);
    GEN p_power_units = zerovec(order);
    GEN current_power, power, exp;

    int i;
    int j;
    for (i=1;i<cyc_order+1;++i) {
        gel(cyc,i) = gpow(p, stoi(3), DEFAULTPREC);
    }
    GEN exponents = my_get_vect(cyc_order-1, cyc);
    //pari_printf("%Ps\n", exponents);
    for (i=1;i<order+1;++i) {
        exp = gel(exponents, i); 
        current_power = algtobasis(K, gen_1);
        
        for (j=1;j<cyc_order+1;++j) {
            if (itos(gel(exp, j))==2) {
                power = nfpow(K, gel(unit_group, j), gel(exp, j));
            }
            else {
                power = nfpow(K, gel(unit_group, j), gneg(gmul(gen_2,gel(exp, j))));
            }
            
            current_power = nfmul(K, current_power, power);
        }
        gel(p_power_units, i) = algtobasis(K, current_power);
    }

    for (i=1;i<glength(p_power_units)+1;++i) {
        pari_printf("power units %d: %Ps\n", i, gel(p_power_units, i));
    }
    
    return p_power_units;
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

int my_is_in_tors_power_units(GEN K, GEN a, GEN p_power_units) {
    int i; 
    int answer = 0;
    for (i=1;i<glength(p_power_units)+1;++i) {
        if (my_QV_equal(a, gel(p_power_units, i))) {
            answer = 1;
            break;
        }
    }
    return answer;
}

int my_is_in_power_units(GEN K, GEN a, GEN p_power_units) {
    int i; 
    int answer = 0;
    for (i=1;i<glength(p_power_units)+1;++i) {
        if (my_QV_equal(a, gel(p_power_units, i))) {
            answer = 1;
            break;
        }
    }
    return answer;
}

// Returns vector of tuples (a,J) with div(a)+pJ = 0.
GEN my_find_Ja_vect(GEN K, GEN J_vect, GEN p, GEN units_mod_p) {
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
    
    
    return Ja_vect;
}



/*
Given J \in CL(K)_p, find t, I such that J=div(t)+(1-sigma)I 
*/
GEN my_find_I (GEN Labs, GEN K, GEN sigma, GEN i_xJ, GEN class_group)
{
    pari_sp av = avma;
    GEN current_I;
    GEN test_vec;

    GEN class_number = bnf_get_no(Labs);
    int clnr = itos(class_number);
    
    GEN Itest;
    int test_found = 0;
    printf("Getting class group of Labs\n");
    //pari_printf("clgp: %Ps\n", bnf_get_clgp(Labs));
    int n;
    printf("Testing if i(J) principal\n");
    test_vec = bnfisprincipal0(Labs, i_xJ, 1);
    
    if (my_QV_equal0(gel(test_vec, 1)))
    {
        printf(ANSI_COLOR_YELLOW "i_x(J) principal!\n\n" ANSI_COLOR_RESET);
        current_I = idealhnf0(Labs, gen_1, NULL);
        test_found = 1;
    }
    

    else {
        printf(ANSI_COLOR_CYAN "i_x(J) not principal!\n\n" ANSI_COLOR_RESET);
        for (n = 1; n < clnr + 1; ++n) {
            printf("%d/%d\n", n, clnr);
            current_I = gel(class_group, n);
            Itest = idealdiv(Labs, i_xJ, my_1MS_ideal(Labs, sigma, current_I));
            test_vec = bnfisprincipal0(Labs, Itest, 1);

            if (my_QV_equal0(gel(test_vec, 1)))
            {
                //printf(ANSI_COLOR_GREEN "FOUND!\n\n" ANSI_COLOR_RESET);
                
                test_found = 1;
                break;
            } 

        }
    }
    
    if (!test_found)
    {
        printf(ANSI_COLOR_RED "No I found in my_find_I\n\n" ANSI_COLOR_RESET);
        exit(0);
    }
    
    GEN ret = mkvec2(gel(test_vec, 2), current_I);
    ret = gerepilecopy(av, ret);
    return ret;
}

GEN my_find_a_vect (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN J_vect, int p_int) {
    pari_sp av = avma;
    GEN a, t, ext_gen;
    int p_rk = glength(J_vect);
    GEN class_group = my_get_clgp (Labs);
    GEN a_vect = zerovec(p_rk);
    int i;
    
    for (i = 1; i < p_rk+1; ++i)
    {
        
        ext_gen = rnfidealup0(Lrel, gel(J_vect, i), 1);
        t = gel(my_find_I(Labs, K, sigma, ext_gen, class_group),1);
        a = algtobasis(K, nfinv(K, rnfeltnorm(Lrel, rnfeltabstorel(Lrel, t))));
        gel(a_vect, i) = a;
    } 
    a_vect = gerepilecopy(av, a_vect);
    return a_vect;
}

// Returns vector of tuples (a,J) with div(a)+pJ = 0.
GEN my_find_Ja_vect_modified(GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN J_vect, GEN units, GEN p) {
    pari_sp av = avma;
    int l = glength(J_vect);
    GEN Ja_vect = zerovec(l);
    
    int i;
    GEN a_vect = my_find_a_vect (Labs, Lrel, K, sigma, J_vect, itos(p));
    for (i=1; i<l+1; ++i) {
        
        gel(Ja_vect, i) = mkvec2(gel(a_vect, i), gel(J_vect, i));
    }

    Ja_vect = gerepilecopy(av, Ja_vect);
    return Ja_vect;
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

GEN my_find_I_new (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN a, GEN i_xJ, GEN class_group, GEN p_power_units, int p_int)
{
    pari_sp av = avma;
    GEN current_I;
    GEN test_vec, t;

    GEN class_number = bnf_get_no(Labs);
    int clnr = itos(class_number);
    
    GEN Itest;
    int test_found = 0;
    printf("Getting class group of Labs\n");
    //pari_printf("clgp: %Ps\n", bnf_get_clgp(Labs));
    int n;
    printf("Testing if i(J) principal\n");
    test_vec = bnfisprincipal0(Labs, i_xJ, 1);
    
    for (n = 1; n < clnr + 1; ++n) {
        printf("%d/%d\n", n, clnr);
        current_I = gel(class_group, n);
        Itest = idealdiv(Labs, i_xJ, my_1MS_ideal(Labs, sigma, current_I));
        test_vec = bnfisprincipal0(Labs, Itest, 1);
        
        t = gel(test_vec, 2);
        
        if (my_QV_equal0(gel(test_vec, 1)))
        {
            // This needs to be modified in order to take into account "mod p"
            pari_printf(ANSI_COLOR_GREEN "Principal: %Ps\n"  ANSI_COLOR_RESET, algtobasis(K, nfmul(K,rnfeltnorm(Lrel, rnfeltabstorel(Lrel, t)),a)));
            // if (my_is_p_power(bnfisunit(K, algtobasis(K, nfmul(K,rnfeltnorm(Lrel, rnfeltabstorel(Lrel, t)),a))), p_int)) {
            //     test_found = 1;
            //     break;
            // }
            if (my_is_in_power_units(K, algtobasis(K, nfmul(K,rnfeltnorm(Lrel, rnfeltabstorel(Lrel, t)),a)), p_power_units)) {
                test_found = 1;
                break;
            }
            
        } 

    }
    
    GEN ret;
    if (!test_found)
    {
        printf(ANSI_COLOR_RED "No I found in my_find_I\n\n" ANSI_COLOR_RESET);
        
        ret = mkvec2(gen_0,gen_0);
        ret = gerepilecopy(av, ret);
       
        return ret;
    }
    
    ret = mkvec2(gel(test_vec, 2), current_I);
    ret = gerepilecopy(av, ret);
    return ret;
}

// find I such that J=div(t)+(1-sigma)I and N(t)=-a
GEN my_find_I_v2 (GEN Labs, GEN Lrel, GEN sigma, GEN K, GEN a, GEN iJ, int p_int) {
    GEN t, div_t, I, T, div_t_J_inv;

    T = rnfisnorminit(K, rnf_get_pol(Lrel), 1);
    //tu = bnf_get_tuU(K);

    // Find an element t with norm -a
    t = rnfeltreltoabs(Lrel, gel(rnfisnorm(T, nfinv(K,a), 0),1));
    
    div_t = idealhnf0(Labs, t, NULL);
    div_t_J_inv = idealdiv(Labs, div_t, iJ);
    
    if (my_SQ_MAT_equal(div_t_J_inv, idealhnf0(Labs, gen_1, NULL))) {
        //printf("triv\n");
        I = idealhnf0(Labs, gen_1, NULL);
        
    }
    else {
        I = my_find_H90_ideal(Labs, Lrel, K, div_t_J_inv, sigma, p_int);
        //printf("non-triv\n");
    }
    //pari_printf("I: %Ps\n\n", I);
    return I;
}

// find I such that div(t)+(1-sigma)I = 0 and N(t)=-unit
GEN my_find_final_I (GEN Labs, GEN Lrel, GEN sigma, GEN K, GEN unit, int p_int) {
    pari_sp av = avma;
    GEN t, div_t, I, T;
    //pari_printf("u: %Ps\n", unit);
    T = rnfisnorminit(K, rnf_get_pol(Lrel), 1);
    //tu = bnf_get_tuU(K);
    t = rnfeltreltoabs(Lrel, gel(rnfisnorm(T, nfinv(K, unit), 0),1));
    div_t = idealhnf0(Labs, t, NULL);
    
    if (my_SQ_MAT_equal(div_t, idealhnf0(Labs, gen_1, NULL))) {
        //printf("triv\n");
        I = idealhnf0(Labs, gen_1, NULL);
        
    }
    else {
        I = idealred(Labs, my_find_H90_ideal(Labs, Lrel, K, div_t, sigma, p_int));
        //printf("non-triv\n");
    }
    
    GEN ret = gerepilecopy(av, I);
    return ret;
}

int my_is_in_vect(GEN a, GEN vect) {
    int check = 0;
    int i;
    for (i=1;i<glength(vect)+1;++i) {
        if (my_QV_equal(a, gel(vect, i))) {
            check = 1;
            break;
        }
    }
    return check;
}

// Under construction
GEN my_find_final_I_vect(GEN Labs, GEN Lrel, GEN sigma, GEN K, GEN class_group, GEN p_power_units, GEN units_mod_p) {
    pari_sp av = avma;
    int l = glength(units_mod_p);
    int pow_l = glength(p_power_units);
    int cl_nr = glength(class_group);
    GEN I_final = zerovec(l), test_vec, I, t_inv, a, unit, test_units = zerovec(pow_l);
    int check;
    
    int i;
    int j;
    for (i=1;i<l+1;++i) {
        check = 0;
        unit = algtobasis(K, gel(units_mod_p, i));
        
        for (j=1;j<pow_l+1;++j) {
            gel(test_units, j) = algtobasis(K, nfmul(K, unit, gel(p_power_units, j)));
        }
        
        for (j=1;j<cl_nr+1;++j) {
            //printf("%d/%d\n", j, cl_nr);
            I = gel(class_group, j);
            test_vec = bnfisprincipal0(Labs, my_1MS_ideal(Labs, sigma, I), 1);
            if (my_QV_equal0(gel(test_vec, 1))) {
                
                t_inv = gel(test_vec, 2);
                a = algtobasis(K, rnfeltnorm(Lrel, rnfeltabstorel(Lrel, t_inv)));
                pari_printf(ANSI_COLOR_GREEN "Principal: %Ps\n"  ANSI_COLOR_RESET, bnfisunit(K, a));
                
                if (my_is_in_vect(a, test_units)) {
                    
                    gel(I_final, i) = I;
                    printf(ANSI_COLOR_GREEN "Found I_final[%d]: \n" ANSI_COLOR_RESET, i);
                    check = 1;
                    break;
                }
                else if (my_is_in_vect(algtobasis(K, nfinv(K, a)), test_units)) {
                    gel(I_final, i) = idealred(Labs, idealinv(Labs, I));
                    printf(ANSI_COLOR_GREEN "Found I_final[%d]: \n" ANSI_COLOR_RESET, i);
                    check = 1;
                    break;
                }
            }
        }
        if (!check) {
            printf(ANSI_COLOR_RED "Did not find I_final[%d]: \n" ANSI_COLOR_RESET, i);
        }
    }
    
    GEN ret = gerepilecopy(av, I_final);
    return ret;
}

GEN my_find_I_vect2 (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN J_vect, GEN Ja_vect, GEN units_mod_2, int p_int) {
    
    int p_rk = glength(J_vect);
    int u_length = glength(units_mod_2);
    GEN I_vect = zerovec(p_rk+u_length);
    GEN I;
    GEN ext_gen;
    
    //pari_printf("%Ps\n", gel(ext_gens, 1));
    
    int i;
    
    for (i = 1; i < p_rk+1; ++i)
    {
        ext_gen = rnfidealup0(Lrel, gel(J_vect, i), 1);
        I = my_find_I_v2(Labs, Lrel, sigma, K, gel(gel(Ja_vect, i),1), ext_gen, p_int);
        
        gel(I_vect, i) = I;
    } 
    
    for (i = p_rk+1; i < p_rk+u_length+1; i++)
    {
        gel(I_vect, i) = my_find_final_I(Labs, Lrel, sigma, K, gel(units_mod_2, i-p_rk), p_int);
        
    }
    
    return I_vect;
}

GEN my_find_I_vect (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN Ja_vect, GEN class_group, GEN p_power_units, int p_int) {
    pari_sp av = avma;
    int p_rk = glength(Ja_vect);
    GEN I_vect = zerovec(p_rk);
    GEN I;
    GEN ext_gen;
    //pari_printf("%Ps\n", gel(ext_gens, 1));
    
    // if (!(glength(class_group)==itos(bnf_get_no(Labs)))) {
    //     printf(ANSI_COLOR_RED "Wrong class group (my_find_I_vect)\n\n" ANSI_COLOR_RESET);
    //     pari_close();
    //     exit(0);
    // }
    int i;
    
    for (i = 1; i < p_rk+1; ++i)
    {
        //ext_gen = idealred(Labs, rnfidealup0(Lrel, gel(J_vect, i), 1));
        ext_gen = rnfidealup0(Lrel, gel(gel(Ja_vect, i),2), 1);
        printf("Lifted J: %d/%d\n", i, p_rk);

        I = gel(my_find_I_new(Labs, Lrel, K, sigma, gel(gel(Ja_vect, i),1), ext_gen, class_group, p_power_units, p_int),2);
        pari_printf("I: %Ps\n", I);
        if (gequal0(I)) {
            
            I = my_find_I_v2(Labs, Lrel, sigma, K, gel(gel(Ja_vect, i),1), ext_gen, p_int);
            gel(I_vect, i) = I;
        }
        else {
            gel(I_vect, i) = I;
        }
        
    } 
    I_vect = gerepilecopy(av, I_vect);
    return I_vect;
}

GEN my_find_I_vect3 (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN Ja_vect, GEN units_mod_p, int p_int, GEN p_power_units) {
    pari_sp av = avma;
    int p_rk = glength(Ja_vect);
    int u_length = glength(units_mod_p);
    GEN I_vect = zerovec(p_rk+u_length);
    //pari_printf("%Ps\n", gel(ext_gens, 1));
    GEN class_group = my_get_clgp (Labs);
    GEN I_first = my_find_I_vect(Labs, Lrel, K, sigma, Ja_vect, class_group, p_power_units, p_int);
    printf("First part of I_vect found\n");
    int i;
    
    for (i = 1; i < p_rk+1; ++i)
    {
        gel(I_vect, i) = gel(I_first, i);
    } 
    GEN I_final = my_find_final_I_vect(Labs, Lrel, sigma, K, class_group, p_power_units, units_mod_p);
    for (i = p_rk+1; i < p_rk+u_length+1; ++i)
    {
        if (gequal0(gel(I_final, i-p_rk))) {
            gel(I_vect, i) = my_find_final_I(Labs, Lrel, sigma, K, gel(units_mod_p, i-p_rk), p_int);
            printf("Found I %d/%d\n", i, p_rk+u_length);
        }
        else {
            gel(I_vect, i) = gel(I_final, i-p_rk);
        }
        
    } 
    I_vect = gerepilecopy(av, I_vect);
    return I_vect;
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
        test_vec = bnfisprincipal0(Labs, Itest, 1);

        if (my_QV_equal0(gel(test_vec, 1)))
        {
            t = gel(test_vec, 2);
            t_rel = rnfeltabstorel(Lrel, t);
            norm_t = rnfeltnorm(Lrel, t_rel);
            diff = nfmul(K, norm_t, a);
            for ( i = 1; i < l+1; i++)
            {
                exponents = bnfisunit(K, nfdiv(K, diff, gel(norms, i)));
                if (l%1000==0)
                {
                    pari_printf("exponents[%d]: %Ps\n", i, exponents);
                }
                
                if (my_is_p_power(exponents, p_int, roots_of_unity_nr))
                {
                    printf("I found\n\n");
                    ret = mkvec2(gel(test_vec, 2), current_I);
                    ret = gerepilecopy(av, ret);
                    return ret;
                }
                
            }
            
                
            
        } 

    }
    printf(ANSI_COLOR_RED "my_find_I_full ended with a problem\n" ANSI_COLOR_RESET);
   
    printf(ANSI_COLOR_RED "No I found in my_find_I\n\n" ANSI_COLOR_RESET);
    pari_close();
    exit(0);

    
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

GEN my_H90_vect (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN Ja_vect, GEN p) {
    pari_sp av = avma;
    int r_rk = glength(Ja_vect), f, j, done = 0;
    GEN I_vect = zerovec(r_rk), a, iJ, F, F_ker_T, t, t_rel, Nt, diff, exp;
   
    int i;
    
    for (i = 1; i < r_rk+1; ++i)
    {
        a = gel(gel(Ja_vect, i), 1);
        iJ = rnfidealup0(Lrel, gel(gel(Ja_vect, i),2), 1);
        F = my_H90(Labs, iJ, sigma);
        /*
            Find F_ker_T --------------------
        */
        for (j = 1; j < f+1; j++)
        {
            t = gel(bnfisprincipal0(Labs, idealdiv(Labs, iJ, gel(F_ker_T, j)), 1), 2);
            t_rel = rnfeltabstorel(Lrel, t);
            Nt = rnfeltnorm(Lrel, t_rel);
            diff = nfmul(K, Nt, a);
            exp = bnfisunit(K, diff);
            /* check if exp lies in the image of the operator associated to N: (O_L)^x -> (O_K)^x */

            if (matsolvemod(my_norm_operator(Labs, Lrel, K, p), 0, exp, 0))
            {
                gel(I_vect, i) = gel(F_ker_T, j);
                done = 1;
                break;
            }
            
        }
        if (!done)
        {
            printf(ANSI_COLOR_RED "my_H90_vect ended with a problem \n" ANSI_COLOR_RESET);
            pari_close();
            exit(0);
        }
        
         
        
    } 
    
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

        pari_printf("norm[%d]: %Ps\n", i, bnfisunit(K, norm));

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

        pari_printf("diff[%d]: %Ps\n", i, bnfisunit(Labs, diff));

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
        image = bnfisunit(Labs, gel(units_group, i));

        pari_printf("%Ps\n", image);

    }
    printf("\nsigma:\n\n");
    for (i = 1; i < l+1; i++)
    {
        image = bnfisunit(Labs, galoisapply(Labs, sigma, gel(units_group, i)));

        pari_printf("%Ps\n", image);

    }
    printf("\nsigma^2:\n\n");
    for (i = 1; i < l+1; i++)
    {
        image = bnfisunit(Labs, galoisapply(Labs, sigma, galoisapply(Labs, sigma, gel(units_group, i))));

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