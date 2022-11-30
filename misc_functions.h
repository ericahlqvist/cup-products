

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
    // printf("Hoooj\n\n");
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
} 

GEN my_1MS_ideal (GEN L, GEN sigma, GEN I) 
{
    GEN output = idealmul(L, I, idealinv(L, galoisapply(L, sigma, I)));
    return output;
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
    // printf("HEJ\n\n");
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
    // output(primes_and_es_in_factorization);
    GEN prime_vect = gel(primes_and_es_in_factorization,1);
    GEN primes_under = my_find_primes_under(LyRel, K, prime_vect);
    GEN primes_and_es_by_primes_under = my_sort_primes_and_es(LyRel, K, primes_and_es_in_factorization, primes_under);
    printf("Factorization for H90 done\n\n");
    

    GEN H90_ideal = idealhnf0(LyAbs, gen_1, NULL);
    GEN new_ideal = idealhnf0(LyAbs, gen_1, NULL);
    int i;
    
    for (i = 1; i < glength(primes_under)+1; i++)
    {
        new_ideal = my_find_H90_ideal_single_prime(LyAbs, LyRel, K, gel(gel(primes_and_es_by_primes_under, i), 1), gel(gel(primes_and_es_by_primes_under, i), 2), sigma_y, p);
        // printf("my_find_H90_ideal_single_prime[%d] done", i);
        H90_ideal = idealmul(LyAbs, H90_ideal, new_ideal);
    }
    H90_ideal = gerepilecopy(av, H90_ideal);
    return H90_ideal;
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
        gel(class_group, n) = idealred(K, current_I); 
    }
    
    class_group = gerepilecopy(av, class_group);
    return class_group;
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

GEN my_find_p_gens (GEN K, GEN p)
{
    GEN my_n = bnf_get_cyc(K);
    int l = glength(my_n);
    GEN p_gens = zerovec(l);
    GEN current_gen;
    
    // int l = glength(my_n);
    GEN my_gens = bnf_get_gen(K);
    pari_printf("my_gens: %Ps\n\n", my_gens);
    int i;
    for (i = 1; i < l+1; ++i)
    {
        // pari_printf("my_n, my_n/p, gen^my_n: %Ps, %Ps, %Ps\n\n", gel(my_n, i), gdiv(gel(my_n, i),p), gel(bnfisprincipal0(K, idealpow(K, gel(my_gens, i), gel(my_n, i)), 1), 1));
        current_gen = idealpow(K,gel(my_gens, i), gdiv(gel(my_n, i),p));
        //pari_printf("current_gen: %Ps\n\n", current_gen);
        gel(p_gens, i) = idealred(K, current_gen);
    }
    return p_gens;
}

/*
Given J \in CL(K)_p, find t, I such that J=div(t)+(1-sigma)I 
*/
GEN my_find_I (GEN Labs, GEN K, GEN sigma, GEN i_xJ)
{
    pari_sp av = avma;
    GEN current_I;
    GEN test_vec;

    GEN class_number = bnf_get_no(Labs);
    int clnr = itos(class_number);
    
    GEN Itest;
    int test_found = 0;
    GEN class_group = my_get_clgp (Labs);
    int n;
    
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
            current_I = gel(class_group, n);
            Itest = idealmul(Labs, i_xJ, my_1MS_ideal(Labs, sigma, current_I));
            test_vec = bnfisprincipal0(Labs, Itest, 1);
            
            if (my_QV_equal0(gel(test_vec, 1)))
            {
                printf(ANSI_COLOR_GREEN "FOUND!\n\n" ANSI_COLOR_RESET);
                
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
    
    GEN ret = mkvec2(gel(test_vec, 2), idealinv(Labs, current_I));
    ret = gerepilecopy(av, ret);
    return ret;
}

GEN my_find_I_vect (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN J_vect) {
    int nr_comp = glength(J_vect);
    GEN I_vect = zerovec(nr_comp);
    GEN I;
    int class_number = itos(bnf_get_no(Labs));
    GEN ext_gen;
    //pari_printf("%Ps\n", gel(ext_gens, 1));
    
    int i;
    
    for (i = 1; i < nr_comp+1; ++i)
    {
        ext_gen = rnfidealup0(Lrel, gel(J_vect, i), 1);
        if (class_number==0) {
            I = idealhnf(Labs, gen_1);
        }
        else {
            I = my_find_I(Labs, K, sigma, ext_gen);
        }
        
        gel(I_vect, i) = gel(I, 2);
    } 
    return I_vect;
}