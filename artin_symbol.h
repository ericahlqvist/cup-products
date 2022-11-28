GEN my_p_Artin_symbol(GEN Labs, GEN Lrel, GEN K, GEN K_factorization, GEN p, GEN Gal_rel, GEN p_exp) {
    pari_sp av = avma;
    GEN p_Artin_symbol;
    // Find prime below to determine size of residue field
    GEN prime_below = gel(gel(gel(K_factorization, 1), 1), 1);
    // Define the prime from factorization
    GEN prime = idealhnf0(K, gel(gel(gel(K_factorization, 1), 1), 1), gel(gel(gel(K_factorization, 1), 1), 2));
    // Lift the prime
    GEN prime_lift = rnfidealup0(Lrel, prime, 1);
    GEN factorization = idealfactor(Labs, prime_lift);
    GEN exp;

    printf("Fact of lift: ");
    output(factorization);
    printf("\n\n");
    if (glength(gel(factorization, 1))==itos(p))
    {
        printf(ANSI_COLOR_GREEN "Split\n\n" ANSI_COLOR_RESET);
        p_Artin_symbol = gen_0;
        return p_Artin_symbol;
    }
    else {
        exp = prime_below;
        pari_printf("Norm p: %Ps\n\nExp: %Ps\n\n", idealnorm(K, prime), exp);
        printf(ANSI_COLOR_YELLOW "Inert\n\n" ANSI_COLOR_RESET);
    }
     //idealnorm(K, prime);
    pari_printf("Norm p: %Ps\n\n", exp);

    // pari_printf("K_factorization: %Ps\n\n", K_factorization);
    // pari_printf("Factorization: %Ps\n\n", factorization);
    GEN prime_lift_1 = idealhnf0(Labs, gel(gel(gel(factorization, 1), 1), 1), gel(gel(gel(factorization, 1), 1), 2));
    
    GEN inertia_index = gel(gel(gel(K_factorization, 1), 1), 4); // Always 1
    pari_printf("Inertia: %Ps\n\n", inertia_index);
    
    GEN sigma;
    
    int l = 2*itos(p);
    GEN generator = zerocol(l);
    gel(generator, 2) = gen_1;
    int test = 0;
    GEN test_vec;
    int i;
    GEN elem1, elem2;
    
    GEN prinit = nfmodprinit(Labs, gel(gel(idealfactor(Labs, prime_lift_1), 1), 1));
    // printf("prinit: ");
    // output(prinit);
    // printf("\n");
    // printf("mod prime: ");
    // output(gel(gel(idealfactor(Labs, prime_lift_1), 1), 1));
    // printf("\n");

    // pari_printf("exp: %Ps\n\n", exp);

    for (i = 1; i < glength(Gal_rel)+1; i++)
    {
        sigma = gel(Gal_rel, i);
        
        elem1 = galoisapply(Labs, sigma, generator);
        elem2 = nfpow(Labs, generator, exp);
        //prinit = gel(gel(idealfactor(Labs, prime_lift_1), 1), 1);
        //nf_to_Fq_init
        // Page 297 in User's guide to the pari lib
        
        
        test_vec = nfsub(Labs,elem1,elem2);
        // printf("test_vec before quot: ");
        // output(test_vec);
        // printf("\n\n");

        test_vec = nf_to_Fq(bnf_get_nf(Labs), test_vec, prinit);
       
        // pari_printf("elem1: %Ps\nelem2: %Ps\n\n", elem1, elem2);
        // printf("test_vec: ");
        // output(test_vec);
        // printf("\n\n");

        if (gequal0(test_vec))
        {
            // Gal_rel ordered as [sigma, sigma^2, ..., sigma^p=id]
            p_Artin_symbol = stoi(itos(p_exp)*i);
            //printf("p_exp: %ld, i: %d\n\n", itos(p_exp), i);
            test = 1;
            break;
        } 
    }
    
    if (!test)
    {
        printf(ANSI_COLOR_RED "ERROR - no galois elem for Artin \n\n" ANSI_COLOR_RESET);
        pari_close();
        exit(0);
    }
    

    return gerepilecopy(av, p_Artin_symbol);
}

GEN my_Artin_symbol (GEN Labs, GEN Lrel, GEN K, GEN I_K, int p, GEN sigma) {
    pari_sp av = avma;
    if (my_QV_equal0(gel(bnfisprincipal0(K, I_K, 1),1)))
    {
        printf(ANSI_COLOR_YELLOW "Principal\n\n" ANSI_COLOR_RESET);
        return gen_0;
    }
    
    int Artin_symbol = 0;
    GEN factorization = idealfactor(K, I_K);
    pari_printf("factorization: %Ps\n\n", factorization);
    GEN primes_and_es_in_factorization = my_find_primes_in_factorization(K, factorization);
    pari_printf("primes_and_es: %Ps\n\n", primes_and_es_in_factorization);
    GEN prime_vect = gel(primes_and_es_in_factorization,1);
    pari_printf("Prime_vect: %Ps\n\n", prime_vect);
    
    GEN e_vect = gel(primes_and_es_in_factorization,2);
    int p_Artin_symbol;
    GEN prime, p_exp;
    GEN Gal_rel = my_Gal_rel(Labs, Lrel, K, sigma, p);
    
    int i;
    for (i = 1; i < glength(prime_vect)+1; i++)
    {
        prime = idealfactor(K, gel(prime_vect, i));
        pari_printf("Prime -> p-Artin: %Ps\n\n", prime);
        p_exp = gel(e_vect, i);
        
        if (itos(gel(e_vect, i))%p == 0 || my_QV_equal0(gel(bnfisprincipal0(K, gel(prime_vect, i), 1),1)))
        {
            
            p_Artin_symbol = 0;
            // printf(ANSI_COLOR_YELLOW "%ld\n" ANSI_COLOR_RESET, itos(gel(e_vect, i)));
            printf(ANSI_COLOR_YELLOW "Trivial\n\n" ANSI_COLOR_RESET);
            printf("p_Artin_symbol: %d\n\n", 0);
        }
        else {
            // printf(ANSI_COLOR_YELLOW "%ld\n" ANSI_COLOR_RESET, itos(gel(e_vect, i)));
            printf(ANSI_COLOR_GREEN "Non-trivial\n\n" ANSI_COLOR_RESET);
            p_Artin_symbol = itos(my_p_Artin_symbol(Labs, Lrel, K, prime, stoi(p), Gal_rel, p_exp));
            printf("p_Artin_symbol: %d\n\n", p_Artin_symbol);
            
        }
        Artin_symbol = (p_Artin_symbol+Artin_symbol);
        // printf("Artin_symbol: %d\n\n", Artin_symbol);
    }
    
    GEN ret = stoi(smodis(stoi(Artin_symbol), p));
    ret = gerepilecopy(av, ret);
    return ret;
}


