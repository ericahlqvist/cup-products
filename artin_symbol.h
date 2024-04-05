/*-------------- The Artin Symbol ---------------
Let L/K be a Galois extension and let p be a prime in O_K which is unramified in L (in our case this holds for all primes in O_K since L/K is unramified). Choose a prime q in O_L lying above p. Then there is a unique element sigma_p in G_q:={sigma in Gal(L/K) : sigma(q)=q} such that sigma = (-)^N(p) when acting on the residue field k(q). 
(1) If p is split in L, then G_q = {1} and k(q) = k(p) which means that (-)^N(p) = (-) = id and hence sigma_p = 1 (so 0 when identified with Z/pZ). 
(2) Otherwise, since we assume L/K to be unramified of prime degree, the only other option is that p is inert in L and hence G_q=Gal(L/K). 
The element sigma_p is independent of the choice of q since we assume L/K to be abelian. 

Since k(q)^x is cyclic, it is enough to check that sigma_p(g) = g^N(p) for a generator g of k(q)^x.
*/


GEN my_p_Artin_symbol(GEN Labs, GEN Lrel, GEN K, GEN K_factorization, GEN p, GEN Gal_rel, GEN p_exp) {
    pari_sp av = avma;
    GEN p_Artin_symbol, exp, sigma, y=pol_x(fetch_user_var("y"));;
    
    // Define the prime from factorization
    GEN prime = idealhnf0(K, gel(gel(gel(K_factorization, 1), 1), 1), gel(gel(gel(K_factorization, 1), 1), 2));

    // Lift the prime. This gives a fractional ideal.
    GEN prime_lift = rnfidealup0(Lrel, prime, 1);

    // Factorize the fractional ideal
    GEN factorization = idealfactor(Labs, prime_lift);
    
    // Check if the prime is split
    if (glength(gel(factorization, 1))==itos(p))
    {
        //printf(ANSI_COLOR_GREEN "Split\n\n" ANSI_COLOR_RESET);
        p_Artin_symbol = gen_0;
        return p_Artin_symbol;
    }
    // Otherwise it most be inert
    else {
        // Determine the size of the residue field: N(p)
        exp = idealnorm(K, prime); 
        
        //printf(ANSI_COLOR_YELLOW "Inert\n\n" ANSI_COLOR_RESET);
    }

    // Define the prime from the factorization
    GEN prime_lift_1 = idealhnf0(Labs, gel(gel(gel(factorization, 1), 1), 1), gel(gel(gel(factorization, 1), 1), 2));
    
    // GEN inertia_index = gel(gel(gel(K_factorization, 1), 1), 4); // Always 1
    // pari_printf("Inertia: %Ps\n\n", inertia_index);
    
    // Define a generator for k(q) (Make sure the choose this such that we get a generator for k(q)^x).
    int l = nf_get_degree(bnf_get_nf(Labs));
    
    GEN generator;

    int test = 0, i;
    GEN test_vec, elem1, elem2;

    
    GEN prinit = nfmodprinit(Labs, gel(gel(idealfactor(Labs, prime_lift_1), 1), 1));
    
    // Define the generator
    generator = algtobasis(Labs, y);

    for (i = 1; i < glength(Gal_rel)+1; i++)
    {
        sigma = gel(Gal_rel, i);
        //pari_printf("sigma: %Ps\n\nGenerator: %Ps\n", sigma, generator);
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
        //printf(ANSI_COLOR_YELLOW "Principal\n\n" ANSI_COLOR_RESET);
        return gen_0;
    }
    
    int Artin_symbol = 0;

    // Factorize the fractional ideal into primes
    GEN factorization = idealfactor(K, I_K);
    //pari_printf("factorization: %Ps\n\n", factorization);
    GEN primes_and_es_in_factorization = my_find_primes_in_factorization(K, factorization);
    //pari_printf("primes_and_es: %Ps\n\n", primes_and_es_in_factorization);
    GEN prime_vect = gel(primes_and_es_in_factorization,1);
    //pari_printf("Prime_vect: %Ps\n\n", prime_vect);
    
    GEN e_vect = gel(primes_and_es_in_factorization,2);
    int p_Artin_symbol;
    GEN prime, p_exp;
    
    GEN Gal_rel = my_Gal_rel(Labs, Lrel, K, sigma, p);
    
    int i;
    for (i = 1; i < glength(prime_vect)+1; i++)
    {
        // The prime in "prid" format
        prime = idealfactor(K, gel(prime_vect, i));
        //pari_printf("Prime -> p-Artin: %Ps\n\n", prime);
        p_exp = gel(e_vect, i);
        
        if (itos(gel(e_vect, i))%p == 0 || my_QV_equal0(gel(bnfisprincipal0(K, gel(prime_vect, i), 1),1)))
        {
            
            p_Artin_symbol = 0;
            // printf(ANSI_COLOR_YELLOW "%ld\n" ANSI_COLOR_RESET, itos(gel(e_vect, i)));
            //printf(ANSI_COLOR_YELLOW "Trivial\n\n" ANSI_COLOR_RESET);
            //printf("p_Artin_symbol: %d\n\n", 0);
        }
        else {
            // printf(ANSI_COLOR_YELLOW "%ld\n" ANSI_COLOR_RESET, itos(gel(e_vect, i)));
            //printf(ANSI_COLOR_GREEN "Non-trivial\n\n" ANSI_COLOR_RESET);

            // Compute the Artin symbol for the prime 
            p_Artin_symbol = itos(my_p_Artin_symbol(Labs, Lrel, K, prime, stoi(p), Gal_rel, p_exp));
            //printf("p_Artin_symbol: %d\n\n", p_Artin_symbol);
            
        }
        Artin_symbol = (p_Artin_symbol+Artin_symbol);
        // printf("Artin_symbol: %d\n\n", Artin_symbol);
    }
    
    GEN ret = stoi(smodis(stoi(Artin_symbol), p));
    ret = gerepilecopy(av, ret);
    return ret;
}


