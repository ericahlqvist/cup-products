import psycopg2

with psycopg2.connect("dbname=cup-db") as conn:

    # Open a cursor to perform database operations
    with conn.cursor() as cur:

        # Query the database and obtain data as Python objects.
        cur.execute("ALTER TABLE p_rank_deg_2_4_16 ADD p_rank int;")
        cur.execute("SELECT Polynomial, cyc FROM p_rank_deg_2_4_16;")
        
        cur.fetchone()
        cur_2 = conn.cursor()

        # You can use `cur.fetchmany()`, `cur.fetchall()` to return a list
        # of several records, or even iterate on the cursor
        for record in cur:
            pol = record[0]
            cyc = record[1]
            l = cyc.split(',')
            p_rank = len(l)
            
            cur_2.execute("UPDATE p_rank_deg_2_4_16 SET p_rank = "+str(p_rank)+" WHERE Polynomial=\'"+pol+"\';")
        # Make the changes to the database persistent
    conn.commit()