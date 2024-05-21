import psycopg2

with psycopg2.connect("dbname=cup-db") as conn:

    # Open a cursor to perform database operations
    with conn.cursor() as cur:

        # Query the database and obtain data as Python objects.
        cur.execute("ALTER TABLE p_rank_deg_2_4_2 ADD is_inf int;")
        cur.execute("SELECT parent_pol, is_inf FROM p_rank_deg_2_4_4;")
        cur.fetchone()
        cur_2 = conn.cursor()

        # You can use `cur.fetchmany()`, `cur.fetchall()` to return a list
        # of several records, or even iterate on the cursor
        for record in cur:
            pol = record[0]
            is_inf = record[1]
            
            if is_inf == 1:
                cur_2.execute("UPDATE p_rank_deg_2_4_2 SET is_inf = 1 WHERE Polynomial=\'"+pol+"\';")
        # Make the changes to the database persistent
    conn.commit()