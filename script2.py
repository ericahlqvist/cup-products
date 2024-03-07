

import subprocess
import psycopg2

# Connect to an existing database
executable_path = './main-pol-sta'
p = "2" # A prime 
with psycopg2.connect("dbname=cup-db") as conn:

    # Open a cursor to perform database operations
    with conn.cursor() as cur:

        # Query the database and obtain data as Python objects.
        cur.execute("SELECT Polynomial FROM p_rank_deg_2_4_4")
        cur.fetchone()
        
        for record in cur:
            #print(record[0])
            pol = record[0]
            
            program_arguments = [p, pol]
            # Construct the command by combining the executable path and arguments
            command = [executable_path]+program_arguments
            
            # Run the C program with arguments and capture its output
            result = subprocess.run(command, stdout=subprocess.PIPE, text=True)
            if result.returncode !=0:
                cur.execute("UPDATE p_rank_deg_2_4_4 SET Cup_rk = "+str(result.returncode)+" WHERE Polynomial=\'"+pol+"\';")

        # Make the changes to the database persistent
        conn.commit()








