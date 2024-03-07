import psycopg2
import subprocess, threading

class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            
            self.process = subprocess.Popen(self.cmd, shell=True)
            self.process.communicate()
            #print('Thread finished')

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            #print('Terminating process')
            self.process.terminate()
            thread.join()
            return self.process.returncode
# Connect to an existing database
with psycopg2.connect("dbname=cup-db") as conn:

    # Open a cursor to perform database operations
    with conn.cursor() as cur:

        # Query the database and obtain data as Python objects.
        cur.execute("SELECT Polynomial FROM p_rank_deg_2_4_8 LIMIT 1000")
        cur.fetchone()
        

        # You can use `cur.fetchmany()`, `cur.fetchall()` to return a list
        # of several records, or even iterate on the cursor
        for record in cur:
            #print(record[0])
            p = "2" # A prime 
            pol = record[0]
            #print("./main-script-sta "+p+" "+pol.replace(" ", ""))
            command = Command("./main-script-sta "+p+" "+pol.replace(" ", ""))
            code = command.run(timeout=300)

        # Make the changes to the database persistent
        conn.commit()