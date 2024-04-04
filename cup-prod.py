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


# Replace 'my_program' with the actual name of your compiled C program
executable_path = './main-pol-sta'



with psycopg2.connect("dbname=cup-db") as conn:

    # Open a cursor to perform database operations
    with conn.cursor() as cur:

        # Query the database and obtain data as Python objects.
        cur.execute("SELECT Polynomial, ROW_NUMBER() OVER () AS row_number FROM p_rank_deg_2_4_8 WHERE cup_rk IS NULL;")
        cur.fetchone()
        cur2 = conn.cursor()
        for row in cur:
            #print(record[0])
            pol = row[0]
            print('Pol: ', pol)
            print('Count:', row[1])
            program_arguments = ['2', pol]
            # Construct the command by combining the executable path and arguments
            command = [executable_path]+program_arguments
            
            # Run the C program with arguments and capture its output
            result = subprocess.run(command, stdout=subprocess.PIPE, text=True)
            #print(result.returncode)
            if result.returncode != 111:
                #print("UPDATE p_rank_deg_2_4_4 SET Cup_rk = "+str(result.returncode)+" WHERE Polynomial=\'"+pol+"\';")
                print("cup rank: ", str(result.returncode))
                cur2.execute("UPDATE p_rank_deg_2_4_8 SET Cup_rk = "+str(result.returncode)+" WHERE Polynomial=\'"+pol+"\';")
            else:
                print("--------------CRASH--------------\n")
                cur2.execute("UPDATE p_rank_deg_2_4_8 SET Cup_rk = \'-1000\' WHERE Polynomial=\'"+pol+"\';")
            conn.commit()
        # Make the changes to the database persistent
        

