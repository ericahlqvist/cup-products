

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


p = "2" # A prime 
open_file = "p_2_[2,2,2,2]_disc_8_mod_16.txt" # A file

# Dmod = "3" # 3, 7, 4, 8
# mod = ""
# if (Dmod == "3" or Dmod == "7"):
#     mod = "8"
# else:
#     mod = "16"

# open_file = ""

# if p == "2":
#     open_file = "p_"+p+"_cyc_2_2_disc_"+Dmod+"_mod_"+mod+".txt"
# else:
#     open_file = "p_"+p+"_disc_"+Dmod+"_mod_"+mod+".txt"

file = open("data/discriminants/"+open_file)
lines = file.readlines()

for line in lines[0:100]:
    my_str = ''.join(map(str, line))
    pol = "s^2+"+my_str[1:]
    command = Command("./main-script-sta "+p+" "+pol)
    code = command.run(timeout=300)
    
    # if code == -15:
    #     res_file = open("output/"+p+"_"+Dmod+"mod"+mod+".txt", "a")
    #     print(my_str)
    #     res_file.write("{\"p\": \""+p+"\", \"D\": \""+my_str.strip()+"\", \"Z-rk\": \"-\", \"K-cyc\": \"-\", \"Lx-cyc\": \"-\", \"Ly-cyc\": \"-\", \"ZM\": \"-\"},\n")
    
file.close()
