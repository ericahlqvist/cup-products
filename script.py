

import subprocess, threading, json

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


p = "5" # A prime 

with open('output/5_4mod16.json', 'r') as file:
    data = json.load(file)

for item in data[:1000]:
    pol = "s^2+"+str(item['D'])[1:]
    command = Command("./main-script-sta "+p+" "+pol)
    code = command.run(timeout=300)
    ##print(code.stdout)
    # if code == -15:
    #     res_file = open("output/"+p+"_"+Dmod+"mod"+mod+".txt", "a")
    #     print(my_str)
    #     res_file.write("{\"p\": \""+p+"\", \"D\": \""+my_str.strip()+"\", \"Z-rk\": \"-\", \"K-cyc\": \"-\", \"Lx-cyc\": \"-\", \"Ly-cyc\": \"-\", \"ZM\": \"-\"},\n")
    
file.close()
