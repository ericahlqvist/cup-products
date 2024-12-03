# cup-products

This repo contains a program which computes cup products in étale cohomology of the ring of integers of a number field, together with a tiny bit of data. The program is written in C using the library [PARI](http://pari.math.u-bordeaux.fr/) which makes it very fast for number fields of low degree. 

The code is not very well documented and lots of the files and functions are old and no longer used. If you want to use the program you should probably just ask me first how it works.   

To run the main program: ./main-pol p "pol(s)"
where you replace pol(s) by a defining polynomial of the number field K in the variable s, and replace p by a prime number dividing the class number. 