import psycopg2
from gappy import gap
import json

gap.LoadPackage("anupq")

gap.eval("F:=FreeGroup( \"a\", \"b\", \"c\", \"d\");")


gap.eval("a:=F.1;")
gap.eval("b:=F.2;")
gap.eval("c:=F.3;")
gap.eval("d:=F.4;")

gap.eval("a_b:=(a^-1)*(b^-1)*a*b;")
gap.eval("a_c:=(a^-1)*(c^-1)*a*c;")
gap.eval("a_d:=(a^-1)*(d^-1)*a*d;")
gap.eval("b_c:=(b^-1)*(c^-1)*b*c;")
gap.eval("b_d:=(b^-1)*(d^-1)*b*d;")
gap.eval("c_d:=(c^-1)*(d^-1)*c*d;")

gap.eval("a2:=a^2;")
gap.eval("b2:=b^2;")
gap.eval("c2:=c^2;")
gap.eval("d2:=d^2;")

with psycopg2.connect("dbname=cup-db") as conn:

    # Open a cursor to perform database operations
    with conn.cursor() as cur:

        # Query the database and obtain data as Python objects.
        cur.execute("ALTER TABLE p_rank_deg_2_4_2 ADD size_automorphisms int;")
        cur.execute("SELECT polynomial, relations FROM p_rank_deg_2_4_2;")
        cur.fetchone()
        cur_2 = conn.cursor()

        # You can use `cur.fetchmany()`, `cur.fetchall()` to return a list
        # of several records, or even iterate on the cursor
        for record in cur:
            pol = record[0]
            relations = record[1]
            #print(relations)
            arg="rels:="+relations+";"
            #print(arg)
            rels = gap.eval(arg)
            #print(rels)
            G = gap.eval("G:=F/rels;")
            hom = gap.eval("hom:=EpimorphismPGroup(G,2,2);")
            Q_F = gap.eval("Q_F:=Image(hom);")
            size = gap.eval("Size(AutomorphismGroup(Q_F));")
            print("pol: ", pol, "      SIZE: ", size)
            cur_2.execute("UPDATE p_rank_deg_2_4_2 SET size_automorphisms = "+str(size)+" WHERE Polynomial=\'"+pol+"\';")
        # Make the changes to the database persistent
    conn.commit()

    

    
    

