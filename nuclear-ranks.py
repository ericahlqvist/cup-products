from gappy import gap
import json

gap.LoadPackage("anupq")

F = gap.FreeGroup( "a", "b", "c", "d")

a=F.GeneratorsOfGroup()[0]
b=F.GeneratorsOfGroup()[1]
c=F.GeneratorsOfGroup()[2]
d=F.GeneratorsOfGroup()[3]

a_b=gap((a^b)*(a^-1))
a_c=gap((a^c)*(a^-1))
a_d=gap((a^d)*(a^-1))
b_c=gap((b^c)*(b^-1))
b_d=gap((b^d)*(b^-1))
c_d=gap((c^d)*(c^-1))

a2=gap(a^2)
b2=gap(b^2)
c2=gap(c^2)
d2=gap(d^2)

f = open("gap.txt", "a")
file = open('data/discriminants/relators-[2,2,2,2]_4_mod_16.json')
data = json.load(file)

for i in data[:10]:
    
    s=gap(i['relations'])
    f.write(str(s)+"\n")
    rels=gap.EvalString(s)
    G=gap(F/rels)

    hom=gap.EpimorphismPGroup(G,2,2)
    Q_F=gap.Image(hom)
    f.write(str(gap.NuclearRank(Q_F))+"\n")

file.close()
f.close()