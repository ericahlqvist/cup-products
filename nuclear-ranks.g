
F := FreeGroup( "a", "b", "c", "d" );;

a:=F.1;;
b:=F.2;;
c:=F.3;;
d:=F.4;;

a_b:=a^b*a^-1;;
a_c:=a^c*a^-1;;
a_d:=a^d*a^-1;;
b_c:=b^c*b^-1;;
b_d:=b^d*b^-1;;
c_d:=c^d*c^-1;;

json_file := InputTextFile("data/discriminants/relators-[2,2,2,2]_4_mod_16.json");;
IsStream(json_file);
r := JsonStreamToGap( json_file );;

r[1];