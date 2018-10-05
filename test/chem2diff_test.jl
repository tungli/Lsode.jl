kon = rand()
koff = rand()

c1 = chem2diff("""
               $(kon),$(koff): 2 A -> B
          """)
c2 = chem2diff("""
               $(kon),$(koff): 2A -> B
               """)

c3 = chem2diff("""
               $(kon),$(koff):A +A -> B
               """)

n = 2
z = 1.0*ones(n)
dz = zeros(n)

dz[1] = 2*(-kon*z[1]^2 + koff*z[2])
dz[2] = kon*z[1]^2 - koff*z[2]

#This works only because of the fixed usage of `y`,`dy`
y = ones(n) 
dy = zeros(n)

eval(c1)
@test all(dy .≈ dz)

eval(c2)
@test all(dy .≈ dz)

eval(c3)
@test all(dy .≈ dz)

