#These test/demonstrate the basic funcionality

@diff_eq F begin
    dy[1] = -y[1]
end
@diff_eq_jac J begin
    pd[1,1] = -1.0
end

#From the book Numerical recipes
@diff_eq F2 begin
    dy[1] = 998*y[1] + 1998*y[2]  
    dy[2] = -999*y[1] - 1999*y[2]
end
@diff_eq_jac J2 begin
    pd[1,1] = 998
    pd[1,2] = 1998
    pd[2,1] = -999
    pd[2,2] = -1999
end

y0 = [1.0]
t = linspace(0,1,10)

y02 = [1.0,0.0]

s_exp = exp.(-t)

@show s_exp2 = hcat( [2*exp.(-t) - exp.(-1000*t),
                -exp.(-t) + exp.(-1000*t)]... )

tol = 1.0e-2

s1 = squeeze(hcat(ode(F,t,y0,verbose=false)...)',2)   #jacobian generated, stiff
@test sqrt(sum(abs2, s1 - s_exp)) < tol
s2 = squeeze(hcat(ode(F,J,t,y0,verbose=false)...)',2)  #jacobian supplied
@test sqrt(sum(abs2, s2 - s_exp)) < tol
s3 = squeeze(hcat(ode(F,t,y0,stiff=false,verbose=false)...)',2)  #jacobian ignored, not stiff
@test sqrt(sum(abs2, s3 - s_exp)) < tol
s4 = squeeze(hcat(ode(F,J,t,y0,[1.0e-4],verbose=false)...)',2)  #increased absolute tolerance
@test sqrt(sum(abs2, s2 - s_exp)) > sqrt(sum(abs2, s4 - s_exp))
s5 = squeeze(hcat(ode(F,J,t,y0,[1.0e-4],1.0e-4,verbose=false)...)',2)  #increased relative tolerance
@test sqrt(sum(abs2, s5 - s_exp)) < sqrt(sum(abs2, s4 - s_exp))

#Stiffs
s6 = hcat(ode(F2,t,y02,verbose=false)...)'   #jacobian generated, stiff
@test sqrt(sum(abs2, s6 - s_exp2)) < tol
s7 = hcat(ode(F2,J2,t,y02,verbose=false)...)'  #jacobian supplied
@test sqrt(sum(abs2, s7 - s_exp2)) < tol
s8 = hcat(ode(F2,t,y02,stiff=false,verbose=false)...)'  #jacobian ignored, not stiff
@test sqrt(sum(abs2, s8 - s_exp2)) < tol



