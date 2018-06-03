#Macros for ODE system creation

"""
`@diff_eq f_name sys_eq`

Creates a function named `f_name` which sets the derivatives vector of an ODE system.

For now the user is limited to using only the `dy`, `y` and `t` symbols.

Example:
```
@diff_eq F begin
dy[1] = -y[1]
end
```
"""
macro diff_eq(f_name,sys_eq)
    a = quote
        $(f_name) = cfunction((pn::Ref{Int64}, pt::Ref{Float64},py::Ptr{Float64}, pdy::Ptr{Float64}) -> begin
                                  n = unsafe_load(pn)
                                  t = unsafe_load(pt)
                                  y = unsafe_wrap(Array, py, n)
                                  dy = unsafe_wrap(Array, pdy, n)

                                  $(sys_eq)
                                  return nothing
                              end,
                              Void,(Ptr{Int64},Ptr{Float64},Ptr{Float64},Ptr{Float64}))
    end
    esc(a)
end


"""
`@diff_eq_jac jac_name jac_eqs`

Creates a function named `jac_name` which sets the values of the Jacobian matrix associated with an ODE system.

For now the user is limited to using only the `pd` (Jacobian matrix), `y` and `t` symbols.

Example:
```
@diff_eq_jac J begin
pd[1,1] = -1.0
end
```
"""
macro diff_eq_jac(jac_name,jac_eqs)
    #pd(i,j) = df(i)/dy(j)
    a = quote
        $(jac_name) = cfunction((pn::Ref{Int64}, pt::Ref{Float64},py::Ptr{Float64}, 
                                 pml::Ptr{Int64}, pmu::Ptr{Int64},ppd::Ptr{Float64},
                                 pnrpd::Ptr{Int64}) -> begin
                                n = unsafe_load(pn)
                                t = unsafe_load(pt)
                                y = unsafe_wrap(Array, py, n)
                                pd = unsafe_wrap(Array, ppd, (n,n))
                                ml = unsafe_load(pml)
                                mu = unsafe_load(pmu)
                                nrpd = unsafe_load(pnrpd)

                                $(jac_eqs)
                                return nothing

                            end
                            ,Void,(Ptr{Int64},Ptr{Float64},Ptr{Float64},
                                   Ptr{Int64},Ptr{Int64},Ptr{Float64},Ptr{Int64}))

    end
    esc(a)
end



