module Lsode
#This is a wrap of the DLSODE function from the odepack FORTRAN library.
#I created this mainly to try out the fortran-calling functionality in Julia.
#For solving differential equations in Julia check out the DifferentialEquations package.


#TODO Extend the `IOPT` functionality parameter
#TODO Copy the ITOL behaviour from the odepack.f
#TODO Make macros more flexible -- e.g. other symbol use

export ode,@diff_eq,@diff_eq_jac

include("macros.jl")

const istate_dict = Dict{Int64,String}(
                                       1 =>"Nothing was done, because TOUT was equal to T.",
                                       2 =>"DLSODE was successful (otherwise, negative).",
                                       -1 =>"Excess work done on this call (perhaps wrong MF).",
                                       -2 =>"Excess accuracy requested (tolerances too small).",
                                       -3 =>"Illegal input detected (see printed message).",
                                       -4 =>"Repeated error test failures (check all inputs).",
                                       -5 =>"Repeated convergence failures (perhaps bad Jacobian supplied or wrong choice of MF or tolerances).",
                                       -6 =>"Error weight became zero during problem solution component i vanished, and ATOL or atol[i] = 0.0).")


"""
`dlsode!(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, [verbose=true])`

Call DLSODE from the FORTRAN 'odepack' library. Updates the variables `Y`, `T`, `ISTATE`
"""
function dlsode!(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,ISTATE, IOPT,
                 RWORK, LRW, IWORK, LIW, JAC, MF;verbose=true)

    ccall((:dlsode_, "libodepack"),Void,
          (Ptr{Void},Ref{Int64},Ptr{Array{Float64,1}},Ref{Float64},Ref{Float64},Ref{Int64},
           Ref{Float64},Ptr{Array{Float64,1}},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Float64},
           Ref{Int64},Ref{Int64},Ref{Int64},Ptr{Void},Ref{Int64}),
          F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT,
          RWORK, LRW, IWORK, LIW, JAC, MF)
    verbose ? println("$(istate_dict[ISTATE.x])") : nothing
end

"""
`ode(c_ode_system::Ptr{Void}, [c_jacobian::Ptr{Void}], t_vec::AbstractVector, y0::AbstractVector, atol::AbstractVector=1.0e-3*ones(length(y0)),rtol::Real = 1.0e-3,[verbose=true], [tol_mode = 2])`

ODE solver. The system and the (optional) Jacobian functions should be created with `@diff_eq` and `@diff_eq_jac` macros.
Outputs a vector of solutions at `t_vec` times. 
"""
function ode(c_ode_system::Ptr{Void}, c_jacobian::Ptr{Void},
             t_vec::AbstractVector, y0::AbstractVector,
             atol::AbstractVector=1.0e-3*ones(length(y0)),
             rtol::Real = 1.0e-3;
             verbose=true, tol_mode = 2)

    y = deepcopy(y0)
    t_vec_out = deepcopy(t_vec)
    y_vec = Vector{Vector{Float64}}(length(t_vec))
    n = length(y0)
    @assert length(atol) == n "
        Length of absolute tolerance vector ($(length(atol))) does not match the length of initial conditions vector ($(n))."
    itol = tol_mode
    itask = 1 #Can this be anything else?
    iopt = 0 #TODO: optional arguments
    istate = Ref{Int64}(1)
    mf = 21
    lrw = Int64(22 + 9*n + n^2)
    liw = Int64(20 + n)
    rwork = zeros(Vector{Float64}(lrw))
    iwork = zeros(Vector{Int64}(liw))

    t = Ref{Float64}(t_vec[1])

    for (i,tout) in enumerate(t_vec)
        dlsode!(c_ode_system, n, y, t,
                tout, itol, rtol, atol, itask, istate, iopt,
                rwork, lrw, iwork, liw, c_jacobian, mf,verbose=verbose)
        y_vec[i] = deepcopy(y)
    end

    y_vec
end

function ode(c_ode_system::Ptr{Void},
             t_vec::AbstractVector, y0::AbstractVector,
             atol::AbstractVector=1.0e-3*ones(length(y0)),
             rtol::Real=1.0e-3;
             verbose=true, stiff=true, tol_mode=2)

    y = deepcopy(y0)
    y_vec = Vector{Vector{Float64}}(length(t_vec))
    n = length(y0)
    @assert length(atol) == n
    itol = tol_mode
    itask = 1 #Can this be anything else?
    iopt = 0 #TODO: optional arguments
    istate = Ref{Int64}(1)
    if stiff
        mf = 22
        lrw = 22 + 9*n + n^2
        liw = 20 + n
    else
        mf = 10
        lrw = 20 + 16*n
        liw = 20
    end

    rwork = zeros(Vector{Float64}(lrw))
    iwork = zeros(Vector{Int64}(liw))

    t = Ref{Float64}(t_vec[1])
    function fake_jac_calc(x1,x2,x3,x4,x5,x6,x7)
        return nothing
    end
    fake_jac = cfunction(fake_jac_calc,Void,(Ptr{Int64},Ptr{Float64},Ptr{Float64},
                                             Ptr{Int64},Ptr{Int64},Ptr{Float64},Ptr{Int64}))

    for (i,tout) in enumerate(t_vec)
        dlsode!(c_ode_system, n, y, t,
                tout, itol, rtol, atol, itask, istate, iopt,
                rwork, lrw, iwork, liw, fake_jac, mf,verbose=verbose)
        y_vec[i] = deepcopy(y)
    end

    y_vec
end

end # module
