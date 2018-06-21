"""
`function chem2diff(s::String)`

From chemical equations get differential equations in the form needed for `@diff_eq`

Example:
```
s = "1.0,2.0: 2A -> B"
chem2diff(s)
```
"""
function chem2diff(s::String)
    #TODO Jacobian
    reactions = split(strip(s),'\n')
    for r in reactions
        @assert sum([i == ':' for i in r]) == 1 "Missing ':' separator for rates in $(r)"
    end
    species = unique(vcat(get_species.(reactions)...))
    diff_eq = ["dy[$(i)] = " for i in 1:length(species) ]
    info("\n"*join([ "$s <---> dy[$i]"  for (i,s) in enumerate(species) ],'\n'))

    for (i,s) in enumerate(species)
        for r in reactions
            ai, bi = get_stoch(r,s)
            c = bi - ai
            if c != 0.0
                kon, koff = get_rates(r)
                diff_eq[i] *= "$(c)*("
                lhs,rhs = "",""
                for (j,t) in enumerate(species)
                    aj,bj = get_stoch(r,t)
                    if aj == 0.0
                    elseif aj == 1.0
                        lhs *= "*y[$(j)]"
                    else
                        lhs *= "*y[$(j)]^$(aj)"
                    end
                    if bj == 0.0
                    elseif bj == 1.0
                        rhs *= "*y[$(j)]"
                    else
                        rhs *= "*y[$(j)]^$(bj)"
                    end
                end
                if !isempty(lhs)
                    diff_eq[i] *= "$(kon)"*lhs
                else
                    diff_eq[i] *= "0"
                end
                if !isempty(rhs)
                    diff_eq[i] *= "-$(koff)"*rhs
                else
                    diff_eq[i] *= "0"
                end
                diff_eq[i] *= ") + "
            end
        end
        diff_eq[i] = diff_eq[i][1:end-3]*"   #y[$(i)] <--> $(s)"
    end
    join("#".*reactions,'\n')*"\n"*join(diff_eq,'\n')
end

function get_rates(reaction::AbstractString)
    rates = split(reaction,':')[1]
    a = sum([ i == ',' for i in rates ])
    @assert any(a .== [0,1]) "Incorrect number ($(a)) of rates in $(reaction)"
    if a == 0
        warn("Reaction $(reaction) has only one rate. Assuming this is the forward rate...")
        return (parse(Float64,rates),0.0)
    else
        rates = split(rates,',')
        return (parse(Float64,rates[1]),parse(Float64,rates[2]))
    end    
end

function get_lhs(reaction::AbstractString)
    c = split(reaction,':')[2]
    lhs = strip(split(c,"->")[1])
end

function get_rhs(reaction::AbstractString)
    c = split(reaction,':')[2]
    rhs = strip(split(c,"->")[2])
end

function sep_coefficient(s::AbstractString)
    a = '1'; i = 0; coef = ""
    while isnumber(a) || a == '.'
        i += 1
        a = s[i]
        coef *= string(a)
    end
    if isempty(coef[1:end-1])
        return 1.0,strip(s)
    end
    coef = coef[1:end-1]
    c = parse(Float64,coef)
    sp = strip(replace(s,coef,""))
    c,sp
end

function get_species(reaction::AbstractString)
    c = split(reaction,':')[2]
    a = replace(c,"->"," + ")
    s = [ String(sep_coefficient(i)[2]) for i in strip.(split(a,'+')) ]
end

function get_stoch(reaction::AbstractString, species::AbstractString)
    a = get_lhs(reaction)
    st_lhs = 0.0
    for i in strip.(split(a,'+'))
        c,sp = sep_coefficient(i)
        if sp == species
            st_lhs += c
        end
    end
    a = get_rhs(reaction)
    st_rhs = 0.0
    for i in strip.(split(a,'+'))
        c,sp = sep_coefficient(i)
        if sp == species
            st_rhs += c
        end
    end
    return st_lhs,st_rhs
end

