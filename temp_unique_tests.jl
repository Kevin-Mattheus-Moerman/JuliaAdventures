using StaticArrays

function unique_ind(x) 
    T = eltype(x)    
    xUni = Vector{T}()        
    indUni = Vector{Int64}()
    indRec = zeros(Int64,length(x))
    for i ∈ eachindex(x) #Loop over entries                
        ind=findfirst(isequal(x[i]),xUni)
        if isnothing(ind) #If not found yet               
            push!(xUni, x[i]) #Add current value
            push!(indUni, i) #Add current index
            indRec[i]=length(xUni)
        else            
            indRec[i]=ind
        end
    end
    return xUni, indUni, indRec        
end

function cunique_ind(x)         
    xUni = Vector{eltype(x)}() #Initialize unique component vector        
    indUni = Vector{Int64}() #Initialize unique index vector
    indRec = zeros(Int64,length(x)) #Allocate reconstruction index vector
    c = Vector{Int64}() #Initialize count vector
    for i ∈ eachindex(x) #Loop over entries                
        ind=findfirst(isequal(x[i]),xUni)
        if isnothing(ind) #New entry
            push!(xUni, x[i]) #Add current entry 
            push!(indUni, i) #Add current index
            indRec[i]=length(xUni) #Add to reconstruction index vector (last index)
            push!(c, 1) #Add entry to count vector, start at 1
        else #Entry was encountered before           
            indRec[i]=ind #Keep index pointing to previously found entry
            c[ind]+=1 #Increment counter for previously found entry
        end
    end
    return xUni, indUni, indRec, c        
end

function cunique(x)     
    xUni = Vector{eltype(x)}()        
    c = Vector{Int64}()
    for i ∈ eachindex(x) #Loop over entries                
        ind=findfirst(isequal(x[i]),xUni)
        if isnothing(ind) #If not found yet               
            push!(xUni, x[i]) #Add current value
            push!(c, 1) #Add current count            
        else            
            c[ind]+=1
        end
    end
    return xUni, c        
end

function funique(itr)     
    T = eltype(itr)    
    seen = Set{T}()
    for x in itr
        if !in(x, seen)
            push!(seen, x)            
        end
    end
    return convert.(T,seen)
end

function getUniqueOccurance(F)     
    Fs=sort(F,dims=2)              
    d = Dict{Int64,Bool}()    
    for i=1:1:size(Fs,1)
        d[i]=true
    end
    for i=1:1:size(Fs,1)       
        # now remove if found elsewhere
        for k in d #Candidate indices other than i
            if i!=k.first && isequal(Fs[i,:],Fs[j,:]) #If we found another               
                #Remove both i and j from list  
                delete!(d,i)
                break             
            end
        end            
    end
    return collect(keys(d))
end

function isNonSharedFace1(F)
    # NB! works only when faces are shared an even number of times!
    Ft = sort.(SVector{4}.(eachrow(F)))
    d = Dict{eltype(Ft),Int64}()
    for i ∈ eachindex(Ft)
        if haskey(d, Ft[i])
            delete!(d, Ft[i])                    
        else
            d[Ft[i]] = i
        end
    end
    return collect(values(d))
end

function isNonSharedFace2(F)
    # NB! works only when faces are shared an even number of times!   
    d = Dict{Vector{Int64},Int64}()
    for i=1:size(F,1)
        fs=sort(F[i,:])
        if haskey(d, fs)
            delete!(d, fs)                    
        else
            d[fs] = i
        end
    end
    return collect(values(d))
end

n=10
m=100000*n
# x=rand(collect(1:n),m)
x=[3 2 2 1 1 2 3 3 5 0 100 3 3 3 3]
# a=unique(x)

#= F=[1 2 3 4; 
   4 5 6 7;
   4 3 2 1;
   5 6 7 8;
   8 7 6 5;
   9 10 11 12]
F=[repeat(rand(1:50,50,4),10,1); rand(1:50,100000,4)].*100 =#

# unique(F)