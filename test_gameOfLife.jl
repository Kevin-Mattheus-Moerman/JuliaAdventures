using GLMakie

GLMakie.activate!(inline=false) # To avoid plotting in plotpane as per: https://github.com/MakieOrg/Makie.jl/issues/2956

np=75

nMax=500

A=zeros(Int64,np,np)

function randInd(np)
    r1=round.(Int64,(np-1).*rand(1)).+1
    r2=round.(Int64,(np-1).*rand(1)).+1
    return r1[1],r2[1]
end

function makeKiteRightUp(A,iStart,jStart)
    if (iStart+2)<size(A,1) && (jStart+1)<size(A,2) && (jStart-1)>1
    A[iStart,jStart]=1
    A[iStart+1,jStart+1]=1
    A[iStart+2,jStart+1]=1
    A[iStart+2,jStart]=1
    A[iStart+2,jStart-1]=1 
    end   
    return A
end

function makeKiteLeftDown(A,iStart,jStart)
    if (iStart+2)<size(A,1) && (jStart+2)<size(A,2)
    A[iStart,jStart]=1
    A[iStart+1,jStart]=1
    A[iStart+2,jStart]=1
    A[iStart+2,jStart+1]=1
    A[iStart+1,jStart+2]=1    
    end
    return A
end

function makeKiteLeftDown(A,iStart,jStart)
    if (iStart+2)<size(A,1) && (jStart+2)<size(A,2)
    A[iStart,jStart]=1
    A[iStart+1,jStart]=1
    A[iStart+2,jStart]=1
    A[iStart+2,jStart+1]=1
    A[iStart+1,jStart+2]=1    
    end
    return A
end

function makePantsLeft(A,i,j)
    if (i+2)<size(A,1) && (j+2)<size(A,2)
    A[i,j]=1   #Origin
    A[i+1,j]=1
    A[i+2,j]=1
    A[i+2,j+1]=1
    A[i+2,j+2]=1
    A[i+1,j+2]=1
    A[i,j+2]=1
    end
    return A
end



# A=makeKiteRightUp(A,3,3)
# A=makeKiteLeftDown(A,10,10)

function randMakeKiteLeftDown(B,np,n)
    for q=1:1:n
        r1,r2=randInd(np)
        B=makeKiteLeftDown(B,r1,r2)
    end
    return B
end

function randmakePantsLeft(B,np,n)
for q=1:1:n
    r1,r2=randInd(np)
    B=makePantsLeft(B,r1,r2)
end
return B
end


A=randMakeKiteLeftDown(A,np,10)

function testShape(A,i,j)
    # T-shape
    # A[i,j]=1   #Origin
    # A[i+1,j]=1
    # A[i+1,j-1]=1
    # A[i+1,j+1]=1

    # #The K
    # A[i  ,j  ]=1
    # A[i-1,j-1]=1
    # A[i-1,j+1]=1
    # A[i+1,j-1]=1
    # A[i+1,j+1]=1
    # A[i-1,j  ]=1

    # #The H
    # A[i  ,j  ]=1
    # A[i-1,j-1]=1
    # A[i-1,j+1]=1
    # A[i+1,j-1]=1
    # A[i+1,j+1]=1
    # A[i-1,j  ]=1
    # A[i+1,j  ]=1

    # #The 0
    # A[i  ,j  ]=1
    # A[i+1,j  ]=1
    # A[i+2,j  ]=1
    # A[i+1,j+1]=1
    # A[i+2,j+1]=1
    # A[i  ,j+2]=1
    # A[i+1,j+2]=1
    # A[i+2,j+2]=1    

    # Love heart
    A[i,j]=1
    A[i-1,j]=1
    A[i-1,j-1]=1
    A[i-1,j+1]=1
    A[i,j+4]=1
    A[i-1,j+4]=1
    A[i-1,j+3]=1
    A[i-1,j+5]=1
    A[i-2,j]=1
    A[i-2,j+1]=1
    A[i-2,j+2]=1
    A[i-2,j+3]=1
    A[i-2,j+4]=1
    A[i-3,j+1]=1
    A[i-3,j+2]=1
    A[i-3,j+3]=1
    A[i-4,j+2]=1

    return A
end

im=round(Int64,np/2)
jm=round(Int64,np/2)

A=testShape(A,im,jm)

function gridShape(A,s)
for i=1:s:size(A,1)
    for j=1:s:size(A,2)
        A=makeKiteLeftDown(A,i,j)
    end
end
return A
end

#A=gridShape(A,8)

# A=randmakePantsLeft(A,np,25)
# A=randMakeKiteLeftDown(A,np,150)

function gameOfLifeIter(A)    
    B=zeros(Int64,size(A))

    for i ∈ 1:size(B,1)
        for j ∈ 1:size(B,2)                     
            
            if i<size(A,1)
                nu=A[i+1,j]
            else
                nu=0
            end
            
            if i>1
                nb=A[i-1,j]
            else
                nb=0
            end

            if j<size(A,2)
                nr=A[i,j+1]
            else
                nr=0
            end
            
            if j>1
                nl=A[i,j-1]
            else
                nl=0
            end

            if i<size(A,1) && j>1
                nul=A[i+1,j-1]
            else
                nul=0
            end

            if i<size(A,1) && j<size(A,2)
                nur=A[i+1,j+1]
            else
                nur=0
            end

            if i>1 && j>1
                nbl=A[i-1,j-1]
            else
                nbl=0
            end

            if i>1 && j<size(A,2)
                nbr=A[i-1,j+1]
            else
                nbr=0
            end

            s=nu+nb+nr+nl+nul+nur+nbl+nbr

            #Die, under population
            if A[i,j]==1 && s<2 
                B[i,j]=0
            end

            #Stay alive
            if A[i,j]==1 && s==2 || s==3
                B[i,j]=1             
            end            

            #Die, over population
            if A[i,j]==1 && s>3
                B[i,j]=0
            end

            #Come alive
            if A[i,j]==0 && s==3
                B[i,j]=1
            end

        end
    end
    return B
end

# show(stdout, "text/plain", A)


function gameOfLife(A,n)
    for q ∈ 1:n
        A=gameOfLifeIter(A)
    end
    return A
end


fig = Figure(resolution=(800,800))
sl_step = Slider(fig[2, 1], range = 0:1:nMax, startvalue = 0)


titleString = lift(sl_step.value) do stepIndex
    "Step: "*string(stepIndex)
end

M = lift(sl_step.value) do stepIndex
    gameOfLife(collect(A'),stepIndex)
end

ax  = Axis(fig[1,1][1,1],aspect = DataAspect(),title = titleString,yreversed = true)

heatmap!(ax, M)
Colorbar(fig[1,1][1,2])
fig