using Printf

function grid3D(x,y,z)

    X = [ i for i ∈ x          , j ∈ 1:length(y), k ∈ 1:length(z) ]
    Y = [ j for i ∈ 1:length(x), j ∈ y          , k ∈ 1:length(z) ]
    Z = [ k for i ∈ 1:length(x), j ∈ 1:length(y), k ∈ z           ]

    return X, Y, Z
end

function occursOnce(virtualFaceIndices)
    L=zeros(Bool,size(virtualFaceIndices))
    for q=1:1:length(virtualFaceIndices)
        L[q]=count(virtualFaceIndices.==virtualFaceIndices[q])==1
    end
    return L
end

function getBoundaryFaces(F)
        
    #Select boundary faces only
    Fs=sort(F,dims=2) #Sorted version of F for uniqueness test

    sizVirtual=maximum(max,F) #-minimum(min,F)+1
    virtualFaceIndices=sub2indn(sizVirtual.*ones(1,size(F,2)),Fs)

    #virtualFaceIndices_uni=unique(virtualFaceIndices)    
    # indFacesUnique=indexin(virtualFaceIndices_uni,virtualFaceIndices)
    # F=F[indFacesUnique,:] #Unique faces
    
    #Remove non-boundary faces (internal faces are shared by >1 element)
    #c=[count(x->x==q,virtualFaceIndices)==1 for q in virtualFaceIndices]
    L=occursOnce(virtualFaceIndices)

    return L

end

function sub2indn(siz,A)
    
    numDim = length(siz)
    k = cumprod([siz[i] for i in 1:length(siz)],dims=1)
    
    if numDim < 2
         error("Invalid size specified. The number of dimensions should be equal or larger than 2")
    end
    
    if any(size(A).==0)
         A=[];
    end
    
    if !isempty(A)    
        
        if size(A,2) != numDim
             error("The specified array size and number of subscript index columns do not match")
        end
    
        #Verify subscripts are within range
        (m,indMax)=findmax(A,dims=1)
        if any(A.<1) || any(m.>siz)
            error("Index out of range")
        end
            
        ind=A[:,1];
        for q=2:1:numDim
            ind = ind .+ (A[:,q].-1).*k[q-1]
        end    

    else
         ind=[]
    end

    return ind
end

function sub2ind(siz,I,J,K)
    
    numDim = length(siz);
    k = cumprod([siz[i] for i in 1:length(siz)],dims=1)
    
    if numDim < 2
         error("Invalid size specified. The number of dimensions should be equal or larger than 2")
    end
    
    if isempty(I) || isempty(J) || isempty(K)
        ind=[]
    else        
    
        #Verify subscripts are within range        
        if any(I.<1) || any(I.>siz[1])
            error("Row index out of range")
        end
           
        if any(J.<1) || any(J.>siz[2])
            error("Column index out of range")
        end

        if any(K.<1) || any(K.>siz[3])
            error("Slice index out of range")
        end

        ind = I .+ (J.-1).*k[1].+ (K.-1).*k[2]
        
    end

    return ind
end

function ind2subn(siz,ind)

    numDim = length(siz);
    k = cumprod([siz[i] for i in 1:length(siz)],dims=1)

    if numDim < 2
        error("Invalid size specified. The number of dimensions should be equal or larger than 2")
    end
   
    #Verify subscripts are within range    
    if any(ind.<1) || any(ind.>prod(siz))
        error("Index out of range")
    end

    A=zeros(Int64,length(ind),numDim) #Initializing output array
    for q=numDim:-1:1 #For all dimensions
        if q==1 #First 1st dimension
            A[:,1]=rem.(ind.-1,k[1]).+1;                
        else
            p=rem.(ind.-1,k[q-1]) .+ 1 # "previous"
            A[:,q]=(ind.-p)./k[q-1] .+ 1 # Current        
            ind=p #Set indices as "previous"
        end    
    end
    return A
end

function ind2subn(siz,ind)

    numDim = length(siz);
  
    k = cumprod([siz[i] for i in 1:length(siz)],dims=1)

    if numDim < 2
        error("Invalid size specified. The number of dimensions should be equal or larger than 2")
   end
   
    #Verify subscripts are within range    
    if any(ind.<1) || any(ind.>prod(siz))
        error("Index out of range")
    end

    A=zeros(Int64,length(ind),numDim) #Initializing output array
    for q=numDim:-1:1 #For all dimensions
        if q==1 #First 1st dimension
            A[:,1]=rem.(ind.-1,k[1]).+1;                
        else
            p=rem.(ind.-1,k[q-1]) .+ 1 # "previous"
            A[:,q]=(ind.-p)./k[q-1] .+ 1 # Current        
            ind=p #Set indices as "previous"
        end    
    end

return A
end

function ind2sub(siz,ind)

    numDim = length(siz);
    k = cumprod([siz[i] for i in 1:length(siz)],dims=1)

    if numDim < 2
        error("Invalid size specified. The number of dimensions should be equal or larger than 2")
    end
   
    #Verify subscripts are within range    
    if any(ind.<1) || any(ind.>prod(siz))
        error("Index out of range")
    end

    A=zeros(Int64,length(ind),numDim) #Initializing output array
    for q=numDim:-1:1 #For all dimensions
        if q==1 #First 1st dimension
            A[:,1]=rem.(ind.-1,k[1]).+1;                
        else
            p=rem.(ind.-1,k[q-1]) .+ 1 # "previous"
            A[:,q]=(ind.-p)./k[q-1] .+ 1 # Current        
            ind=p #Set indices as "previous"
        end
    end          

return A[:,1],A[:,2],A[:,3]
end


function ind2faces(siz,ijk)
 
    Fi=[ijk[:,1]    ijk[:,1].+1 ijk[:,1].+1 ijk[:,1]   ; #Top
        ijk[:,1]    ijk[:,1].+1 ijk[:,1].+1 ijk[:,1]   ; #Bottom
        ijk[:,1]    ijk[:,1].+1 ijk[:,1].+1 ijk[:,1]   ; #Front
        ijk[:,1]    ijk[:,1].+1 ijk[:,1].+1 ijk[:,1]   ; #Back
        ijk[:,1]    ijk[:,1]    ijk[:,1]    ijk[:,1]   ; #Side 1
        ijk[:,1].+1 ijk[:,1].+1 ijk[:,1].+1 ijk[:,1].+1; #Side 2
        ] 

    Fj=[ijk[:,2].+1 ijk[:,2].+1 ijk[:,2]    ijk[:,2]   ; #Top
        ijk[:,2]    ijk[:,2]    ijk[:,2].+1 ijk[:,2].+1; #Bottom
        ijk[:,2]    ijk[:,2]    ijk[:,2]    ijk[:,2]   ; #Front
        ijk[:,2].+1 ijk[:,2].+1 ijk[:,2].+1 ijk[:,2].+1; #Back
        ijk[:,2]    ijk[:,2].+1 ijk[:,2].+1 ijk[:,2]   ; #Side 1
        ijk[:,2].+1 ijk[:,2]    ijk[:,2]    ijk[:,2].+1; #Side 2
        ]

    Fk=[ijk[:,3]    ijk[:,3]    ijk[:,3]    ijk[:,3]   ; #Top
        ijk[:,3].+1 ijk[:,3].+1 ijk[:,3].+1 ijk[:,3].+1; #Bottom
        ijk[:,3]    ijk[:,3]    ijk[:,3].+1 ijk[:,3].+1; #Front
        ijk[:,3].+1 ijk[:,3].+1 ijk[:,3]    ijk[:,3]   ; #Back
        ijk[:,3].+1 ijk[:,3].+1 ijk[:,3]    ijk[:,3]   ; #Side 1
        ijk[:,3].+1 ijk[:,3].+1 ijk[:,3]    ijk[:,3]   ; #Side 2
        ]
    
    F = sub2ind(siz.+1,Fi,Fj,Fk)
    n=size(ijk,1)
    N = [repeat([-1  0  0 ],n,1); 
         repeat([ 0  0  1 ],n,1);
         repeat([-1  0  0 ],n,1);
         repeat([ 1  0  0 ],n,1);
         repeat([ 0 -1  0 ],n,1);
         repeat([ 0  1  0 ],n,1);
         ]
    return F,N    
end

function toSTL(F,V,N,fileName)
    fileHandle = open(fileName, "w")

    #Start solid section
    write(fileHandle, "solid part \n")

    #Loop over faces
    for qf=1:1:size(F,1) #For all faces
        write(fileHandle, "  facet normal " *  @sprintf("%d",N[qf,1]) * " " * @sprintf("%d",N[qf,2]) * " " * @sprintf("%d",N[qf,3]) * "\n")

        write(fileHandle, "    outer loop \n")
        for qv=1:1:size(F,2) #For all face vertices
            #@sprintf("%.12e",pi)
            write(fileHandle,"      vertex  " * @sprintf("%.12e",V[F[qf,qv],1]) * " " * @sprintf("%.12e",V[F[qf,qv],2]) * " " * @sprintf("%.12e",V[F[qf,qv],3]) * "\n")    
        end
        write(fileHandle, "    endloop \n")
        write(fileHandle, "  endfacet \n")
    end
    write(fileHandle, "endsolid part \n")
    close(fileHandle)
end

function im2patch(imageSize,indFound)
    ijkFound=ind2subn(imageSize,indFound)

    #Composed raw faces array
    F,N=ind2faces(imageSize,ijkFound)

    logicBoundaryFaces=getBoundaryFaces(F)
    F=F[logicBoundaryFaces,:]
    N=N[logicBoundaryFaces,:]

    #Check which nodes are used
    indUsed=unique(F)
 
    #Fix face indices in anticipation of a reduced coordinate set
    indFix1=1:1:length(indUsed)
    indFix2=zeros(Int64,maximum(max,F),1)
    indFix2[indUsed]=indFix1
    F=indFix2[F]

    #Create coordinate array
    I,J,K=ind2sub(imageSize.+1,indUsed)

    Iv=convert.(Float64,I).-0.5
    Jv=convert.(Float64,J).-0.5
    Kv=convert.(Float64,K).-0.5

    V=[Iv Jv Kv] #Coordinates

    return F,V,N
end

function toOBJ(F,V,N,fileName)

    splitName = splitext(fileName)    
    fileNameNoExt = splitName[1] 
    fileNameNoExt = splitdir(fileNameNoExt)
    pathName = fileNameNoExt[1] 
    objName = fileNameNoExt[2]

    createUseMTL=false

    # Setup OBJ file for writing
    fileHandleOBJ = open(fileName, "w") #OBJ file handle

    # Comment line
    write(fileHandleOBJ, "# Created using Julia code \n")

    if createUseMTL
        # Setup MTL file for writing    
        fileNameMTL = joinpath(pathName,objName*".mtl")        
        
        # mtllib line
        write(fileHandleOBJ, "mtllib " * objName * ".mtl \n")
    end
        
    # Start object line 
    write(fileHandleOBJ, "o " * objName * " \n")

    # Define vertex field
    for q=1:1:size(V,1) #For all vertices
        write(fileHandleOBJ,"v " * @sprintf("%.12e",V[q,1]) * " " * @sprintf("%.12e",V[q,2]) * " " * @sprintf("%.12e",V[q,3]) * "\n")    
    end

    # Define normal vector field
    for q=1:1:size(N,1) #For all normal vectors
        write(fileHandleOBJ,"vn " * @sprintf("%.12e",N[q,1]) * " " * @sprintf("%.12e",N[q,2]) * " " * @sprintf("%.12e",N[q,3]) * "\n")    
    end
    
    # Solid? line
    write(fileHandleOBJ, "s 1 \n")

    if createUseMTL
        # MTL spec line        
        write(fileHandleOBJ, "usemtl " * objName * " \n")
    end

    # Define face field
    for q=1:1:size(F,1) #For all faces
        textToWrite="f "
        for qf=1:1:size(F,2)
            textToWrite*= @sprintf("%d",F[q,qf]) * "//" * @sprintf("%d",F[q,qf]) * " "           
        end
        write(fileHandleOBJ,textToWrite * "\n")    
    end
    close(fileHandleOBJ)

    # Create MTL file
    if createUseMTL
        fileHandleMTL = open(fileNameMTL, "w") 
        write(fileHandleMTL, "newmtl " * objName * " \n")
        write(fileHandleMTL, "Ns 0.000000 \n")
        write(fileHandleMTL, "Ka 1.000000 1.000000 1.000000 \n")
        write(fileHandleMTL, "Ks 0.000000 0.000000 0.000000 \n")
        write(fileHandleMTL, "Ke 0.000000 0.000000 0.000000 \n")
        write(fileHandleMTL, "Ni 1.450000 \n")
        write(fileHandleMTL, "illum 1 \n")
        write(fileHandleMTL, "Ni 1.450000 \n")
        write(fileHandleMTL, "map_Kd " * objName * ".png" * " \n")
        write(fileHandleMTL, "map_d " * objName * ".png" * " \n")
        close(fileHandleMTL)
    end
end

################## 


s=0.2
r=1;
X,Y,Z= grid3D(-r:s:r,-r:s:r,-r:s:r)

M = sqrt.(X.^2 .+ Y.^2 .+ Z.^2)

imageSize=size(M)

indFound=[sub2ind(imageSize,i[1],i[2],i[3]) for i ∈ findall(M.<=(r+s/100))]

#ijkFound=ind2subn(imageSize,indFound);

F,V,N = im2patch(imageSize,indFound)

FT=[F[:,1] F[:,2] F[:,3]; F[:,3] F[:,4] F[:,1]] #Convert to triangles
NT=[N;N]

toSTL(FT,V,NT,"/home/kevin/Desktop/temp.stl")

fileName = "/home/kevin/Desktop/temp.obj"
toOBJ(F,V,N,fileName)
