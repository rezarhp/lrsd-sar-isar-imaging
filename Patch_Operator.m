function Fout=Patch_Operator(f,Nx,Ny,Mx,My,Ns,rt,ct)
    global indexes
    indexes=zeros(Ns,2);
    F=zeros(Mx*My,Ns);
    f_matrix=reshape(f,[Ny Nx]);
    ind=0;
    for r=rt
        rs=r:r+My-1;
        for c=ct
            patch=f_matrix(rs,c:c+Mx-1);
            ind=ind+1;
            indexes(ind,:)=[r,c];
            F(:,ind)=patch(:);
        end
    end
    Fout=F(:);
end