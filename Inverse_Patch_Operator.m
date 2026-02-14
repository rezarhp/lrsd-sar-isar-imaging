function f=Inverse_Patch_Operator(Fin,Nx,Ny,Mx,My,Ns)
    global indexes
    F=reshape(Fin,Mx*My,Ns);
    f_matrix=zeros(Ny,Nx);
    f_nominator=zeros(Ny,Nx);
    f_denominator=zeros(Ny,Nx);
    Iones=ones(My,Mx);
    for ind=1:Ns
        patch=reshape(F(:,ind),[My Mx]);
        r=indexes(ind,1);
        c=indexes(ind,2);
        rs=r:r+My-1;
        cs=c:c+Mx-1;
        f_nominator(rs,cs)=f_nominator(rs,cs)+patch;
        f_denominator(rs,cs)=f_denominator(rs,cs)+Iones;
    end
    for i=1:Ny
        for j=1:Nx
            if f_denominator(i,j)==0, continue; end
            f_matrix(i,j)=f_nominator(i,j)/f_denominator(i,j);
        end
    end
    f=f_matrix(:);
end