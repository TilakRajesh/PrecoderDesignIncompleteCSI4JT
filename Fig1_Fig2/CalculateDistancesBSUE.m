function Dij = CalculateDistancesBSUE(BSs, UEs, rmax)
Dij=zeros(UEs,BSs);
if (exist('Opcion.mat','file') == 0) 
    pos = 27; % cluster center

     XY_or=[-2 4;-1 4;0 4;1 4;
        -2 3;-1 3;0 3;1 3;2 3;
        -3 2;-2 2;-1 2;0 2;1 2;2 2;3 2;
        -3 1;-2 1;-1 1;0 1;1 1;2 1;3 1;
        -3 0;-2 0;-1 0;0 0;1 0;2 0;3 0;4 0;
        -3 -1;-2 -1;-1 -1;0 -1;1 -1;2 -1;3 -1;
        -3 -2;-2 -2;-1 -2;0 -2;1 -2;2 -2;3 -2;
        -2 -3;-1 -3;0 -3;1 -3;2 -3;
        -2 -4;-1 -4;0 -4;1 -4];

    BS = [-4 0;2 4;2 -4];
    % Translate the stations to the 1st quadrant with offset [4 4]
    BS = BS + 4*ones(size(BS));
    XY_or = XY_or + 4*ones(size(XY_or));

    h = (sqrt(3)/2)*rmax;
    Ax = rmax/4;
    Ay = h/4;
    AAx = Ax/4;
    AAy = Ay/4;
    BS(:,1) = Ax*BS(:,1);
    BS(:,2) = Ay*BS(:,2);
    XY_or(:,1)=Ax*XY_or(:,1);
    XY_or(:,2)=Ay*XY_or(:,2);
    
    Usuarios=pos*ones(1,UEs);
    
    c = rsmak('circle');
    Opcion=fnplt(fncmb(c,diag([AAx,AAy]))).';
    save Opcion.mat Opcion Usuarios XY_or BS
else
    load Opcion.mat
end


for ue=1:UEs
    in=ceil(size(Opcion,1)*rand);
    Dist_user=XY_or(Usuarios(ue),:)+Opcion(in,:);
    for bs=1:BSs
        Dij(ue,bs)=sqrt(((BS(bs,1)-Dist_user(1))^2)+((BS(bs,2)-Dist_user(2))^2));
    end
end

end
