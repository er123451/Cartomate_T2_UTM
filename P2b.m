%-------------Archivo a cargar-------------
file = "sal-P2a";
%------------------------------------------

fid = fopen("sal-P2b","w");
data = load(file);
a = data(1,1);
invf = data(2,1);
long0 = data(1,2);
k0 = data(2,2);
[m,n] = size(data);
data1 = data(3:m,1:2);
data2 = data(1:2,1:2);
fprintf(fid,'%.6f\t%.6f\n',data2(1,:));
fprintf(fid,'%.6f\t%.6f\n',data2(2,:));
for i = 1:(m-3)
    [long,lat] = commands2.GKinv(data1(i,1),data1(i,2),long0,k0,0,0,a,invf);
    fprintf(fid,'%.6f\t%.6f\n',long,lat);
end

fclose(fid);
