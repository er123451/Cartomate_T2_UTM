%-------------Archivo a cargar-------------
file = "datos-practica-2.txt";
%------------------------------------------

fid = fopen("sal-P2a-02785633B","w");
data = load(file);
a = data(1,1);
invf = data(2,1);
long0 = data(1,2);
k0 = data(2,2);
[m,n] = size(data);
data1 = data(3:m,:);
data2 = data(1:2,:);
fprintf(fid,'%%Ernesto Hontecillas Molina, DNI: 02785633B\n');
fprintf(fid,'%f\t%f\t0\t0\n',data2(1,:));
fprintf(fid,'%.9f\t%f\t0\t0\n',data2(2,:));
for i = 1:(m-3)
    [x,y,convm,k] = commands2.GKdir(data1(i,1),data1(i,2),long0,k0,0,0,a,invf);
    fprintf(fid,'%.4f\t%.4f\t%.10f\t%.10f\n',x,y,convm,k);
end
fclose(fid);
