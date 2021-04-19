fid = fopen("1.txt","w");
data = load("datos-practica-2.txt");
a = data(1,1);
invf = data(2,1);
long0 = data(1,2);
k0 = data(2,2);
[m,n] = size(data);
data1 = data(3:m,:);
for i = 1:(m-3)
    [x,y] = commands2.GKdir(data1(i,1),data1(i,2),long0,k0,0,0,a,invf);
    fprintf(fid,'\n%.4f\t%.4f\n',x,y);
end
fclose(fid);

