close all;
NX=60;
NZ=28;
NGPS=11;
load('NezsorGPS.dat');
for i=1:NZ
    
    x(NX*(i-1)+1:NX*(i))=1:NX;
    y(NX*(i-1)+1:NX*(i))=0.;
    z(NX*(i-1)+1:NX*(i))=i;
end

for i=1:NGPS
    figure
    %quiver(x',z',zeros(1,NX*NZ)',ones(1,NX*NZ)')
    quiver(x',z',1000*NezsorGPS(1+(i-1)*NX*NZ:(i)*NX*NZ,2),1000*NezsorGPS(1+(i-1)*NX*NZ:(i)*NX*NZ,1))
    %quiver3(x',y,z',NezsorGPS(1+(i-1)*NX*NZ:(i)*NX*NZ,1),NezsorGPS(1+(i-1)*NX*NZ:(i)*NX*NZ,2),NezsorGPS(1+(i-1)*NX*NZ:(i)*NX*NZ,3))
end