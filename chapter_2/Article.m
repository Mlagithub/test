clear
clc

load 10-1.txt
load 10-2.txt
load 10-3.txt
load 10-4.txt
load 10-5.txt

x1 = X10_1(:,1);
y1 = X10_1(:,2);
x2 = X10_2(:,1);
y2 = X10_2(:,2);
x3 = X10_3(:,1);
y3 = X10_3(:,2);
x4 = X10_4(:,1);
y4 = X10_4(:,2);

if (X10_5 == 1)
    load 6mm.csv
    x5 = X6mm(:,1);
    y5 = X6mm(:,2);
elseif (X10_5 == 2)
    load 8mm.csv
    x5 = X8mm(:,1);
    y5 = X8mm(:,2);
elseif (X10_5 == 3)
    load 10mm.csv
    x5 = X10mm(:,1);
    y5 = X10mm(:,2);
else
    load V.csv
    x5 = V(:,1);
    y5 = V(:,2);
end 

plot(x1,y1,'sr',x2,y2,'dg',x3,y3,'xb',x4,y4,'--m',x5,y5,'pk')

delete('E:\Fortran\Projects\Article\Article\10-*.txt')