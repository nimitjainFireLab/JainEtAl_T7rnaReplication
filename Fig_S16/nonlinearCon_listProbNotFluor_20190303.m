clc; 
clear all;

a=[0.5688073394
0.6751269036
0.6549707602
0.6570048309
0.6558441558
0.6798245614
0.6501128668
0.652173913
0.606442577
0.6067961165
0.642578125
0.6460176991
0.671875
0.6092124814
0.6075533662
0.6044303797];

c=0.01:0.01:10;
d=zeros(length(c),1);
for i=1:length(c)
    d(i,1)=exp(-c(i))+c(i)*exp(-c(i));
end

plot(c,d,'r-') %monotonicity here for positive c is a good sign--means that a lambda can be found--use nonlinear equation solver for this (with constraint that lambda must be >=0)

b=zeros(length(a),1);
resVec=zeros(length(a),1);
for i=1:length(a)
    [b(i,1),resVec(i,1)]=lsqnonlin(@(x)(exp(-x))+(x*exp(-x))-a(i,1),1,0);
end

exp(-b)+b.*exp(-b)
a

