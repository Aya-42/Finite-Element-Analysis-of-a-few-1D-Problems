clc
close all 
clear all

E=29e6;  

L=[]; %Add values

A=[];

I=[];

N=8; %number of elements
DOF=2;
nodes=7;
n=nodes*DOF;

wc1=200;
wc2=350;
wl1=150;
wl2=150;
External=1500;
%----------------------------


%Local K set up
for i=1:N
    K(:,:,i)=2*E*I(i)/L(i)^3[6 3*L(i) -6 3*L(i)
                             3*L(i) 2*L(i)^2 -3*L(i) L(i)^2
                            -6 -3*L(i) 6 -3*L(i)
                            3*L(i) L(i)^2 -3*L(i) 2*L(i)^2]
end

%Global K set up
KG1(n,n)=0; 
KG2(n,n)=0; 
KG3(n,n)=0;
KG4(n,n)=0; 
KG5(n,n)=0;
KG6(n,n)=0;  
KG7(n,n)=0; 
KG8(n,n)=0;

% ELEMENT 1   %left is local position/ right is global position
P(1)=1; 
P(2)=2;
P(3)=3;
P(4)=4;

for i=1:DOF*2
for j=1:DOF*2
KG1(P(i),P(j))=Ke(i,j,1); 
end   
end 

% ELEMENT 2   %left is local position/ right is global position
P(1)=4; 
P(2)=5;
P(3)=6;
P(4)=10;

for i=1:DOF*2
for j=1:DOF*2
KG2(P(i),P(j))=Ke(i,j,2); 
end   
end 

% ELEMENT 3    %left is local position/ right is global position
P(1)=4;
P(2)=5;
P(3)=6;
P(4)=7;
for i=1:DOF*2
for j=1:DOF*2 
KG3(P(i),P(j))=Ke(i,j,3); 
end   
end 

% ELEMENT 4   %left is local position/ right is global position
P(1)=4;
P(2)=5;
P(3)=6;
P(4)=16;
for i=1:DOF*2
for j=1:DOF*2
KG4(P(i),P(j))=Ke(i,j,4); 
end   
end 

% ELEMENT 5  %left is local position/ right is global position
P(1)=10;   
P(2)=11;
P(3)=12;
P(4)=13;
for i=1:DOF*2
for j=1:DOF*2
KG5(P(i),P(j))=Ke(i,j,5); 
end   
end 

% ELEMENT 6 %left is local position/ right is global position
P(1)=16;    
P(2)=17;
P(3)=18;
P(4)=19;
for i=1:DOF*2
for j=1:DOF*2 
KG6(P(i),P(j))=Ke(i,j,6); 
end   
end 

% ELEMENT 7 %left is local position/ right is global position
P(1)=7;     
P(2)=8;
P(3)=9;
P(4)=16;
for i=1:DOF*2
for j=1:DOF*2 
KG7(P(i),P(j))=Ke(i,j,7); 
end   
end 

% ELEMENT 8 %left is local position/ right is global position
P(1)=7;    
P(2)=8;
P(3)=9;
P(4)=10;

for i=1:DOF*2
for j=1:DOF*2 
KG8(P(i),P(j))=Ke(i,j,8); 
end   
end 

KG=KG1+KG2+KG3+KG4+KG5+KG6+KG7+KG8;
Kf=KG;


%Global Force Matrix 
for i=1:n
    F1(i)=0;
    F2(i)=0;
    F3(i)=0;
    F4(i)=0;
    F5(i)=0;
    F6(i)=0;
    F7(i)=0;
    F8(i)=0;
end

% ELEMENT 8 %left is local position/ right is global position
P(1)=7;   
P(2)=8;
P(3)=9;
P(4)=10;
%PVM
F8(P(1))=0;
F8(P(2))=-wc1*L(8)/2-3*wl1*L(8)/20;
F8(P(3))=-wc1*L(8)^2/12-wl1*L(8)^2/30;
F8(P(4))=0;

% ELEMENT 7 %left is local position/ right is global position
P(1)=7;    
P(2)=8;
P(3)=9;
P(4)=16;
%PVM
F7(P(1))=0;
F7(P(2))=-wc1*L(7)/2-3*wl1*L(7)/20;
F7(P(3))=wc1*L(7)^2/12+wl1*L(7)^2/30;
F7(P(4))=0;

% ELEMENT 5     %left is local position/ right is global position
P(1)=10;    
P(2)=11;
P(3)=12;
P(4)=13;
%PVM
F5(P(1))=0;
F5(P(2))=-wc2*L(5)/2-3*wl2*L(5)/20;
F5(P(3))=-wc2*L(5)^2/12-wl2*L(5)^2/30;
F5(P(4))=0;

% ELEMENT 6   %left is local position/ right is global position
P(1)=16; 
P(2)=17;
P(3)=18;
P(4)=19;
%PVM
F6(P(1))=0;
F6(P(2))=-wc2*L(6)/2-3*wl2*L(6)/20;
F6(P(3))=wc2*L(6)^2/12+wl2*L(6)^2/30;
F6(P(4))=0;

F=F1+F2+F3+F4+F5+F6+F7+F8;
Ff=F;

%BCs
KG(2,:)=[];
KG(:,2)=[];
KG(1,:)=[];
KG(:,1)=[];

F(2)=[];
F(1)=[];

%Solution
U=KG\F';
U=[0,0,U']'

%Other information
R=Kf*U-Ff'






