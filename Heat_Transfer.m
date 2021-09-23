clc
clear all
close all

T1=200;
Tf=50;

q1=70;  %thermal conductivity 
q2=40;
q3=20;

L1=2e-2;
L2=2.5e-2;
L3=4e-2;

h=10;
A=1;
n=4; %nodes
sum=0;

%Local and Global Stiffness Matrices 
k1=q1*A/L1*[1 -1;-1 1];
k2=q2*A/L2*[1 -1;-1 1];
k3=q3*A/L3*[1 -1;-1 1]+h*A*[0 0;0 1];

K1=zeros(n);
K1(1:2,1:2)=k1(:,:);
K2=zeros(n);
K2(2:3,2:3)=k2(:,:);
K3=zeros(n);
K3(3:4,3:4)=k3(:,:);

K=K1+K2+K3

%Forces
f1=[0;0];
f2=[0;0];
f3=h*Tf*A*[0;1];

F1(4)=0;
F1(:,1:2)=f1(:,:);
F2(4)=0;
F2(:,2:3)=f2(:,:);
F3(4)=0;
F3(:,3:4)=f3(:,:);

F=F1+F2+F3

%BC's and Temperature 
F(1)=T1;
K(1,:)=0;
K(1,1)=1;

%Solution
T=K\F'


