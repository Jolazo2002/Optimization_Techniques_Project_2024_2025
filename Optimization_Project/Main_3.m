clearvars; % clear workspace
clc; % clear command window
close all;

%% Constants

a=[1.25 1.25 1.25 1.25 1.25 1.5 1.5 1.5 1.5 1.5 1 1 1 1 1 1 1];

c=[54.13 21.56 34.08 49.19 33.03 21.84 29.96 24.87 47.24 33.97 26.89 32.76 39.98 37.12 53.83 61.65 59.73 115];

%V=85;

%% Problem construction

f=@(x) Obj_func(x,a,c);

x0=20*ones(1,17);
x0=[x0 100];

%% Constaints

Aeq=zeros(9,18);
b=zeros(9,1);

% edge A

Aeq(1,1:4)=1;
Aeq(1,18)=-1;

% edge B

Aeq(2,1)=1;

Aeq(2,[5 6])=-1;

% edge C

Aeq(3,2)=1;

Aeq(3,[7 8])=-1;

% edge D

Aeq(4,[3 8 9])=1;

Aeq(4,[11 12 13])=-1;

% edge E

Aeq(5,4)=1;

Aeq(5,[9 10])=-1;

% edge F

Aeq(6,[5 14])=1;

Aeq(6,16)=-1;

% edge G

Aeq(7,[6 7 13])=1;

Aeq(7,[14 15])=-1;

% edge H

Aeq(8,[10 11])=1;

Aeq(8,17)=-1;

% edge I

Aeq(9,[12 15 16 17])=1;
Aeq(9,18)=-1;


% lower/upper bounds

lb=zeros(1,17);
lb=[lb 85];
ub=c;

[sol,val]=fmincon(f,x0,[],[],Aeq,b,lb,ub);