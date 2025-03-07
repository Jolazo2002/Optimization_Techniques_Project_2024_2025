clearvars; % clear workspace
clc; % clear command window
close all;

%% Constants

a=[1.25 1.25 1.25 1.25 1.25 1.5 1.5 1.5 1.5 1.5 1 1 1 1 1 1 1];

c=[54.13 21.56 34.08 49.19 33.03 21.84 29.96 24.87 47.24 33.97 26.89 32.76 39.98 37.12 53.83 61.65 59.73];

V=100;

%% Problem construction

f=@(x) Obj_func(x,a,c);


%% Constaints

Aeq=zeros(9,17);
b=zeros(9,1);

% edge A

Aeq(1,1:4)=1;
b(1)=V;

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
b(9)=V;




%% Genetic Algorithm

Psize=200;
Pd=0.7;
Pm=0.1;

init=20*ones(1,17);


n=normrnd(0,2,[Psize-1 1]);

n=[0;n];

population=init+n;


i=0;
loop=150;

lfit=zeros(1,loop);

minfit=100000; %this will be used to track the best fitness

mingene=zeros(1,17); %this will be used to track the gene that gives the best fitness


while i<loop

    i=i+1;
    fitness=zeros(Psize,1);
    
    % Projection
    for j=1:Psize
        
        while norm(Aeq*population(j,:)'-b)>10^(-3) || min(population(j,:))<0 || max(population(j,:)-c)>0
            population(j,:)=population(j,:)-(Aeq'*pinv(Aeq*Aeq')*(Aeq*population(j,:)'-b))';

            population(j,:)=min(max(population(j,:), 0), c);
        end
        
        
        
        fitness(j)=f(population(j,:));
        

    end

    %10 tournament selection
    k=10;

    Parents=zeros(Psize,1);

    for j=1:Psize
        perm=randperm(Psize);

        selected=perm(1:k);

        fitselect=fitness(selected);

        [~,ind]=min(fitselect);

        Parents(j)=selected(ind);



    end

    Parents=reshape(Parents,Psize/2,2);

    randcrossp=rand(size(Parents,1),1);


    %deciding which parents are going to breed
    crosscoeff=rand();

    crosscheck= randcrossp < Pd;


    %breeding
    offsprings1=(crosscoeff)*population(Parents(:,1),:)+(1-crosscoeff)*population(Parents(:,2),:);

    
    offsprings2=(crosscoeff)*population(Parents(:,2),:)+(1-crosscoeff)*population(Parents(:,1),:);


    offsprings=[offsprings1(crosscheck == 1,:);offsprings2(crosscheck == 1,:)];

    temp=[population;offsprings];
    
    %deciding which offsprings are going to mutate
    randmut=rand(size(temp,1),17);

    mutchech=double(randmut < Pm);

    N=(normrnd(0,1,size(mutchech)).*mutchech);

    temp=temp+N;

    %Recalculating the fitness for the old and new samples
    fitness2=zeros(length(N),1);
    
    % Projection
    for j=1:length(temp)
        
        while norm(Aeq*temp(j,:)'-b)>10^(-3) || min(temp(j,:))<0 || max(temp(j,:)-c)>0
            temp(j,:)=temp(j,:)-(Aeq'*pinv(Aeq*Aeq')*(Aeq*temp(j,:)'-b))';

            temp(j,:)=min(max(temp(j,:), 0), c);

   
        end
        
        
        fitness2(j)=f(temp(j,:));
        

    end

    %sorting and keeping only the most fitted samples
    [t,sortind]=sort(fitness2,'ascend');

    population=temp(sortind(1:Psize),:);

    lfit(i)=t(1);

    if t(1)<minfit
        minfit=t(1);
        mingene=population(1,:);

    end




end