clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Patricia Amado Caballero
% Email: patricia.amado@uva.es
% Date: 2025-08-27

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path('Data',path)

parameters=["p1" "p2" "mean11" "mean12" "mean21" "mean22"  "sigma11" "sigma12" "sigma21" "sigma22" "rho1" "rho2"];

load('MediumOcclusionMap_samplePatient.mat')

[X,Y]  =meshgrid((1:1:100),(1:1:45));
data=map_patient;
titulo='Oclussion Map  ';
data=double(data);
X0=getX0(data,X,Y);
[solution,finalCurve] =   gmm2D(data,X,Y,X0);

%Visualization
figure;mesh(X,Y, data)
colormap(jet)
colorbar
title([titulo],'Interpreter','none');
hold on
scatter(solution(3),solution(4),100,RGBorange,'o','filled')
hold on
scatter(solution(5),solution(6),100,RGBPurple,'o','filled')
hold off   

save('GMM_SamplePatient.mat'  ,'solution','finalCurve','parameters','X0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [solution,finalCurve] =   gmm2D(data,X,Y,X0)

parameters.X =X;
parameters.Y =Y;
parameters.data =data;


Aeq = [1 1 0 0 0 0 0 0 0 0 0 0];
Beq = 1;
LB  = [0 0  0 0 0 0 0 0 0 0 -1 -1]';
UB  = [1 1 100 45 100 45 100 45 100 45 1  1]';


options = optimoptions('fmincon','Display','off');
solution = fmincon(@(x) myFunction(x,parameters), X0,[],[], Aeq, Beq, LB, UB,[],options);

finalCurve  =solution(1)*normpdf(X,solution(3),solution(7)).*normpdf(Y,solution(4)+solution(11) * solution(8)/solution(7)*(X-solution(3)),solution(8)*sqrt(1-solution(11)^2))+...
             solution(2)*normpdf(X,solution(5),solution(9)).*normpdf(Y,solution(6)+solution(12) * solution(10)/solution(9)*(X-solution(5)),solution(10)*sqrt(1-solution(12)^2));
         
end

function output = myFunction(x, parameters);
  data= parameters.data;
    X=parameters.X;
    Y=parameters.Y;
    buffer = x(1)*normpdf(X,x(3),x(7)).*normpdf(Y,x(4)+x(11) * x(8)/x(7)*(X-x(3)),x(8)*sqrt(1-x(11)^2))+...
             x(2)*normpdf(X,x(5),x(9)).*normpdf(Y,x(6)+x(12) * x(10)/x(9)*(X-x(5)),x(10)*sqrt(1-x(12)^2));
         
    output = sum(abs(buffer(:)-data(:)));
end

function X0 = getX0(data, X, Y, index);

      maxY=double(max(data));
      aux=islocalmax(maxY);
       indAux=find(aux==1);
       localMax=sort(maxY(indAux));
       maxX=max(data,[],2);
       [coordy1 coordx1]=find(data==localMax(length(localMax)));
        max1=localMax(length(localMax));
        
       if(length(localMax)>1)
           [coordy2 coordx2]=find(data==localMax(length(localMax)-1));
           max2=localMax(length(localMax)-1);
           total=max1+max2;
            p1=max1/total;
       else
           localMax2=sort(maxY);
            [coordy2 coordx2]=find(data==localMax2(length(localMax2)-1));
            max2=localMax2(length(localMax2)-1);
            p1=1-(max2/(max1*100));      
       end

    mean11=X(coordy1,coordx1); 
    mean12=Y(coordy1,coordx1); 
    mean21=X(coordy2,coordx2); 
    mean22=Y(coordy2,coordx2); 
    p2=1-p1;
    
    X0=[p1 p2 mean11 mean12 mean21 mean22  1 1 1 1   0.5 0.5];
end

