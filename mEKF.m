%-----------------------------------------
%Title: An optimised EKF
%Author: Dylan Ilsley
%Date: August 2019
%----------------------------------------

%% Variables
numState=3;
N=50;%number of measurements
stdM=0.01;%the deviation of the measurement
stdP=0;%currently not being used in simulation
%not going to worry about passing Q or R for now
filename=[num2str(N),'_1S.xlsx']
%% Mazzoni
e=1;%tolerance level


%% Creating the EKF param


f=@(x)[x(1)+x(3)*cos(x(2)*pi/180);x(2)+1.5;x(3)+2;];%TEST CASE 1
%f=@(x)[x(1)+1.5;x(2)+1.5;x(3)+2;];%Mazzo test

%% Prepping measurements and true values ahead of time
X=simulateEvent(f,numState,N);%simulating true values

Z=getMeasurements(X,stdM,numState);
%disp(Z(:,1:3));

f=@(x,dt)[x(1)+dt*x(3)*cos(x(2)*pi/180);x(2)+1.5*dt;x(3)+2*dt;];%used for predictions %final x is the time
%f=@(x,dt)[x(1)+dt*1.5;x(2)+1.5*dt;x(3)+2*dt;]; %FOR
%MAZZO
dt=single(1);
%s=single([0.1 45 1.667  X(1,end)]);%initial state ESTIMATE %putting the timestamp on
s=Z(1,:);
%-----------------------------------------------------
%PLOTTING THE SIMULATION VERSUS THE INITIAL READINGS
%plot(Z(:,1))
%hold on
%plot(X(:,1))
%hold off
%-----------------------------------------------



%% parallel set up
%Z=gpuArray(Z);%sending to gpu after having been created

init=s;%;gpuArray(s);
syms x [1 numState] 
a=f(x,dt);
A=jacobian(a,x);
A=matlabFunction(A);

%  X=gather(X);
%  Z=gather(Z);
%  init=gather(init);%converting back to cpu

%% Live simulation 
tic;
spmd 
    switch labindex
        case 1 %EKF WORKER
             sV=ekfTwo(f,init,numState,stdM,dt,A,x);%set up should only require the initial state estimate and not any measurements
             
        case 2 %MEASUREMENT READINGS
            
            for i=2:N
                %pause(0.001);
                labSend(Z(i,:),1,1)%lab send is sending just one?
                
            end
            labSend(Z(i,:),1,0);%telling ekf that there won't be any more sending
            %labSendReceive(1,1,"",2);
       
            %going to use a shared array for update purposes
       
            
    end
    
end
toc
%% Plotting of figure
%plotly_path = fullfile(pwd, 'plotly');
%addpath(genpath(plotly_path));
%plotlysetup('dylanIlsley', 'y4dVRQrV5gJwwh3IDItN');
%figure();

sV=sV{1};%gets what sV is in worker 1
numVar=1;
lineWidth=8;
return;
figure()
subplot(3,1,1);
plot(sV(1:end-1,end),sV(1:end-1,numVar),'-o')
hold on;
plot(X(:,numVar))
plot(Z(:,numVar))
xlabel('time (years)');
ylabel('Population 1');
hold off;

subplot(3,1,2);
plot(sV(1:end-1,end),sV(1:end-1,numVar+1),'-o')
hold on;
plot(X(:,numVar+1))
plot(Z(:,numVar+1))
xlabel('time (years)');
ylabel('Population 2');
hold off;

subplot(3,1,3);
plot(sV(1:end-1,end),sV(1:end-1,numVar+2),'-o')
hold on;
plot(X(:,numVar+2))
plot(Z(:,numVar+2))
xlabel('time (years)');
ylabel('Population 3');
hold off;


%% Writing to Excel tables
sV=gather(sV);
X=gather(X);
Z=gather(Z);

saveas(gcf,'epsFig.eps','epsc') 

writematrix(sV,filename,'Sheet',1);
writematrix(X,filename,'Sheet',2);
writematrix(Z,filename,'Sheet',3);

%% Writing to plotly

%plotly(gather(sV(:,end)),gather(sV(:,numVar+1)));
%fig=fig2plotly(gcf);

%saveplotlyfig(fig, 'EKF_Example.png')
%% Live simulation (NEW)
% 
% p=gcp();
% 
% fun=parfeval(p,@ekfTwo,2,f,x,2); 
% [values,P]=fetchOutputs(fun)
% cancelFutures = onCleanup(@() cancel(fun));%making sure it isn't runnning at end of program
%% pre-simulation functions
%function to simulate the true values
function X = simulateEvent(f,numState,N)

X=zeros(N, numState+1,'single','gpuArray');%always 6 from the time
dtTime=1;
time=datenum(dtTime);

%want in a double format
X(1,:)=([rand(1,numState),time]);
%X(1,:)=[[0.1 45 1.667 ],time]; %initial state of the bug


 for t=2:N %simulating one minute of travel
%      prevV=X(t-1,3);
%      prevT=X(t-1,2);
%      prevD=X(t-1,1);
%      d=prevV*cos(prevT*pi/180);
    dtTime=dtTime+seconds(1);
    time=datenum(dtTime);%been one more minute
     
      XnT=X(t-1,1:end-1);%X no Time
     XnT=f(XnT);
     
     X(t,:)=[XnT',time];%[prevD+d prevT prevV];
    
    % X(t,:)=f(X(t-1,:));
 end

end %function
%function to create the observed measurements including the uncertainty of
%the measurement
function Z = getMeasurements(X,std,numState)%
    Z=zeros(size(X,1), numState+1,'single','gpuArray');
    Z=[X(:,1:numState)+std.*randn(size(X,1),numState), X(:,numState+1:end)];%need to fix this since right now it isn't truly random
end

function data = func(x, dt)
        if nargin<2
         dt = 1;
        end
        r=0.5;
        data=r.*x.*(1-x);
        data=data.';
end


%% WTEC







    

