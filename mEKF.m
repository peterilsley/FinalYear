%-----------------------------------------
%Title: An optimised EKF
%Author: Dylan Ilsley
%Date: August 2019
%----------------------------------------

%% Variables
numState=3;
N=10;%number of measurements
stdM=0.4;%the deviation of the measurement
stdP=0;%currently not being used in simulation
%not going to worry about passing Q or R for now

%% Mazzoni
e=1;%tolerance level


%% Creating the EKF param



f=@(x)[x(1)+x(3)*cos(x(2)*pi/180);x(2)+5;x(3);];%for simulations, doesn't use dt

%% Prepping measurements and true values ahead of time
X=simulateEvent(f,numState,N);%simulating true values
X(:,1:3)
Z=getMeasurements(X,stdM,numState);
%disp(Z(:,1:3));

f=@(x,dt)[x(1)+dt*x(3)*cos(x(2)*pi/180);x(2)+5*dt;x(3)];%used for predictions %final x is the time
dt=single(1);
s=single([0.1 45 1.667 X(1,end)]);%initial state ESTIMATE %putting the timestamp on

%-----------------------------------------------------
%PLOTTING THE SIMULATION VERSUS THE INITIAL READINGS
%plot(Z(:,1))
%hold on
%plot(X(:,1))
%hold off
%-----------------------------------------------



%% parallel set up
Z=gpuArray(Z);%sending to gpu after having been created

x=gpuArray(s);
%% Live simulation (OLD)
tic;
spmd 
    switch labindex
        case 1 %EKF WORKER
             sV=ekfTwo(f,x,numState,stdM,dt);%set up should only require the initial state estimate and not any measurements
             
        case 2 %MEASUREMENT READINGS
            
            for i=2:N
                %pause(0.1);
                labSend(Z(i,:),1,1)%lab send is sending just one?
                
            end
            labSend(Z(i,:),1,0);%telling ekf that there won't be any more sending
            %labSendReceive(1,1,"",2);
       
            %going to use a shared array for update purposes
       
            
    end
    
end
toc

%plotly_path = fullfile(pwd, 'plotly');
%addpath(genpath(plotly_path));
%plotlysetup('dylanIlsley', 'y4dVRQrV5gJwwh3IDItN');
figure();
sV=sV{1};%gets what sV is in worker 1
numVar=1;
lineWidth=8;

subplot(3,1,1);
plot(sV(:,end),sV(:,numVar),'-o')
hold on;
plot(X(:,numVar))
plot(Z(:,numVar))
xlabel('time (s)');
ylabel('Distance (m)');
hold off;

subplot(3,1,2);
plot(sV(:,end),sV(:,numVar+1),'-o')
hold on;
plot(X(:,numVar+1))
plot(Z(:,numVar+1))
xlabel('time (s)');
ylabel('Velocity (m/s)');
hold off;

subplot(3,1,3);
plot(sV(:,end),sV(:,numVar+1),'-o')
hold on;
plot(X(:,numVar+1))
plot(Z(:,numVar+1))
xlabel('time (s)');
ylabel('Angle (Theta)');
hold off;


plotly(gather(sV(:,end)),gather(sV(:,numVar+1)));
fig=fig2plotly(gcf);

saveplotlyfig(fig, 'EKF_Example.png')
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

X=single(zeros(N, numState+1));%always 6 from the time
dtTime=1;
time=datenum(dtTime);

%want in a double format
X(1,:)=[[0.1 45 1.667],time]; %initial state of the bug


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
    Z=zeros(size(X,1), numState+1);
    Z=[X(:,1:3)+std.*randn(size(X,1),3), X(:,4:end)];%need to fix this since right now it isn't truly random
end



%% Live measurement function



    

