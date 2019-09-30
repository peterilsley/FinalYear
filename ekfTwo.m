%% EkfTwo
%------------------------------------
%
%
%
%-------------------------------------

function [x_,P_] = ekfTwo(fstate,init,N,stdM,dt)
%% Constants
I=eye(N,'gpuArray');
%% Variables specified
x_=single(init)
t=0;
q=0.1;
P_= single(eye(N,'gpuArray')); 

P=eye(N,'gpuArray');
Q=q^2*eye(N);%NOT USED RIGHT NOW
R=stdM^2;
H=gpuArray([1 0 0]);%no unit transformation
%% Mazzoni
e=0.04;
 B=0.8;
 prevA=zeros(N,'gpuArray');%create these on the GPU
 nowA=zeros(N,'gpuArray');
%% Precalculations
x=[];
syms x [1 N] 

    a=fstate(x,dt);
    A=jacobian(a,x);%the transition matrix %assumes a linear dt
x1='';
x2='';
x3='';%Need to find a way to not need these
check=1/dt;
odt=dt;
prevT=0;
%% Main operation
while 1
     if (t<1)%only want to predict one ahead for now
         
        if (t+dt>1 && t<1)
            odt=dt;
         dt=1-t  ;
        end
        t=t+dt;
        tic;
        predict;
        prevT=toc;%save how long it took
         
         %dt=odt;
        
       
        
     else 
        [dataAvail,srcWkrIdx,tag] = labProbe;
        if (dataAvail==1)%checking if new measurement has been sent
            [data,srcWkrIdx,tag] = labReceive;
        
            t=0;
            if (tag==0)%checking for done flag
               break;
             elseif (tag==1)
                z=data(1:end-1);
                correct;
            end
       
        end
     end
end 

% labSendReceive(2,2,"Done",1);
    
%%% Project the error covariance ahead(ALSO LOOPED)

%% Predict function
%----------------------------------
%Description: This function is in charge of doing the predicition for the
%             next time state of the function
%
%Comment: This has been done as a nested function so that it can have access to the
%         EKF variables and so that no passing needs to be done between the two
%         while making the main code look cleaner
%
%-----------------------------------
function predict
    %disp("Prediction Phase")
    height=size(x_,1);
    
   x_(height+1,:)=[fstate(x_(height,1:end-1),dt); x_(height,end)+dt]';%time mustn't be included in the prediction 
    
    xN=x_(height+1,:)-x_(height,:);%find better notation
    
    A1=subs(A,x,xN(1:end-1));%ignore the time
    A1=single(vpa(A1));
    
    
    P_=A1*P*A1.';
    %% Mazzoni section
     
    mazzoni;
    %eps is the minimum possible step that can be performed in matlab due to
    %floating point difference
    
function mazzoni() 
    prevA=nowA; %probably slowed down by all of these
    nowA=A1;
    %prevA=single(vpa(subs(A,x,x_(height,:))));
    %nowA=subs(A,x,x_(height+1,:));%need to tie dt into the jacobian
     
    E=dt^2/2*((nowA-prevA)/(3*dt)-prevA.^2/6);%*f(u);%what is f???
    E=single(vpa(E));
    
    C = bsxfun(@rdivide, E, x_(height+1,1:end-1));%%will use gpu for computing
   
    EMax=max(abs(C(:)));%%returns the maximum element in the E matrix
%     if (prevT~=0 && t~=1)%MY STUFF FOR REAL-TIME
%         timeMes=1;
%         CalcTime=t+((timeMes-t)/prevT-1)*dt;
%         alpha=timeMes/CalcTime;
%         e=alpha^2*e;
%     end
    
    Qcoeff=B*sqrt(e./EMax);
    dt=gather(dt*Qcoeff)%why???????
    
    
   %classUnderlying(dt)
   
end 
end

%% Correct function
function correct
   % disp("Correction phase")
    
    %H is typically found by taking the jacobian of the measurement
    %equation at point x
    %compute the kalman gain;
    
   % K=P_*H'*inv(H*P_*H'+R);
   L=chol(H*P_*H'+R);
   height=size(x_,1);
   
   %K=P_*H'*inv(H*P_*H'+R);
   
  U=L\eye(1);%need to check if eye(1) is fine at some point
   A1=L'\U;
  K=P_*H'*A1;
    %update estimate with measurement
    x_(height,1:end-1)=x_(height,1:end-1)+K'.*((z-x_(height,1:end-1)));%musn't update the timestamp
   % x_(height,:)=x_(height,:)+(K'.*(z-x_(height,:)));%x being the current estimate
    
    %Update the error covariance
 
    P =(I-K*H)*P_;
   
    
end



end


