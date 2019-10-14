%% EkfTwo
%------------------------------------
%
%
%
%-------------------------------------

function [x_,P_] = ekfTwo(fstate,init,N,stdM,dt,A,x)
%% Constants
I=eye(N,'single','gpuArray');
%% Variables specified
x_=init;
t=0;
q=0.1;
P_= eye(N,'single','gpuArray'); 

P=eye(N,'single','gpuArray');
%Q=q^2*eye(N);%NOT USED RIGHT NOW
R=stdM^2;
%H=gpuArray([1 0 0 ]);%no unit transformation
H=ones(1,N,'single','gpuArray');
%% Mazzoni
e=0.0000000003;
 B=0.8;
 prevA=zeros(N,'gpuArray');%create these on the GPU
 nowA=zeros(N,'gpuArray');
 %a=fstate(x,dt);
 %A=jacobian(a,x)
 
 %% MY MAZZONI
 dtLog=ones(1,5,'single','gpuArray');
 dtCount=1;
 tMin=0.0001; 
 WCET=0;
 first=true;
 times=zeros(50,1,'single','gpuArray');
%% Precalculations

%syms(sym('x', [1 N]))=''


%the transition matrix %assumes a linear dt
%might just have to do the jacobian regularly?

check=1/dt;
odt=dt;
prevT=0;
count=0;
prediction=0;
timeMes=0.08;
%% CPU conversion

% I=gather(I);
% P_= gather(P_); 
% P=gather(P);
% H=gather(H);

%% Main operation
height=size(x_,1);
startInterval=tic;

while 1
     if (t<timeMes)%only want to predict one ahead for now
         
        if (t+dt>timeMes && t<timeMes)
            odt=dt;
         dt=timeMes-t  ;
        end
        t=t+dt;
        startPred=tic;
        %%DIAGNOSTIC PURPOSES
        %----------------------------------
        %f=@()predict;
       % predictTime=gputimeit(f)
        %------------------------------------
        predict;
        height=height+1;
        x_(height,:)=prediction;
        prevT=toc(startPred);%save how long it took to predict
        
         first=false;
      
       %dt=odt;
        
     else 
        [dataAvail,srcWkrIdx,tag] = labProbe;
        if (dataAvail==1)%checking if new measurement has been sent
            len=toc(startInterval);
            [data,srcWkrIdx,tag] = labReceive;
            count=count+1;
            times(count)=len;
            t=0;
            if (tag==0)%checking for done flag
               disp(mean(times))
               disp(max(times))
               disp(sum(times(:)>0.1))
               break;
             elseif (tag==1)
               
                z=data(1:end-1);
                
         
                correct;
                startInterval=tic;
               %g=@()correct;
                %correctTime=gputimeit(g)
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
    
    previous=x_(height,1:end-1);

    prediction=[fstate(previous,dt); x_(height,end)+dt]';%gpu computations are faster when done like this
   
    xN=prediction(1:end-1)-previous;
    xD = mat2cell(xN,1,ones(1,numel(xN)));
    %
    
    A1=A(xD{2:end});
   %A1=subs(A,x,xN);%ignore the time
    
   % A1=single(vpa(A1))
    
    P_=A1*P*A1.';
    
    %% Mazzoni section
     
  mazzoni;
    %eps is the minimum possible step that can be performed in matlab due to
    %floating point difference
    
function mazzoni() 
    prevA=nowA;
    nowA=A1;
    %prevA=single(vpa(subs(A,x,x_(height,:))));
    %nowA=subs(A,x,x_(height+1,:));%need to tie dt into the jacobian
     
    E=dt^2/2*((nowA-prevA)/(3*dt)-prevA.^2/6);
    %E=single(vpa(E));
    
    C = bsxfun(@rdivide, E, prediction(1:end-1));%%will use gpu for computing
   
    EMax=max(abs(C(:)));%%returns the maximum element in the E matrix
    if (t~=timeMes && first==false)%MY STUFF FOR REAL-TIME
        timeTaken=toc(startInterval);
        %dtLog(dtCount)=dt;
        %dtCount=dtCount+1;
        %if (dtCount>5)
         %   dtCount=1;
        %end
   
        if (prevT>WCET)
            WCET=prevT
        end
       
       CalcTime=timeTaken+((timeMes-timeTaken)/WCET-1)*dt;
        
       %maxTime=true;
       if (CalcTime<timeMes)
            alpha=timeMes/CalcTime;
            e=alpha^2*e;
       end
        
     end
    
    Qcoeff=B*sqrt(e./EMax);
    dt=dt*Qcoeff;
    if (dt<tMin)
        dt=tMin;
    end
    
    
   %classUnderlying(dt)
   
end 
end

%% Correct function
function correct
   % disp("Correction phase")
    
    %H is typically found by taking the jacobian of the measurement
    %equation at point x
    %compute the kalman gain;
    PHt=P_*H';
    %% Normal Kalman calc
    %K=PHt*inv(H*PHt+R);
    
    R=H*PHt+R;
   L=chol(R);
   U=L\eye(1);%need to check if eye(1) is fine at some point
   R1=L'\U;
   K=PHt*R1;
   
   

  
    %update estimate with measurement
    currentEstimate=x_(height,1:end-1);
    x_(height,1:end-1)=currentEstimate+K'.*((z-currentEstimate));%musn't update the timestamp
   % x_(height,:)=x_(height,:)+(K'.*(z-x_(height,:)));%x being the current estimate
    
    %Update the error covariance
 
    P =(I-K*H)*P_;
 
    
end



end


