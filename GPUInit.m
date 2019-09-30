n=3;      %number of state
q=0.1;    %std of process 
r=0.1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement  
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
h=@(x)x(1);                               % measurement equation
s=[0;0;1];                                % initial state
x=s+q*randn(3,1); %initial state          % initial state with noise
P = eye(n);                               % initial state covraiance
N=100;                                     % total dynamic steps
xV = zeros(n,N, 'gpuArray');          %estmate        % allocate memory
sV = zeros(n,N, 'gpuArray');          %actual
zV = zeros(1,N, 'gpuArray');

%--------------------------------------------
%Transferring to GPU - ONLY NEEDED VARIABLES
%-------------------------------------------
Q=gpuArray(Q);
R=gpuArray(R);
x=gpuArray(x);
P=gpuArray(P);

%-------------------------------------------

for k=1:N %this has been to predetermine all the measurements and true states so that a parfor can be used below
    s = f(s) + q*randn(3,1); 
    z = h(s) + r*randn; 
    zV(k)  = z; 
    sV(:,k)= s; 
end
tic();


for k=1:N 
  [x, P] = GPUekf(f,x,P,h,zV(k),Q,R);            % ekf 
  xV(:,k) = x;                            % save estimate
end
toc()
 for k=1:3                                 % plot results
  subplot(3,1,k)
  plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
 end