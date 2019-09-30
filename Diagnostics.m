N=1000;
times=zeros(5,2);
H=rand(N);
P_=rand(N);
R=rand(N);
Sigma=[10 10; 5 5];
df=2;
W=wishrnd(Sigma,df)
for i=1:5
    f=@()( cho(W));
  
   times(i,2)=gputimeit(f,1);
  f=@()(in(W));
   times(i,1)=timeit(f);
   
end
plot(times(:,1))
hold on
plot(times(:,2))
hold off
times


function A=cho(W)
    for i=1:10
        L=chol(W);
        U=L\eye(1);
        A=L'\U;
    end
end

function A=in(W)
     for i=1:10
       A=inv(W);
    end

end
