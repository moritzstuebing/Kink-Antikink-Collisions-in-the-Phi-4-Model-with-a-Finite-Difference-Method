clc, clear all, clearvars;

%variables

v = input('Give an initial velocity:'); %Initial velocity input
a = 7; %starting distance from origin
gamma = 1/((1-v^2)^(1/2));

L = 150; %Length of x-axis either side of origin
N = 30000; %Number of spatial points
hx = (2*L)/N;
x = linspace(-L,L,N); %Vector for x

tim = 200; %Time limit
K = 25000; %Number of time points
ht = tim/K;
t = linspace(0,tim,K); %Vector for t

%Defining u0
for i=1:N
  if x(i) < 0
    u0(i) = tanh(x(i) + a);
  elseif (x(i) >= 0)
    u0(i) = -tanh(x(i) - a);
  end
end

%Defining u1
for i=1:N
  if x(i) < 0
    u1(i) = -v*gamma*(sech(gamma*(x(i) + a)))^2;
  elseif (x(i) >= 0)
    u1(i) = -v*gamma*(sech(gamma*(x(i) - a)))^2;
  end
end

%Matrix to store solutions
U = zeros(N,K);
U(:,1) = u0;
U(:,2) = u0 + ht*u1;

%Boundary conditions
U(1,:) = -1;
U(N,:) = -1;

%Main loop
for k = 2:K-1
    for j = 2:N-1
        U(j, k+1) = ((ht^2)/(hx^2))*(U(j+1,k) - 2*U(j,k) + U(j-1,k)) - (ht^2)*(-2*U(j,k) + 2*U(j,k).^3) + 2*U(j,k) - U(j,k-1);
    end
end
