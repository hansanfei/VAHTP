clear all;
close all;
clc;

%% Problem dimension
M = 128;                % Row number of the dictionary matrix
N = M* 2;               % Column number of the dictionary matrix
L = 5;                  % Number of measurement vectors
K = 30;                 % Number of nonzero rows (i.e. source number) in the solution matrix
Phi = randn(M,N);
Phi = Phi./(ones(M,1)*sqrt(sum(Phi.^2)));
beta = ones(K,1)*0.9;    
% Temporal correlation of each nonzero row in the solution matrix.
% Generate the K nonzero rows, each row being an AR(1) process. All the AR(1)
% processes have different AR coefficients, which are randomly chosen from [0.7,1)
nonzeroW=zeros(K,L*100);

MaxNbIter=100;
tolRes=1e-5;

%     nonzeroW(:,1) = randn(K,1);
nonzeroW(:,1) = ones(K,1);
for ii = 2 : L*100
    nonzeroW(:,ii) = beta .* nonzeroW(:,ii-1) + sqrt(1-beta).*(ones(K,1).*randn(K,1));
end
nonzeroW = nonzeroW(:,end-L+1:end);   % Ensure the AR processes are stable

% Normalize each row
nonzeroW = nonzeroW./( sqrt(sum(nonzeroW.^2,2)) * ones(1,L) );

% Locations of nonzero rows are randomly chosen
ind = randperm(N);
indice = ind(1:K);
Wgen = zeros(N,L);
Wgen(indice,:) = nonzeroW;
% Noiseless signal

signal = Phi * Wgen;

%% VAHTP algorithm
[X_VAHTP] = VAHTP(signal,Phi,MaxNbIter,tolRes);
%    [X_mhtp,~,~] = GMHTP(signal,Phi,MaxNbIter,tolRes,Eps);
Mse_mhtp = norm(X_VAHTP - Wgen,'fro')/norm(Wgen,'fro');
fprintf('Error is %f\n',Mse_mhtp )

x_coordinate=1:N;

%% figures of the first and third snapshots of the recovered signal, you can try other snapshots

figure;
stem(x_coordinate,Wgen(:,1),'kO','LineWidth',2)
hold on;
% stem(x_coordinate,z_tailL1,'r','LineWidth',2)
stem(x_coordinate,X_VAHTP(:,1),'r*','LineWidth',1)
title('1th snapshots')
legend('real','recovery by VAHTP');

figure;
stem(x_coordinate,Wgen(:,3),'kO','LineWidth',2)
hold on;
% stem(x_coordinate,z_tailL1,'r','LineWidth',2)
stem(x_coordinate,X_VAHTP(:,3),'r*','LineWidth',1)
title('3th snapshots')
legend('real','recovery by VAHTP');


