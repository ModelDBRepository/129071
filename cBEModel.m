function [T,Y,mun,sigman,validPoints]=cBEModel(b,E)

% Internal parameters
probability_leak_tolerance = 1e-6; % Maximum amount of probability that we allow to leave the system of differential equations.
P=1-E;                             % Power
tspan = [0 25];
dim=1000;                          % Number of p(n,t)'s used, n=1:dim 

% Initial condition
y0 = zeros(1,dim);
y0(1) = 1;

% Terminal segment count
ngamma=1:dim;

% Define transition rate matrix ...
Rho = b.*sparse(-diag(ngamma.^P)+diag((ngamma(1:end-1)).^P,-1));
Rho(dim,dim) = 0;                                               % Gather all lost probility in remainder
Sparsity=spones(Rho);

% ...  and solve the problem using ode15s
options = odeset('Jacobian',Rho,'JPattern',Sparsity);
[T,Y]=ode15s(@RhoModel,tspan,y0, options);

% Calculate mean and variance of the number of terminal segments 
mun=ngamma*Y';
sigman=sqrt(ngamma.^2*Y'-mun.^2);

% Determine points for which we are within the probability leak tolerance
validPoints= (Y(:,end)<probability_leak_tolerance );

    function dydt=RhoModel(t,y) %#ok<INUSL>
        dydt = Rho*y;
    end

end 