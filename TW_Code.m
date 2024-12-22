% Number of grid points and initial domain size
N = 10000; L = 1000;

% Spatial step size
dX = L / (N-1);

% Initial spatial domain
x = linspace(0,L,N);

% Initial condition
u0 = (1+tanh(10*L*(0.02-x/L)))/2;
uInit = [u0'; ones(N,1); zeros(N,1)];

% Number of iterations/lengths to remesh the Eulerian-Lagrangian domains
numIter = 1500;
iterLength = 3;
outputsPerIter = 2;

%wave speed c
epsilon = 0.001;
c = 30;

%This is the linear domain growth rate, wich we relate to the wave speed c.
%NB: This is because L is the length of our Lagrangian domain.
r=c/L;

%Smoothed heaviside H(u-epsilon)
G = @(u)(1+tanh(10*L^2*(u-epsilon)))/2;

% Growth function S, kinetics f, and diffusion rate D
S = @(t,u)r./(1+r*t) + 0*u;
%f = @(t,u)u.*(1-u).*G(u);
f = @(t,u)u.*(1-u);
D = 1;

% Indices of variables
uI = 1:N;
muI = N+1:2*N;
xI = 2*N+1:3*N;

us = reshape(uInit(uI), 1, []);
xs = reshape(cumtrapz(ones(N,1)), 1, []) * dX;
ts = [0];
tCurrent = 0;

% Construct a sparse Jacobian
Z = sparse(zeros(N)); I = speye(N); e = ones(N,1);
lap = spdiags([e,-2*e,e],[1,0,-1],N,N);
JPattern = [lap, lap, Z; I, I, Z; I, Z, Z];
odeParams = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'JPattern', JPattern, 'MaxStep', 1);


%% Solve and iterate.
tic
for iter = 1 : numIter
    % Display iteration count.
    fprintf(['\r',repmat(' ',1,length(num2str(numIter))*2 + 3)])
    fprintf(['\r',num2str(iter),' / ',num2str(numIter)])

    % Times at which to generate the solution during this iteration.
    T = linspace(tCurrent, tCurrent+iterLength, outputsPerIter);

    % RHS of the system. Note dX changes each iteration.
    F = @(t,U)[f(t,U(uI))-S(t,U(uI)).*U(uI)+D*(1/dX)^2*Lap(U(uI),U(muI));...
        S(t,U(uI)).*U(muI); S(t,U(uI))];

    % Solve the system using a stiff solver and low tolerances.
    sol = ode15s(F, T, uInit, odeParams);
    % Evaluate the solution on T.
    U = deval(sol,T)';

    % Extract the solution.
    u = U(:,uI);
    mu = U(:,muI);

    % Form the Eulerian domain.
    x = cumtrapz(exp(U(:,xI)),2)*dX;

    % Store the output.
    us = [us; u(2:end,:)];
    xs = [xs; x(2:end,:)];
    ts = [ts; T(2:end)'];
    tCurrent = T(end);

    % Prepare for the next iteration.
    % Compute the current domain length.
    Ln = xs(end,end);
    % Compute a uniform Lagrangian domain and dX.
    X = linspace(0, Ln, N);
    dX = Ln / (N-1);


    % Construct the initial condition from the final timepoint.
    u0 = reshape(interp1(x(end,:), u(end,:), X), [], 1);
    uInit = [u0; ones(N,1); zeros(N,1)];

    %plot(u0); drawnow;

    I = find(U(1,uI)<0.5,1,'first');
    
    % Display the total time elapsed, along with the current domain length, location of the front, and far-field value of u.

    fprintf(['. Cumulative time taken: ',num2str(toc,4), '; Domain length: ', num2str(Ln,4), '; Midpoint: ', num2str(X(I),4), '; u(0.9*L) =  ', num2str(U(1,uI(round(0.9*end))),4)])


if U(1,uI(end))>0.9
    fprintf(['. Wave hit the rightmost boundary; ending simulation early.'])
    break;
end

end
fprintf('\n')

centred = 0;
plot_TW

% 1D Laplacian incorporating local volume form mu.
% This looks like: D*(D*u_x)_x, where D=1/mu.
% Division by dX does not occur here!
function u = Lap(u,mu)
% u(i-1)-2*u(i)+u(i+1)
D = 1./mu; D2 = D.^2;
u = [(D2(1)+D(1)*D(2)).*(u(2)-u(1));...
    D2(2:end-1).*(u(1:end-2)-2*u(2:end-1)+u(3:end))+D(2:end-1).*(D(3:end).*(u(3:end)-u(2:end-1)) +D(1:end-2).*(u(1:end-2)-u(2:end-1)));...
    (D2(end)+D(1)*D(2)).*(u(end-1)-u(end))]./2;
end

