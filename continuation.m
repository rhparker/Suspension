%% continuation code for wave speed c

% grid for x
L = 20;
N = 512;

% starting parameters
par.c = 0;
par.b = -2;

% method configuration
%  Fourier

config.method = 'Fourier';
config.BC = 'periodic';
config.degree = 3;
xout = linspace(-L, L, N+1)';
xout = xout(1:end-1);
D = D_fourier(N, L, config.degree);

% % Chebyshev
% config.method = 'Chebyshev';
% config.degree = 3;
% config.Dirichlet = 'LR';
% % config.Neumann = 'LR';
% config.num_Dirichlet = 2;
% [D, xout] = D_cheb(N, L, config.degree, config);

u = lefton(xout, 0, 20, par);
fval = Bfamily(xout, u, par, D, config);


%% secant continuation code in parameter c

% number of iterations
iterations = 400;

% continuation parameters
contPar.numContSteps    = iterations;
contPar.Name            = 'c';      % continuation parameter continued in

contPar.ds          = 0.03;       % continuation step size: should always be positive!
contPar.initial_ds  = 0.03;       % initial step: sign determines direction of continuation

% system parameters
% initial condition is our exact solution from above
u0 = u;

%% find two initial points on the bifurcation curve

% first point
options = optimset('Display','iter','Algorithm','levenberg-marquardt','MaxIter',500,'Jacobian','on');
options.TolFun = 1e-12;
options.TolX = 1e-12;

[u1,fval,exitflag,output,jacobian1]  = fsolve( @(u) Bfamily_int(xout,u,par,D,config), u0, options);

v0 = [u1; getfield(par,contPar.Name)];  % v0 is the first point with parameter name

parameter = getfield(par,contPar.Name); % Set a few facility vectors
normu     = norm(u1);
contdata  = v0;

% second point
par = setfield(par,contPar.Name,getfield(par,contPar.Name)+contPar.initial_ds); % increase par by ds

[u2,fval,exitflag,output,jacobian1]  = fsolve( @(u) Bfamily_int(xout,u,par,D,config), u1, options);
    
v1 = [u2; getfield(par,contPar.Name)]; % v1 is the second point with parameter name

parameter = [parameter getfield(par,contPar.Name)];
normu     = [normu     norm(u2)];
contdata  = [contdata  v1];

%% Continuation code
% At each continuation step
for index = 1:contPar.numContSteps

  % Predictor
  v = v1 + (v1 - v0)/norm(v1 - v0, 2) * contPar.ds;
  
  disp(['Predictor = ',num2str(v(end))]);
  % Call fsolve predictor/corrector function
  
  [v,res,exitflag,output,jacobian2] = fsolve(@(v) FixedPointSecantPredictorCorrector(xout,v,v1,v0,@Bfamily_int,par,D,config,contPar),v,options); 
 
  disp(['Step = ',int2str(index),' Parameter = ',num2str(v(end)), ' Norm(residual,inf) = ',num2str(norm(res,inf))]);
  % Update output parameters
  parameter = [parameter v(end)];
%  normu     = [normu norm(v(1:end-1))];
  normu     = [normu norm(v(1:end-1))];
  contdata  = [contdata v];

  plot(parameter,normu,'-o');   % plot the bifurcation diagram on the fly
  xlabel(contPar.Name);
  ylabel('Norm of the solution');drawnow;

  % Prepare for the next continuation step
  v0 = v1;
  v1 = v;
  
  % Save intermediate data so we don't lose it
  if mod(index, 500) == 0
      uc = contdata;
    save uc_out x uc config;
  end
end

% save data at end
uc = contdata;
save uc_out xout uc config;





