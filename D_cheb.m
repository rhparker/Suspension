%% Chebyshev differentiation matrices
% for interval [-L, L]

function [D, x] = D_cheb(N, L, num, config)

% generate Chebyshev matrices for num derivatives
[x, DM] = chebdif(N, num);

% size of final matrices; may be smaller if we have
% Dirichlet BCs
final_size = N;

% % deal with Dirichlet BCs by zeroing out
% if isfield(config, 'Dirichlet')
%     % recall x for Chebyshev is written backwards
%     if strcmp(config.Dirichlet, 'R') 
%         DM(1, 1:end, :) = zeros(1, N, num);
%         DM(1:end, 1, :) = zeros(N, 1, num);
%     elseif strcmp(config.Dirichlet, 'L')
%         DM(end, 1:end, :) = zeros(1, N, num);
%         DM(1:end, end, :) = zeros(N, 1, num);
%     elseif strcmp(config.Dirichlet, 'LR')
%         DM(1, 1:end, :) = zeros(1, N, num);
%         DM(1:end, 1, :) = zeros(N, 1, num);
%         DM(end, 1:end, :) = zeros(1, N, num);
%         DM(1:end, end, :) = zeros(N, 1, num);
%     end
% end

% deal with Dirichlet BCs by removing points
if isfield(config, 'Dirichlet')
    % R side; remove first point
    % recall x for Chebyshev is written backwards
    if strcmp(config.Dirichlet, 'R') 
        DM = DM(2:end,2:end,:);
        x = x(2:end);
        final_size = final_size - 1;
    % L side; remove last point
    elseif strcmp(config.Dirichlet, 'L')
        DM = DM(1:end-1,1:end-1,:);
        x = x(1:end-1);
        final_size = final_size - 1;
    % both sides; remove first and last point
    elseif strcmp(config.Dirichlet, 'LR')
        DM = DM(2:end-1,2:end-1,:);
        x = x(2:end-1);
        final_size = final_size - 2;
    end
end

% scale for length
if isfield(config,'half_line')
    % scale factor for [0, L]
    scale = 2/L;
    x = (L/2).*(x+1);
else
    % scale factor for [-L, L]
    scale = 1/L;
    x = L.*x;
end

D = zeros([final_size, final_size, num]);

% scale matrices appropriately
for index = 1:num
    M = DM(:,:,index);
    % scale to interval [-L, L]
    M = M * (scale^index);
    % put into array D
    D(:,:,index) = M;
end

% deal with Neumann BCs, if any

if isfield(config, 'Neumann')
    % Neumann on R
    if strcmp(config.Neumann,'R')
        M1 = diag(L - x);
        M2 = diag(1./(L - x));
        D1 = M1*D(:,:,1)*M2 - M2;
        D2 = ( M1*D(:,:,2) - 2*D(:,:,1) )*M2;
        D3 = ( M1*D(:,:,3) - 3*D(:,:,2) )*M2;
        D(:,:,1) = D1;
        D(:,:,2) = D2;
        D(:,:,3) = D3;
    end
    
    % Neumann on R
    if strcmp(config.Neumann,'L')
        M1 = diag(L + x);
        M2 = diag(1./(L + x));
        D1 = M1*D(:,:,1)*M2 + M2;
        D2 = ( M1*D(:,:,2) + 2*D(:,:,1) )*M2;
        D3 = ( M1*D(:,:,3) + 3*D(:,:,2) )*M2;
        D(:,:,1) = D1;
        D(:,:,2) = D2;
        D(:,:,3) = D3;
    end
    
    % Neumann on both sides
    if strcmp(config.Neumann,'LR')
        % interpolating polynomial to be of form
        % p(x) = (L-x)(L+x) q(x)
        M1 = diag(L^2 - x.^2);
        M2 = diag( 1./ (L^2 - x.^2) );
        M3 = diag( x./ (L^2 - x.^2) );
        M4 = diag( x ); 
        D1 = M1*D(:,:,1)*M2 -  2*M3;
        D2 = (M1*D(:,:,2) -  4*M4*D(:,:,1) -  2)*M2;
        D3 = (M1*D(:,:,3) -  6*M4*D(:,:,2) -  6*D(:,:,1))*M2;
    end
end


end