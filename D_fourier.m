%% Fourier differentiation matrices
% grid points: N
% interval: [-L, L]
% number of derivatives: num

function D = D_fourier(N, L, num)

% initialize array
D = zeros([N, N, num]);

% scale factor for [-L, L]
scale = pi/L;

for index = 1:num
    % generate Fourier matrix
    [~,M]  = fourdif(N,index);
    % scale to interval [-L, L]
    M = M * (scale^index);
    % put into array D
    D(:,:,index) = M;
end

end