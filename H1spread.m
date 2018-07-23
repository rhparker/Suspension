function sp = H1spread(x,U,D)

norms = H1norm(x,U,D);

diffs = abs( norms(2:end) - norms(1:end-1) );

sp = var(diffs);

end