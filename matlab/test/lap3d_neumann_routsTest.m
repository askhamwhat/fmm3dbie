%
% This file tests the Helmholtz impedance problem
%
%
run ../startup.m
S = surfer.sphere(6, 1, 2, 11);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;
[~, npts] = size(srcvals);

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;

rhs = lap3d.kern(src_info, S, 'sprime');

[sig] = lap3d.neumann.solver(S, rhs, eps);

targ_info = [];
targ_info.r = xyz_out;

pot = lap3d.neumann.eval(S, sig, eps, targ_info);
pot_ex = lap3d.kern(src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));

