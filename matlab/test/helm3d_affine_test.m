Sinput = geometries.ellipsoid([1,2,3],[2,2,2],[1,2,3]);
mat = diag([1,1/2,1/3]);
shift = -mat*[1,2,3].';
 [S] = affine_transf(Sinput,mat,shift);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;

zk = 1.1;
rep_pars = [1.0, 0.0];
ndeg = 1;

jn = sqrt(pi/2/zk)*besselj(ndeg+0.5,zk);
hn = sqrt(pi/2/zk)*besselh(ndeg+0.5,1,zk);

rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
rhs = S.r(3,:)./rr;

rhs = rhs(:);
eps = 1e-7;


p = helm3d.dirichlet.eval(S,rhs,S,eps,zk,rep_pars);

p_ex = rhs*jn*hn*1j*zk;

err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential=%d\n',err1);

%% Now test eval routine with precomputed quadrature corrections
opts_quad = [];
opts_quad.format='rsc';
Q = helm3d.dirichlet.get_quadrature_correction(S, ...
   eps,zk,rep_pars,S,opts_quad);
opts = [];
opts.precomp_quadrature = Q;
p = helm3d.dirichlet.eval(S,rhs,S,eps,zk,rep_pars,S,opts);
err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential with precomp corr=%d\n',err1);


%% Now test the solver + kernel evaluation routines

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zk,src_info,S,'s');

rep_pars = [-1j*zk, 1];
sig = helm3d.dirichlet.solver(S,rhs,eps,zk,rep_pars);

targ_info = [];
targ_info.r = xyz_out;

dat = helm3d.kern(zk,S,targ_info,'c',rep_pars);

pot = dat*(sig.*wts);
pot_ex = helm3d.kern(zk,src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));



%%%%%%%%%%%%%%%%%%%% now test translate and scale


Sinput = geometries.ellipsoid([1,2,3],[2,2,2],[1,2,3]);
sfs = [1,1/2,1/3];
S = scale(Sinput,sfs);
shift = -sfs.'.*[1,2,3].';
 [S] = translate(S,shift);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;

zk = 1.1;
rep_pars = [1.0, 0.0];
ndeg = 1;

jn = sqrt(pi/2/zk)*besselj(ndeg+0.5,zk);
hn = sqrt(pi/2/zk)*besselh(ndeg+0.5,1,zk);

rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
rhs = S.r(3,:)./rr;

rhs = rhs(:);
eps = 1e-7;


p = helm3d.dirichlet.eval(S,rhs,S,eps,zk,rep_pars);

p_ex = rhs*jn*hn*1j*zk;

err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential=%d\n',err1);

%% Now test eval routine with precomputed quadrature corrections
opts_quad = [];
opts_quad.format='rsc';
Q = helm3d.dirichlet.get_quadrature_correction(S, ...
   eps,zk,rep_pars,S,opts_quad);
opts = [];
opts.precomp_quadrature = Q;
p = helm3d.dirichlet.eval(S,rhs,S,eps,zk,rep_pars,S,opts);
err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential with precomp corr=%d\n',err1);


%% Now test the solver + kernel evaluation routines

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zk,src_info,S,'s');

rep_pars = [-1j*zk, 1];
sig = helm3d.dirichlet.solver(S,rhs,eps,zk,rep_pars);

targ_info = [];
targ_info.r = xyz_out;

dat = helm3d.kern(zk,S,targ_info,'c',rep_pars);

pot = dat*(sig.*wts);
pot_ex = helm3d.kern(zk,src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));
