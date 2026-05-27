      implicit real *8 (a-h,o-z) 
      implicit integer *8 (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      
      real *8 ts(2), rres
      real *8, allocatable :: rfacs(:,:), errs(:)
      character *100 fname
      integer *8 ipars(2)
      integer *8, allocatable :: row_ptr(:),col_ind(:)
      integer *8, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      integer *8, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer *8, allocatable :: ixyzso(:),nfars(:)

      integer *8, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3),xyz_in(3),stracmat(3,3),smat(3,3), dmat(3,3)
      real *8 xyz_src(3),xyz_targ(3)
      real *8 velgrad(3,3), vel(3), pre, tractemp(3)
      real *8 sigout(3), uin(3), uintest(3), dpars(2), st1(3), du1(3)
      real *8 udir(3), uneu(3,10), uavecomp(3), uavetest(3)
      real *8 st2(3), du2(3), uconst(3)
      real *8 v(3), omega(3), r0(3), udiff(3,10), udiff2(3,10)      
      real *8, allocatable :: uval(:,:), tracval(:,:), soln(:,:)
      complex * 16 zpars
      integer *8 int8_0,int8_3,int8_9
      character *100, igeom
      real *8 c0(3)
      integer *8 nuv(2)
      
      
      call prini(6,13)

      int8_0 = 0
      int8_1 = 1
      int8_3 = 3
      int8_9 = 9
      done = 1
      pi = atan(done)*4

      norder = 4
c
c     patch type, iptype0 = 1, for triangles with RV nodes
c     = 11, for quadrangles with GL nodes
c     = 12, for quadrangles with Cheb nodes      
c     
      iptype0 = 11

      igeom = 'sphere'
      if (trim(igeom).eq.'sphere') then
         a = 1
         na = 3
         c0(1) = 0
         c0(2) = 0
         c0(3) = 0
         call get_sphere_npat_mem(a, na, c0, norder, iptype0, 
     1        npatches, npts)

         allocate(srcvals(12,npts), srccoefs(9,npts))
         allocate(ixyzs(npatches+1), iptype(npatches))
         allocate(norders(npatches))

         call get_sphere_npat(a, na, c0, norder, iptype0,
     1        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

         xyz_out(1) = 3.17d0
         xyz_out(2) = -0.03d0
         xyz_out(3) = 3.15d0

         xyz_in(1) = 0.17d0
         xyz_in(2) = 0.23d0
         xyz_in(3) = -0.11d0
      elseif (trim(igeom).eq.'stellarator') then

         nuv(1) = 10
         nuv(2) = 30
         
         call get_stellarator_npat_mem(nuv, norder, iptype0, 
     1        npatches, npts)

         allocate(srcvals(12,npts), srccoefs(9,npts))
         allocate(ixyzs(npatches+1), iptype(npatches))
         allocate(norders(npatches))

         call get_stellarator_npat(nuv, norder, iptype0,
     1        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)


         xyz_in(1) = -4.501d0
         xyz_in(2) = 1.7d-3
         xyz_in(3) = 0.00001d0

         xyz_out(1) = -3.5d0
         xyz_out(2) = 3.1d0
         xyz_out(3) = 20.1d0
      endif

      ifinout = 0
      
      if(ifinout.eq.0) then
        xyz_src(1) = xyz_out(1) 
        xyz_src(2) = xyz_out(2) 
        xyz_src(3) = xyz_out(3)

        xyz_targ(1) = xyz_in(1)
        xyz_targ(2) = xyz_in(2)
        xyz_targ(3) = xyz_in(3)
      else
        xyz_src(1) = xyz_in(1) 
        xyz_src(2) = xyz_in(2) 
        xyz_src(3) = xyz_in(3)

        xyz_targ(1) = xyz_out(1)
        xyz_targ(2) = xyz_out(2)
        xyz_targ(3) = xyz_out(3)
      endif

      ifplot = 0

      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      allocate(uval(3,npts),tracval(3,npts))

      sigout(1) = 1.1d0
      sigout(2) = -0.27d0
      sigout(3) = .31d0
      
      do i=1,npts
         call st3d_slp_vec(xyz_src,int8_3,srcvals(1,i),
     1        int8_0,dpars,int8_0,zpars,int8_0,ipars,smat)
         uval(1,i) = smat(1,1)*sigout(1) + smat(1,2)*sigout(2)
     1        + smat(1,3)*sigout(3)
         uval(2,i) = smat(2,1)*sigout(1) + smat(2,2)*sigout(2)
     1        + smat(2,3)*sigout(3)
         uval(3,i) = smat(3,1)*sigout(1) + smat(3,2)*sigout(2)
     1        + smat(3,3)*sigout(3)
      enddo



      call st3d_slp_vec(xyz_src,int8_3,xyz_targ,int8_0,
     1     dpars,int8_0,zpars,int8_0,ipars,smat)
      
       uintest(1) = smat(1,1)*sigout(1) + smat(1,2)*sigout(2)
     1        + smat(1,3)*sigout(3)
       uintest(2) = smat(2,1)*sigout(1) + smat(2,2)*sigout(2)
     1        + smat(2,3)*sigout(3)
       uintest(3) = smat(3,1)*sigout(1) + smat(3,2)*sigout(2)
     1        + smat(3,3)*sigout(3)


      npt1 = 1
      allocate(ipatch_id(npt1),uvs_targ(2,npt1))
      do i=1,npt1
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
      
c
c     solve dirichlet (velocity) problem
c

      alpha = 1
      beta = 1
      
      dpars(1) = alpha
      dpars(2) = beta
      numit = 200

      allocate(soln(3,npts),errs(numit+1))

      eps_gmres = 1d-7
      eps = 1d-7

      write(*,*) npts, ' points...'
      
      call stok_comb_vel_solver(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,dpars,numit,
     2     ifinout,uval,eps_gmres,niter,errs,rres,soln)


      call prin2('gmres errs *',errs,niter)
      call prin2('gmres rres *',rres,1)

      ndt_in = 3
      nt_in = 1
      call stok_comb_vel_eval(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndt_in,nt_in,xyz_targ,
     2     ipatch_id,uvs_targ,eps,dpars,soln,udir)

      udir(1) = udir(1) 
      udir(2) = udir(2) 
      udir(3) = udir(3) 

      sum = 0
      sumrel = 0

      do i = 1,3
         sum = sum + (udir(i)-uintest(i))**2
         sumrel = sumrel + (uintest(i))**2
      enddo
      
      call prin2('rel err in velocity *',sqrt(sum/sumrel),1)

c
c  solve scattering problem
c
      do i=1,npts
        uval(1,i) = -1
        uval(2,i) = 0
        uval(3,i) = 0
        soln(1,i) = 0
        soln(2,i) = 0
        soln(3,i) = 0
      enddo
      
      call stok_comb_vel_solver(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,dpars,numit,
     2     ifinout,uval,eps_gmres,niter,errs,rres,soln)


      call prin2('gmres errs *',errs,niter)
      call prin2('gmres rres *',rres,1)

      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,soln,'stell-stok-scat-soln-ref1.vtk','a')

      
      stop
      end
