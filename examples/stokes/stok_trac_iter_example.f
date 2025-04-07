      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      
      real *8 ts(2), rres
      real *8, allocatable :: rfacs(:,:), errs(:)
      character *100 fname
      integer ipars(2)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ixyzso(:),nfars(:)

      integer, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3,10),xyz_in(3,10),stracmat(3,3),smat(3,3),
     1     dmat(3,3)
      real *8 xyz_src(3,10),xyz_targ(3,10)
      real *8 velgrad(3,3), vel(3), pre, tractemp(3)
      real *8 sigsrc(3,10),uin(3),utargtest(3,10),dpars(2),st1(3),du1(3)
      real *8 usol(3,10), uneu(3,10), uavecomp(3), uavetest(3)
      real *8 st2(3), du2(3), uconst(3)
      real *8 v(3), omega(3), r0(3), udiff(3,10), udiff2(3,10)      
      real *8, allocatable :: uval(:,:), tracval(:,:), soln(:,:)
      real *8 :: xtrans(3), wcross(3), vtmp(3)
      complex * 16 zpars

      call prini(6,13)

      done = 1
      pi = atan(done)*4

      nsrctarg = 8

c
c       select geometry type
c       igeomtype = 1 => sphere
c       igeomtype = 2 => stellarator
c 
      igeomtype = 1
      if(igeomtype.eq.1) ipars(1) = 1
      if(igeomtype.eq.2) ipars(1) = 5*2

      if(igeomtype.eq.1) then
        npatches = 12*(4**ipars(1))
      endif
      if(igeomtype.eq.2) then
        ipars(2) = ipars(1)*3
        npatches = 2*ipars(1)*ipars(2)
      endif


      if(igeomtype.eq.1) then

         xyz_out(1,1) = 3.17d0
         xyz_out(2,1) = -0.03d0
         xyz_out(3,1) = 3.15d0

         xyz_in(1,1) = 0.17d0
         xyz_in(2,1) = 0.23d0
         xyz_in(3,1) = -0.11d0

         do j = 2,nsrctarg
            dr = 4 + cos(5d0*j)
            th = 20d0*j
            ph = 10d0*j
            xyz_out(1,j) = dr*sin(ph)*cos(th)
            xyz_out(2,j) = dr*sin(ph)*sin(th)
            xyz_out(3,j) = dr*cos(ph)

            dr = 0.25 + 0.125*cos(3d0*j)
            th = 19d0*j
            ph = 13d0*j
            xyz_in(1,j) = dr*sin(ph)*cos(th)
            xyz_in(2,j) = dr*sin(ph)*sin(th)
            xyz_in(3,j) = dr*cos(ph)
         enddo
      endif

      if(igeomtype.eq.2) then
        xyz_in(1,1) = -4.501d0
        xyz_in(2,1) = 1.7d-3
        xyz_in(3,1) = 0.00001d0

        xyz_out(1,1) = -3.5d0
        xyz_out(2,1) = 3.1d0
        xyz_out(3,1) = 20.1d0
      endif
      ifinout = 0

      if(ifinout.eq.0) then

         do j = 1,nsrctarg
            xyz_src(1,j) = xyz_out(1,j) 
            xyz_src(2,j) = xyz_out(2,j) 
            xyz_src(3,j) = xyz_out(3,j)

            xyz_targ(1,j) = xyz_in(1,j)
            xyz_targ(2,j) = xyz_in(2,j)
            xyz_targ(3,j) = xyz_in(3,j)
         enddo
      else
         do j = 1,nsrctarg
            xyz_src(1,j) = xyz_in(1,j) 
            xyz_src(2,j) = xyz_in(2,j) 
            xyz_src(3,j) = xyz_in(3,j)

            xyz_targ(1,j) = xyz_out(1,j)
            xyz_targ(2,j) = xyz_out(2,j)
            xyz_targ(3,j) = xyz_out(3,j)
         enddo
      endif

      norder = 8
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(targs(3,npts))
      ifplot = 0

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
      allocate(ixyzso(npatches+1),nfars(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      allocate(uval(3,npts),tracval(3,npts))

      sigsrc(1,1) = 1.1d0
      sigsrc(2,1) = -0.27d0
      sigsrc(3,1) = .31d0
      do j = 2,nsrctarg
         sigsrc(1,j) = cos(j*17d0)
         sigsrc(2,j) = sin(j*21d0)
         sigsrc(3,j) = cos(51d0*j)
      enddo

c     boundary data 
      
      do i=1,npts
         tracval(1,i) = 0
         tracval(2,i) = 0
         tracval(3,i) = 0         
         do j = 1,nsrctarg
            call st3d_strac_vec(9,xyz_src(1,j),3,srcvals(1,i),0,dpars,
     1           0,zpars,0,ipars,smat)
            tracval(1,i) = smat(1,1)*sigsrc(1,j) + smat(1,2)*sigsrc(2,j)
     1           + smat(1,3)*sigsrc(3,j) + tracval(1,i)
            tracval(2,i) = smat(2,1)*sigsrc(1,j) + smat(2,2)*sigsrc(2,j)
     1           + smat(2,3)*sigsrc(3,j) + tracval(2,i)
            tracval(3,i) = smat(3,1)*sigsrc(1,j) + smat(3,2)*sigsrc(2,j)
     1           + smat(3,3)*sigsrc(3,j) + tracval(3,i)
         enddo
      enddo

c
c     exact solution
c
      
      do i = 1,nsrctarg
         utargtest(1,i) = 0
         utargtest(2,i) = 0
         utargtest(3,i) = 0         
         do j = 1,nsrctarg
            call st3d_slp_vec(9,xyz_src(1,j),3,xyz_targ(1,i),
     1           0,dpars,0,zpars,0,ipars,smat)
            utargtest(1,i) = smat(1,1)*sigsrc(1,j)+smat(1,2)*sigsrc(2,j)
     1           + smat(1,3)*sigsrc(3,j) + utargtest(1,i)
            utargtest(2,i) = smat(2,1)*sigsrc(1,j)+smat(2,2)*sigsrc(2,j)
     1           + smat(2,3)*sigsrc(3,j) + utargtest(2,i)
            utargtest(3,i) = smat(3,1)*sigsrc(1,j)+smat(3,2)*sigsrc(2,j)
     1           + smat(3,3)*sigsrc(3,j) + utargtest(3,i)
         enddo
      enddo


      npt1 = 1
      allocate(ipatch_id(npt1),uvs_targ(2,npt1))
      do i=1,npt1
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
      
c
c     solve Neumann (traction) problem 
c

      numit = 200

      allocate(soln(3,npts),errs(numit+1))

      eps_gmres = 1d-6
      eps = 1d-6
      
      call stok_s_trac_solver(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,numit,
     2     ifinout,tracval,eps_gmres,niter,errs,rres,soln)
      
      
      call prin2('gmres errs *',errs,niter)
      call prin2('gmres rres *',rres,1)
      
      ndt_targ = 3
      dpars(1) = 1
      dpars(2) = 0
      call stok_comb_vel_eval(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndt_targ,nsrctarg,xyz_targ,
     2     ipatch_id,uvs_targ,eps,dpars,soln,usol)


c     
c     find a rigid body motion based on first 3 points of 
c     exact (solutions only determined up to a rigid body
c     rotation)
c     

      nfit = 3
      call find_rigid_body_3D(utargtest,usol,xyz_targ,nfit,
     1     xtrans,wcross)

      do j = 1,nsrctarg
         call eval_rigid_body_3D(xyz_targ(1,j),xtrans,wcross,vtmp)
         do i = 1,3
            usol(i,j) = usol(i,j) + vtmp(i)
         enddo
      enddo
      
      sum = 0
      sumrel = 0

      do j = 1,nsrctarg
         sum0 = 0
         do i = 1,3
            sum0 = sum0 + (usol(i,j)-utargtest(i,j))**2
            sumrel = sumrel + (utargtest(i,j))**2
         enddo
c         write(*,*) j, dsqrt(sum0), sqrt(sumrel)
         sum = sum + sum0
      enddo
      
      call prin2('rel err in velocity *',sqrt(sum/sumrel),1)

      stop
      end

      subroutine eval_rigid_body_3D(pt,x,w,v)
      implicit none
      real *8 :: pt(3), x(3), w(3), v(3)

      v(1) = x(1) + w(2)*pt(3) - w(3)*pt(2)
      v(2) = x(2) - w(1)*pt(3) + w(3)*pt(1)
      v(3) = x(3) + w(1)*pt(2) - w(2)*pt(1)
      
      return
      end

      subroutine find_rigid_body_3D(us,vs,pts,npts,x,w)
      implicit none
      real *8 :: us(3,npts), vs(3,npts), pts(3,npts), x(3), w(3)
      integer :: npts
c     local
      real *8, allocatable :: diffs(:), mat(:,:), work(:)
      integer :: i,j, m, n, nrhs, lda, ldb, info, lwork

      lwork = 3*npts*6
      allocate(diffs(3*npts),mat(3*npts,6),work(lwork))
      
      do j = 1,npts
         do i = 1,3
            diffs(i+(j-1)*3) = us(i,j) - vs(i,j)
         enddo
      enddo

      do j = 1,6
         do i = 1,npts*3
            mat(i,j) = 0
         enddo
      enddo

c     x gets added to both

      do j = 1,npts
      
         mat(1+3*(j-1),1) = 1
         mat(2+3*(j-1),2) = 1
         mat(3+3*(j-1),3) = 1

c     coefs of w \times pts(:,1) and pts(:,2)
      
         mat(1+3*(j-1),5) = pts(3,j)
         mat(1+3*(j-1),6) = -pts(2,j)      
         mat(2+3*(j-1),4) = -pts(3,j)
         mat(2+3*(j-1),6) = pts(1,j) 
         mat(3+3*(j-1),4) = pts(2,j)
         mat(3+3*(j-1),5) = -pts(1,j) 
      enddo
      
c     solve
      
      m = 3*npts
      n = 6
      nrhs = 1
      lda = m
      ldb = m
      info = 0
      call dgels('N',m,n,nrhs,mat,lda,diffs,ldb,work,lwork,info)
      
      x(1) = diffs(1)
      x(2) = diffs(2)
      x(3) = diffs(3)
      w(1) = diffs(4)
      w(2) = diffs(5)
      w(3) = diffs(6)

      write(*,*) 'info ', info
      
      return
      end
      
      

      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end

