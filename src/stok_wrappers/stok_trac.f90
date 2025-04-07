!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! quadrature generation, layer potential evaluation, and solver
! for the traction problem with the basic single layer representation
!

subroutine getnearquad_stok_s_trac(npatches, norders, &
     ixyzs, iptype, npts, srccoefs, srcvals, eps, &
     iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
  !
  !  This subroutine generates the near field quadrature
  !  for the representation:
  !
  !  u = S_{stok}[\sigma]  (1) 
  !
  !  and returns quantities related to evaluating the traction of u on surface
  !  at the surface discretization nodes.
  !
  !  If values at other points on the surface is desired then the 
  !  data should be reinterpolated to those points from the
  !  discretization nodes
  !
  !
  !  On imposing the boundary condition, we get the following operator
  !
  !  u = -\pm \sigma/2 +  S'_{stok}[\sigma]
  !
  !  The quadrature is computed by the following strategy
  !  targets within a sphere of radius rfac0*rs
  !  of a patch centroid is handled using adaptive integration
  !  where rs is the radius of the bounding sphere
  !  for the patch
  !  
  !  All other targets in the near field are handled via
  !  oversampled quadrature
  !
  !  The recommended parameter for rfac0 is 1.25d0
  !  
  !  Input arguments:
  !    - npatches: integer
  !        number of patches
  !    - norders: integer(npatches)
  !        order of discretization on each patch 
  !    - ixyzs: integer(npatches+1)
  !        ixyzs(i) denotes the starting location in srccoefs,
  !        and srcvals array corresponding to patch i
  !    - iptype: integer(npatches)
  !        type of patch
  !        iptype = 1, triangular patch discretized using RV nodes
  !        iptype = 11, quadrangular patch discretized with GL nodes
  !        iptype = 12, quadrangular patch discretized with Chebyshev 
  !                     nodes
  !    - npts: integer
  !        total number of discretization points on the boundary
  !    - srccoefs: real *8 (9,npts)
  !        basis expansion coefficients of xyz, dxyz/du,
  !        and dxyz/dv on each patch. 
  !        For each point 
  !          * srccoefs(1:3,i) is xyz info
  !          * srccoefs(4:6,i) is dxyz/du info
  !          * srccoefs(7:9,i) is dxyz/dv info
  !    - srcvals: real *8 (12,npts)
  !        xyz(u,v) and derivative info sampled at the 
  !        discretization nodes on the surface
  !          * srcvals(1:3,i) - xyz info
  !          * srcvals(4:6,i) - dxyz/du info
  !          * srcvals(7:9,i) - dxyz/dv info
  !          * srcvals(10:12,i) - normals info
  !    - eps: real *8
  !        precision requested
  !    - iquadtype: integer
  !        quadrature type
  !          * iquadtype = 1, use ggq for self + adaptive integration
  !            for rest
  !    - nnz: integer
  !        number of source patch-> target interactions in the near field
  !    - row_ptr: integer(npts+1)
  !        row_ptr(i) is the pointer
  !        to col_ind array where list of relevant source patches
  !        for target i start
  !    - col_ind: integer (nnz)
  !        list of source patches relevant for all targets, sorted
  !        by the target number
  !    - iquad: integer(nnz+1)
  !        location in wnear_ij array where quadrature for col_ind(i)
  !        starts for a single kernel. In this case the different kernels
  !        are matrix entries are located at (m-1)*nquad+iquad(i), where
  !        m is the kernel number
  !    - rfac0: real *8
  !        radius parameter for switching to predetermined quadarature
  !        rule        
  !    - nquad: integer
  !        number of near field entries corresponding to each source 
  !        target pair
  !
  !  Output arguments
  !    - wnear: real *8(6,nquad)
  !        The desired near field quadrature for
  !        S'_{stok}. Since the tensor is symmetric, quadratures are only
  !        generated for the upper half of the tensor
  !        * wnear(1,:) - stores the quadratures for (1,1) entry of the 
  !                       tensor 
  !        * wnear(2,:) - stores the quadratures for (1,2) entry of the 
  !                       tensor 
  !        * wnear(3,:) - stores the quadratures for (1,3) entry of the 
  !                       tensor 
  !        * wnear(4,:) - stores the quadratures for (2,2) entry of the 
  !                       tensor 
  !        * wnear(5,:) - stores the quadratures for (2,3) entry of the 
  !                       tensor 
  !        * wnear(6,:) - stores the quadratures for (3,3) entry of the 
  !                       tensor 
  !
  implicit none 
  integer, intent(in) :: npatches, norders(npatches), npts, nquad
  integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
  real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
  real *8, intent(in) :: rfac0
  integer, intent(in) :: iquadtype
  integer, intent(in) :: nnz
  integer, intent(in) :: row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
  real *8, intent(out) :: wnear(6,nquad)

  integer :: ndtarg, ntarg
  integer, allocatable :: ipatch_id(:)
  real *8, allocatable :: uvs_targ(:,:), wneartmp(:)
  integer i
  integer ndd, ndz, ndi
  real *8 :: dpars
  complex *16 :: zpars
  integer :: ipars
  procedure (), pointer :: fker
  external st3d_strac
  integer :: ijloc(2,6), ii, ipv


  ndtarg = 12
  ntarg = npts
  allocate(ipatch_id(npts), uvs_targ(2,npts))

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
  do i = 1,npts
     ipatch_id(i) = -1
     uvs_targ(1,i) = 0
     uvs_targ(2,i) = 0
  enddo
  !$OMP END PARALLEL DO

  call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
       ipatch_id, uvs_targ)
  !
  !        initialize the appropriate kernel function
  !

  ndd = 0
  ndi = 2
  ndz = 0
  fker => st3d_strac
  ipv = 1

  ijloc(1,1) = 1
  ijloc(2,1) = 1
  ijloc(1,2) = 1
  ijloc(2,2) = 2
  ijloc(1,3) = 1
  ijloc(2,3) = 3
  ijloc(1,4) = 2
  ijloc(2,4) = 2
  ijloc(1,5) = 2
  ijloc(2,5) = 3
  ijloc(1,6) = 3
  ijloc(2,6) = 3

  allocate(wneartmp(nquad))


  if(iquadtype.eq.1) then

     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)          
     do ii = 1,nquad
        wneartmp(ii) = 0
     enddo
     !$OMP END PARALLEL DO

     do i = 1,6   
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
             iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
             ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
             zpars, ndi, ijloc(1,i), nnz, row_ptr, col_ind, iquad, &
             rfac0, nquad, wneartmp)

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)          
        do ii = 1,nquad
           wnear(i,ii) = wneartmp(ii)
        enddo
        !$OMP END PARALLEL DO
     enddo
  endif


  return
end subroutine getnearquad_stok_s_trac

!
!
subroutine lpcomp_stok_s_trac_addsub(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
     ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
     novers, nptso, ixyzso, srcover, whtsover, lwork, work, ndim, &
     sigma, trac)
  !
  !  This subroutine evaluates the *traction* data corresponding to
  !  the following integral representation:
  !  
  !  u = S_{stok}[\sigma]
  !
  !  where the near field is precomputed and stored
  !  in the row sparse compressed format.
  !
  ! In particular, this routine evaluates
  !
  ! S'_{stok}[\sigma](x) + beta*n(x) \int n(y)\cdot digma dS
  !
  ! The nullspace correction (integral with beta coefficient) is optional.
  !        
  !  The fmm is used to accelerate the far-field and 
  !  near-field interactions are handled via precomputed quadrature
  !
  !  Using add and subtract - no need to call tree and set fmm parameters
  !  can directly call existing fmm library
  !
  !  Input arguments:
  !    - npatches: integer
  !        number of patches
  !    - norders: integer(npatches)
  !        order of discretization on each patch 
  !    - ixyzs: integer(npatches+1)
  !        ixyzs(i) denotes the starting location in srccoefs,
  !        and srcvals array corresponding to patch i
  !    - iptype: integer(npatches)
  !        type of patch
  !        iptype = 1, triangular patch discretized using RV nodes
  !        iptype = 11, quadrangular patch discretized with GL nodes
  !        iptype = 12, quadrangular patch discretized with Chebyshev nodes
  !    - npts: integer
  !        total number of discretization points on the boundary
  !    - srccoefs: real *8 (9,npts)
  !        basis expansion coefficients of xyz, dxyz/du,
  !        and dxyz/dv on each patch. 
  !        For each point 
  !          * srccoefs(1:3,i) is xyz info
  !          * srccoefs(4:6,i) is dxyz/du info
  !          * srccoefs(7:9,i) is dxyz/dv info
  !    - srcvals: real *8 (12,npts)
  !        xyz(u,v) and derivative info sampled at the 
  !        discretization nodes on the surface
  !          * srcvals(1:3,i) - xyz info
  !          * srcvals(4:6,i) - dxyz/du info
  !          * srcvals(7:9,i) - dxyz/dv info
  !          * srcvals(10:12,i) - normals info
  !    - eps: real *8
  !        precision requested
  !    - ndd: integer
  !        number of real parameters defining the kernel/
  !        integral representation. should have ndd=0 or 1
  !    - dpars: real *8(ndd)
  !        if ndd=1, dpars(1) is the parameter beta above
  !    - ndz: integer
  !        number of complex parameters defining the kernel/
  !        integral representation (unused in this routine)
  !    - zpars: complex *16(ndz)
  !        complex parameters defining the kernel/
  !        integral represnetation (unused in this routine)
  !    - ndi: integer
  !        number of integer parameters defining the kernel/
  !        integral representation (unused in this routine)
  !    - ipars: integer(ndi)
  !        integer parameters defining the kernel/
  !        integral represnetation (unused in this routine)
  !    - nnz: integer
  !        number of source patch-> target interactions in the near field
  !    - row_ptr: integer(npts+1)
  !        row_ptr(i) is the pointer
  !        to col_ind array where list of relevant source patches
  !        for target i start
  !    - col_ind: integer (nnz)
  !        list of source patches relevant for all targets, sorted
  !        by the target number
  !    - iquad: integer(nnz+1)
  !        location in wnear_ij array where quadrature for col_ind(i)
  !        starts for a single kernel. In this case the different kernels
  !        are matrix entries are located at (m-1)*nquad+iquad(i), where
  !        m is the kernel number
  !    - nquad: integer
  !        number of near field entries corresponding to each source target
  !        pair
  !    - nker: integer
  !        number of kernels in quadrature correction, must be 6
  !    - wnear: real *8(nker, nquad)
  !        Precomputed near field quadrature for
  !        \alpha S_{stok} + \beta D_{stok}. Since the Stokes 
  !        tensor is symmetric, quadratures are only
  !        generated for the lower half of the tensor
  !        * wnear(1,:) - stores the quadratures for (1,1) entry of the 
  !                       tensor 
  !        * wnear(2,:) - stores the quadratures for (1,2) entry of the 
  !                       tensor 
  !        * wnear(3,:) - stores the quadratures for (1,3) entry of the 
  !                       tensor 
  !        * wnear(4,:) - stores the quadratures for (2,2) entry of the 
  !                       tensor 
  !        * wnear(5,:) - stores the quadratures for (2,3) entry of the 
  !                       tensor 
  !        * wnear(6,:) - stores the quadratures for (3,3) entry of the 
  !                       tensor 
  !    - novers: integer(npatches)
  !        order of discretization for oversampled sources and
  !        density
  !    - ixyzso: integer(npatches+1)
  !        ixyzso(i) denotes the starting location in srcover,
  !        corresponding to patch i
  !    - nptso: integer
  !        total number of oversampled points
  !    - srcover: real *8 (12,nptso)
  !        oversampled set of source information
  !    - whtsover: real *8 (nptso)
  !        smooth quadrature weights at oversampled nodes
  !    - lwork: integer
  !        size of work array (must be npts)
  !    - work: real *8(lwork)
  !        work array
  !        * work(1:npts) stores the smooth quadrature weights for
  !          integration on the discretized surface
  !    - ndim: integer
  !        number of densities per point on the surface,
  !        must be 3 for this routine
  !    - sigma: real *8(3,npts)
  !        The density sigma above                                        
  !
  !  Output arguments:
  !    - trac: real *8(3,npts)
  !        traction corresponding to representation
  !
  
  
  implicit none
  integer, intent(in) :: npatches, npts
  integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches)
  real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)

  real *8, intent(in) :: eps
      
  integer, intent(in) :: ndd, ndz, ndi
  real *8, intent(in) :: dpars(ndd)
  complex *16, intent(in) :: zpars(ndz)
  integer, intent(in) :: ipars(ndi)

  integer, intent(in) :: nnz, nquad
  integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
  integer, intent(in) :: iquad(nnz+1)

  integer, intent(in) :: nker
  real *8, intent(in) :: wnear(nker,nquad)

  integer, intent(in) :: nptso
  integer, intent(in) :: ixyzso(npatches+1), novers(npatches)
  real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)

  integer, intent(in) :: lwork
  real *8, intent(in) :: work(lwork)

  integer, intent(in) :: ndim

  real *8, intent(in) :: sigma(ndim,npts)

  real *8, intent(out) :: trac(ndim,npts)
  
  integer norder,npols,nover,npolso
  real *8, allocatable :: potsort(:)

  real *8, allocatable :: sources(:,:), targvals(:,:)
  real *8, allocatable :: sigmaover(:,:)
  real *8, allocatable :: stoklet(:,:), strslet(:,:)
  real *8, allocatable :: pottarg(:,:), strsvec(:,:)
  real *8, allocatable :: gradtarg(:,:,:), presstarg(:)  
  integer ns, nt
  real *8 alpha, beta
  integer ifstoklet, ifstrslet
  integer ifppreg, ifppregtarg
  real *8 tmp(10), val
  real *8 over4pi

  real *8 w11, w12, w13, w21, w22, w23, w31, w32, w33
  real *8 sig1, sig2, sig3
  real *8 dn1, dn2, dn3

  integer i, j, jpatch, jquadstart, jstart

  real *8 pottmp
  real *8, allocatable :: sttmp2(:,:), strstmp2(:,:), strsvec2(:,:)

  integer nmax
  real *8 timeinfo(10), t1, t2, omp_get_wtime

  real *8, allocatable :: srctmp2(:,:), tracdirect(:,:)
  real *8 thresh, ra
  real *8 rr, rmin
  integer nss, ii, l, npover

  integer nd, ntarg0, ntarg
  integer ier, iper

  integer ndtmp, ndsigma, ndfmm

  real *8 ttot, done, pi

  real *8 potv(3), pv, gradv(3,3)
  real *8 sigma_sum, smat(3,3)

  parameter (nd=1,ntarg0=1)
  data over4pi/0.07957747154594767d0/
      
  ns = nptso
  done = 1
  pi = atan(done)*4

  ntarg = npts
  
  !
  !    estimate max number of sources in near field of 
  !    any target
  !
  
  
  nmax = 0
  call get_near_corr_max(ntarg, row_ptr, nnz, col_ind, npatches, &
       ixyzso, nmax)
  
  allocate(srctmp2(3,nmax), sttmp2(3,nmax), strstmp2(3,nmax))
  allocate(strsvec2(3,nmax))
  
  ifppreg = 0
  ifppregtarg = 3
  allocate(sources(3,ns), targvals(3,ntarg))
  allocate(stoklet(3,ns))
  allocate(sigmaover(3,ns))
  allocate(pottarg(3,ntarg),gradtarg(3,3,ntarg),presstarg(ntarg))

  ! 
  !       oversample density
  !
  
  ndsigma = 3
  call oversample_fun_surf(ndsigma, npatches, norders, ixyzs, &
       iptype, npts, sigma, novers, ixyzso, ns, sigmaover)


  !
  !      set relevant parameters for the fmm
  !
  
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
  do i=1,ns
     sources(1,i) = srcover(1,i)
     sources(2,i) = srcover(2,i)
     sources(3,i) = srcover(3,i)
     
     stoklet(1,i) = sigmaover(1,i)*whtsover(i)*over4pi
     stoklet(2,i) = sigmaover(2,i)*whtsover(i)*over4pi
     stoklet(3,i) = sigmaover(3,i)*whtsover(i)*over4pi
  enddo
  !$OMP END PARALLEL DO      
      
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
  do i = 1,ntarg
     targvals(1,i) = srcvals(1,i)
     targvals(2,i) = srcvals(2,i)
     targvals(3,i) = srcvals(3,i)
  enddo
  !$OMP END PARALLEL DO      


  ifstoklet = 1
  ifstrslet = 0

  iper = 0
  ier = 0
  
  !
  !
  !       call the fmm
  !
  
  call cpu_time(t1)
  !$      t1 = omp_get_wtime()
  ndfmm = 1
  call stfmm3d(ndfmm, eps, ns, sources, ifstoklet, stoklet, &
       ifstrslet, tmp, tmp, ifppreg, tmp, tmp, tmp, ntarg, &
       targvals, ifppregtarg, pottarg, presstarg, gradtarg, ier)
  call cpu_time(t2)
  !$      t2 = omp_get_wtime()

  timeinfo(1) = t2-t1
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
  do i = 1,ntarg
     dn1 = srcvals(10,i)
     dn2 = srcvals(11,i)
     dn3 = srcvals(12,i)     
     trac(1,i) = -dn1*presstarg(i) + &
          (dn1*gradtarg(1,1,i)+dn2*gradtarg(2,1,i)+dn3*gradtarg(3,1,i)) + &
          (dn1*gradtarg(1,1,i)+dn2*gradtarg(1,2,i)+dn3*gradtarg(1,3,i))
     trac(2,i) = -dn2*presstarg(i) + &
          (dn1*gradtarg(1,2,i)+dn2*gradtarg(2,2,i)+dn3*gradtarg(3,2,i)) + &
          (dn1*gradtarg(2,1,i)+dn2*gradtarg(2,2,i)+dn3*gradtarg(2,3,i))
     trac(3,i) = -dn3*presstarg(i) + &
          (dn1*gradtarg(1,3,i)+dn2*gradtarg(2,3,i)+dn3*gradtarg(3,3,i)) + &
          (dn1*gradtarg(3,1,i)+dn2*gradtarg(3,2,i)+dn3*gradtarg(3,3,i))
  enddo
  !$OMP END PARALLEL DO 

  
  !
  !        compute threshold for ignoring local computation
  !
  ndtmp = 3 
  call get_fmm_thresh(ndtmp, ns, sources, ndtmp, ntarg, targvals, &
       thresh)
  !
  !       add in precomputed quadrature
  !

  call cpu_time(t1)
  !$      t1 = omp_get_wtime()

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j, jpatch, jquadstart) &
  !$OMP PRIVATE(jstart, npols, l, sig1, sig2, sig3, w11, w12, w13, w21) &
  !$OMP PRIVATE(w22, w23, w31, w32, w33)
  do i = 1,ntarg
     do j = row_ptr(i),row_ptr(i+1)-1
        jpatch = col_ind(j)
        npols = ixyzs(jpatch+1) - ixyzs(jpatch)
        jquadstart = iquad(j)
        jstart = ixyzs(jpatch) 
        do l=1,npols
           sig1 = sigma(1,jstart+l-1)
           sig2 = sigma(2,jstart+l-1)
           sig3 = sigma(3,jstart+l-1)

           w11 = wnear(1,jquadstart+l-1)
           w12 = wnear(2,jquadstart+l-1)
           w13 = wnear(3,jquadstart+l-1)

           w21 = w12 
           w22 = wnear(4,jquadstart+l-1)
           w23 = wnear(5,jquadstart+l-1)

           w31 = w13 
           w32 = w23 
           w33 = wnear(6,jquadstart+l-1)


           trac(1,i) = trac(1,i) + w11*sig1 + w12*sig2 + w13*sig3 
           trac(2,i) = trac(2,i) + w21*sig1 + w22*sig2 + w23*sig3 
           trac(3,i) = trac(3,i) + w31*sig1 + w32*sig2 + w33*sig3 
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, srctmp2) &
  !$OMP PRIVATE(sttmp2, strstmp2, strsvec2, nss, l, jstart, ii, npover) &
  !$OMP PRIVATE(potv, pv, gradv, smat)
  do i = 1,ntarg
     nss = 0
     do j = row_ptr(i), row_ptr(i+1)-1
        jpatch = col_ind(j)
        do l = ixyzso(jpatch),ixyzso(jpatch+1)-1
           nss = nss + 1
           srctmp2(1,nss) = srcover(1,l)
           srctmp2(2,nss) = srcover(2,l)
           srctmp2(3,nss) = srcover(3,l)

           sttmp2(1,nss) = stoklet(1,l)
           sttmp2(2,nss) = stoklet(2,l)
           sttmp2(3,nss) = stoklet(3,l)
        enddo
     enddo

     potv(1:3) = 0
     pv = 0
     gradv(1:3,1:3) = 0

     call st3ddirectstokg(nd, srctmp2, sttmp2, &
          nss, targvals(1,i), ntarg0, potv, pv, &
          gradv, thresh)

     dn1 = srcvals(10,i)
     dn2 = srcvals(11,i)
     dn3 = srcvals(12,i)     
     
     trac(1,i) = trac(1,i) + dn1*pv - &
          (dn1*gradv(1,1)+dn2*gradv(2,1)+dn3*gradv(3,1)) - &
          (dn1*gradv(1,1)+dn2*gradv(1,2)+dn3*gradv(1,3))
     trac(2,i) = trac(2,i) + dn2*pv - &
          (dn1*gradv(1,2)+dn2*gradv(2,2)+dn3*gradv(3,2)) - &
          (dn1*gradv(2,1)+dn2*gradv(2,2)+dn3*gradv(2,3))
     trac(3,i) = trac(3,i) + dn3*pv - &
          (dn1*gradv(1,3)+dn2*gradv(2,3)+dn3*gradv(3,3)) - &
          (dn1*gradv(3,1)+dn2*gradv(3,2)+dn3*gradv(3,3))

  enddo

  call cpu_time(t2)
  !$      t2 = omp_get_wtime()

  ! add \int n_t \times n_s sigma nullspace correction

  if (ndd .eq. 1) then
     if (abs(dpars(1)) .gt. 0d0) then
        sigma_sum = 0
        
        do i = 1,npts
           sigma_sum = sigma_sum + work(i)*(srcvals(10,i)*sigma(1,i) + &
                srcvals(11,i)*sigma(2,i) + srcvals(12,i)*sigma(3,i))
        enddo

        sigma_sum = sigma_sum*dpars(1)
        do i = 1,npts
           trac(1,i) = trac(1,i) + srcvals(10,i)*sigma_sum
           trac(2,i) = trac(2,i) + srcvals(11,i)*sigma_sum
           trac(3,i) = trac(3,i) + srcvals(12,i)*sigma_sum
        enddo
     endif
  endif
  
  return
end subroutine lpcomp_stok_s_trac_addsub
!
!
!        
subroutine stok_s_trac_solver(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, eps, numit, ifinout, &
     rhs, eps_gmres, niter, errs, rres, soln)

  !
  !  This subroutine solves the Stokes traction problem
  !  on the exterior/interior of an object where the potential
  !  is represented using a single layer integral representation.
  !
  !  This subroutine is the simple interface as opposed to the
  !  _solver_guru routine which is called after initialization
  !  in this routine.
  !
  !  Representation:
  !    u = S_{stok}[\sigma]
  !
  !  Boundary condition:
  !    trac[u] = f
  !
  !  The linear system is solved iteratively using GMRES
  !  until a relative residual of eps_gmres is reached
  !
  !
  !  Input arguments:
  !    - npatches: integer
  !        number of patches
  !    - norders: integer(npatches)
  !        order of discretization on each patch 
  !    - ixyzs: integer(npatches+1)
  !        ixyzs(i) denotes the starting location in srccoefs,
  !        and srcvals array corresponding to patch i
  !    - iptype: integer(npatches)
  !        type of patch
  !        iptype = 1, triangular patch discretized using RV nodes
  !        iptype = 11, quadrangular patch discretized with GL nodes
  !        iptype = 12, quadrangular patch discretized with Chebyshev nodes
  !    - npts: integer
  !        total number of discretization points on the boundary
  !    - srccoefs: real *8 (9,npts)
  !        basis expansion coefficients of xyz, dxyz/du,
  !        and dxyz/dv on each patch. 
  !        For each point 
  !          * srccoefs(1:3,i) is xyz info
  !          * srccoefs(4:6,i) is dxyz/du info
  !          * srccoefs(7:9,i) is dxyz/dv info
  !    - srcvals: real *8 (12,npts)
  !        xyz(u,v) and derivative info sampled at the 
  !        discretization nodes on the surface
  !          * srcvals(1:3,i) - xyz info
  !          * srcvals(4:6,i) - dxyz/du info
  !          * srcvals(7:9,i) - dxyz/dv info
  !          * srcvals(10:12,i) - normals info
  !    - eps: real *8
  !        precision requested for computing quadrature and fmm
  !        tolerance
  !    - numit: integer
  !        max number of gmres iterations
  !    - ifinout: integer
  !        ifinout = 0, interior problem
  !        ifinout = 1, exterior problem
  !    - rhs: real *8(3,npts)
  !        velocity data
  !    - eps_gmres: real *8
  !        gmres tolerance requested
  !      
  !
  !  Output arguments:
  !    - niter: integer
  !        number of gmres iterations required for relative residual
  !        to converge to eps_gmres
  !    - errs:  real *8 (numit+1)
  !        relative residual as a function of iteration
  !        number (only errs(1:niter) is populated))
  !    - rres: real *8
  !        relative residual for computed solution
  !    - soln: real *8(3,npts)
  !        density which solves the velocity problem \sigma
  !				 
  implicit none
  integer, intent(in) :: npatches, npts
  integer, intent(in) :: ifinout
  integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches)
  real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
  real *8, intent(in) :: eps, eps_gmres
  real *8, intent(in) :: rhs(3,npts)
  integer, intent(in) :: numit
  real *8, intent(out) :: soln(3,npts)
  real *8, intent(out) :: errs(numit+1)
  real *8, intent(out) :: rres
  integer, intent(out) :: niter

  integer norder, npols
  real *8, allocatable :: targs(:,:)
  integer, allocatable :: ipatch_id(:)
  real *8, allocatable :: uvs_targ(:,:)
  integer ndtarg, ntarg



  integer nover, npolso, nptso
  integer nnz,nquad
  integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
  real *8, allocatable :: wnear(:,:)

  real *8, allocatable :: srcover(:,:), wover(:)
  integer, allocatable :: ixyzso(:), novers(:)

  real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

  integer i, j, jpatch, jquadstart, jstart

  integer ipars
  complex *16 zpars
  real *8 timeinfo(10), t1, t2, omp_get_wtime


  real *8 ttot, done, pi
  real *8 rfac, rfac0
  integer iptype_avg, norder_avg
  integer ikerorder, iquadtype, npts_over

  !
  !
  !       gmres variables
  !
  real *8 did, dtmp
  complex *16 ztmp
  integer nker


  done = 1
  pi = atan(done)*4
  !
  !
  !        setup targets as on surface discretization points
  ! 
  ndtarg = 3
  ntarg = npts
  allocate(targs(ndtarg,npts), uvs_targ(2,ntarg), ipatch_id(ntarg))

  !$OMP PARALLEL DO DEFAULT(SHARED)
  do i=1,ntarg
     targs(1,i) = srcvals(1,i)
     targs(2,i) = srcvals(2,i)
     targs(3,i) = srcvals(3,i)
     ipatch_id(i) = -1
     uvs_targ(1,i) = 0
     uvs_targ(2,i) = 0
  enddo
  !$OMP END PARALLEL DO   


  !
  !    initialize patch_id and uv_targ for on surface targets
  !
  call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
       ipatch_id, uvs_targ)

  !
  !
  !        this might need fixing
  !
  iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
  norder_avg = floor(sum(norders)/(npatches+0.0d0))

  call get_rfacs(norder_avg, iptype_avg, rfac, rfac0)

  allocate(cms(3,npatches), rads(npatches), rad_near(npatches))
  
  call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, &
       srccoefs, cms, rads)

  !$OMP PARALLEL DO DEFAULT(SHARED) 
  do i=1,npatches
     rad_near(i) = rads(i)*rfac
  enddo
  !$OMP END PARALLEL DO      

  !
  !    find near quadrature correction interactions
  !
  print *, "entering find near mem"
  call findnearmem(cms, npatches, rad_near, ndtarg, targs, npts, nnz)

  allocate(row_ptr(npts+1), col_ind(nnz))

  call findnear(cms, npatches, rad_near, ndtarg, targs, npts, &
       row_ptr, col_ind)

  allocate(iquad(nnz+1)) 
  call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind, &
       iquad)

  ikerorder = 0

  !
  !    estimate oversampling for far-field, and oversample geometry
  !

  allocate(novers(npatches), ixyzso(npatches+1))

  print *, "beginning far order estimation"

  ztmp = 0

  call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
       rads, npts, srccoefs, ndtarg, npts, targs, ikerorder, ztmp, &
       nnz, row_ptr, col_ind, rfac, novers, ixyzso)

  npts_over = ixyzso(npatches+1) - 1

  allocate(srcover(12,npts_over), wover(npts_over))

  call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
       srccoefs, srcvals, novers, ixyzso, npts_over, srcover)

  call get_qwts(npatches, novers, ixyzso, iptype, npts_over, &
       srcover, wover)

  !
  !   compute near quadrature correction
  !
  nquad = iquad(nnz+1) - 1
  nker = 6
  allocate(wnear(nker,nquad))

  !$OMP PARALLEL DO DEFAULT(SHARED)      
  do i=1,nquad
     wnear(1,i) = 0
     wnear(2,i) = 0
     wnear(3,i) = 0
     wnear(4,i) = 0
     wnear(5,i) = 0
     wnear(6,i) = 0
  enddo
  !$OMP END PARALLEL DO    


  iquadtype = 1

  print *, "starting to generate near quadrature"
  call cpu_time(t1)
  !$      t1 = omp_get_wtime()      

  call getnearquad_stok_s_trac(npatches, norders, &
       ixyzs, iptype, npts, srccoefs, srcvals, eps, iquadtype, &
       nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
  call cpu_time(t2)
  !$      t2 = omp_get_wtime()     
  
  print *, "done generating near quadrature, now starting gmres"
  
  call stok_s_trac_solver_guru(npatches, norders, ixyzs, &
       iptype, npts, srccoefs, srcvals, eps, numit, ifinout, &
       rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
       npts_over, ixyzso, srcover, wover, eps_gmres, niter, &
       errs, rres, soln)

  return
end subroutine stok_s_trac_solver
!
!

subroutine stok_s_trac_solver_guru(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, eps, numit, ifinout, &
     rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
     nptso, ixyzso, srcover, whtsover, eps_gmres, niter, &
     errs, rres, soln)
  !
  !
  !  This subroutine solves the Stokes traction problem
  !  on the exterior of an object where the potential
  !  is represented using the single layer representation.
  !
  !
  !  Representation:
  !    u = S_{stok}[\sigma]
  !
  !  Boundary condition:
  !    trac[u] = f
  !
  !  The linear system is solved iteratively using GMRES
  !  until a relative residual of eps_gmres is reached
  !
  !
  !  Input arguments:
  !    - npatches: integer
  !        number of patches
  !    - norders: integer(npatches)
  !        order of discretization on each patch 
  !    - ixyzs: integer(npatches+1)
  !        ixyzs(i) denotes the starting location in srccoefs,
  !        and srcvals array corresponding to patch i
  !    - iptype: integer(npatches)
  !        type of patch
  !        iptype = 1, triangular patch discretized using RV nodes
  !        iptype = 11, quadrangular patch discretized with GL nodes
  !        iptype = 12, quadrangular patch discretized with Chebyshev nodes
  !    - npts: integer
  !        total number of discretization points on the boundary
  !    - srccoefs: real *8 (9,npts)
  !        basis expansion coefficients of xyz, dxyz/du,
  !        and dxyz/dv on each patch. 
  !        For each point 
  !          * srccoefs(1:3,i) is xyz info
  !          * srccoefs(4:6,i) is dxyz/du info
  !          * srccoefs(7:9,i) is dxyz/dv info
  !    - srcvals: real *8 (12,npts)
  !        xyz(u,v) and derivative info sampled at the 
  !        discretization nodes on the surface
  !          * srcvals(1:3,i) - xyz info
  !          * srcvals(4:6,i) - dxyz/du info
  !          * srcvals(7:9,i) - dxyz/dv info
  !          * srcvals(10:12,i) - normals info
  !    - eps: real *8
  !        precision requested for computing quadrature and fmm
  !        tolerance
  !    - numit: integer
  !        max number of gmres iterations
  !    - ifinout: integer
  !        ifinout = 0, interior problem
  !        ifinout = 1, exterior problem
  !    - rhs: real *8(3,npts)
  !        velocity data
  !    - nnz: integer
  !        number of source patch-> target interactions in the near field
  !    - row_ptr: integer(npts+1)
  !        row_ptr(i) is the pointer
  !        to col_ind array where list of relevant source patches
  !        for target i start
  !    - col_ind: integer (nnz)
  !        list of source patches relevant for all targets, sorted
  !        by the target number
  !    - iquad: integer(nnz+1)
  !        location in wnear_ij array where quadrature for col_ind(i)
  !        starts for a single kernel. In this case the different kernels
  !        are matrix entries are located at (m-1)*nquad+iquad(i), where
  !        m is the kernel number
  !    - nquad: integer
  !        number of near field entries corresponding to each source target
  !        pair
  !    - nker: integer
  !        number of kernels in quadrature correction, must be 6
  !    - wnear: real *8(nker, nquad)
  !        precomputed quadrature corrections 
  !    - novers: integer(npatches)
  !        order of discretization for oversampled sources and
  !        density
  !    - ixyzso: integer(npatches+1)
  !        ixyzso(i) denotes the starting location in srcover,
  !        corresponding to patch i
  !    - nptso: integer
  !        total number of oversampled points
  !    - srcover: real *8 (12,nptso)
  !        oversampled set of source information
  !    - whtsover: real *8 (nptso)
  !        smooth quadrature weights at oversampled nodes
  !    - eps_gmres: real *8
  !        gmres tolerance requested
  !
  !  output
  !    - niter: integer
  !        number of gmres iterations required for relative residual
  !        to converge to eps_gmres
  !    - errs:  real *8 (numit+1)
  !        relative residual as a function of iteration
  !        number (only errs(1:niter) is populated))
  !    - rres: real *8
  !        relative residual for computed solution
  !    - soln: real *8(3,npts)
  !        density which solves the neumann problem \sigma
  !				 
  
  implicit none
  integer, intent(in) :: npatches, npts
  integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches)
  real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
  real *8, intent(in) :: eps, eps_gmres
  integer, intent(in) :: ifinout
  real *8, intent(in) :: rhs(3,npts)
  integer, intent(in) :: numit

  real *8, intent(out) :: errs(numit+1)
  real *8, intent(out) :: rres
  integer, intent(out) :: niter

  integer, intent(in) :: nnz, nquad
  integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
  integer, intent(in) :: iquad(nnz+1)

  integer, intent(in) :: nker
  real *8, intent(in) :: wnear(nker,nquad)
  
  integer, intent(in) :: nptso
  integer, intent(in) :: novers(npatches), ixyzso(npatches+1)
  real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)


  real *8, intent(out) :: soln(3,npts)

  real *8 did
  real *8 dpars, pi, done, rsurf
  integer ndd_use

  procedure (), pointer :: fker
  external lpcomp_stok_s_trac_addsub

  integer ndd, ndi, ndz, lwork, ndim
  real *8 work
  integer ipars, nkertmp
  complex *16 zpars

  integer ndtarg
  complex *16 ima
  data ima/(0.0d0,1.0d0)/
  real *8, allocatable :: wts(:)

  integer i

  did = (-1)**(ifinout)/2d0
  fker => lpcomp_stok_s_trac_addsub

  ndd = 1
  ndi = 0
  ndz = 0

  lwork = 0
  ndim = 3
  allocate(wts(npts))

  call get_qwts(npatches, norders, ixyzs, iptype, npts, &
       srcvals, wts)
  rsurf = 0
  !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rsurf)      
  do i=1,npts
     rsurf = rsurf + wts(i)
  enddo
  !$OMP END PARALLEL DO      
  done = 1.0d0
  pi = atan(done)*4.0d0

  dpars = -2*pi/rsurf
  
  call dgmres_guru(npatches, norders, ixyzs, &
       iptype, npts, srccoefs, srcvals, wts, &
       eps, ndd_use, dpars, ndz, zpars, ndi, ipars, &
       nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
       nptso, ixyzso, srcover, whtsover, npts, wts, &
       ndim, fker, did, rhs, numit, eps_gmres, niter, errs, &
       rres, soln)
  
  return
end subroutine stok_s_trac_solver_guru
