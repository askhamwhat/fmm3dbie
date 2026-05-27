c
c     This file has the following user callable routines
c
c     ?getnearquad_guru - guru interface for computing near
c              field quadrature
c        
c         z - complex
c         d - double precision
c
c     These routines handle:
c     - selecting GGQ or adaptive for self interactions
c     - vector-valued kernels
c     - interpolating/passing second derivative information
c       to kernels
c     - interpolating/passing data to kernels
c     - selecting other parameters and obtaining diagnostics
c      
c     Notes for developers:
c     - complex kernel calculations are handled by the
c     real-valued routine, setting nker = 2*nker0, where nker0
c     is the size of the complex vector-valued kernel.
c
c     TODO:
c     - implement vector-valued adaptive self 
c     - timing experiments to find any places where scalar
c     kernels suffer 
c     - debug quadparams_adap. deal with large outlier orders?
c     - chase through srccoefs 
c     - point legacy routines here (after timing)
c
      
      subroutine dgetnearquad_guru(npatches,norders,
     1     ixyzs,iptype,npts,isd,ndsc,ndsv,srccoefs,srcvals,
     2     ndtarg,ntarg,targvals,ipatch_id,uvs_targ,
     3     ifcustomdens,nordersdens,ldrdens,idenstype,eps,ipv,
     4     fker,nker,ndd,dpars,ndz,zpars,ndi,ipars,
     5     liopts,iopts,ldopts,dopts,nnz,row_ptr,col_ind,iquad,
     6     nquad,wnear,linfo,info)
c
c------------------------
c  This subroutine generates the near field quadrature
c  for a given kernel which is assumed to be
c  a compact/principal value integral operator
c  where the near field is specified by the user 
c  in row sparse compressed format.
c
c  The quadrature is computed by the following strategy
c  targets within a sphere of radius rfac0*rs
c  of a chunk centroid is handled using adaptive integration
c  where rs is the radius of the bounding sphere
c  for the patch
c  
c  All other targets in the near field are handled via
c  oversampled quadrature
c
c
c  Input arguments:
c
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization on each patch
c    - ixyzs: integer *8(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer *8(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quadrangular patch discretized using GL nodes,
c                        full degree polynomials
c        *  iptype = 12, quadrangular patch discretized using cheb nodes,
c                        full degree polynomials
c    - npts: integer *8
c        total number of discretization points on the boundary
c    - isd: integer *8
c        * isd = 0, no second der info in srccoefs
c        * isd = 1, second der info in srccoefs, pass to kern 
c    - ndsc: integer *8
c        leading dimension of srccoefs. >= 9 if isd=0 and >= 18 if isd=1
c    - ndsv: integer *8
c        leading dimension of srcvals. >= 12.
c    - srccoefs: double precision (ndsc,npts)
c        basis expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. For each patch
c        * srccoefs(1:3,i) is xyz info
c        * srccoefs(4:6,i) is dxyz/du info
c        * srccoefs(7:9,i) is dxyz/dv info
c        if (isd .eq. 1) 
c        * srccoefs(10:12,i) is d2xyz/du2 info
c        * srccoefs(13:15,i) is d2xyz/duv info
c        * srccoefs(16:18,i) is d2xyz/dv2 info
c        --> data in entries after 9 (isd = 0) or 18 (isd = 1) are
c        interpolated and passed as part of srcinfo to kernel
c        in entries after 12 (isd = 0) or 30 (isd = 1). 
c        --> unreasonably large ndsc will slow things down
c     - srcvals: double precision (ndsv,npts)
c        xyz(u,v) and derivative info sampled at the 
c        discretization nodes on the surface
c        * srcvals(1:3,i) - xyz info
c        * srcvals(4:6,i) - dxyz/du info
c        * srcvals(7:9,i) - dxyz/dv info
c        * srcvals(10:12,i) - normals info
c        it is not required that second order info be stored in
c        srcvals to use isd .eq. 1
c    - ndtarg: integer *8
c        leading dimension of target array
c    - ntarg: integer *8
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - ipatch_id: integer *8(ntarg)
c        patch id of target, patch_id = -1, if target off-surface
c    - uvs_targ: double precision (2,ntarg)
c        local uv coordinates on patch if target on surface
c    - ifcustomdens: integer *8 
c        flag. ifcustomdens = 0, ignores the next three inputs
c        ifcustomdens = 1, specify density basis type per patch
c    - nordersdens: integer *8(npatches)
c        nordersdens(i) is an integer parameter in the basis 
c        for the density on patch i, usually the order of the basis 
c    - ldrdens: integer *8(npatches+1)       
c        ldrdens(i+1)-ldrdens(i) is the number of basis functions 
c        on patch i. this determines the number of integrals 
c        to be computed for this patch.
c     - eps: double precision
c        precision requested
c    - ipv: integer *8
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators
c        * ipv = 2, for hypersingular operators
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c        val is real*8(nker)
c    - nker: integer *8
c        size of output of fker
c    - ndd: integer *8
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer *8
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer *8
c        number of integer *8 parameters
c    - ipars: integer *8(ndi)
c        integer *8 parameters
c    - liopts: integer *8
c        length of integer options array 
c    - iopts: integer *8(liopts)
c        For experts. Set various options for quadrature generation.  
c        parameters are listed below with [default values]. any input
c        with an index later than liopts is ignored. liopts=0 gives all
c        default. 
c        * iopts(1) = iselftype [1]. iselftype = 1 (ggq)
c                                    iselftype = 2 (adaptive)
c        * iopts(2) = intype [2]. see documentation of dtria/dquadints
c        * iopts(3) = istrat [2]. see documentation of dtria/dquadints
c        * iopts(4) = ifp [0]. for istrat 1 or 3. see documentation of
c                        dtria/dquadints
c    - ldopts: integer *8
c        length of real options array 
c    - dopts: real *8(ldopts)
c        For experts. Set various options for quadrature generation.  
c        parameters are listed below with [default values]. any input
c        with an index later than ldopts is ignored. ldopts=0 gives all
c        default. 
c        * dopts(1) = rfac0 [1.25d0]. radius parameter for near field
c    - nnz: integer *8
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer *8(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer *8 (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer *8(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - nquad: integer *8
c     number of entries in wnear
c     - linfo: integer *8
c     length of info array. 
c
c  Output parameters:
c
c    - wnear: double precision(nquad)
c     near field quadrature corrections
c     - info: integer *8(linfo)
c        information about run 
c----------------------------------------------------               
c
      
      implicit none
      integer *8, intent(in) :: ndi,ndd,ndz,ipv
      integer *8, intent(in) :: ipars(ndi)
      integer *8, intent(in) :: ndtarg, nker, isd, ndsc, ndsv
      integer *8, intent(in) :: npatches,norders(npatches),npts
      integer *8, intent(in) :: liopts, iopts(liopts)
      integer *8, intent(in) :: ldopts, nquad
      real *8, intent(in) :: dopts(ldopts)
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(ndsc,npts),srcvals(ndsv,npts),eps
      integer *8, intent(in) :: ntarg,ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: nnz, linfo
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(out) :: wnear(nker,nquad)
      integer *8, intent(out) :: info(linfo)
      integer *8, intent(in) :: ifcustomdens,nordersdens(npatches)
      integer *8, intent(in) :: ldrdens(npatches+1),idenstype(npatches)

      integer *8 ntrimax
      real *8, allocatable :: cms(:,:),rads(:)
      real *8, allocatable :: targ_near(:,:),targ_far(:,:)
      integer *8, allocatable :: iind_near(:),iind_far(:)
      real *8, allocatable :: umatr(:,:),vmatr(:,:),uvs(:,:),wts(:)

c
c        temporary variables
c
      integer *8, allocatable :: col_ptr(:),row_ind(:),iper(:)
      real *8, allocatable :: sints_n(:,:,:),svtmp_n(:,:,:)
      real *8, allocatable :: sints_f(:,:,:),svtmp_f(:,:,:)
      real *8, allocatable :: svtmp2(:,:,:)

      real *8, allocatable :: xyztarg2(:,:)

      real *8, allocatable :: qnodes_tri(:,:),qwts_tri(:)
      real *8, allocatable :: qnodes_quad(:,:),qwts_quad(:)
      real *8 ra

      real *8 ff1,ff2,cra1,cra2
      integer *8 nlev, nqorder_f, nqordert, nqorderq
      real *8 rfac0,rsc,rr,tmp(3),epsp
      real *8 done,dzero,eps_adap, rfac,rn1
      integer *8 norder_avg
      integer *8 ipoly, iselftype 
      character *1 ttype
      integer *8 int8_1, int8_11
      integer *8 i, ifmetric, ifp, ii, ii2, iii, iiif
      integer *8 intype, ipatch, iqstart, istrat, itarg2
      integer *8 itargptr,j,jpatch,jpt,jtarg,l,n2,ncols
      integer *8 nnodes, istart, norder, npmax, npols
      integer *8 npts_f_quad, npts_f_tri, nquad_f, ntarg2, ntarg2m
      integer *8 ntarg_f, ntarg_n, ntest0, ntri_f, ier
      integer *8 ndens

      external fker

      int8_1 = 1
      int8_11 = 11
      done = 1
      dzero = 0

      ier = 0
      if (((isd .eq. 1) .and. (ndsc .lt. 18)) .or.
     1     (ndsc .lt. 9)) then
         ier = 1
         if (linfo .gt. 0) info(1) = 1
         return
      endif

c     
c     unpack options
c

c defaults       

c     ggq for self
      iselftype = 1 
c     use xiao gimbutas nodes      
      intype = 2
c     use adaptive integration (see documentation of ctriaints)     
      istrat = 2
c     relevant parameter for istrat =1/3      
      ifp = 0
c     based on experiments
      rfac0 = 1.25d0

      if (liopts .gt. 0) iselftype = iopts(1) 
      if (liopts .gt. 1) intype = iopts(2) 
      if (liopts .gt. 2) istrat = iopts(3)
      if (liopts .gt. 3) ifp = iopts(4)

      if (ldopts .gt. 0) rfac0 = dopts(1)
      
c     
c     
c     transpose the row sparse compressed format
c     to the column sparse compressed format needed for the
c     the adaptive integration routines
c     
      allocate(row_ind(nnz),col_ptr(npatches+1),iper(nnz))

      call rsc_to_csc(npatches,ntarg,nnz,row_ptr,col_ind,col_ptr,
     1     row_ind,iper)


      allocate(cms(3,npatches),rads(npatches))

c     
c     find chunk radius and centroid again
c     

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,
     1     srccoefs,cms,rads)

c
c     TODO: return to this. return total kernel calls in info?
c     
c     suppress computation of 
c     some relevant metrics for analyzing adaptive
c     quadrature performance
c     
      ifmetric = 0
      rn1 = 0
      n2 = 0
      
      npmax = 3000
      norder_avg = floor(sum(norders)/(npatches+0.0d0)+0.49d0)
      if(norder_avg.ge.8) npmax = 6000

c     
c     Near quadrature params
c     - nqorder is the order of quadrature nodes
c     used on each triangle for handling the adaptive
c     integration. Optimal choice depends on order of
c     discretization and accuracy requested
c     
c     Far quadrature params
c     - nlev is the number of refinements of the standard simplex/
c     square     
c     - nqorder_f is the order of XG nodes on each of the smaller
c     triangles/quads.
c
c     nlev and nqorder_f depend on eps,norder_avg,
c     

c     TODO: eps_adap gets chosen by this routine. could in theory
c     depend on quad vs tri? probably not important
c     
      eps_adap = eps
      
c     
c     get triangle parameters
c
      
      call get_quadparams_adap0(eps,int8_1,norder_avg,rfac0,
     1     nqordert,eps_adap,nlev,nqorder_f)

      call triasymq_pts(nqorder_f,nnodes)
      
      ntri_f = 4**nlev
      npts_f_tri = ntri_f*nnodes
      allocate(qnodes_tri(2,npts_f_tri),qwts_tri(npts_f_tri))

      call gen_xg_unif_nodes_tri(nlev,nqorder_f,nnodes,npts_f_tri,
     1     qnodes_tri,qwts_tri)

c     
c     get quad parameters
c     
      
      call get_quadparams_adap0(eps,int8_11,norder_avg,rfac0,
     1     nqorderq,eps_adap,nlev, nqorder_f)
      
      call squarearbq_pts(nqorder_f,nnodes)
      
      nquad_f = 4**nlev
      npts_f_quad = nquad_f*nnodes
      allocate(qnodes_quad(2,npts_f_quad),qwts_quad(npts_f_quad))

      call gen_xg_unif_nodes_quad(nlev,nqorder_f,nnodes,npts_f_quad,
     1     qnodes_quad,qwts_quad)

c     number of source patches to be processed per batch
c     ntest0 must be 1 now that we are allowing different bases
c     per patch
      ntest0 = 1
      itargptr = 1

c     this defines the notion of well-separated boxes for
c     initial tree in istrat 1 or 3.
c      
c     TODO: consider tuning rfac based on nqorder,etc
      rfac = 1.0d0     

      
C$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC)
C$OMP$PRIVATE(ipatch,ntarg2,xyztarg2,sints_n,sints_f,svtmp_n,svtmp_f)
C$OMP$PRIVATE(ii,ii2,i,jpatch,svtmp2,iiif,l,ntarg2m)
C$OMP$PRIVATE(j,iii,istart,itarg2,iqstart,jpt,jtarg)
C$OMP$PRIVATE(targ_near,targ_far,iind_near,iind_far,rr)
C$OMP$PRIVATE(ntarg_f,ntarg_n,npols,norder)
C$OMP$PRIVATE(uvs,umatr,vmatr,wts,ndens)
C$OMP$PRIVATE(rsc,tmp,epsp,ipoly,ttype,ncols)
      do ipatch=1,npatches

         npols = ixyzs(ipatch+1)-ixyzs(ipatch)
         norder = norders(ipatch)
         if (ifcustomdens .eq. 1) then
            norderdens = nordersdens(ipatch)
            ndens = ldrdens(ipatch+1)-ldrdens(ipatch)
         else
            norderdens = norder
            ndens = npols
         endif

         allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
         allocate(wts(npols))
         if(iptype(ipatch).eq.11) ipoly = 0
         if(iptype(ipatch).eq.12) ipoly = 1
         ttype = "F"

         call get_disc_exps(norder,npols,iptype(ipatch),uvs,
     1        umatr,vmatr,wts)

c     
c     estimate rescaling of epsp needed to account for the scale of the
c     patch 
c     

         ii = ixyzs(ipatch)
         call cross_prod3d(srcvals(4,ii),srcvals(7,ii),tmp)
         rsc = (tmp(1)**2 + tmp(2)**2 + tmp(3)**2)*wts(1)**2
         do i=2,npols
            call cross_prod3d(srcvals(4,ii+i-1),srcvals(7,ii+i-1),tmp)
            rr = (tmp(1)**2 + tmp(2)**2 + tmp(3)**2)*wts(i)**2
            if(rr.lt.rsc) rsc = rr
         enddo
         epsp = eps_adap*rsc**0.25d0

         ntarg2m = col_ptr(ipatch+1)-col_ptr(ipatch) 
         allocate(xyztarg2(ndtarg,ntarg2m),svtmp2(nker,npols,ntarg2m))
         allocate(targ_near(ndtarg,ntarg2m),iind_near(ntarg2m))
         allocate(targ_far(ndtarg,ntarg2m),iind_far(ntarg2m))

         ii = 1
         ii2 = 1
         do i=col_ptr(ipatch),col_ptr(ipatch+1)-1
            jtarg = row_ind(i)
            jpatch = ipatch_id(jtarg)

            if(ipatch.ne.jpatch) then
               do iiif=1,ndtarg
                  xyztarg2(iiif,ii2) = targvals(iiif,jtarg)
               enddo
               ii2 = ii2 + 1
            endif
         enddo
         ntarg2 = ii2-1

c     
c     split near far targs 
c     
         rr = rfac0*rads(ipatch)
         call get_near_far_split_pt(ndtarg,ntarg2,xyztarg2,rr,
     1        cms(1,ipatch),ntarg_f,ntarg_n,targ_far,targ_near,
     2        iind_far,iind_near)

         allocate(sints_n(nker,npols,ntarg_n))
         allocate(svtmp_n(nker,npols,ntarg_n))
         
         allocate(sints_f(nker,npols,ntarg_f))
         allocate(svtmp_f(nker,npols,ntarg_f))
         
         istart = ixyzs(ipatch)

c     
c     fill out near part of layer potential
c     

         if(iptype(ipatch).eq.1) 
     1        call dtriaints_vec(epsp,istrat,intype,ntest0,norder,npols,
     2        isd,ndsc,srccoefs(1,istart),ndtarg,ntarg_n,targ_near,ifp,
     3        xyztarg2,itargptr,ntarg_n,norder,npols,fker,nker,
     4        ndd,dpars,ndz,zpars,ndi,ipars,nqordert,npmax,rfac,sints_n,
     5        ifmetric,rn1,n2)

         if(iptype(ipatch).eq.11.or.iptype(ipatch).eq.12) 
     1        call dquadints_vec(epsp,istrat,intype,ntest0,norder,ipoly,
     2        ttype,npols,isd,ndsc,srccoefs(1,istart),ndtarg,ntarg_n,
     3        targ_near,ifp,xyztarg2,itargptr,ntarg_n,norder,npols,
     4        fker,nker,ndd,dpars,ndz,zpars,ndi,ipars,nqorderq,npmax,
     5        rfac,sints_n,ifmetric,rn1,n2)


         ncols = nker*ntarg_n

         if (nker.gt.1) then
           call permute_12_3d(sints_n,svtmp_n,nker,npols,ntarg_n)
           call dgemm_guru('t','n',npols,ncols,npols,done,umatr,npols,
     1          svtmp_n,npols,dzero,sints_n,npols)
           call permute_12_3d(sints_n,svtmp_n,npols,nker,ntarg_n) 
         else
           call dgemm_guru('t','n',npols,ncols,npols,done,umatr,npols,
     1          sints_n,npols,dzero,svtmp_n,npols)
         endif
c     
c     fill out far part of layer potential
c     
         if(iptype(ipatch).eq.1) 
     1        call dtriaints_wnodes_vec(ntest0,norder,npols,isd,ndsc,
     2        srccoefs(1,istart),ndtarg,ntarg_f,targ_far,itargptr,
     3        ntarg_f,norder,npols,fker,nker,ndd,dpars,ndz,zpars,
     4        ndi,ipars,npts_f_tri,qnodes_tri,qwts_tri,sints_f)
         if(iptype(ipatch).eq.11.or.iptype(ipatch).eq.12) 
     1        call dquadints_wnodes_vec(ntest0,norder,ipoly,ttype,npols,
     2        isd,ndsc,srccoefs(1,istart),ndtarg,ntarg_f,targ_far,
     3        itargptr,ntarg_f,norder,npols,fker,nker,ndd,dpars,
     4        ndz,zpars,ndi,ipars,npts_f_quad,qnodes_quad,qwts_quad,
     5        sints_f)


         ncols = nker*ntarg_f
         if (nker.gt.1) then
           call permute_12_3d(sints_f,svtmp_f,nker,npols,ntarg_f)
           call dgemm_guru('t','n',npols,ncols,npols,done,umatr,npols,
     1          svtmp_f,npols,dzero,sints_f,npols)
           call permute_12_3d(sints_f,svtmp_f,npols,nker,ntarg_f) 
         else
           call dgemm_guru('t','n',npols,ncols,npols,done,umatr,npols,
     1          sints_f,npols,dzero,svtmp_f,npols)
         endif
c     
c     combine svtmp_f, svtmp_n to fill out svtmp2
c     
         do i=1,ntarg_f

            ii = iind_far(i)
            do j=1,npols
               do l = 1,nker
                  svtmp2(l,j,ii) = svtmp_f(l,j,i)
               enddo
            enddo
         enddo

         do i=1,ntarg_n
            ii = iind_near(i)
            do j=1,npols
               do l = 1,nker
                  svtmp2(l,j,ii) = svtmp_n(l,j,i)
               enddo
            enddo
         enddo

c     
c     
c     fill out relevant sections of wnear, call self quad for
c     targets on same patch
c     
         itarg2 = 1
         do i=col_ptr(ipatch),col_ptr(ipatch+1)-1
            jtarg = row_ind(i)
            jpatch = ipatch_id(jtarg)
            iqstart = iquad(iper(i))-1

            if(ipatch.eq.jpatch) then
               if (iselftype .eq. 1)
     $              call dget_ggq_self_quad_pt_vec(ipv,norder,npols,
     1              uvs_targ(1,jtarg),iptype(ipatch),
     2              umatr,isd,ndsc,srccoefs(1,istart),
     3              ndtarg,targvals(1,jtarg),fker,nker,ndd,dpars,
     4              ndz,zpars,ndi,ipars,wnear(1,iqstart+1))
            endif

            if(ipatch.ne.jpatch) then
               do l=1,npols
                  do j = 1,nker
                     wnear(j,iqstart+l) = svtmp2(j,l,itarg2)
                  enddo
               enddo
               itarg2 = itarg2 + 1
            endif
         enddo

         deallocate(xyztarg2,svtmp2,svtmp_f,svtmp_n,sints_f,sints_n)
         deallocate(targ_near,targ_far,iind_near,iind_far)
         deallocate(uvs,umatr,vmatr,wts)
      enddo
C$OMP END PARALLEL DO      

      return
      end
c
c
      subroutine dget_ggq_self_quad_pt_vec(ipv,norder,npols,uvs,iptype,
     1     umat,isd,ndsc,srccoefs,ndtarg,targ,fker,nker,ndd,dpars,
     2     ndz,zpars,ndi,ipars,dquad)
c
c
c------------------
c  This subroutine evaluates the integral of an integral
c  operator at a target on the interior of a triangluar patch 
c  using the generalized Gaussian quadratures developed by
c  Gimbutas in https://github.com/zgimbutas/selfquad3d
c
c  The quadrature currently cannot handle targets on the boundary
c  of the triangle
c
c  Input arguments:
c    - ipv: integer *8
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators
c        * ipv = 2, for hypersingular operators
c    - isd: integer *8
c        * isd = 0, no second der info in srccoefs
c        * isd = 1, second der info in srccoefs, pass to kern eval 
c    - ndsc: integer *8
c        leading dimension of srccoefs. 9 if isd=0 and 18 if isd=1   
c    - norder: integer *8
c        order of patch discretization
c    - npols: integer *8
c        number of discretization nodes on patch
c    - uvs: double precision(2)
c        local u,v coordinates of target
c    - iptype: integer *8(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quadrangular patch discretized using GL nodes,
c                        full degree polynomials
c        *  iptype = 12, quadrangular patch discretized using cheb nodes,
c                        full degree polynomials
c    - umat: double precision(npols,npols)
c        values to coefficient matrix
c    - srccoefs: double precision (ndsc,npts)
c        basis expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. For each patch
c        * srccoefs(1:3,i) is xyz info
c        * srccoefs(4:6,i) is dxyz/du info
c        * srccoefs(7:9,i) is dxyz/dv info
c        if (isd .eq. 1) 
c        * srccoefs(10:12,i) is d2xyz/du2 info
c        * srccoefs(13:15,i) is d2xyz/duv info
c        * srccoefs(16:18,i) is d2xyz/dv2 info
c    - ndtarg: integer *8
c        leading dimension of target array
c    - targ: double precision (ndtarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c        val is real *8(nker) 
c     - nker: integer *8
c        size of fker output 
c     - ndd: integer *8
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer *8
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer *8
c        number of integer *8 parameters
c    - ipars: integer *8(ndi)
c        integer *8 parameters
c
c  Output parameters:
c    
c    - dquad: double precision(nker,npols)
c        near quadrature correction
c    
c
c


      implicit real *8(a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: norder,npols,ndtarg,nker,isd,ndsc
      real *8, intent(in) :: srccoefs(ndsc,npols)
      real *8, intent(in) :: targ(ndtarg)
      real *8, intent(in) :: uvs(2),umat(npols,npols)
      integer *8, intent(in) :: ndi,ndd,ndz
      integer *8, intent(in) :: ipars(ndi)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      real *8, intent(out) :: dquad(nker,npols)

      real *8 srcvals(30)
      real *8, allocatable :: sigvalstmp(:)
      real *8, allocatable :: srctmp(:,:)
      real *8, allocatable :: qwts(:),sigvals(:,:)
      real *8, allocatable :: xs(:),ys(:),ws(:)
      real *8 uv(2),verts(2,4)
      real *8 alpha,beta
      real *8 fval
      real *8, allocatable :: fint(:,:),fvals(:,:)
      character *1 transa,transb
      real *8 done,dzero 
      integer *8 n9,n1,lda

      external fker
      
      done = 1
      dzero = 0

      nmax = 10000
      allocate(ws(nmax),xs(nmax),ys(nmax))

      if(iptype.eq.1) then
        nv = 3
        verts(1,1) = 0
        verts(2,1) = 0
        verts(1,2) = 1
        verts(2,2) = 0
        verts(1,3) = 0
        verts(2,3) = 1
      endif
      if(iptype.eq.11.or.iptype.eq.12) then
         nv = 4
         verts(1,1) = -1
         verts(2,1) = -1

         verts(1,2) = 1
         verts(2,2) = -1

         verts(1,3) = 1
         verts(2,3) = 1
         
         verts(1,4) = -1
         verts(2,4) = 1
      endif
      allocate(fint(nker,npols),sigvalstmp(npols))

c
c       compute all source info at target point on patch
c
      call get_basis_pols(uvs,norder,npols,iptype,sigvalstmp)

      alpha = 1.0d0
      beta = 0.0d0
      lda = ndsc
      n1 = 1
      call dgemv_guru('n',lda,npols,alpha,srccoefs,lda,sigvalstmp,n1,
     1     beta,srcvals,n1)
      

      do j=1,npols
         do l = 1,nker
            fint(l,j) = 0
         enddo
      enddo
      ns = 0
      call self_quadrature(norder, ipv, verts, nv, uvs(1), uvs(2),
     1    srcvals(4), ns, xs, ys, ws, ier)

      if (isd .eq. 1) then
         ndata = ndsc-9
         ndsv=12 + ndata
      else
         ndata = ndsc-18
         ndsv=30 + ndata 
      endif
      allocate(srctmp(ndsv,ns),qwts(ns),fvals(nker,ns))
      allocate(sigvals(npols,ns))

      do i=1,ns
        uv(1) = xs(i)
        uv(2) = ys(i)
        call get_basis_pols(uv,norder,npols,iptype,sigvals(1,i))
      enddo

      transa = 'N'
      transb = 'N'
      alpha = 1
      beta = 0
      lda = ndsc
      ldb = npols
      ldc = ndsv

      call dgemm_guru(transa,transb,ndsc,ns,npols,alpha,
     1     srccoefs,lda,sigvals,ldb,beta,srctmp,ldc)
c
      call get_srcvals_auxinfo(ns,ws,isd,ndata,ndsv,srctmp,
     1     done,qwts)
      
      do i=1,ns
         call fker(srctmp(1,i),ndtarg,targ,ndd,dpars,
     2        ndz,zpars,ndi,ipars,fvals(1,i))
      enddo
      do i=1,ns
         qtmp = qwts(i)
         do j=1,npols
            sigvals(j,i) = sigvals(j,i)*qtmp
         enddo
      enddo

      call dgemm_guru('n','t',nker,npols,ns,done,fvals,nker,
     1     sigvals,npols,dzero,fint,nker)

      deallocate(srctmp,qwts,sigvals,fvals)

      call dgemm_guru('n','n',nker,npols,npols,done,fint,nker,
     1     umat,npols,dzero,dquad,nker)

      return
      end
c     
c     
c
      
      subroutine dcopytocolsbyinds(m,n,a,b,inds)
c     b(i,inds(j)) <- a(i,j)
c
      implicit none
      real *8 :: a(m,n), b(m,*)
      integer *8 :: m,n,inds(*)
c     local
      integer *8 :: i, j, jj

      do j = 1,n
         jj = inds(j)
         do i = 1,m
            b(i,jj) = a(i,j)
         enddo
      enddo

      return
      end

      subroutine dcopy_seq(a,b,m)
c     b(i) <- a(i)
      implicit none
      real *8 :: a(*), b(*)
      integer *8 :: m
c     local
      integer *8 :: i

      do i = 1,m
         b(i) = a(i)
      enddo

      return
      end
      
c

