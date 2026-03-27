

subroutine dgetnearquad_kernelmenu(npatches,norders, &
     ixyzs,iptype,npts,isd,ndsc,ndsv,srccoefs,srcvals, &
     ndtarg,ntarg,targvals,ipatch_id,uvs_targ,eps,ipv, &
     ckerfam,ckername,nker,ndd,dpars,ndz,zpars,ndi,ipars, &
     liopts,iopts,ldopts,dopts,nnz,row_ptr,col_ind,iquad, &
     nquad,wnear,linfo,info)
  
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
  character *40, intent(in) :: ckerfam, ckername

  procedure (), pointer :: fker

  !lap
  external l3d_comb, l3d_slp, l3d_dlp, l3d_sprime, l3d_qlp
  external l3d_sgrad_vec, l3d_spp_sum_dp
  
  ! stok
  external st3d_slp_vec, st3d_dlp_vec, st3d_slp, st3d_dlp
  external st3d_comb_vec, st3d_comb_vec6, st3d_comb
  external st3d_strac_vec, st3d_strac, st3d_sprime_vec
  external st3d_sprime, st3d_sdiv_vec
  
  integer *8 :: nker0

  info(1) = 0

  !
  !        initialize the appropriate kernel function
  !

  nker0 = 0
  if (trim(ckerfam) .eq. 'l3d') then
     ! laplace kernels 
     if (trim(ckername) .eq. 'comb') then
        fker => l3d_comb
        nker0 = 1
     elseif (trim(ckername) .eq. 'slp') then
        fker => l3d_slp
        nker0 = 1
     elseif (trim(ckername) .eq. 'dlp') then
        fker => l3d_dlp
        nker0 = 1
     elseif (trim(ckername) .eq. 'sprime') then
        fker => l3d_sprime
        nker0 = 1
     elseif (trim(ckername) .eq. 'qlp') then
        fker => l3d_sprime
        nker0 = 1
     elseif (trim(ckername) .eq. 'sgrad_vec') then
        fker => l3d_sgrad_vec
        nker0 = 3
     elseif (trim(ckername) .eq. 'spp_sum_dp') then
        fker => l3d_spp_sum_dp
        nker0 = 1
     else
        info(1) = 2048
     end if
  elseif (trim(ckerfam) .eq. 'st3d') then
     if (trim(ckername) .eq. 'comb') then
        fker => st3d_comb
        nker0 = 1
     elseif (trim(ckername) .eq. 'comb_vec') then
        fker => st3d_comb_vec
        nker0 = 9
     elseif (trim(ckername) .eq. 'comb_vec6') then
        fker => st3d_comb_vec6
        nker0 = 6
     elseif (trim(ckername) .eq. 'slp') then
        fker => st3d_slp
        nker0 = 1
     elseif (trim(ckername) .eq. 'slp_vec') then
        fker => st3d_slp_vec
        nker0 = 9
     elseif (trim(ckername) .eq. 'slp_vec6') then
        fker => st3d_slp_vec
        nker0 = 6
     elseif (trim(ckername) .eq. 'dlp') then
        fker => st3d_dlp
        nker0 = 1
     elseif (trim(ckername) .eq. 'dlp_vec') then
        fker => st3d_dlp_vec
        nker0 = 1
     elseif (trim(ckername) .eq. 'dlp_vec6') then
        fker => st3d_dlp_vec6
        nker0 = 1
     elseif (trim(ckername) .eq. 'sprime') then
        fker => st3d_sprime
        nker0 = 1
     elseif (trim(ckername) .eq. 'sprime_vec') then
        fker => st3d_sprime_vec
        nker0 = 9
     elseif (trim(ckername) .eq. 'sprime_vec6') then
        fker => st3d_sprime_vec
        nker0 = 6
     elseif (trim(ckername) .eq. 'sdiv_vec') then
        fker => st3d_sprime_vec
        nker0 = 3
     else
        info(1) = 2048
     endif
     
  else
     ! kernel family name not recognized
     info(1) = 1024
  end if

  if (nker0 .ne. nker) info(1) = 4096
  
  if (info(1) .ne. 0) return

  call dgetnearquad_guru(npatches,norders, &
       ixyzs,iptype,npts,isd,ndsc,ndsv,srccoefs,srcvals, &
       ndtarg,ntarg,targvals,ipatch_id,uvs_targ,eps,ipv, &
       fker,nker,ndd,dpars,ndz,zpars,ndi,ipars, &
       liopts,iopts,ldopts,dopts,nnz,row_ptr,col_ind,iquad, &
       nquad,wnear,linfo,info)
  
end subroutine dgetnearquad_kernelmenu
