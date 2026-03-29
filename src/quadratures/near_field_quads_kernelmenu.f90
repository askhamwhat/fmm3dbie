! The kernel menu routine does dispatch to dgetnearquad_guru
!
! The aims of this routine are:
! (1) replace redundant code in wrappers
! (2) reduce redundant code for MATLAB interfaces 
!
! Warning: this routine requires that the user correctly
! specifies ipv and nker. For complex kernels, nker should be
! 2 times the length of the complex kernel output
!
! The convention follows the subroutine naming
! convention. e.g. most subroutines have names of the
! form
!
! l3d_slp
!
! the kernel family name is the initial string
! ckerfam = 'l3d*'
! the kernel name is the remainder without the leading _
! ckername = 'slp*'
!
! THE ASTERISKS ARE REQUIRED TO OVERCOME DIFFICULTIES IN
! PASSING AND COMPARING STRINGS
!
! if the kernel name has more underscores they are retained
!
! e.g. st3d_slp_vec6 has
!
! ckerfam = 'st3d*' and ckername = 'slp_vec6*'
!

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

  character *40 :: trim_ast

  procedure (), pointer :: fker

  !!! built-ins 
  
  ! lap 3d
  external l3d_comb, l3d_slp, l3d_dlp, l3d_sprime, l3d_qlp
  external l3d_sgrad_vec, l3d_spp_sum_dp
  
  ! stok 3d
  external st3d_slp_vec, st3d_dlp_vec, st3d_slp, st3d_dlp
  external st3d_comb_vec, st3d_comb_vec6, st3d_comb
  external st3d_strac_vec, st3d_strac, st3d_sprime_vec
  external st3d_sprime, st3d_sdiv_vec
  external st3d_dlp_vec6, st3d_slp_vec6, st3d_sprime_vec6
  
  ! helm 3d
  external h3d_comb, h3d_slp, h3d_dlp, h3d_sprime, h3d_qlp
  external h3d_sgrad_vec, h3d_dprime_diff, h3d_slp_diff
  external h3d_dlp_diff, h3d_sprime_diff

  info(1) = 0

  ! initialize the appropriate kernel function

  if (trim(trim_ast(ckerfam)) .eq. 'l3d') then
     ! laplace kernels 
     if (trim(trim_ast(ckername)) .eq. 'comb') then
        fker => l3d_comb
     elseif (trim(trim_ast(ckername)) .eq. 'slp') then
        fker => l3d_slp
     elseif (trim(trim_ast(ckername)) .eq. 'dlp') then
        fker => l3d_dlp
     elseif (trim(trim_ast(ckername)) .eq. 'sprime') then
        fker => l3d_sprime
     elseif (trim(trim_ast(ckername)) .eq. 'qlp') then
        fker => l3d_sprime
     elseif (trim(trim_ast(ckername)) .eq. 'sgrad_vec') then
        fker => l3d_sgrad_vec
     elseif (trim(trim_ast(ckername)) .eq. 'spp_sum_dp') then
        fker => l3d_spp_sum_dp
     else
        info(1) = 2048
     end if
  elseif (trim(trim_ast(ckerfam)) .eq. 'st3d') then
     ! stokes kernels 
     if (trim(trim_ast(ckername)) .eq. 'comb') then
        fker => st3d_comb
     elseif (trim(trim_ast(ckername)) .eq. 'comb_vec') then
        fker => st3d_comb_vec
     elseif (trim(trim_ast(ckername)) .eq. 'comb_vec6') then
        fker => st3d_comb_vec6
     elseif (trim(trim_ast(ckername)) .eq. 'slp') then
        fker => st3d_slp
     elseif (trim(trim_ast(ckername)) .eq. 'slp_vec') then
        fker => st3d_slp_vec
     elseif (trim(trim_ast(ckername)) .eq. 'slp_vec6') then
        fker => st3d_slp_vec
     elseif (trim(trim_ast(ckername)) .eq. 'dlp') then
        fker => st3d_dlp
     elseif (trim(trim_ast(ckername)) .eq. 'dlp_vec') then
        fker => st3d_dlp_vec
     elseif (trim(trim_ast(ckername)) .eq. 'dlp_vec6') then
        fker => st3d_dlp_vec6
     elseif (trim(trim_ast(ckername)) .eq. 'sprime') then
        fker => st3d_sprime
     elseif (trim(trim_ast(ckername)) .eq. 'sprime_vec') then
        fker => st3d_sprime_vec
     elseif (trim(trim_ast(ckername)) .eq. 'sprime_vec6') then
        fker => st3d_sprime_vec6
     elseif (trim(trim_ast(ckername)) .eq. 'sdiv_vec') then
        fker => st3d_sprime_vec
     else
        info(1) = 2048
     endif
  elseif (trim(trim_ast(ckerfam)) .eq. 'h3d') then
     ! Helmholtz kernels 
     if (trim(trim_ast(ckername)) .eq. 'comb') then
        fker => h3d_comb
     elseif (trim(trim_ast(ckername)) .eq. 'slp') then
        fker => h3d_slp
     elseif (trim(trim_ast(ckername)) .eq. 'dlp') then
        fker => h3d_dlp
     elseif (trim(trim_ast(ckername)) .eq. 'sprime') then
        fker => h3d_sprime
     elseif (trim(trim_ast(ckername)) .eq. 'qlp') then
        fker => h3d_sprime
     elseif (trim(trim_ast(ckername)) .eq. 'sgrad_vec') then
        fker => h3d_sgrad_vec
     elseif (trim(trim_ast(ckername)) .eq. 'dprime_diff') then
        fker => h3d_dprime_diff
     elseif (trim(trim_ast(ckername)) .eq. 'slp_diff') then
        fker => h3d_slp_diff
     elseif (trim(trim_ast(ckername)) .eq. 'dlp_diff') then
        fker => h3d_dlp_diff
     elseif (trim(trim_ast(ckername)) .eq. 'sprime_diff') then
        fker => h3d_sprime_diff
     else
        info(1) = 2048
     end if
  else
     ! kernel family name not recognized
     info(1) = 1024
  end if

  if (info(1) .ne. 0) return

  call dgetnearquad_guru(npatches,norders, &
       ixyzs,iptype,npts,isd,ndsc,ndsv,srccoefs,srcvals, &
       ndtarg,ntarg,targvals,ipatch_id,uvs_targ,eps,ipv, &
       fker,nker,ndd,dpars,ndz,zpars,ndi,ipars, &
       liopts,iopts,ldopts,dopts,nnz,row_ptr,col_ind,iquad, &
       nquad,wnear,linfo,info)
  
end subroutine dgetnearquad_kernelmenu

function trim_ast(input_str) result(output_str)
  ! Explicitly defining 40 characters as requested
  character(len=40), intent(in) :: input_str
  character(len=40)             :: output_str
  integer                       :: ast_pos
  
  ! INDEX finds the position of the first occurrence of '*'
  ast_pos = index(input_str, '*')
  
  if (ast_pos > 0) then
     ! If an asterisk is found, take the substring right up to it.
     ! Fortran automatically pads the remaining characters with spaces.
     output_str = input_str(1:ast_pos-1)
  else
     ! If no asterisk is found, return the original string.
     output_str = input_str
  end if
  
end function trim_ast
