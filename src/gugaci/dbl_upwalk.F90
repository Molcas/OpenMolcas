!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine dbl_upwalk()

use gugaci_global, only: jpad_upwei, jroute_sys, lsm_inn, mxnode, ng_sm, norb_dbl, norb_dz, norb_frz, ns_sm, nu_ad
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iw, lri, lrj, lsmi, lsmid, lsmij, lsmit, lsmj, no_d, no_t, node

nu_ad(:) = 0
jpad_upwei(:) = 0

if (norb_dbl == 1) then

  ! v(1),d(2-9),s(18-25)           for s=0
  ! v(1),d(2-9),s(18-25),d'(26-33)   for s<>0

  mxnode = 17+ng_sm
  lri = norb_frz+1
  lsmi = lsm_inn(lri)
  lsmid = Mul(lsmi,ns_sm)
  ! for node_v
  nu_ad(1) = 1
  jpad_upwei(1) = 1
  ! for node_d
  nu_ad(1+lsmid) = 1+lsmid
  jpad_upwei(1+lsmid) = 1
  ! for node_s
  nu_ad(17+ns_sm) = 17+ns_sm
  jpad_upwei(17+ns_sm) = 1

  if (jroute_sys == 1) return
  mxnode = 25+ng_sm
  ! for node_d'
  nu_ad(25+lsmid) = 25+lsmid
  jpad_upwei(25+lsmid) = 1

else

  nu_ad(1) = 1
  jpad_upwei(1) = 1
  if (norb_dbl == 0) then
    mxnode = 1
    return
  end if
  do lri=norb_frz+1,norb_dz
    lsmi = lsm_inn(lri)
    lsmid = Mul(lsmi,ns_sm)
    no_d = lsmid+1
    jpad_upwei(no_d) = jpad_upwei(no_d)+1
    do lrj=lri+1,norb_dz
      lsmj = lsm_inn(lrj)
      lsmij = Mul(lsmi,lsmj)
      lsmit = Mul(lsmij,ns_sm)
      no_t = lsmit+9
      jpad_upwei(no_t) = jpad_upwei(no_t)+1
    end do
  end do
  ! v(1),d(2-9),t(10-17),s(18-25),d'(26-33),t'(34-41)
  select case (jroute_sys)
    case default ! (1)
      mxnode = 25                     !v,d,t,s
      jpad_upwei(18:25) = jpad_upwei(10:17)
      jpad_upwei(17+ns_sm) = jpad_upwei(17+ns_sm)+norb_dbl
    case (2)
      mxnode = 25+8
      jpad_upwei(18:25) = jpad_upwei(10:17)+jpad_upwei(10:17)
      jpad_upwei(17+ns_sm) = jpad_upwei(17+ns_sm)+norb_dbl
      jpad_upwei(26:33) = jpad_upwei(2:9)
    case (3)
      mxnode = 25+8+8
      jpad_upwei(18:25) = jpad_upwei(10:17)+jpad_upwei(10:17)
      jpad_upwei(17+ns_sm) = jpad_upwei(17+ns_sm)+norb_dbl
      jpad_upwei(26:33) = jpad_upwei(2:9)
      jpad_upwei(34:41) = jpad_upwei(10:17)
  end select

  do node=2,mxnode
    iw = jpad_upwei(node)
    if (iw == 0) cycle
    nu_ad(node) = node
  end do

end if

return

end subroutine dbl_upwalk

subroutine ext_downwalk()

use gugaci_global, only: iseg_downwei, ng_sm, nlsm_ext, norb_ext, nu_ae
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: im, imi, imij, imj, iwmij(8)

nu_ae(:) = 0
nu_ae(1) = 1
do im=1,ng_sm
  nu_ae(1+im) = 1+im
  nu_ae(9+im) = 9+im
  nu_ae(17+im) = 17+im
end do

iwmij = 0
iseg_downwei(nu_ae(1)) = 1
do imi=1,ng_sm
  iseg_downwei(nu_ae(1+imi)) = nlsm_ext(imi)
  do imj=imi,ng_sm
    imij = Mul(imi,imj)
    if (imij /= 1) then
      iwmij(imij) = iwmij(imij)+nlsm_ext(imi)*nlsm_ext(imj)
      cycle
    end if
    iwmij(1) = iwmij(1)+nlsm_ext(imi)*(nlsm_ext(imi)-1)/2
  end do
end do
do im=1,ng_sm
  iseg_downwei(nu_ae(9+im)) = iwmij(im)
  iseg_downwei(nu_ae(17+im)) = iwmij(im)
end do
iseg_downwei(nu_ae(18)) = iseg_downwei(nu_ae(18))+norb_ext

return

end subroutine ext_downwalk

subroutine readdrt(ludrt)

use gugaci_global, only: ja, jb, jd, jj, jm, js, jt, jv, kk, no, norb_inn
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ludrt
integer(kind=iwp) :: id, idisk, idx(2)

idisk = 0
! number of nodes
call idafile(ludrt,2,idx,2,idisk)
idisk = idx(2)
call idafile(ludrt,2,idx,1,idisk)
id = idx(1)
call idafile(ludrt,2,ja,id,idisk)
call idafile(ludrt,2,jb,id,idisk)
call idafile(ludrt,2,jm,id,idisk)
call idafile(ludrt,2,jj,4*(id+1),idisk)
call idafile(ludrt,2,kk,1+id,idisk)
no(:) = 0
call idafile(ludrt,2,no,norb_inn+2,idisk)
call idafile(ludrt,2,idx,1,idisk)
jv = idx(1)
call idafile(ludrt,2,jd,8,idisk)
call idafile(ludrt,2,jt,8,idisk)
call idafile(ludrt,2,js,8,idisk)

return

end subroutine readdrt

! juv,just(nost,nost),jud(nost)
! |   \  1           |
! | d,dd,s(i=i)      |
! |     \ s,t,tt(i<j)|
! |      \       1 2 |     deal with inner of dbl_space
! |ss(i>j)\          |
! |  2 1   \         |
subroutine dbl_downwalk()

use gugaci_global, only: iseg_downwei, iseg_sta, jud, just, lsm_inn, ng_sm, norb_dbl, norb_dz, norb_frz, ns_sm
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: im, ismi, ismij, ismj, lr0, lri, lrj, nnd, nns, nnt !, lsml(10,10)       !to del

if (norb_dbl == 0) then
  !----------- norb_dbl=0 ----------------------------------------------
  do im=1,ng_sm
    nnd = iseg_sta(1+im)
    nnt = iseg_sta(9+im)
    nns = iseg_sta(17+im)
    do lri=norb_dz,norb_frz+1,-1
      ismi = lsm_inn(lri)
      if (ismi /= im) cycle
      jud(lri) = nnd
      nnd = nnd+iseg_downwei(1+im)
    end do
    do lrj=norb_dz,norb_frz+1,-1
      ismj = lsm_inn(lrj)
      do lri=lrj,1,-1
        ismi = lsm_inn(lri)
        ismij = Mul(ismi,ismj)
        if (ismij /= im) cycle
        just(lri,lrj) = nns
        nns = nns+iseg_downwei(17+im)
        if (lri == lrj) cycle
        just(lrj,lri) = nnt
        nnt = nnt+iseg_downwei(9+im)
      end do
    end do
  end do
end if
!----------- norb_dbl=0 ------------------------------------------------
!----------- norb_dbl<>0 -----------------------------------------------
do im=1,ng_sm
  nnd = 0
  nns = 0
  do lri=norb_frz+1,norb_dz
    ismi = Mul(lsm_inn(lri),ns_sm)
    if (ismi /= im) cycle
    jud(lri) = nnd
    nnd = nnd+1
  end do
  do lri=norb_frz+1,norb_dz-1
    ismi = Mul(lsm_inn(lri),ns_sm)
    do lrj=lri+1,norb_dz      !tmp
      ismj = lsm_inn(lrj)
      ismij = Mul(ismi,ismj)
      if (ismij /= im) cycle
      just(lri,lrj) = nns
      nns = nns+1
    end do
  end do
  if (im == ns_sm) then
    do lr0=norb_frz+1,norb_dz
      just(lr0,lr0) = nns
      nns = nns+1
    end do
  end if
  do lri=norb_frz+1,norb_dz-1
    ismi = Mul(lsm_inn(lri),ns_sm)
    do lrj=lri+1,norb_dz      !tmp
      ismj = lsm_inn(lrj)
      ismij = Mul(ismi,ismj)
      if (ismij /= im) cycle
      just(lrj,lri) = nns
      nns = nns+1
    end do
  end do
end do
!write(u6,*) 'ns_sm',ns_sm
!write(u6,*) just(1:4,1:4)
!do i=norb_frz+1,norb_dz                                     !to del
!  lmi = lsm_inn(i)                                          !to del
!  do j=norb_frz+1,norb_dz                                   !to del
!    lmj = lsm_inn(j)                                        !to del
!    lsml(i,j) = Mul(lmi,lmj)                            !to del
!  end do                                                    !to del
!end do                                                      !to del
!write(nf2,*) '   jud ...'                                   !to del
!write(nf2,'(2x,12i8)') (jud(lr),lr=norb_frz+1,norb_dz)      !to del
!write(nf2,*) '   just ...'                                  !to del
!do i=norb_frz+1,norb_dz                                     !to del
!  write(nf2,'(2x,8i3)') (lsml(i,lr),lr=norb_frz+1,norb_dz)  !to del
!  write(nf2,'(2x,8i3)') (just(i,lr),lr=norb_frz+1,norb_dz)  !to del
!  write(nf2,*)
!end do                                                      !to del
!----------- norb_dbl<>0 -----------------------------------------------

return

end subroutine dbl_downwalk
