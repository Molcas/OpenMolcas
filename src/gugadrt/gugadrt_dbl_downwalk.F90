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

subroutine gugadrt_dbl_downwalk()
!     juv,just(nost,nost),jud(nost)
!     |  \  1         |
!     | d,dd,s(i=i)   |
!     |    \ s,t,tt(i<j)|
!     |     \       1 2 |     deal with inner of dbl_space
!     |ss(i>j)\       |
!     |  2 1  \       |

use gugadrt_global, only: iseg_sta, iseg_downwei, lsm_inn, max_innorb, ng_sm, norb_dbl, norb_dz, norb_frz, ns_sm
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: im, ismi, ismij, ismj, lr0, lri, lrj, nnd, nns, nnt
integer(kind=iwp), allocatable :: jud(:), just(:,:)

call mma_allocate(jud,max_innorb,label='jud')
call mma_allocate(just,max_innorb,max_innorb,label='just')

if (norb_dbl == 0) then
  !----------- norb_dbl=0 ------------------------------------------------
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
    do lrj=lri+1,norb_dz !tmp
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
    do lrj=lri+1,norb_dz !tmp
      ismj = lsm_inn(lrj)
      ismij = Mul(ismi,ismj)
      if (ismij /= im) cycle
      just(lrj,lri) = nns
      nns = nns+1
    end do
  end do
end do

call mma_deallocate(jud)
call mma_deallocate(just)

return

end subroutine gugadrt_dbl_downwalk
