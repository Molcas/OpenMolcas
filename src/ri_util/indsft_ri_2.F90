!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine IndSft_RI_2(iCmp,iShell,jBas,lBas,Shijij,iAO,iAOst,ijkl,SOint,nSOint,nSOs,TInt,nTInt,iSO2Ind,iOffA)
!***********************************************************************
!  object: to sift and index the SO integrals.                         *
!                                                                      *
!          the indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          april '90                                                   *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use Symmetry_Info, only: Mul, nIrrep
use sort_data, only: nSkip
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Constants, only: Zero, One
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), iShell(4), jBas, lBas, iAO(4), iAOst(4), ijkl, nSOint, nSOs, nTInt, iSO2Ind(nSOs), &
                                 iOffA(4,0:7)
logical(kind=iwp), intent(in) :: Shijij
real(kind=wp), intent(in) :: SOint(ijkl,nSOint)
real(kind=wp), intent(inout) :: TInt(nTInt)
integer(kind=iwp) :: i1, i12, i2, i3, i34, i4, ij, iOffA_, iOffB_, iSO, ix, j, j1, j12, j2, j3, j4, jSO, jSOj, jSym(0:7), k12, &
                     k34, kSO, lSO, lSOl, lSym(0:7), memSO2, mm_, mx, nijkl, nn
logical(kind=iwp) :: qijij
#ifdef _DEBUGPRINT_
#include "print.fh"
integer(kind=iwp) :: iprint, irout
real(kind=wp) :: r1, r2, tr1 = Zero, tr2 = Zero
real(kind=wp), external :: ddot_
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
k12 = 0
k34 = 0
#ifdef _DEBUGPRINT_
irout = 39
iprint = nprint(irout)
if (iPrint >= 49) then
  r1 = DDot_(ijkl*nSOInt,SOInt,1,[One],0)
  r2 = DDot_(ijkl*nSOInt,SOInt,1,SOInt,1)
  tr1 = tr1+r1
  tr2 = tr2+r2
  write(u6,*) ' Sum=',r1,tr1
  write(u6,*) ' Dot=',r2,tr2
  call RecPrt(' in indsft:SOint ',' ',SOint,ijkl,nSOint)
end if
#endif
memSO2 = 0

! quadruple loop over elements of the basis functions angular
! description. loops are reduced to just produce unique SO integrals
! observe that we will walk through the memory in AOint in a
! sequential way.

i1 = 1
i3 = 1
j1 = 0
j3 = 0
do i2=1,iCmp(2)
  do j=0,nIrrep-1
    ix = 0
    if (iAOtSO(iAO(2)+i2,j) > 0) ix = 2**j
    jSym(j) = ix
  end do
  if (iShell(2) > iShell(1)) then
    i12 = iCmp(2)*(i1-1)+i2
  else
    i12 = iCmp(1)*(i2-1)+i1
  end if
  do i4=1,iCmp(4)
    do j=0,nIrrep-1
      ix = 0
      if (iAOtSO(iAO(4)+i4,j) > 0) ix = 2**j
      lSym(j) = ix
    end do
    if (iShell(4) > iShell(3)) then
      i34 = iCmp(4)*(i3-1)+i4
    else
      i34 = iCmp(3)*(i4-1)+i3
    end if
    if (Shijij .and. (i34 > i12)) cycle
    qijij = Shijij .and. (i12 == i34)
    !write(u6,*) 'i1,i2,i3,i4=',i1,i2,i3,i4

    ! loop over Irreps which are spanned by the basis function.
    ! again, the loop structure is restricted to ensure unique
    ! integrals.

    do j2=0,nIrrep-1
      if (jSym(j2) == 0) cycle
      j12 = Mul(j1+1,j2+1)-1
      if (qijij) then
        if (iShell(1) > iShell(2)) then
          k12 = nIrrep*j1+j2+1
        else
          k12 = nIrrep*j2+j1+1
        end if
      end if

      iOffA_ = iOffA(1,j2)
      iOffB_ = iOffA(3,j2)
      if (j2 /= 0) iOffB_ = iOffB_+1
      mm_ = iOffA(4,j2)
      nn = mm_-iOffA(2,j2)
      mx = nTri_Elem(nn)

      j4 = Mul(j12+1,j3+1)-1
      if (lSym(j4) == 0) cycle
      if (qijij) then
        if (iShell(3) > iShell(4)) then
          k34 = nIrrep*j3+j4+1
        else
          k34 = nIrrep*j4+j3+1
        end if
        if (k34 > k12) cycle
      end if

      memSO2 = memSO2+1
      if ((nSkip(j2+1)+nSkip(j4+1)) /= 0) cycle

      ! Compute absolute starting SO index
      jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
      lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)

      nijkl = 0
      do lSOl=lSO,lSO+lBas-1
        do jSOj=jSO,jSO+jBas-1
          nijkl = nijkl+1
          iSO = jSOj-nBas(j2)
          kSO = lSOl-nBas(j4)

          iSO = iSO2Ind(iSO+iOffB_)+nn
          ij = iTri(iSO,kSO)-mx+iOffA_
          TInt(ij) = SOint(nijkl,memSO2)

        end do
      end do

    end do

  end do
end do

return

end subroutine IndSft_RI_2
