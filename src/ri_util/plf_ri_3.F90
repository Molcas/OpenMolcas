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

subroutine PLF_RI_3(AOint,ijkl,jCmp,kCmp,lCmp,iShell,iAO,iAOst,jBas,kBas,lBas,kOp,TInt,nTInt,iOff,iShlSO,nBasSh,iSOShl,nSO,nShell, &
                    nSym,iSSOff)
!***********************************************************************
!                                                                      *
!  object: to sift and index the petite list format integrals.         *
!                                                                      *
!          the indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          May '90                                                     *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri
use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ijkl, jCmp, kCmp, lCmp, iShell(4), iAO(4), iAOst(4), jBas, kBas, lBas, kOp(4), nTInt, iOff(3), &
                                 nSO, iShlSO(nSO), nShell, nSym, nBasSh(0:nSym-1,nShell), iSOShl(nSO), iSSOff
real(kind=wp), intent(in) :: AOint(ijkl,jCmp,kCmp,lCmp)
real(kind=wp), intent(_OUT_) :: TInt(nTInt)
integer(kind=iwp) :: i2, i3, i4, iAux, iC, iD, iOff1, iShC, iSOs(4), jSOj, kl, kl_B, kSOk, lCmp_Max, lSOl, n3C, nC, nijkl
logical(kind=iwp) :: Shkl

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
Shkl = iShell(3) == iShell(4)
iOff1 = nBas(0)
n3C = iOff(3)
if (iShell(4) > iShell(3)) then
  write(u6,*) 'iShell(4) > iShell(3)'
  call Abend()
end if

do i2=1,jCmp
  iSOs(2) = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
  do i3=1,kCmp
    iSOs(3) = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
    lCmp_Max = lCmp
    if (Shkl) lCmp_Max = i3
    do i4=1,lCmp_Max
      iSOs(4) = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
      !                                                                *
      !*****************************************************************
      !                                                                *
      if (Shkl .and. (i3 == i4)) then
        !                                                              *
        !***************************************************************
        !                                                              *
        nijkl = 0
        do lSOl=iSOs(4),iSOs(4)+lBas-1
          iD = iShlSO(lSOl)                  ! Relative index
          do kSOk=iSOs(3),iSOs(3)+kBas-1
            iC = iShlSO(kSOk)
            iShC = iSOShl(kSOk)
            nC = nBasSh(0,iShC)

            kl = iTri(iC,iD)+iSSOff

            do jSOj=iSOs(2),iSOs(2)+jBas-1
              nijkl = nijkl+1
              if (lSOl > kSOk) cycle
              iAux = jSOj-iOff1

              kl_B = (iAux-1)*n3C+kl
              TInt(kl_B) = AOint(nijkl,i2,i3,i4)

            end do
          end do
        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      else
        !                                                              *
        !***************************************************************
        !                                                              *
        nijkl = 0
        do lSOl=iSOs(4),iSOs(4)+lBas-1
          iD = iShlSO(lSOl)
          do kSOk=iSOs(3),iSOs(3)+kBas-1
            iC = iShlSO(kSOk)
            iShC = iSOShl(kSOk)
            nC = nBasSh(0,iShC)

            if (Shkl) then
              kl = iTri(iC,iD)
            else
              kl = (iD-1)*nC+iC
            end if
            kl = kl+iSSOff

            do jSOj=iSOs(2),iSOs(2)+jBas-1
              iAux = jSOj-iOff1
              nijkl = nijkl+1

              kl_B = (iAux-1)*n3C+kl
              TInt(kl_B) = AOint(nijkl,i2,i3,i4)

            end do
          end do
        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine PLF_RI_3
