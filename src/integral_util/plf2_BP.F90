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

!#define _DEBUGPRINT_
subroutine PLF2_BP(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
!***********************************************************************
!                                                                      *
!  object: to sift and index the petite list format integrals.         *
!                                                                      *
!          the indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          May '90                                                     *
!***********************************************************************

use SOAO_Info, only: iAOtSO
use Definitions, only: wp, iwp
use Breit, only: D_tensor, PAO
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ijkl, iCmp, jCmp, kCmp, lCmp, iAO(4), iAOst(4), iBas, jBas, kBas, lBas, kOp(4)
real(kind=wp), intent(in) :: AOint(ijkl,6,iCmp,jCmp,kCmp,lCmp)
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, iAOj, iAOk, iAOl, iAOSti, iAOStj, iAOStk, iAOStl, iSO, iSOi, &
                     iSOs(4), jSO, jSOj, kSO, kSOk, lSO, lSOl, nijkl
real(kind=wp) :: Prod_ij

! quadruple loop over elements of the basis functions angular
! description. loops are reduced to just produce unique SO integrals
! observe that we will walk through the memory in AOint in a
! sequential way.


If (.Not.Associated(PAO)) Call Abend()

iAOsti = iAOst(1)
iAOstj = iAOst(2)
iAOstk = iAOst(3)
iAOstl = iAOst(4)
iAOi = iAO(1)
iAOj = iAO(2)
iAOk = iAO(3)
iAOl = iAO(4)

do i1=1,iCmp
  iSOs(1) = iAOtSO(iAOi+i1,kOp(1))+iAOsti
  do i2=1,jCmp
    iSOs(2) = iAOtSO(iAOj+i2,kOp(2))+iAOstj
    do i3=1,kCmp
      iSOs(3) = iAOtSO(iAOk+i3,kOp(3))+iAOstk
      do i4=1,lCmp
        iSOs(4) = iAOtSO(iAOl+i4,kOp(4))+iAOstl

        iSO = iSOs(1)
        jSO = iSOs(2)
        kSO = iSOs(3)
        lSO = iSOs(4)

        nijkl = 0
        do lSOl=lSO,lSO+lBas-1
          do kSOk=kSO,kSO+kBas-1
            do jSOj=jSO,jSO+jBas-1
              do iSOi=iSO,iSO+iBas-1
                nijkl = nijkl+1
                Prod_ij = AOint(nijkl,1,i1,i2,i3,i4)*PAO(nijkl,i1,i2,i3,i4)
                D_tensor(1,1) = D_tensor(1,1) + Prod_ij
                Prod_ij = AOint(nijkl,2,i1,i2,i3,i4)*PAO(nijkl,i1,i2,i3,i4)
                D_tensor(2,1) = D_tensor(2,1) + Prod_ij
                Prod_ij = AOint(nijkl,3,i1,i2,i3,i4)*PAO(nijkl,i1,i2,i3,i4)
                D_tensor(3,1) = D_tensor(3,1) + Prod_ij
                Prod_ij = AOint(nijkl,4,i1,i2,i3,i4)*PAO(nijkl,i1,i2,i3,i4)
                D_tensor(2,2) = D_tensor(2,2) + Prod_ij
                Prod_ij = AOint(nijkl,5,i1,i2,i3,i4)*PAO(nijkl,i1,i2,i3,i4)
                D_tensor(2,3) = D_tensor(2,3) + Prod_ij
                Prod_ij = AOint(nijkl,6,i1,i2,i3,i4)*PAO(nijkl,i1,i2,i3,i4)
                D_tensor(3,3) = D_tensor(3,3) + Prod_ij


              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do

end subroutine PLF2_BP
