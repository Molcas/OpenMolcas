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
subroutine PLF2(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
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

use Index_Functions, only: iTri
use SOAO_Info, only: iAOtSO
use k2_arrays, only: Sew_Scr
use lw_Info, only: lwInt, lwSqn, lwSyB
use Gateway_Info, only: ThrInt
use sort_data, only: DimSyB, lSll
use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ijkl, iCmp, jCmp, kCmp, lCmp, iAO(4), iAOst(4), iBas, jBas, kBas, lBas, kOp(4)
real(kind=wp), intent(in) :: AOint(ijkl,iCmp,jCmp,kCmp,lCmp)
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, iAOj, iAOk, iAOl, iAOSti, iAOStj, iAOStk, iAOStl, iBin, ijklCmp, iSO, iSOi, iSOij, &
                     iSOkl, iSOs(4), jSO, jSOj, kSO, kSOk, lSO, lSOl, mij, nij, nijkl, nUt
real(kind=wp) :: A_Int
#ifdef _DEBUGPRINT_
real(kind=wp) :: r1, r2
real(kind=wp), external :: DDot_

r1 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,[One],0)
r2 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
write(u6,*) ' Sum=',r1
write(u6,*) ' Dot=',r2
call RecPrt(' In Plf2: AOInt',' ',AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)
#endif

! Allocate space to store integrals together with their
! Symmetry batch and sequence number.
! To avoid conflicts in using memory this is done in the
! subroutine PSOAO

nUt = -1
nij = DimSyB(1,1)
mij = lSll(1)/nij
!write(u6,*) 'nij,mij=',nij,mij

! quadruple loop over elements of the basis functions angular
! description. loops are reduced to just produce unique SO integrals
! observe that we will walk through the memory in AOint in a
! sequential way.

iAOsti = iAOst(1)
iAOstj = iAOst(2)
iAOstk = iAOst(3)
iAOstl = iAOst(4)
iAOi = iAO(1)
iAOj = iAO(2)
iAOk = iAO(3)
iAOl = iAO(4)

ijklCmp = iCmp*jCmp*kCmp*lCmp
Sew_Scr(lwSyB:lwSyB+ijkl*2*ijklCmp-1) = One

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
            iSOkl = iTri(kSOk,lSOl)
            do jSOj=jSO,jSO+jBas-1
              do iSOi=iSO,iSO+iBas-1
                nijkl = nijkl+1
                A_Int = AOint(nijkl,i1,i2,i3,i4)
                if (abs(A_Int) < ThrInt) cycle
                iSOij = iTri(iSOi,jSOj)

                !write(u6,*) 'iSOij,iSOkl=',iSOij,iSOkl

                nUt = nUt+1
                Sew_Scr(lwInt+nUt) = A_int
                iBin = (iSOkl-1)/mij
                !write(u6,*) 'iBin=',iBin+1
                Sew_Scr(lwSyB+nUt) = real(iBin+1,kind=wp)
                Sew_Scr(lwSqN+nUt) = real((iSOkl-1-iBin*mij)*nij+iSOij,kind=wp)
                !write(u6,*) 'iSq=',Sew_Scr(lwSqN+nUt)

                if (iSOij /= iSOkl) then
                  nUt = nUt+1
                  Sew_Scr(lwInt+nUt) = A_int
                  iBin = (iSOij-1)/mij
                  !write(u6,*) 'iBin=',iBin+1
                  Sew_Scr(lwSyB+nUt) = real(iBin+1,kind=wp)
                  Sew_Scr(lwSqN+nUt) = real((iSOij-1-iBin*mij)*nij+iSOkl,kind=wp)
                  !write(u6,*) 'iSq=',Sew_Scr(lwSqN+nUt)
                end if

              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do

! pass the integral to phase 1 of the bin sorting algorithm

call SORT1A(nUt+1,Sew_Scr(lwInt),Sew_Scr(lwSqN),Sew_Scr(lwSyB))
nUt = 0

return

end subroutine PLF2
