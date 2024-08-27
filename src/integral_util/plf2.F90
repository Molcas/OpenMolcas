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
subroutine PLF2(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
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
use k2_arrays, only: Sew_Scr
use lw_Info, only: lwSyB, lwInt, lwSqn
use Gateway_Info, only: ThrInt
use sort_data, only: DimSyB, lSll
use Constants, only: One

implicit none
integer ijkl, iCmp, jCmp, kCmp, lCmp, iBas, jBas, kBas, lBas
real*8 AOint(ijkl,iCmp,jCmp,kCmp,lCmp)
integer iShell(4), iAO(4), kOp(4), iAOst(4), iSOs(4)
integer i, j, iTri, nUt, iAOSti, iAOStj, iAOStk, iAOStl, iAOi, iAOj, iAOk, iAOl, i1, i2, i3, i4, nij, mij, iSO, jSO, kSO, lSO, &
        iSOi, jSOj, kSOk, lSOl, nijkl, iSOij, iSOkl, ijklCmp, iBin
real*8 AInt
#ifdef _DEBUGPRINT_
real*8 r1, r2
real*8, external :: DDot_
#endif
! Statement function
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

#ifdef _DEBUGPRINT_
r1 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,[One],0)
r2 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
write(6,*) ' Sum=',r1
write(6,*) ' Dot=',r2
call RecPrt(' In Plf2: AOInt',' ',AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)
#endif

! Allocate space to store integrals together with their
! Symmetry batch and sequence number.
! To avoid conflicts in using memory this is done in the
! subroutine PSOAO

nUt = -1
nij = DimSyB(1,1)
mij = lSll(1)/nij
!write(6,*) 'nij,mij=',nij,mij

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
call DCopy_(ijkl*2*ijklCmp,[One],0,Sew_Scr(lwSyB),1)

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
                AInt = AOint(nijkl,i1,i2,i3,i4)
                if (abs(AInt) < ThrInt) Go To 420
                iSOij = iTri(iSOi,jSOj)

                !write(6,*) 'iSOij,iSOkl=',iSOij,iSOkl

                nUt = nUt+1
                Sew_Scr(lwInt+nUt) = Aint
                iBin = (iSOkl-1)/mij
                !write(6,*) 'iBin=',iBin+1
                Sew_Scr(lwSyB+nUt) = dble(iBin+1)
                Sew_Scr(lwSqN+nUt) = dble((iSOkl-1-iBin*mij)*nij+iSOij)
                !write(6,*) 'iSq=',Sew_Scr(lwSqN+nUt)

                if (iSOij /= iSOkl) then
                  nUt = nUt+1
                  Sew_Scr(lwInt+nUt) = Aint
                  iBin = (iSOij-1)/mij
                  !write(6,*) 'iBin=',iBin+1
                  Sew_Scr(lwSyB+nUt) = dble(iBin+1)
                  Sew_Scr(lwSqN+nUt) = dble((iSOij-1-iBin*mij)*nij+iSOkl)
                  !write(6,*) 'iSq=',Sew_Scr(lwSqN+nUt)
                end if

420             continue
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
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(iShell)

end subroutine PLF2
