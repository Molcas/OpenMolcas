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

subroutine PLF_RICD(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,TInt,nTInt,mTInt,iTOff,iOffij,iOffkl)
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
use SOAO_Info, only: iAOtSO
use Definitions, only: wp, iwp
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Constants, only: One
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ijkl, iCmp, jCmp, kCmp, lCmp, iAO(4), iAOst(4), iBas, jBas, kBas, lBas, kOp(4), nTInt, mTInt, &
                                 iTOff, iOffij, iOffkl
real(kind=wp), intent(in) :: AOint(ijkl,iCmp,jCmp,kCmp,lCmp)
real(kind=wp), intent(_OUT_) :: TInt(nTInt,mTInt)
#include "ibas_ricd.fh"
integer(kind=iwp) :: i1, i2, i3, i4, iAOi, iAOj, iAOk, iAOl, iAOsti, iAOstj, iAOstk, iAOstl, ijSOij, iSO, iSOi, iSOij, iSOkl, &
                     iSOs(4), jSO, jSOj, klSOkl, kSO, kSOk, lSO, lSOl, nijkl
#ifdef _DEBUGPRINT_
real(kind=wp) :: r1, r2
real(kind=wp), external :: ddot_
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
r1 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,[One],0)
r2 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
write(u6,*) ' Sum=',r1
write(u6,*) ' Dot=',r2
call RecPrt(' In PLF_RICD: AOInt',' ',AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)
#endif

! quadruple loop over elements of the basis functions angular
! description. loops are reduced to just produce unique SO integrals
! observe that we will walk through the memory in AOint in a
! sequential way.

iAOsti = iAOst(1)
iAOstj = iAOst(2)
iAOstk = iAOst(3)
iAOstl = iAOst(4)
!write(u6,*) 'iAOsti,iAOstj,iAOstk,iAOstl=',iAOsti,iAOstj,iAOstk,iAOstl
iAOi = iAO(1)
iAOj = iAO(2)
iAOk = iAO(3)
iAOl = iAO(4)
!write(u6,*) 'iAOs=',iAO
!write(u6,*) 'kOps=',kOp
!write(u6,*) 'iTOff,iOffij,iOffkl=',iTOff,iOffij,iOffkl
!write(u6,*) 'iBas,jBas,kBas,lBas=',iBas,jBas,kBas,lBas

! The writing of the integrals here are shell blocked.

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

        !write(u6,*)
        !write(u6,*) 'i1,i2,i3,i4,iSOs=',i1,i2,i3,i4,iSOs
        !write(u6,*) 'iBas,jBas,kBas,lBas=',iBas,jBas,kBas,lBas

        nijkl = 0
        do lSOl=lSO,lSO+lBas-1
          do kSOk=kSO,kSO+kBas-1
            if (iAO(3) == iAO(4)) then
              iSOkl = iTri(kSOk,lSOl)+iOffkl
            else
              iSOkl = (kSOk-1)*lCmp*lBas_+lSOl+iOffkl
            end if
            do jSOj=jSO,jSO+jBas-1
              do iSOi=iSO,iSO+iBas-1
                nijkl = nijkl+1
                if (iAO(1) == iAO(2)) then
                  iSOij = iTri(iSOi,jSOj)+iOffij
                else
                  iSOij = (iSOi-1)*jCmp*jBas_+jSOj+iOffij
                end if

                !write(u6,*) 'iSOij,iSOkl,AOint=',iSOij,iSOkl,AOint(nijkl,i1,i2,i3,i4)
                ijSOij = max(iSOij,iSOkl)-iTOff
                klSOkl = min(iSOij,iSOkl)
                TInt(klSOkl,ijSOij) = AOint(nijkl,i1,i2,i3,i4)

              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('TInt','(45G8.2)',TInt,nTInt,mTInt)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine PLF_RICD
