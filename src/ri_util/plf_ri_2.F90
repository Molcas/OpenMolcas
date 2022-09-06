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
! Copyright (C) 1990,2005, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine PLF_RI_2(AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,TInt,nTInt,iSO2Ind,iOffA,nSOs)
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
!          Modified to 2-center RI June '05                            *
!                                                                      *
!***********************************************************************

use SOAO_Info, only: iAOtSO
use Basis_Info, only: nBas

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 AOint(ijkl,jCmp,lCmp), TInt(nTInt)
integer iShell(4), iAO(4), kOp(4), iAOst(4), iSO2Ind(nSOs), iOffA(4)
logical Shijij
! Statement function
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
irout = 109
iPrint = nPrint(irout)
iPrint = 99
if (iPrint >= 49) then
  r1 = DDot_(ijkl*jCmp*lCmp,AOInt,1,[One],0)
  r2 = DDot_(ijkl*jCmp*lCmp,AOInt,1,AOInt,1)
  write(6,*) ' Sum=',r1
  write(6,*) ' Dot=',r2
end if
if (iPrint >= 99) call RecPrt(' In Plf_RI_2: AOInt',' ',AOInt,ijkl,jCmp*lCmp)
#endif

iAOstj = iAOst(2)
iAOstl = iAOst(4)
iAOj = iAO(2)
iAOl = iAO(4)
iOff = nBas(0)
iOffA_ = iOffA(1)
mm_ = iOffA(4)
nn = mm_-iOffA(2)
mx = nn*(nn+1)/2

#ifdef _DEBUGPRINT_
write(6,*) 'nn,mx=',nn,mx
write(6,*) 'iOff=',nn,mx
write(6,*) 'lBas,jBas=',lBas,jBas
write(6,*) 'lCmp,jCmp=',lCmp,jCmp
#endif

do i2=1,jCmp
  jSO = iAOtSO(iAOj+i2,kOp(2))+iAOstj
  do i4=1,lCmp
    lSO = iAOtSO(iAOl+i4,kOp(4))+iAOstl

    nijkl = 0
    do lSOl=lSO,lSO+lBas-1
      kSO = lSOl-iOff

      do jSOj=jSO,jSO+jBas-1

        iSO = jSOj-iOff
        nijkl = nijkl+1
        AInt = AOint(nijkl,i2,i4)

        iSO = iSO2Ind(iSO)+nn
        ij = iTri(iSO,kSO)-mx+iOffA_
        TInt(ij) = AInt

      end do
    end do

  end do
end do

#ifdef _WARNING_WORKAROUND_
if (.false.) then
  call Unused_integer(iCmp)
  call Unused_integer(kCmp)
  call Unused_integer_array(iShell)
  call Unused_logical(Shijij)
  call Unused_integer(iBas)
  call Unused_integer(kBas)
end if
#endif

end subroutine PLF_RI_2
