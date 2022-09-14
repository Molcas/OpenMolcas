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

subroutine PLF_RI_2(AOint,ijkl,jCmp,lCmp,iAO,iAOst,jBas,lBas,kOp,TInt,nTInt,iSO2Ind,iOffA,nSOs)
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
!          Modified to 2-center RI June '05                            *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use SOAO_Info, only: iAOtSO
use Basis_Info, only: nBas
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ijkl, jCmp, lCmp, iAO(4), iAOst(4), jBas, lBas, kOp(4), nTInt, nSOs, iSO2Ind(nSOs), iOffA(4)
real(kind=wp), intent(in) :: AOint(ijkl,jCmp,lCmp)
real(kind=wp), intent(_OUT_) :: TInt(nTInt)
integer(kind=iwp) :: i2, i4, iAOj, iAOl, iAOstj, iAOstl, ij, iOff, iOffA_, iSO, jSO, jSOj, kSO, lSO, lSOl, mm_, mx, nijkl, nn

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
  write(u6,*) ' Sum=',r1
  write(u6,*) ' Dot=',r2
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
mx = nTri_Elem(nn)

#ifdef _DEBUGPRINT_
write(u6,*) 'nn,mx=',nn,mx
write(u6,*) 'iOff=',nn,mx
write(u6,*) 'lBas,jBas=',lBas,jBas
write(u6,*) 'lCmp,jCmp=',lCmp,jCmp
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

        iSO = iSO2Ind(iSO)+nn
        ij = iTri(iSO,kSO)-mx+iOffA_
        TInt(ij) = AOint(nijkl,i2,i4)

      end do
    end do

  end do
end do

end subroutine PLF_RI_2
