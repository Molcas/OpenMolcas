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

subroutine PLF_Cho_Diag(TInt,nInt,AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp)
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
!                                                                      *
!***********************************************************************

use SOAO_Info, only: iAOtSO
use ChoArr, only: iSOSHl, iShlSO, nBstSh
use Constants

implicit real*8(A-H,O-Z)
#include "cholesky.fh"
#include "print.fh"
real*8 AOint(ijkl,iCmp,jCmp,kCmp,lCmp), TInt(nInt)
integer iShell(4), iAO(4), kOp(4), iAOst(4), iSOs(4)
logical Shijij
external ddot_
! Statement function
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

irout = 109
jprint = nprint(irout)
if (jPrint >= 49) then
  r1 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,[One],0)
  r2 = DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
  write(6,*) ' Sum=',r1
  write(6,*) ' Dot=',r2
end if
if (jPrint >= 99) call RecPrt(' In Plf_CD: AOInt',' ',AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)

! Allocate space to store integrals to gether with their
! Symmetry batch and sequence number.
! To avoid conflicts in using memory this is done in the
! subroutine PSOAO

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
                iSOij = iTri(iSOi,jSOj)

                if (iSOij == iSOkl) then
                  ISHLI = ISOSHL(ISOI)
                  ISHLJ = ISOSHL(JSOJ)
                  NUMI = NBSTSH(ISHLI)
                  NUMJ = NBSTSH(ISHLJ)
                  if ((ISHLI == ISHLJ) .and. (ISHLI == SHA)) then
                    KIJ = ITRI(ISHLSO(ISOI),ISHLSO(JSOJ))
                  else
                    if ((ISHLI == SHA) .and. (ISHLJ == SHB)) then
                      KIJ = NUMI*(ISHLSO(JSOJ)-1)+ISHLSO(ISOI)
                    else if ((ISHLJ == SHA) .and. (ISHLI == SHB)) then
                      KIJ = NUMJ*(ISHLSO(ISOI)-1)+ISHLSO(JSOJ)
                    else
                      call CHO_QUIT('Integral error',104)
                      KIJ = -999999
                    end if
                  end if
                  TInt(KIJ) = AOint(nijkl,i1,i2,i3,i4)
                end if

              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(iShell)
  call Unused_logical(Shijij)
end if

end subroutine PLF_Cho_Diag
