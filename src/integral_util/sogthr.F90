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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine SOGthr(SOInt,iBas,jBas,nSOInt,PrpInt,nPrp,lOper,iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)
!***********************************************************************
!                                                                      *
! Object: to gather elements, from the Fock or 1st order density matrix*
!         in SO-basis, which are associated with a shell pair.         *
!         OBSERVE that the matrix is folded to triangular form, i.e.   *
!         the off diagonal elements has twice there value. This in     *
!         order to reduce the summation to i>=j. However, this will    *
!         not be used for the diagonal blocks. Hence, these elements   *
!         will only be assigned half the value from the original       *
!         matrix.                                                      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use SOAO_Info, only: iAOtSO
use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Constants

implicit none
integer iBas, jBas, nSOInt, lOper, iCmp, jCmp, iShell, jShell, iAO, jAO, nPrp
real*8 SOInt(iBas*jBas,nSOInt), PrpInt(nPrp)
logical AeqB
integer, external :: iPntSO
integer i, j, iTri, lSO, j1, i1, j2, j12, jjMx, i2, iSO1, iSO2, indAO1, indAO2, ipij, Indi, Indj, iPnt, nRow
real*8 Fact
! Statement function
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt(' In SOGthr: PrpInt',' ',PrpInt,1,nPrp)
write(6,*) ' iAO, jAO=',iAO,jAO
write(6,*) ' iShell, jShell=',iShell,jShell
#endif
lSO = 0
do j1=0,nIrrep-1
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle

    ! Gather the SO's from lower rectangular blocks and triangular
    ! diagonal blocks.

    do j2=0,j1
      j12 = ieor(j1,j2)
      if (iand(lOper,2**j12) == 0) Go To 300
      jjMx = jCmp
      if ((iShell == jShell) .and. (j1 == j2)) jjMx = i1
      do i2=1,jjMx
        if (iAOtSO(jAO+i2,j2) < 0) cycle
        lSO = lSO+1
        iSO1 = iAOtSO(iAO+i1,j1)
        iSO2 = iAOtSO(jAO+i2,j2)

        iPnt = iPntSO(j1,j2,lOper,nbas)
        !write(6,*) iSO1,iSO2,iPnt,i1,j1,i2,j2
        do indAO1=0,iBas-1
          ! Diagonal block (only unique elements).
          do indAO2=0,jBas-1
            ipij = indAO2*iBas+indAO1+1
            Indi = iSO1+indAO1
            Indj = iSO2+indAO2

            ! Move one element and unfold.

            Fact = Half
            if (Indi == Indj) Fact = One

            if (j1 == j2) then
              SOInt(ipij,lSO) = Fact*PrpInt(iPnt+iTri(Indi,Indj))
            else
              nRow = nBas(j1)
              SOInt(ipij,lSO) = Fact*PrpInt(iPnt+nRow*(Indj-1)+Indi)
            end if

          end do
        end do

      end do
300   continue
    end do

  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt(' In SOGthr: SOInt',' ',SOInt,iBas*jBas,nSOInt)
#endif

return
! Avoid unused argument warnings
if (.false.) call Unused_logical(AeqB)

end subroutine SOGthr
