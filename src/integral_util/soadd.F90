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
subroutine SOAdd(SOInt,iBas,jBas,nSOInt,PrpInt,nPrp,lOper,iCmp,jCmp,iShell,jShell,iAO,jAO)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '91                                              *
!***********************************************************************

use Index_Functions, only: iTri
use SOAO_Info, only: iAOtSO
use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iBas, jBas, nSOInt, nPrp, lOper, iCmp, jCmp, iShell, jShell, iAO, jAO
real(kind=wp), intent(in) :: SOInt(iBas*jBas,nSOInt)
real(kind=wp), intent(inout) :: PrpInt(nPrp)
integer(kind=iwp) :: i1, i2, IndAO1, IndAO2, Indij, ip, iPnt, iSO, iSO1, iSO2, j1, j12, j2, jSO, lSO, nRow
integer(kind=iwp), external :: iPntSO

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt(' In SOAdd:SOInt',' ',SOInt,iBas*jBas,nSOInt)
#endif

lSO = 0
do j1=0,nIrrep-1
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle

    ! Scatter the SO's onto lower rectangular blocks and triangular
    ! diagonal blocks.

    do j2=0,j1
      j12 = ieor(j1,j2)
      if (.not. btest(lOper,j12)) cycle

      do i2=1,jCmp
        if (iAOtSO(jAO+i2,j2) < 0) cycle
        if ((iShell == jShell) .and. (j1 == j2) .and. (i1 < i2)) cycle

        lSO = lSO+1
        iSO1 = iAOtSO(iAO+i1,j1)
        iSO2 = iAOtSO(jAO+i2,j2)

        iPnt = iPntSO(j1,j2,lOper,nbas)
        do indAO1=1,iBas
          do indAO2=1,jBas
            ip = (indAO2-1)*iBas+indAO1

            ! Move one electron integral.

            iSO = iSO1+IndAO1-1
            jSO = iSO2+IndAO2-1

            ! Diagonal block. Store only unique elements
            if ((j1 == j2) .and. (iSO1 == iSO2) .and. (iSO < jSO)) cycle

            if (j1 == j2) then
              ! Diagonal symmetry block
              Indij = iPnt+iTri(iSO,jSO)
            else
              ! Off-diagonal symmetry block j1>j2
              nRow = nBas(j1)
              Indij = iPnt+nRow*(jSO-1)*nRow+iSO
            end if

            PrpInt(Indij) = PrpInt(Indij)+SOInt(ip,lSO)

          end do
        end do

      end do
    end do

  end do
end do

return

end subroutine SOAdd
