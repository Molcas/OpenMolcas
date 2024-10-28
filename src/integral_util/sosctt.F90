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
subroutine SOSctt(SOInt,iBas,jBas,nSOInt,PrpInt,nPrp,lOper,iCmp,jCmp,iShell,jShell,iAO,jAO,rHrmt)
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
real(kind=wp) , intent(in):: SOInt(iBas*jBas,nSOInt), rHrmt
real(kind=wp) , intent(inout):: PrpInt(nPrp)
integer(kind=iwp) :: i1, i2, indAO1, indAO2, Indi, Indj, ip, iPnt, iSO1, iSO2, j1, j12, j2, jBsMax, jjMx, lSO, nRow
integer(kind=iwp), external :: iPntSO

#ifdef _DEBUGPRINT_
call RecPrt(' In SOSctt:SOInt',' ',SOInt,iBas*jBas,nSOInt)
call RecPrt(' In SOSctt:PrpInt',' ',PrpInt,1,nPrp)
#endif

lSO = 0
do j1=0,nIrrep-1
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle

    ! Scatter the SO's onto lower rectangular blocks and triangular
    ! diagonal blocks.

    do j2=0,nIrrep-1
      j12 = ieor(j1,j2)
      if (.not. btest(lOper,j12)) cycle
      jjMx = jCmp
      if ((iShell == jShell) .and. (j1 == j2)) jjMx = i1
      do i2=1,jjMx
        if (iAOtSO(jAO+i2,j2) < 0) cycle
        lSO = lSO+1
        iSO1 = iAOtSO(iAO+i1,j1)
        iSO2 = iAOtSO(jAO+i2,j2)
        !write(u6,*) iSO1,iAO,i1,j1,iSO2,jAO,i2,j2

        iPnt = iPntSO(max(j1,j2),min(j1,j2),lOper,nbas)
        do indAO1=0,iBas-1
          ! Diagonal block. Store only unique elements
          jBsMax = jBas-1
          if ((j1 == j2) .and. (iSO1 == iSO2)) jBsMax = indAO1
          do indAO2=0,jBsMax
            ip = indAO2*iBas+indAO1+1

            ! Move one-electron integral.

            if (j1 == j2) then
              ! Diagonal symmetry block
              Indi = iSO1+indAO1
              Indj = iSO2+indAO2
              PrpInt(iPnt+iTri(Indi,Indj)) = SOInt(ip,lSO)
            else
              ! Off-diagonal symmetry block j1>j2
              if (j1 > j2) then
                Indi = iSO1+indAO1
                Indj = iSO2+indAO2
                nRow = nBas(j1)
                PrpInt(iPnt+nRow*(Indj-1)+Indi) = SOInt(ip,lSO)
              else
                Indj = iSO1+indAO1
                Indi = iSO2+indAO2
                nRow = nBas(j2)
                PrpInt(iPnt+nRow*(Indj-1)+Indi) = rHrmt*SOInt(ip,lSO)
              end if
            end if

          end do
        end do

      end do
    end do

  end do
end do

return

end subroutine SOSctt
