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
subroutine SOGthr(SOInt,iBas,jBas,nSOInt,PrpInt,nPrp,lOper,iCmp,jCmp,iShell,jShell,iAO,jAO)
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

use Index_Functions, only: iTri
use SOAO_Info, only: iAOtSO
use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Constants, only: One, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iBas, jBas, nSOInt, nPrp, lOper, iCmp, jCmp, iShell, jShell, iAO, jAO
real(kind=wp), intent(out) :: SOInt(iBas*jBas,nSOInt)
real(kind=wp), intent(in) :: PrpInt(nPrp)
integer(kind=iwp) :: i1, i2, indAO1, indAO2, Indi, Indj, ipij, iPnt, iSO1, iSO2, j1, j12, j2, jjMx, lSO, nRow
real(kind=wp) :: Fact
integer(kind=iwp), external :: iPntSO

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt(' In SOGthr: PrpInt',' ',PrpInt,1,nPrp)
write(u6,*) ' iAO, jAO=',iAO,jAO
write(u6,*) ' iShell, jShell=',iShell,jShell
#endif
lSO = 0
do j1=0,nIrrep-1
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle

    ! Gather the SO's from lower rectangular blocks and triangular
    ! diagonal blocks.

    do j2=0,j1
      j12 = ieor(j1,j2)
      if (.not. btest(lOper,j12)) cycle
      jjMx = jCmp
      if ((iShell == jShell) .and. (j1 == j2)) jjMx = i1
      do i2=1,jjMx
        if (iAOtSO(jAO+i2,j2) < 0) cycle
        lSO = lSO+1
        iSO1 = iAOtSO(iAO+i1,j1)
        iSO2 = iAOtSO(jAO+i2,j2)

        iPnt = iPntSO(j1,j2,lOper,nbas)
        !write(u6,*) iSO1,iSO2,iPnt,i1,j1,i2,j2
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
    end do

  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt(' In SOGthr: SOInt',' ',SOInt,iBas*jBas,nSOInt)
#endif

return

end subroutine SOGthr
