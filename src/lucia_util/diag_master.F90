!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine diag_master(nTU,TU,nTUVX,TUVX)
! To do in this subroutine:
!
! - Make sure all calling parameters are accounted for
! - Make sure mxntts is set in previous subroutine
! - Make sure nsmst is set in previous subroutine
!
! Set up the diagonal for the CI calculation

use lucia_data, only: INT1, nIrrep, NGAS, NGSSH
use CandS, only: ISSM
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTU, nTUVX
real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
integer(kind=iwp) :: IADD, iSym, ITU, MTU, NAT, NT, NU

INT1(:) = Zero
MTU = 0
ITU = 0
IADD = 0
do iSym=1,nIrrep
  NAT = sum(NGSSH(iSym,1:nGAS))
  if (NAT == 0) cycle
  do NT=1,NAT
    MTU = MTU+IADD
    do NU=1,NT
      MTU = MTU+1
      ITU = ITU+1
      INT1(ITU) = TU(MTU)
    end do
  end do
  IADD = IADD+NAT
end do

call GASCI(ISSM,1,nTUVX,TUVX)

end subroutine diag_master
