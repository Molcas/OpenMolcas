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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GetGrad_ER(Functional,GradNorm,R,CMO,nBasis,nOrb2Loc,Timing)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: compute ER Functional and its gradient norm.
!          The R matrix is computed as R(i,j) = (ij|jj).
!
! Note: symmetry is NOT allowed (but is not tested!).

use Data_structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use Constants, only: Zero, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBasis, nOrb2Loc
real(kind=wp), intent(out) :: Functional, GradNorm, R(nOrb2Loc,nOrb2Loc)
real(kind=wp), intent(in) :: CMO(nBasis,nOrb2Loc)
logical(kind=iwp), intent(in) :: Timing
integer(kind=iwp), parameter :: nSym = 1
integer(kind=iwp) :: i, irc, j, nOcc(nSym)
character(len=80) :: Txt
type(DSBA_Type) :: CMOt(1)
character(len=*), parameter :: SecNam = 'GetGrad_ER'

! Initialization.
! ---------------

Functional = Zero
GradNorm = Zero
if ((nOrb2Loc < 1) .or. (nBasis < 1)) return

! Transpose CMO (only the part to be localised).
! ----------------------------------------------

call Allocate_DT(CMOt(1),[nOrb2Loc],[nBasis],nSym)
do i=1,nOrb2Loc
  CMOt(1)%SB(1)%A2(i,:) = CMO(:,i)
end do

! Compute R.
! ----------

nOcc(1) = nOrb2Loc
irc = -1
call Cho_Get_Rij(irc,CMOt,nOcc,R,Timing)
if (irc /= 0) then
  write(Txt,'(A,I6)') 'Cho_Get_Rij returned',irc
  call SysAbendMsg(SecNam,'Calculation of ER gradient failed:',Txt)
end if
call Deallocate_DT(CMOt(1))

! Compute gradient norm and functional.
! -------------------------------------

do i=1,nOrb2Loc-1
  Functional = Functional+R(i,i)
  do j=i+1,nOrb2Loc
    GradNorm = GradNorm+(R(i,j)-R(j,i))**2
  end do
end do
Functional = Functional+R(nOrb2Loc,nOrb2Loc)
GradNorm = Four*sqrt(GradNorm)

end subroutine GetGrad_ER
