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

use Data_structures, only: DSBA_Type
use Data_structures, only: Allocate_DSBA, Deallocate_DSBA
implicit real*8(a-h,o-z)
real*8 R(nOrb2Loc,nOrb2Loc), CMO(nBasis,nOrb2Loc)
logical Timing
character(LEN=10), parameter :: SecNam = 'GetGrad_ER'
character*80 Txt
integer, parameter :: nSym = 1
integer nOcc(nSym)
type(DSBA_Type) CMOt

! Initialization.
! ---------------

Functional = 0.0d0
GradNorm = 0.0d0
if ((nOrb2Loc < 1) .or. (nBasis < 1)) return

! Transpose CMO (only the part to be localised).
! ----------------------------------------------

call Allocate_DSBA(CMOt,[nOrb2Loc],[nBasis],nSym)
do i=1,nOrb2Loc
  CMOt%SB(1)%A2(i,:) = CMO(:,i)
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
call Deallocate_DSBA(CMOt)

! Compute gradient norm and functional.
! -------------------------------------

do i=1,nOrb2Loc-1
  Functional = Functional+R(i,i)
  do j=i+1,nOrb2Loc
    GradNorm = GradNorm+(R(i,j)-R(j,i))**2
  end do
end do
Functional = Functional+R(nOrb2Loc,nOrb2Loc)
GradNorm = 4.0d0*sqrt(GradNorm)

end subroutine GetGrad_ER
