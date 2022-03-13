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

subroutine MoldenDump(iC,CooRef,nP,nA,nC)

use Constants, only: Angstrom
use Definitions, only: wp, iwp, u6

implicit none
#include "maxi.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iC(3), nP, nA, nC
real(kind=wp) :: CooRef(MxCen,3)
integer(kind=iwp) :: ind, iP, jC, kk
real(kind=wp) :: Coo(MxCen,3) !IFG
logical(kind=iwp) :: ValidOrNot

! Clarifying words.

write(u6,*)
write(u6,*)
write(u6,*) '   * Coordinates given in form for Molden *'
write(u6,*)
write(u6,*) ' Put everything within the lines in a separate file and view with Molden.'
write(u6,*) ' Observe that the identity of molecules that are not valid water molecules is unknown.'
write(u6,*)
write(u6,*) '------------------------------------------------------------------------------------'

! Print total number of particles.

write(u6,*) '  Substitue this line with number of atoms.'
write(u6,*)
do iP=1,nP
  ind = nC*(iP-1)
  do jC=1,nC
    Coo(jC,1) = Work(iC(1)+ind+jC-1)
    Coo(jC,2) = Work(iC(2)+ind+jC-1)
    Coo(jC,3) = Work(iC(3)+ind+jC-1)
  end do
  call IsItValid(Coo,CooRef,ValidOrNot)
  if (.not. ValidOrNot) then
    do jC=1,nC
      write(u6,92) 'C  ',(Angstrom*Coo(jC,kk),kk=1,3)
    end do
  else
    write(u6,92) 'O  ',(Angstrom*Coo(1,kk),kk=1,3)
    write(u6,92) 'H  ',(Angstrom*Coo(2,kk),kk=1,3)
    write(u6,92) 'H  ',(Angstrom*Coo(3,kk),kk=1,3)
  end if
end do
write(u6,*)
write(u6,*) '------------------------------------------------------------------------------------'

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nA)

92 format(A,3(F10.6))

end subroutine MoldenDump
