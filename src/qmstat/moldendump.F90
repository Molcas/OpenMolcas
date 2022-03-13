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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "WrkSpc.fh"
#include "constants.fh"
#include "real.fh"
parameter(AuAng=1d10*CONST_BOHR_RADIUS_IN_SI_)
dimension CooRef(MxCen,3), Coo(MxCen,3)
dimension iC(3)
logical ValidOrNot

! Clarifying words.

write(6,*)
write(6,*)
write(6,*) '   * Coordinates given in form for Molden *'
write(6,*)
write(6,*) ' Put everything within the lines in a separate file and view with Molden.'
write(6,*) ' Observe that the identity of molecules that are not valid water molecules is unknown.'
write(6,*)
write(6,*) '------------------------------------------------------------------------------------'

! Print total number of particles.

write(6,*) '  Substitue this line with number of atoms.'
write(6,*)
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
      write(6,92) 'C  ',(AuAng*Coo(jC,kk),kk=1,3)
    end do
  else
    write(6,92) 'O  ',(AuAng*Coo(1,kk),kk=1,3)
    write(6,92) 'H  ',(AuAng*Coo(2,kk),kk=1,3)
    write(6,92) 'H  ',(AuAng*Coo(3,kk),kk=1,3)
  end if
end do
write(6,*)
write(6,*) '------------------------------------------------------------------------------------'

! Formats

92 format(A,3(F10.6))

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nA)

end subroutine MoldenDump
