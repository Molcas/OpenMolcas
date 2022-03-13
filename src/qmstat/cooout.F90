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

subroutine Cooout(Head,Cordst,nPart,nCent)

implicit real*8(a-h,o-z)
#include "maxi.fh"
dimension Cordst(MxCen*MxPut,3)
character Head*200

write(6,*)
write(6,*)
write(6,'(A)') Head
kaunter = 0
do i=1,nPart
  write(6,*) 'Molecule ',i
  do j=1,nCent
    kaunter = kaunter+1
    write(6,*) (Cordst(kaunter,ii),ii=1,3)
  end do
end do

return

end subroutine Cooout
