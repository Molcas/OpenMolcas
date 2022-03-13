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

subroutine BornMayerBK(iQ_Atoms,BoMaH,BoMaO)
! With the Brdarski-Karlstrom scheme, construct the Born-Mayer parameters.

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
dimension BoMaH(MxAt), BoMaO(MxAt)
dimension rBdi(MxCen), rBdiQ(MxAt)
parameter(cjhr=0.1734d0)

! The solvent part.

do i=1,2
  rdi2 = 0
  do j=1,3
    rdi2 = rdi2+quadi(j,i)
  end do
  rBdi(i) = sqrt(rdi2/charDi(i))
end do

! The solute part.

do i=1,iQ_Atoms
  rdi2 = 0
  do j=1,3
    rdi2 = rdi2+QuadiQ(j,i)
  end do
  rBdiQ(i) = sqrt(rdi2/charDiQ(i))
end do

! Put together.

do i=1,iQ_Atoms
  BoMaH(i) = 1/(cjhr*(RBdiQ(i)+rbdi(1)))
  BoMaO(i) = 1/(cjhr*(RBdiQ(i)+rbdi(2)))
  if (iPrint >= 8) then
    write(6,*) '   Born-Mayer parameters.'
    write(6,'(A,i2,A,2(f12.4))') '    Atom ',i,' (H/O):',BoMaH(i),BoMaO(i)
  end if
end do

return

end subroutine BornMayerBK
