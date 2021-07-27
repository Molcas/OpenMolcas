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

subroutine PrGrad(Label,Grad,nGrad,Names,iPrint)
!***********************************************************************
!                                                                      *
! Object: to print set gradient with respect to the symmetrical dis-   *
!         placements.                                                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use Symmetry_Info, only: lIrrep

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
real*8 Grad(nGrad)
real*8 CGrad(3,MxAtom)
character CNames(MxAtom)*(LENIN5)
character Label*(*), Names(nGrad)*(LENIN6)
character Namei*(LENIN5)

write(6,*)
call Banner(Label,1,len(Label)+30)
write(6,*)
!if (iPrint == 4) then
if (.true.) then
  call TrGrd_Alaska_(CGrad,CNames,Grad,nGrad,iCen)
  write(6,'(1x,A,A)') ' Irreducible representation: ',lIrrep(0)
  write(6,'(1x,90A     )') ('-',i=1,90)
  write(6,'(7x,3(23x,A))') 'X','Y','Z'
  write(6,'(1x,90A     )') ('-',i=1,90)
  do iGrad=1,iCen
    TempX = CGrad(1,iGrad)
    TempY = CGrad(2,iGrad)
    TempZ = CGrad(3,iGrad)
    Namei = CNames(iGrad)
    write(6,'(2X,A,3X,3ES24.14)') Namei,TempX,TempY,TempZ
  end do
  write(6,'(1x,90A     )') ('-',i=1,90)
else

  ! Modified by Luca De Vico november 2005 Teokem
  ! I need to print the full gradient vector
  ! to use it for constrained optimizations
  ! with TRANSVERSE option in SLAPAF MODULE

  !mGrad = min(21,nGrad)

  mGrad = nGrad
  write(6,'(15x,A,A)') ' Irreducible representation: ',lIrrep(0)
  write(6,*)
  do iGrad=1,mGrad
    Temp = Grad(iGrad)
    !if (abs(Temp) < 1.0D-15) Temp = Zero
    write(6,'(16X,A,15X,ES15.7)') Names(iGrad),Temp
  end do

  !if (nGrad > 21) then
  !  write(6,*)
  !  write(6,*) '   ... list is truncated ...'
  !  write(6,*)
  !end if

  ! End of modifications

end if
write(6,*)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(iPrint)

end subroutine PrGrad
