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

subroutine PrGrad(Label,Grad,nGrad,Names)
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
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(in) :: Grad(nGrad)
character(len=LenIn6), intent(in) :: Names(nGrad)
integer(kind=iwp) :: iCen, iGrad, mGrad
real(kind=wp) :: CGrad(3,MxAtom), Temp, TempX, TempY, TempZ
character(len=LenIn5) :: CNames(MxAtom), Namei

write(u6,*)
call Banner(Label,1,len(Label)+30)
write(u6,*)
if (.true.) then
  call TrGrd_Alaska(CGrad,CNames,Grad,nGrad,iCen)
  write(u6,'(1x,A,A)') ' Irreducible representation: ',lIrrep(0)
  write(u6,'(1x,A)') repeat('-',90)
  write(u6,'(7x,3(23x,A))') 'X','Y','Z'
  write(u6,'(1x,A)') repeat('-',90)
  do iGrad=1,iCen
    TempX = CGrad(1,iGrad)
    TempY = CGrad(2,iGrad)
    TempZ = CGrad(3,iGrad)
    Namei = CNames(iGrad)
    write(u6,'(2X,A,3X,3ES24.14)') Namei,TempX,TempY,TempZ
  end do
  write(u6,'(1x,A)') repeat('-',90)
else

  ! Modified by Luca De Vico november 2005 Teokem
  ! I need to print the full gradient vector
  ! to use it for constrained optimizations
  ! with TRANSVERSE option in SLAPAF MODULE

  !mGrad = min(21,nGrad)

  mGrad = nGrad
  write(u6,'(15x,A,A)') ' Irreducible representation: ',lIrrep(0)
  write(u6,*)
  do iGrad=1,mGrad
    Temp = Grad(iGrad)
    !if (abs(Temp) < 1.0e-15_wp) Temp = Zero
    write(u6,'(16X,A,15X,ES15.7)') Names(iGrad),Temp
  end do

  !if (nGrad > 21) then
  !  write(u6,*)
  !  write(u6,*) '   ... list is truncated ...'
  !  write(u6,*)
  !end if

  ! End of modifications

end if
write(u6,*)

return

end subroutine PrGrad
