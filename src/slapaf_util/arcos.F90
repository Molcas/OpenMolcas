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

function arCos(Arg)

implicit real*8(a-h,o-z)
real*8 ArCos
character*72 Warning
#include "real.fh"

A = Arg
Delta = 1.0D-12
if (abs(A) > One) then
  write(Warning,3) A
3 format(1X,'Warning argument of aCos= ',1F21.18)
  if (abs(A) < One+Delta) then
    !call WarningMessage(1,Warning)
    A = sign(One,A)
  else
    call WarningMessage(2,Warning)
    call Abend()
  end if
end if

ArCos = acos(A)

return

end function arCos
