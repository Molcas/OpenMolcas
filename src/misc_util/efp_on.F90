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

logical function EFP_On()

#ifdef _EFP_
use EFP_Module

implicit real*8 (a-h,o-z)

call Get_lScalar('EFP',EFP)

EFP_On = lEFP
#else
EFP_On = .false.
#endif

return

end function EFP_On
