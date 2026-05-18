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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************

pure function IPROW(IROW,NQOT,NREM)

  use definitions, only: iwp

  integer(kind=iwp) :: IPROW
  integer(kind=iwp), intent(in) :: IROW, NQOT, NREM
  integer(kind=iwp) :: TMP

  TMP = IROW-NREM*(NQOT+1)
  if (TMP > 0) then
    IPROW = (TMP-1)/NQOT+NREM+1
  else
    IPROW = (IROW-1)/(NQOT+1)+1
  end if

end function IPROW
