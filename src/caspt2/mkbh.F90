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
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MKBH()
! For completeness, even case H has formally S and B
! matrices. This costs nothing, and saves conditional
! looping, etc in the rest of the routines.

use constants, only: Zero
use EQSOLV, only: IDBMAT
use caspt2_global, only: LUSBT
use caspt2_module, only: NSYM, NINDEP
use definitions, only: iwp, wp

implicit none
real(kind=wp) DUM(1)
integer(kind=iwp) ISYM, ICASE, NIN, IDISK

DUM(1) = Zero
do ISYM=1,NSYM
  do ICASE=12,13
    NIN = NINDEP(ISYM,ICASE)
    if (NIN > 0) then
      IDISK = IDBMAT(ISYM,ICASE)
      call DDAFILE(LUSBT,1,DUM,1,IDISK)
    end if
  end do
end do

end subroutine MKBH
