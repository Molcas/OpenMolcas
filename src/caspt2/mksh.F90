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

subroutine MKSH()
! For completeness, even case H has formally S and B
! matrices. This costs nothing, and saves conditional
! looping, etc in the rest  of the routines.

use caspt2_module, only: NINDEP, NSYM
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSMAT
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ICASE, IDISK, ISYM, NIN
real(kind=wp) :: Dum(1)

DUM(1) = One
do ISYM=1,NSYM
  do ICASE=12,13
    NIN = NINDEP(ISYM,ICASE)
    if (NIN > 0) then
      IDISK = IDSMAT(ISYM,ICASE)
      call DDAFILE(LUSBT,1,DUM,1,IDISK)
    end if
  end do
end do

end subroutine MKSH
