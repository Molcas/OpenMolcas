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
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! Written by Steven Vancoillie, May 2011
! A set of subroutines that can handle RHS arrays in either a serial or
! parallel environment, depending on the situation.
!***********************************************************************
! --> when running serially, the RHS arrays are stored on LUSOLV and are
! loaded into the WORK array when needed.
! --> when running in parallel, the RHS arrays are stored on disk as
! disk resident arrays (DRAs) with filename RHS_XX_XX_XX, where XX is a
! number referring to the case, symmetry, and RHS vector respectively,
! and are loaded onto a global array when needed.
!***********************************************************************

subroutine RHS_FPRINT(CTYPE,IVEC)

use caspt2_module, only: NASUP, NINDEP, NISUP, NSYM
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: CTYPE
integer(kind=iwp), intent(in) :: IVEC
integer(kind=iwp) :: ICASE, ISYM, lg_W, NAS, NIN, NIS, NROW
real(kind=wp) :: FP(8)
real(kind=wp), external :: RHS_DDOT

!-SVC: print out DNRM2 of the all RHS components
NROW = 0 ! dummy initialize
do ICASE=1,13
  do ISYM=1,NSYM

    NAS = NASUP(ISYM,ICASE)
    NIN = NINDEP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)

    if (CTYPE == 'C') then
      NROW = NAS
    else if (CTYPE == 'SR') then
      NROW = NIN
    else
      write(u6,'(1X,A)') 'RHS_FPRINT: invalid type: '//CTYPE
      call ABEND()
    end if

    if ((NAS /= 0) .and. (NIN /= 0) .and. (NIS /= 0)) then
      call RHS_ALLO(NROW,NIS,lg_W)
      call RHS_READ(NROW,NIS,lg_W,iCASE,iSYM,iVEC)
      FP(ISYM) = sqrt(RHS_DDOT(NROW,NIS,lg_W,lg_W))
      call RHS_FREE(lg_W)
    else
      FP(ISYM) = Zero
    end if
  end do
  write(u6,'(1X,I2,1X,8F21.14)') ICASE,FP(1:NSYM)
end do

end subroutine RHS_FPRINT
