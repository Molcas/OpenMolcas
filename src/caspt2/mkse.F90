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

subroutine MKSE(DREF,NDREF)
! Set up the matrix SE(t,x)
! Formula used:
!    SE(t,x)=2*dtx - Dtx

use definitions, only: iwp, wp
use constants, only: Two
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSMAT
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NSYM, NINDEP, NASH, NAES

implicit none
integer(kind=iwp), intent(in) :: NDREF
real(kind=wp), intent(in) :: DREF(NDREF)
real(kind=wp), allocatable :: SE(:)
integer(kind=iwp) ISYM, NINP, NINM, NAS, NSE, IT, ITABS, IX, IXABS, ISE, ID, IDISK

do ISYM=1,NSYM
  NINP = NINDEP(ISYM,6)
  if (NINP == 0) cycle
  NINM = NINDEP(ISYM,7)
  NAS = NASH(ISYM)
  NSE = (NAS*(NAS+1))/2
  if (NSE > 0) call mma_allocate(SE,NSE,Label='SE')
  do IT=1,NAS
    ITABS = IT+NAES(ISYM)
    do IX=1,IT
      IXABS = IX+NAES(ISYM)
      ISE = (IT*(IT-1))/2+IX
      ID = (ITABS*(ITABS-1))/2+IXABS
      if (ITABS == IXABS) then
        SE(ISE) = Two-DREF(ID)
      else
        SE(ISE) = -DREF(ID)
      end if
    end do
  end do

  ! Write to disk
  if ((NSE > 0) .and. (NINDEP(ISYM,6) > 0)) then
    IDISK = IDSMAT(ISYM,6)
    call DDAFILE(LUSBT,1,SE,NSE,IDISK)
    if ((NINM > 0) .and. (NINDEP(ISYM,7) > 0)) then
      IDISK = IDSMAT(ISYM,7)
      call DDAFILE(LUSBT,1,SE,NSE,IDISK)
    end if
    call mma_deallocate(SE)
  end if
end do

end subroutine MKSE
