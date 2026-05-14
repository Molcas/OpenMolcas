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
! Copyright (C) 1984,1989-1993, Jeppe Olsen                            *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ICISPS()
! Number of dets and combinations
! per symmetry for each type of internal space
!
! ====================
! Output  XISPSM is calculated
! ====================
!
! Jeppe Olsen, Winter 1991
! Last revision April 1991

use Str_Info, only: NOCTYP, STR
use MCLR_Data, only: IACTI, IASTFI, IBSTFI, IDC, MNR1IC, MNR3IC, MXR1IC, MXR3IC, MXSB, MXSOOB, NAELCI, NBELCI, NICISP, XISPSM
use input_mclr, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IATP, IBTP, ICI, IIDC, ISYM, MXCEXP, MXS, MXSOO, NCOMB
real(kind=wp) :: XNCOMB
integer(kind=iwp), allocatable :: LBLTP(:), LCVST(:)

! Local memory
call mma_allocate(LBLTP,nIrrep,Label='LBLTP')
!if ((IDC == 3) .or. (IDC == 4)) call mma_allocate(LCVST,nIrrep,Label='LCVST')
call mma_allocate(LCVST,nIrrep,Label='LCVST')

! Obtain array giving symmetry of sigma v reflection times string symmetry.
!if ((IDC == 3) .or. (IDC == 4)) call SIGVST(LCVST,nIrrep)

MXSB = 0
MXSOOB = 0
MXCEXP = 0
XISPSM(:,:) = Zero
do ICI=1,NICISP
  do ISYM=1,nIrrep
    IATP = IASTFI(ICI)
    IBTP = IBSTFI(ICI)
    !MS write(u6,*) ' NRASDT : ICI IATP IBTP ',ICI,IATP,IBTP
    if (NAELCI(ICI) == NBELCI(ICI)) then
      IIDC = IDC
    else
      IIDC = 1
    end if
    if (IACTI(ICI) == 1) then
      call ZBLTP(ISYM,nIrrep,IIDC,LBLTP,LCVST)
      call NRASDT(MNR1IC(ICI),MXR1IC(ICI),MNR3IC(ICI),MXR3IC(ICI),ISYM,nIrrep,NOCTYP(IATP),NOCTYP(IBTP),Str(IATP)%EL1, &
                  Str(IBTP)%EL1,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,Str(IATP)%EL3,Str(IBTP)%EL3,NCOMB,XNCOMB,MXS,MXSOO,LBLTP)
      XISPSM(ISYM,ICI) = XNCOMB
      MXSOOB = max(MXSOOB,MXSOO)
      MXSB = max(MXSB,MXS)
      MXCEXP = max(MXCEXP,NCOMB)
    !else
    !  XISPSM(ISYM,ICI) = Zero
    end if
  end do
end do
call mma_deallocate(LCVST)
call mma_deallocate(LBLTP)

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of internal combinations per symmetry'
write(u6,*) ' ==========================================='

do ICI=1,NICISP
  if (IACTI(ICI) == 1) then
    write(u6,*) ' Internal CI space ',ICI
    call WRTMAT(XISPSM(1,ICI),1,nIrrep,1,nIrrep)
  end if
end do
write(u6,*) ' Largest CI space                 ',MXCEXP
write(u6,*) ' Largest symmetry block           ',MXSB
write(u6,*) ' Largest Symmetry-type-type block ',MXSOOB
#endif

end subroutine ICISPS
