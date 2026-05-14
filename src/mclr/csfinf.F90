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

subroutine CsfInf(lSym,iSpin,iSPC,nsym)

use Str_Info, only: CFTP, CNSM, DFTP, DTOC, NELEC, NOCTYP, STR
use MCLR_Data, only: i1, iAnders, IASTFI, IBSTFI, iDC, iRefSM, lConf, llDET, LuCSF2SD, MAXOP, MINOP, MNR1IC, MS2, MXR3IC, NACOB, &
                     NCNATS, NELCI, NORB1, NORB2, NORB3, PSSIGN
use CandS, only: ICSM, ICSPC, ISSM, ISSPC
use input_mclr, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lSym, iSpin, iSPC, nsym
integer(kind=iwp) :: IA, IATP, IBTP, idum(1), ISYM, LLCSF, MNELR1, MXELR3, NCOMB, NEL, NOCTPA, NOCTPB, NOOS
integer(kind=iwp), allocatable :: IOOS1(:), NOOS1(:), SBLTP(:), SIOIO(:)

! Sorry about this  but this is just to tell the program
! that no CSF<->SD coefficents is in core
i1 = -9
iAnders = -9
ICSM = lSym
ISSM = lSym
ICSPC = 1
ISSPC = 1
NEL = NELCI(ISPC)
IATP = IASTFI(ISPC)
IBTP = IBSTFI(ISPC)
NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
MNELR1 = MNR1IC(ISPC)
!MXELR3 = MNR1IC(ISPC)
MXELR3 = MXR3IC(ISPC)
iRefSm = lsym
! Obtain OOS pointer array
call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
call IAIBCM_MCLR(MNR1IC(ISSPC),MXR3IC(ISSPC),NOCTPA,NOCTPB,Str(IATP)%EL1,Str(IATP)%EL3,Str(IBTP)%EL1,Str(IBTP)%EL3,SIOIO)
call mma_allocate(SBLTP,nIrrep,Label='SBLTP')
NOOS = NOCTPA*NOCTPB*nIrrep
call mma_allocate(IOOS1,NOOS,Label='IOOS1')
call mma_allocate(NOOS1,NOOS,Label='NOOS1')
call INTCSF(NACOB,NEL,iSpin,MS2,NORB1,NORB2,NORB3,MNELR1,MXELR3,LLCSF,1,PSSIGN,lconf,lldet)

! Calculate CG COEFFICENTS ETC

call CSDTMT(DFTP,CFTP,DTOC,PSSIGN)

! Calculate the reordering vector and write it to disk

iA = 0
do iSym=1,nSym
  !OOS arrayy
  call ZBLTP(ISYM,nIrrep,IDC,SBLTP,idum)
  call ZOOS(ISYM,SBLTP,nIrrep,SIOIO,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,NOCTPA,NOCTPB,idc,IOOS1,NOOS1,NCOMB,0)
  !EAW call CNFORD(CNSM(1)%ICTS,CNSM(1)%ICONF,iSym,NACOB,DFTP,NCNATS(1,ISYM),NEL,0,IDUM,IDUM,IASTFI(ISPC),IBSTFI(ISPC),IOOS1, &
  !                NORB1,NORB2,NORB3,MNELR1,MXELR3,NELEC(IASTFI(ISPC)),NELEC(IBSTFI(ISPC)),MINOP,MAXOP)
  call CNFORD(CNSM(1)%ICTS,CNSM(1)%ICONF,iSym,NACOB,DFTP,NCNATS(1,ISYM),NEL,0,IDUM,IDUM,IASTFI(ISPC),IBSTFI(ISPC),IOOS1,NORB1, &
              NORB2,NORB3,MNELR1,MXELR3,NELEC(IASTFI(ISPC)),NELEC(IBSTFI(ISPC)),MINOP,MAXOP,PSSIGN)

  call iDAFILE(LUCSF2SD,1,CNSM(1)%ICTS,lldet,iA)
  call iDAFILE(LUCSF2SD,1,CNSM(1)%ICONF,lconf,iA)
end do

call mma_deallocate(CNSM(1)%ICTS)
call mma_deallocate(CNSM(1)%ICONF)
call mma_deallocate(IOOS1)
call mma_deallocate(NOOS1)
call mma_deallocate(SBLTP)
call mma_deallocate(SIOIO)

end subroutine CsfInf
