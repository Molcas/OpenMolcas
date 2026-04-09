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

subroutine PROPER(PROP,ISTATE,JSTATE,TDMZZ,WDMZZ)

use rassi_global_arrays, only: JBNUM
use RASSI_AUX, only: TocM
use Cntrl, only: FnTOM, IRREP, lSym1, lSym2, LuTOM, NPROP, NSTATE, PNAME, PTYPE, ToFile
use Symmetry_Info, only: Mul, nIrrep
use rassi_data, only: NBASF, NBST, NTDMZZ
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: PROP(NSTATE,NSTATE,NPROP), TDMZZ(NTDMZZ), WDMZZ(NTDMZZ)
integer(kind=iwp) :: ISTATE, JSTATE
integer(kind=iwp) :: I, ICALL = 0, iDisk, iDIskSav = 0, IndCall, IOFF(8), iProp, iSy12, iType, J, JOB1, JOB2, Mask, NIP, nScr
character(len=8) :: LABEL
real(kind=wp), allocatable :: IP(:), SCR(:,:)

! COMBINED SYMMETRY OF STATES:
JOB1 = JBNUM(ISTATE)
JOB2 = JBNUM(JSTATE)
LSYM1 = IRREP(JOB1)
LSYM2 = IRREP(JOB2)
ISY12 = MUL(LSYM1,LSYM2)
! THE SYMMETRY CHECK MASK:
MASK = 2**(ISY12-1)
! ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
NIP = 4+(NBST*(NBST+1))/2
call mma_allocate(IP,NIP,Label='IP')
! FIRST SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS OF TDMSCR
call mk_IOFF(IOFF,nIrrep,NBASF,ISY12)
! CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
! AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES
NSCR = (NBST*(NBST+1))/2
call mma_allocate(SCR,nSCR,4,LABEL='SCR')
SCR(:,:) = Zero
call MK_TWDM(nIrrep,TDMZZ,WDMZZ,nTDMZZ,SCR,nSCR,iOFF,NBASF,ISY12)

! AT THIS POINT, THE SYMMETRICALLY AND ANTISYMMETRICALLY FOLDED
! DENSITY MATRICES, AND WE-REDUCED SPIN DENSITY MATRICES, HAVE BEEN
! CALCULATED BEGINNING IN SCR.
! LOOP OVER ALL REQUIRED ONE-ELECTRON OPERATORS:

!-------------------------------------------
!If requested by user, put SCR in an unformatted file for later
!use by another program. (A.Ohrn)
if (ToFile) then
  call DaName(LuToM,FnToM)
  if (iCall == 0) then  !Make room for table-of-contents
    iDisk = 0
    TocM(1:nState*(nState+1)/2) = -1
    call iDaFile(LuToM,1,TocM,nState*(nstate+1)/2,iDisk)
    TocM(1) = iDisk
    iDiskSav = iDisk
    iCall = 1
  end if
  if (iState < jState) then

    ! For the rest of the code to work this cannot be violated.

    write(u6,*) 'Proper: iState < jState'
    call Abend()
  end if
  i = iState
  j = jState
  indCall = i*(i-1)/2+j  !Which call this is
  ToCM(indCall) = iDiskSav
  iDisk = iDiskSav
  !write(u6,*) 'IndCall,iDisk=',IndCall,iDisk
  call dDaFile(LuToM,1,SCR,4*nSCR,iDisk) !The THING.
  iDiskSav = iDisk  !Save diskaddress.
  iDisk = 0
  call iDaFile(LuToM,1,TocM,nState*(nState+1)/2,iDisk)
  !Put table of contents.
  call DaClos(LuToM)
end if
!End of ToFile
!write(u6,*) 'ISTATE,JSTATE=',ISTATE,JSTATE
do IPROP=1,NPROP
  PROP(ISTATE,JSTATE,IPROP) = Zero
  LABEL = PNAME(IPROP)

  call UPCASE(LABEL)

  ! If the user wants the ASD term, it is the same as
  ! the EF2 term without the nuclear contribution
  ! the new magnetic integrals are calculated from X2C
  ! the old spin-dependent part is thus denoted as ASDO
  ! O for old
  if (LABEL(1:4) == 'ASDO') then
    LABEL(1:4) = 'EF2 '
    !write(u6,*) 'EF2---->ASD Here'
  end if
  if (LABEL(1:4) == 'TMOM') cycle

  if (LABEL(1:4) == 'PSOP') then
    LABEL(1:4) = 'PSOI'
    !write(u6,*) 'PSOI---->PSOP Here'
  end if

  if (LABEL(1:6) == 'DMP   ') LABEL(1:6) = 'DMS  1'

  ITYPE = 0
  if (PTYPE(IPROP) == 'HERMSING') ITYPE = 1
  if (PTYPE(IPROP) == 'ANTISING') ITYPE = 2
  if (PTYPE(IPROP) == 'HERMTRIP') ITYPE = 3
  if (PTYPE(IPROP) == 'ANTITRIP') ITYPE = 4
  if (ITYPE == 0) then
    write(u6,*) 'RASSI/PROPER internal error.'
    write(u6,*) 'Erroneous property type.'
    write(u6,*) 'PTYPE(IPROP)=',PTYPE(IPROP)
    call ABEND()
  end if

  call MK_PROP(PROP,IPROP,ISTATE,JSTATE,LABEL,ITYPE,IP,NIP,SCR,NSCR,MASK,ISY12,IOFF)

end do
call mma_deallocate(SCR)
call mma_deallocate(IP)

end subroutine PROPER
