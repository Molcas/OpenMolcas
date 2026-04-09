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
subroutine INTCSF(NACTOB,NACTEL,MULTP,MS2,NORB1,NORB2,NORB3,NEL1MN,NEL3MX,LCSF,NCNSM,PSSIGN,lconf,ldet)
! Initializing routine for CSF-DET expansions of internal space
!
! Set up common block /CSFDIM/
! This gives information about the number of dets,csf's for
! each symmetry
!
! find local memory requirements for CSF routines
! Largest local memory requirements in CNFORD,CSFDET_MCLR is returned in
! LCSF
!
!   DFTP : OPEN SHELL DETERMINANTS OF PROTO TYPE
!   CFTP : BRANCHING DIAGRAMS FOR PROTO TYPES
!   DTOC  : CSF-DET TRANSFORMATION FOR PROTO TYPES
!   CNSM(I)%ICONF : SPACE FOR STORING NCNSM CONFIGURATION EXPANSIONS
!
! If PSSIGN /= 0, spin combinations are used !!
!
! Last modification : Sept 20 : sign and address of dets goes together
!                      in CNSM(:)%ICTS

use Str_Info, only: CFTP, CNSM, DFTP, DTOC
use MCLR_Data, only: MAXOP, MINOP, MS2P, MULTSP, MXPCSM, NCNASM, NCNATS, NCPCNT, NCSASM, NDPCNT, NDTASM, NTYP
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NACTOB, NACTEL, MULTP, MS2, NORB1, NORB2, NORB3, NEL1MN, NEL3MX, NCNSM
integer(kind=iwp), intent(out) :: LCSF, lconf, ldet
real(kind=wp), intent(in) :: PSSIGN
integer(kind=iwp) :: IAEL, IBEL, ICL, ICNSM, IEL1, IEL2, IEL3, ILCNF, ILLCNF, IMSCMB, IOP, IOP1, IOP2, IOP3, IOPEN, ISYM, ITP, &
                     ITYP, IWEYLF, LCNFOR, LCSFDT, LDTOC, LICS, LIDT, LLCONF, MULTS, MXCNSM, MXDT, MXPCTP, MXPTBL, NEL
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ITYPE
#endif
integer(kind=iwp), external :: IBINOM

if (PSSIGN /= Zero) then
  IMSCMB = 1
else
  IMSCMB = 0
end if

! Define parameters in SPINFO

MULTSP = MULTP
MS2P = MS2
MULTS = MULTSP

NEL = NACTEL
! ===============================
! Allowed number of open orbitals
! ===============================
MINOP = abs(MS2)
MAXOP = 0
do IEL1=NEL1MN,2*NORB1
  do IEL3=0,NEL3MX
    IEL2 = NACTEL-IEL1-IEL3
    if (IEL2 < 0) cycle
    IOP1 = min(NORB1,2*NORB1-IEL1,IEL1)
    IOP2 = min(NORB2,2*NORB2-IEL2,IEL2)
    IOP3 = min(NORB3,2*NORB3-IEL3,IEL3)
    IOP = IOP1+IOP2+IOP3
    MAXOP = max(MAXOP,IOP)
  end do
end do
!write(u6,*) ' MAXOP with RAS constraints :',MAXOP
NTYP = MAXOP-MINOP+1

MXPCTP = size(NCPCNT)
if (NTYP > MXPCTP) then
  write(u6,*) '  NUMBER OF CONFIGURATION TYPES TO LARGE'
  write(u6,*) '  CHANGE PARAMETER MXPCTP TO AT LEAST ',NTYP
  write(u6,*) '  CURRENT VALUE OF MXPCTP ',MXPCTP
  write(u6,*) ' MTYP IN LUSPIN TO SMALL'
  call Abend()
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' MINOP MAXOP NTYP ',MINOP,MAXOP,NTYP
#endif
! ===============================================
! Number of sd's and csf's per configuration type
! ===============================================
do ITP=1,NTYP
  IOPEN = MINOP+ITP-1
  IAEL = (IOPEN+MS2)/2
  IBEL = (IOPEN-MS2)/2
  if (IAEL+IBEL == IOPEN) then
    NDPCNT(ITP) = IBINOM(IOPEN,IAEL)
    if ((IMSCMB /= 0) .and. (IOPEN /= 0)) NDPCNT(ITP) = NDPCNT(ITP)/2
    if (IOPEN >= MULTS-1) then
      NCPCNT(ITP) = IWEYLF(IOPEN,MULTS)
    else
      NCPCNT(ITP) = 0
    end if
  else
    NDPCNT(ITP) = 0
    NCPCNT(ITP) = 0
  end if
end do
#ifdef _DEBUGPRINT_
write(u6,'(/A)') ' Information about prototype configurations'
write(u6,'(A)') ' =========================================='
write(u6,'(/A)')
if (IMSCMB == 0) then
  write(u6,'(/A)') ' Combinations = Slater determinants'
else
  write(u6,'(/A)') ' Combinations = Spin combinations'
end if
write(u6,'(/A)') '  Open orbitals   Combinations    CSFs'
do IOPEN=MINOP,MAXOP,2
  ITYPE = IOPEN-MINOP+1
  write(u6,'(5X,I3,10X,I6,7X,I6)') IOPEN,NDPCNT(ITYPE),NCPCNT(ITYPE)
end do
#endif
! ==============================================
! Number of Combinations and CSF's per  symmetry
! ==============================================
call CISIZE(NORB1,NORB2,NORB3,NEL1MN,NEL3MX,NACTEL,MINOP,MAXOP,MXPCTP,MXPCSM,NCNATS,NCNASM,NDTASM,NCSASM,NDPCNT,NCPCNT)

! ===========================================
! Permanent and local memory for csf routines
! ===========================================
!
! memory for CSDTMT arrays.
! Largest block of proto type determinants.
! Largest number of prototype determinants
! All configurations (of any specific symmetry)

LIDT = 0
LICS = 0
LDTOC = sum(NCPCNT(1:NTYP)*NDPCNT(1:NTYP))
MXPTBL = 0
MXDT = max(0,maxval(NDPCNT(1:NTYP)))
do ITP=1,NTYP
  IOPEN = MINOP+ITP-1
  LIDT = LIDT+NDPCNT(ITP)*IOPEN
  LICS = LICS+NCPCNT(ITP)*IOPEN
  MXPTBL = max(NDPCNT(ITP)*IOPEN,MXPTBL)
end do
! local memory for CSFDET_MCLR
LCSFDT = MXPTBL+MAXOP
! local memory for CNFORD
LCNFOR = max(2*NTYP+NACTOB,(MXDT+2)*NACTEL)
! local memory for any routine used in construction of csf basis
LCSF = max(LCSFDT,LCNFOR)
! Memory needed to store ICONF array
LCONF = 0
LDET = 0
ILCNF = 0
do ISYM=1,MXPCSM
  ILLCNF = 0
  LLCONF = 0
  LDET = max(LDET,NDTASM(ISYM))
  do ITYP=1,NTYP
    IOPEN = ITYP+MINOP-1
    ICL = (NEL-IOPEN)/2
    LLCONF = LLCONF+NCNATS(ITYP,ISYM)*(IOPEN+ICL)
    ILLCNF = ILLCNF+NCNATS(ITYP,ISYM)
  end do
  !write(u6,*) ' MEMORY FOR HOLDING CONFS OF SYM... ',ISYM,LLCONF
  LCONF = max(LCONF,LLCONF)
  ILCNF = max(ILCNF,ILLCNF)
end do

! notice the ILCNF number ! yma

#ifdef _DEBUGPRINT_
write(u6,'(/A,I8)') '  Memory for holding largest list of configurations ',LCONF
write(u6,'(/A,I8)') '  Size of largest CI expansion (combinations) ',LDET
write(u6,'(/A,I8)') '  Size of largest CI expansion (confs) ',ILCNF
call xflush(u6) !yma
#endif

! permanent memory for csf proto type arrays

call mma_allocate(DFTP,LIDT,Label='DFTP')
call mma_allocate(CFTP,LICS,Label='CFTP')
call mma_allocate(DTOC,LDTOC,Label='DTOC')

! Permanent arrays for reordering and phases
MXCNSM = size(CNSM)
if (NCNSM > MXCNSM) then
  write(u6,'(A,2I2)') '  TROUBLE IN CSFDIM NCNSM > MXCNSM : NCNSM,MXCNSM',NCNSM,MXCNSM
  write(u6,*) ' CSFDIM : NCNSM  IS GREATER THAN MXCNSM'
  call Abend()
end if
do ICNSM=1,NCNSM
  call mma_allocate(CNSM(ICNSM)%ICONF,LCONF,Label='ICONF')
  call mma_allocate(CNSM(ICNSM)%ICTS,LDET,Label='ICTS')
end do

end subroutine INTCSF
