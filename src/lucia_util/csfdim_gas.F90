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
! Copyright (C) 1994, Jeppe Olsen                                      *
!               2024, Giovanni Li Manni                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine CSFDIM_GAS(IOCCLS,NOCCLS,ISYM)
! Initializing routine for CSF-DET expansions
!
! information about the number of dets,csf's for
! each symmetry. CI space is defined by the NOCCLS
! occupation classes IOCCLS
!
! DETERMINE BASE ADDRESSES
!             DFTP : OPEN SHELL DETERMINANTS OF PROTO TYPE
!             CFTP : BRANCHING DIAGRAMS FOR PROTO TYPES
!             DTOC  : CSF-DET TRANSFORMATION FOR PROTO TYPES
!             ONF_OCC%I(I) : SPACE FOR STORING  NCNSM
!                        CONFIGURATION EXPANSIONS
! (Spin signaled by PSSIGN in CIINFO)
!
! Adapted for GAS calculations and LUCIA, Dec. 2001
! G. Li Manni, June 2024: Scale-up capability for single SD ROHF type calculations

use Data_Structures, only: Allocate_DT
use lucia_data, only: CFTP, CONF_OCC, CONF_REO, DFTP, DTOC, I2ELIMINATED_IN_GAS, I_ELIMINATE_GAS, IB_CONF_OCC, IB_CONF_REO, &
                      IB_SD_FOR_OPEN, IBCONF_ALL_SYM_FOR_OCCLS, IELIMINATED_IN_GAS, IPRCIX, MAXOP, MINOP, MS2, MULTS, MXPORB, &
                      N_2ELIMINATED_GAS, N_ELIMINATED_GAS, NCONF_ALL_SYM, NCONF_PER_OPEN, NCONF_PER_SYM, NCONF_TOT, NCSF_HEXS, &
                      NCSF_PER_SYM, NGAS, NOBPT, NOCOB, NPCMCNF, NPCSCNF, NPDTCNF, NSD_PER_SYM, PSSIGN, REO_PTDT, SDREO, SDREO_I, &
                      Z_PTDT
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NOCCLS, IOCCLS(NGAS,NOCCLS), ISYM
#include "warnings.h"
integer(kind=iwp) :: HEXS_CNF(MXPORB+1), I, IAEL, IALPHA, IB, IBEL, ICL, IDOREO, IDUM, IDUM_ARR(1), IELIM, IGAS, ILCNF, ILLCNF, &
                     INITIALIZE_CONF_COUNTERS, IOPEN, ITP, ITYP, J, JGAS, JOCCLS, LCONF, LDTOC, LENGTH_LIST, LICS, LIDT, LLCONF, &
                     LPTDT, LZ, maxingas(N_ELIMINATED_GAS), maxingas2(N_2ELIMINATED_GAS), MXDT, MXPTBL, NCMB, NCONF_OCCLS, NCSF, &
                     NELEC, NSD, TMP_CNF(MXPORB+1)
integer(kind=iwp), external :: IBINOM, IWEYLF

IDUM = 0
IDUM_ARR(1) = 0

#ifdef _DEBUGPRINT_
write(u6,*) '  PSSIGN : ',PSSIGN
write(u6,*) ' MULTS, MS2 = ',MULTS,MS2
#endif
NELEC = sum(IOCCLS(:,1))

! Define parameters in SPINFO

! Allowed number of open orbitals
MINOP = abs(MS2)
call MAX_OPEN_ORB(MAXOP,IOCCLS,NGAS,NOCCLS,NOBPT)
if (IPRCIX >= 6) write(u6,*) ' MINOP MAXOP ',MINOP,MAXOP

! Number of prototype sd's and csf's per configuration prototype

do IOPEN=MINOP,MAXOP
  ITP = IOPEN+1
  ! Unpaired electrons :
  IAEL = (IOPEN+MS2)/2
  IBEL = (IOPEN-MS2)/2
  if ((IAEL+IBEL == IOPEN) .and. (IAEL-IBEL == MS2) .and. (IAEL >= 0) .and. (IBEL >= 0)) then
    if ((PSSIGN == Zero) .or. (IOPEN == 0)) then
      ! Number of determinants is in general set to number of combinations
      NPDTCNF(ITP) = IBINOM(IOPEN,IAEL)
    else
      NPDTCNF(ITP) = IBINOM(IOPEN,IAEL)/2
    end if
    if (IOPEN >= MULTS-1) then
      NPCSCNF(ITP) = IWEYLF(IOPEN,MULTS)
    else
      NPCSCNF(ITP) = 0
    end if
  else
    NPDTCNF(ITP) = 0
    NPCSCNF(ITP) = 0
  end if
end do
NPCMCNF(MINOP+1:MAXOP+1) = NPDTCNF(MINOP+1:MAXOP+1)

if (IPRCIX >= 5) then
  if (PSSIGN == Zero) then
    write(u6,*) '  (Combinations = Determinants)'
  else
    write(u6,*) '  (Spin combinations in use)'
  end if
  write(u6,'(/A)') ' Information about prototype configurations'
  write(u6,'( A)') ' =========================================='
  write(u6,'(/A)') '  Open orbitals   Combinations    CSFs'
  do IOPEN=MINOP,MAXOP,2
    write(u6,'(5X,I3,10X,I6,7X,I6)') IOPEN,NPCMCNF(IOPEN+1),NPCSCNF(IOPEN+1)
  end do
end if

! Number of Configurations per occupation type

!if (NOCCLS > MXPCSM) then
!  write(u6,*) ' A known bug has reoccurred -- It seems that'
!  write(u6,*) ' the named constant MXPCSM must be increased'
!  write(u6,*) ' from its current value MXPCSM=',MXPCSM
!  write(u6,*) ' to AT LEAST NOCCLS=',NOCCLS
!  write(u6,*) ' This parameter is found in the module'
!  write(u6,*) '  <molcas>/src/lucia_util/lucia_data.F90'
!  write(u6,*) ' Change it. Then ''cd'' to molcas root'
!  write(u6,*) ' directory and give command ''make''.'
!  write(u6,*) ' But this may also be a bug. Please tell the'
!  write(u6,*) ' molcas developers!'
!  call Quit(_RC_INTERNAL_ERROR_)
!end if
!MGD : max occupation in removed GAS spaces
do i=1,N_ELIMINATED_GAS
  iGAS = IELIMINATED_IN_GAS(i)
  maxingas(i) = max(0,maxval(IOCCLS(iGAS,:)))
end do
do i=1,N_2ELIMINATED_GAS
  iGAS = I2ELIMINATED_IN_GAS(i)
  maxingas2(i) = max(0,maxval(IOCCLS(iGAS,:)))
end do
HEXS_CNF(:) = 0
NCONF_PER_OPEN(:,ISYM) = 0

call mma_allocate(IBCONF_ALL_SYM_FOR_OCCLS,NOCCLS,Label='IBCONF_ALL_SYM_FOR_OCCLS')
NCONF_ALL_SYM = 0
do JOCCLS=1,NOCCLS
  if (JOCCLS == 1) then
    INITIALIZE_CONF_COUNTERS = 1
  else
    INITIALIZE_CONF_COUNTERS = 0
  end if

  IBCONF_ALL_SYM_FOR_OCCLS(JOCCLS) = NCONF_ALL_SYM+1

  IDOREO = 0
  TMP_CNF(:) = 0
  call GEN_CONF_FOR_OCCLS(IOCCLS(:,JOCCLS),IDUM,INITIALIZE_CONF_COUNTERS,NGAS,ISYM,MINOP,MAXOP,1,NOCOB,NOBPT,TMP_CNF,NCONF_OCCLS, &
                          IB_CONF_REO,IB_CONF_OCC,IDUM_ARR,IDOREO,IDUM_ARR,NCONF_ALL_SYM,idum_arr,nconf_tot)
  NCONF_PER_OPEN(:,ISYM) = NCONF_PER_OPEN(:,ISYM)+TMP_CNF(:)
  !MGD add to hexs_cnf only if the configuration does not have max occupation
  ! in the selected GAS space
  if (I_ELIMINATE_GAS > 0) then
    ielim = 0
    if ((I_ELIMINATE_GAS == 1) .or. (I_ELIMINATE_GAS == 3)) then
      do j=1,N_ELIMINATED_GAS
        jGAS = IELIMINATED_IN_GAS(j)
        if (IOCCLS(jGAS,JOCCLS) == maxingas(j)) ielim = 1
      end do
    end if
    if (I_ELIMINATE_GAS > 1) then
      do j=1,N_2ELIMINATED_GAS
        jGAS = I2ELIMINATED_IN_GAS(j)
        if (IOCCLS(jGAS,JOCCLS) >= maxingas2(j)-1) ielim = 1
      end do
    end if
    if (ielim == 0) HEXS_CNF(:) = HEXS_CNF(:)+TMP_CNF(:)
  end if
  ! testing
  !write(u6,*) 'nconf_per_open after first call of gen_conf_for_occls'
  !call iwrtma(nconf_per_open,1,4,1,4)

  ! NCONF_ALL_SYM is accumulated, so
  !if (JOCCLS == 1) then
  !  NCONF_ALL_SYM_FOR_OCCLS(JOCCLS) = NCONF_ALL_SYM
  !  IBCONF_ALL_SYM_FOR_OCCLS(JOCCLS) = 1
  !else
  !  ! PAM2009: It was discovered that these two arrays could be overrun.
  !  ! The arrays are declared in lucia_data, and their dimension
  !  ! is MXPCSM, which is set in lucia_data -- both included above.
  !  ! So MXPCSM is now increased from 20 to 40 -- if this is not a final
  !  ! solution remains to be discovered:
  !  ! IFG 2025: It doesn't make sense for MXPCSM to be larger than 8,
  !  !           so if this is still a problem, there must be some bug.
  !  NCONF_ALL_SYM_FOR_OCCLS(JOCCLS) = NCONF_ALL_SYM-NCONF_ALL_SYM_PREV
  !  IBCONF_ALL_SYM_FOR_OCCLS(JOCCLS) = IBCONF_ALL_SYM_FOR_OCCLS(JOCCLS-1)+NCONF_ALL_SYM_FOR_OCCLS(JOCCLS-1)
  !end if
end do
! Number of CSF's in expansion
NCSF = sum(NCONF_PER_OPEN(1:MAXOP+1,ISYM)*NPCSCNF(1:MAXOP+1))
! Number of combinations in expansion
NCMB = sum(NCONF_PER_OPEN(1:MAXOP+1,ISYM)*NPCMCNF(1:MAXOP+1))
! Number of SD's in expansion
NSD = NCMB
if (PSSIGN /= Zero) NSD = NSD+NCMB-NCONF_PER_OPEN(1,ISYM)*NPDTCNF(1)
!MGD
if (I_ELIMINATE_GAS > 0) then
  NCSF_HEXS = sum(HEXS_CNF(1:MAXOP+1)*NPCSCNF(1:MAXOP+1))
else
  NCSF_HEXS = 0
end if

NCSF_PER_SYM(ISYM) = NCSF
NSD_PER_SYM(ISYM) = NCMB
NCONF_PER_SYM(ISYM) = sum(NCONF_PER_OPEN(:,ISYM))
if (IPRCIX >= 5) then
  write(u6,*) ' Number of CSFs  ',NCSF
  write(u6,*) ' Number of SDs   ',NSD
  write(u6,*) ' Number of Confs ',NCONF_PER_SYM(ISYM)
  write(u6,*) ' Number of CMBs  ',NCMB
end if

! Total number of configurations and length of configuration list
!    INFO_CONF_LIST(NCONF_PER_OPEN,MAXOP,NEL,LENGTH_LIST,NCONF_TOT,IB_REO,IB_OCC)
call INFO_CONF_LIST(NCONF_PER_OPEN(:,ISYM),MAXOP,NELEC,LENGTH_LIST,NCONF_TOT,IB_CONF_REO,IB_CONF_OCC)
! Permanent and local memory for csf routines

! memory for CSDTMT arrays.
! Largest block of proto type combinations .
! Largest number of prototype csf's

LIDT = 0
LICS = 0
LDTOC = 0
MXPTBL = 0
MXDT = 0
LCONF = 0
do IOPEN=MINOP,MAXOP
  ITP = IOPEN+1
  LIDT = LIDT+NPCMCNF(ITP)*IOPEN
  LICS = LICS+NPCSCNF(ITP)*IOPEN
  LDTOC = LDTOC+NPCSCNF(ITP)*NPCMCNF(ITP)
  MXDT = max(MXDT,NPCMCNF(ITP))
  MXPTBL = max(NPCMCNF(ITP)*IOPEN,MXPTBL)
end do
! Memory needed to store ICONF array
LCONF = 0
ILCNF = 0

!LDET = NCMB
LLCONF = 0
ILLCNF = 0
do IOPEN=MINOP,MAXOP
  ITYP = IOPEN+1
  ICL = (NELEC-IOPEN)/2
  LLCONF = LLCONF+NCONF_PER_OPEN(ITYP,ISYM)*(IOPEN+ICL)
  ILLCNF = ILLCNF+NCONF_PER_OPEN(ITYP,ISYM)
end do
!write(u6,*) ' MEMORY FOR HOLDING CONFS OF SYM... ',ISYM,LLCONF
LCONF = max(LCONF,LLCONF)
ILCNF = max(ILCNF,ILLCNF)

if (IPRCIX >= 5) then
  write(u6,'(/A,I8)') '  Memory for holding list of configurations ',LCONF
  write(u6,'(/A,I8)') '  Size of CI expansion (combinations)',NCMB
  write(u6,'(/A,I8)') '  Size of CI expansion (confs)',ILCNF
end if

! permanent memory for csf proto type arrays

call mma_allocate(DFTP,max(1,LIDT),Label='DFTP')
call mma_allocate(CFTP,max(1,LICS),Label='CFTP')
call mma_allocate(DTOC,max(1,LDTOC),Label='DTOC')

! PERMANENT ARRAYS FOR
! HOLDING CONFIGURATION EXPANSIONS AND REORDER ARRAYS

! Occupation of configurations
call mma_allocate(CONF_OCC(ISYM)%A,LCONF,Label='CONF_OCC()')
! Reorder array for configurations
call mma_allocate(CONF_REO(ISYM)%A,NCONF_TOT,Label='CONF_REO()')
! Reorder array for determinants, index and sign
call mma_allocate(SDREO_I(ISYM)%A,NCMB,Label='SDREO_I()')
SDREO(1:NCMB) => SDREO_I(ISYM)%A(:)

! Arrays for addressing prototype determinants for each number of orbitals

call Allocate_DT(Z_PTDT,[MINOP+1,MAXOP+1],Label='Z_PTDT')
call Allocate_DT(REO_PTDT,[MINOP+1,MAXOP+1],Label='REO_PTDT')
do IOPEN=MINOP,MAXOP
  ITYP = IOPEN+1

  IALPHA = (IOPEN+MS2)/2
  LZ = IOPEN*IALPHA
  LPTDT = IBINOM(IOPEN,IALPHA)
  call mma_allocate(Z_PTDT(ITYP)%A,LZ,Label='Z_PTDT()')
  call mma_allocate(REO_PTDT(ITYP)%A,LPTDT,Label='REO_PTDT()')
end do

! Array giving first determinant with given number of electrons
! in list of determinants ordered  according to the number of open orbitals

IB_SD_FOR_OPEN(:) = 0
IB = 1
do IOPEN=MINOP,MAXOP
  IB_SD_FOR_OPEN(IOPEN+1) = IB
  if (mod(IOPEN-MS2,2) == 0) IB = IB+NCONF_PER_OPEN(IOPEN+1,ISYM)*NPCMCNF(IOPEN+1)
end do

end subroutine CSFDIM_GAS
