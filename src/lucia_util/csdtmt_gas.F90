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
! Copyright (C) 1992,2001, Jeppe Olsen                                 *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine CSDTMT_GAS()
! Construct in IDTFTP list of proto type combinations in IDFTP
! Construct in ICFTP list of proto type CSF's in ICFTP
! Construct in DTOC matrix expanding proto type CSF's in terms of
! prototype combinations in DTOC
!
! Construct also array for going from lexical address of prototype det
! to address in IDFTP
!
! where MXPDBL is size of largest prototype combination block .
!
! Jeppe Olsen
!
! Changed to Combination form, June 92
! Adapted to LUCIA December 2001

use lucia_data, only: CFTP, DFTP, DTOC, IPRCIX, MAXOP, MINOP, MS2, MULTS, NPCMCNF, NPCSCNF, PSSIGN, REO_PTDT, Z_PTDT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: IALPHA, ICDCBS, ICSBS, IDET, IDTBS, IFLAG, IOPEN, IOPEN_, ITP, LSCR, MAX_DC, NNDET
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NNCM, NNCS
#endif
integer(kind=iwp), allocatable :: iSCR2(:)
real(kind=wp), allocatable :: SCR1(:)

! Size of largest csf-sd block
MAX_DC = max(0,maxval(NPCMCNF(1:MAXOP+1)))
#ifdef _DEBUGPRINT_
write(u6,*) ' Size of largest D to C block ',MAX_DC
#endif
LSCR = max(MAX_DC,MAXOP)
LSCR = MAX_DC*MAXOP+MAXOP
call mma_allocate(SCR1,LSCR,Label='SCR1')

! Set up combinations and upper determinants

if (IPRCIX >= 5) then
  write(u6,*)
  write(u6,*) ' *************************************'
  write(u6,*) ' Generation of proto type determinants'
  write(u6,*) ' *************************************'
  write(u6,*)
end if
! Still tired of stupid compiler warnings
IDTBS = 0
ICSBS = 0
do IOPEN=0,MAXOP
  ITP = IOPEN+1
  if (IPRCIX >= 5) then
    write(u6,*)
    write(u6,'(A,I3,A)') '       Type with ',IOPEN,' open orbitals'
    write(u6,'(A)') '       **********************************'
    write(u6,*)
  end if
  if (ITP == 1) then
    IDTBS = 1
    ICSBS = 1
  else
    IDTBS = IDTBS+(IOPEN-1)*NPCMCNF(ITP-1)
    ICSBS = ICSBS+(IOPEN-1)*NPCSCNF(ITP-1)
  end if

  if (IOPEN /= 0) then
    ! Proto type combinations and branching diagram for
    ! proto type combinations
    if (MS2+1 == MULTS) then
      IFLAG = 2
      call SPNCOM(IOPEN,MS2,NNDET,DFTP(IDTBS),CFTP(ICSBS),IFLAG,PSSIGN)
      !    SPNCOM(NOPEN,MS2,NDET,IABDET,IABUPP,IFLAG,PSSIGN)
    else
      IFLAG = 1
      call SPNCOM(IOPEN,MS2,NNDET,DFTP(IDTBS),CFTP(ICSBS),IFLAG,PSSIGN)
      IFLAG = 3
      call SPNCOM(IOPEN,MULTS-1,NNDET,DFTP(IDTBS),CFTP(ICSBS),IFLAG,PSSIGN)
    end if
  end if
end do
! End of loop over number of open orbitals

! Set up z-matrices for addressing prototype determinants with
! a given number of open orbitals, and for readdressing to
! the order given in IDCNF
! Scr : largest block of 2*NOPEN + (NALPHA+1)*(NOPEN+1)
LSCR = 0
do IOPEN=MINOP,MAXOP
  if (mod(IOPEN-MS2,2) == 0) then
    IALPHA = (IOPEN+MS2)/2
    LSCR = max(2*IOPEN+(IALPHA+1)*(IOPEN+1),LSCR)
  end if
end do
call mma_allocate(iSCR2,LSCR,Label='iSCR2')

IDTBS = 1
!-jwk do IOPEN=0,MAXOP
do IOPEN=MINOP,MAXOP
  ITP = IOPEN+1
  if (ITP == 1) then
    IDTBS = 1
  else
    IDTBS = IDTBS+(IOPEN-1)*NPCMCNF(ITP-1)
  end if

  IALPHA = (IOPEN+MS2)/2
  IDET = NPCMCNF(IOPEN+1)
  !write(u6,*) ' IOPEN, IDET = ',IOPEN,IDET
  IOPEN_ = IOPEN
  call REO_PTDET(IOPEN_,IALPHA,Z_PTDT(ITP)%A,REO_PTDT(ITP)%A,DFTP(IDTBS),IDET,iSCR2)
end do

! matrix expressing csf's in terms of determinants

! Tired of compiler warnings
IDTBS = 0
ICSBS = 0
ICDCBS = 0
do IOPEN=0,MAXOP
  ITP = IOPEN+1
  if (ITP == 1) then
    IDTBS = 1
    ICSBS = 1
    ICDCBS = 1
  else
    IDTBS = IDTBS+(IOPEN-1)*NPCMCNF(ITP-1)
    ICSBS = ICSBS+(IOPEN-1)*NPCSCNF(ITP-1)
    ICDCBS = ICDCBS+NPCMCNF(ITP-1)*NPCSCNF(ITP-1)
  end if
  !if (NPCMCNF(ITP)*NPCSCNF(ITP) == 0) goto 30
  if (IPRCIX >= 5) then
    write(u6,*)
    write(u6,*) ' ***********************************'
    write(u6,*) ' CSF - SD/COMB transformation matrix'
    write(u6,*) ' ***********************************'
    write(u6,'(A)')
    write(u6,'(A,I3,A)') '  Type with ',IOPEN,' open orbitals'
    write(u6,'(A)') '  ************************************'
    write(u6,*)
  end if
  if (IOPEN == 0) then
    DTOC(ICDCBS) = One
  else
    call CSFDET(IOPEN,DFTP(IDTBS),NPCMCNF(ITP),CFTP(ICSBS),NPCSCNF(ITP),DTOC(ICDCBS),SCR1,size(SCR1),PSSIGN)
    !    CSFDET(NOPEN,IDET,NDET,ICSF,NCSF,CDC,WORK,PSSIGN)
  end if
end do
! End of loop over number of open shells

call mma_deallocate(SCR1)
call mma_deallocate(iSCR2)

#ifdef _DEBUGPRINT_
write(u6,*) ' List of CSF-SD transformation matrices'
write(u6,*) ' ======================================'
write(u6,*)
IDTBS = 1
ICSBS = 1
ICDCBS = 1
do IOPEN=0,MAXOP
  ITP = IOPEN+1
  if (ITP == 1) then
    IDTBS = 1
    ICSBS = 1
    ICDCBS = 1
  else
    IDTBS = IDTBS+(IOPEN-1)*NPCMCNF(ITP-1)
    ICSBS = ICSBS+(IOPEN-1)*NPCSCNF(ITP-1)
    ICDCBS = ICDCBS+NPCMCNF(ITP-1)*NPCSCNF(ITP-1)
  end if
  NNCS = NPCSCNF(ITP)
  NNCM = NPCMCNF(ITP)
  if ((NNCS > 0) .and. (NNCM > 0)) then
    write(u6,*) ' Number of open shells : ',IOPEN
    write(u6,*) ' Number of combinations per conf ',NNCM
    write(u6,*) ' Number of CSFs per conf         ',NNCS
    call WRTMAT(DTOC(ICDCBS),NNCM,NNCS,NNCM,NNCS)
  end if
end do
#endif

end subroutine CSDTMT_GAS
