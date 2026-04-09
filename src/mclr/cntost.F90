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
subroutine CNTOST(ICONF,ICTSDT,NAEL,NBEL,IPRODT,IREFSM,NORB,NEL,IGENSG,ISGNA,ISGNB,IAGRP,IBGRP,IOOS,PSSIGN)
! Obtain pointer abs(ICTSDT(I)) giving address of determinant I in
! STRING ordering for determinant I in CSF ordering.
! Going between the two formats can involve a sign change. this is
! stored in the sign of ICTSDT)
! SGNCTS is thus to be multiplied with vector ordered in CSF ordering.
!
! December 1990 : NCNFCN,ICNFOK added
! January 1991  : IGENSG,ISGNA,ISGNB added
! April   1991  : LUCIA version
! September 1993 > Sign and address stored together

use MCLR_Data, only: MINOP, NCNATS, NDPCNT, NTYP
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICONF(*), NAEL, NBEL, IPRODT(*), IREFSM, NORB, IGENSG, ISGNA(*), ISGNB(*), IAGRP, IBGRP, IOOS(*)
integer(kind=iwp), intent(_OUT_) :: ICTSDT(*)
integer(kind=iwp), intent(out) :: NEL
real(kind=wp), intent(in) :: PSSIGN
integer(kind=iwp) :: IABNUM, IC, ICL, ICNBS, ICNBS0, ICNF, IDET, IJKL_NUM, IOCC, IOPEN, IPBAS, IPSFAC, ISGN, ISGNAB, ITYP, JDET, &
                     JDTABS, MXDT
integer(kind=iwp), allocatable :: LDTBL(:), LIA(:), LIB(:), SCR23(:)

! IWORK should at least be of length (MXDT+2)*NEL,
! where MXDT is the largest number of prototype determinants occuring
! in a single block.

NEL = NAEL+NBEL

! Local memory

! Largest number of dets for a given type
MXDT = max(0,maxval(NDPCNT(1:NTYP)))
call mma_allocate(LDTBL,MXDT*NEL,Label='LDTBL')
call mma_allocate(LIA,NAEL,Label='LIA')
call mma_allocate(LIB,NBEL,Label='LIB')
call mma_allocate(SCR23,NEL,Label='SCR23')

! Loop over configurations and generate determinants in compact form

ICNF = 0
JDTABS = 0
IPSFAC = 0 ! Removes a compiler error
ISGNAB = 0 ! Removes a compiler error
ICNBS0 = 0 ! dummy initialize
IPBAS = 0  ! dummy initialize
ijkl_num = 0 !yma counter
do ITYP=1,NTYP
  IDET = NDPCNT(ITYP)
  IOPEN = ITYP+MINOP-1
  ICL = (NEL-IOPEN)/2
  IOCC = IOPEN+ICL
  if (ITYP == 1) then
    ICNBS0 = 1
  else
    ICNBS0 = ICNBS0+NCNATS(ITYP-1,IREFSM)*(NEL+IOPEN-1)/2
  end if
  ! Base for prototype determinants
  if (ITYP == 1) then
    IPBAS = 1
  else
    IPBAS = IPBAS+NDPCNT(ITYP-1)*(IOPEN-1)
  end if
  ! Determinants for this configuration
  do IC=1,NCNATS(ITYP,IREFSM)
    ICNF = ICNF+1
    ICNBS = ICNBS0+(IC-1)*(IOPEN+ICL)
    ! Check orbital occupancy with additional constraints
#   ifdef _DEBUGPRINT_
    write(u6,*) ' IC ICNF ICNBS',IC,ICNF,ICNBS
#   endif
    call CNDET_MCLR(ICONF(ICNBS),IPRODT(IPBAS),IDET,NEL,IOCC,IOPEN,ICL,LDTBL)
    ! Separate determinants into strings and determine string number.
    do JDET=1,IDET
      !write(117,'(1X,I8,1X,A,1X)',advance='no') ITYP,'ITYP'  ! yma
      JDTABS = JDTABS+1
      call DETSTR_MCLR(LDTBL(1+(JDET-1)*NEL),LIA,LIB,NEL,NAEL,NBEL,ISGN,SCR23)
      ijkl_num = ijkl_num+1
      ! Find number (and sign)of this determinant in string ordering
      ICTSDT(JDTABS) = IABNUM(LIA,LIB,IAGRP,IBGRP,IGENSG,ISGNA,ISGNB,ISGNAB,IOOS,NORB,IPSFAC,PSSIGN)
      if (ISGN*ISGNAB*IPSFAC == -1) ICTSDT(JDTABS) = -ICTSDT(JDTABS)
    end do
  end do
end do
call mma_deallocate(SCR23)
call mma_deallocate(LIB)
call mma_deallocate(LIA)
call mma_deallocate(LDTBL)

end subroutine CNTOST
