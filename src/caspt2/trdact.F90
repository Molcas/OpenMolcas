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

subroutine TRDACT(IVEC,JVEC,DTU)
! Add to the active-active block of transition density matrix,
!    D(t,u) = Add <IVEC| E(t,u) |JVEC> = <0| W1T E(t,u) W2 |0>
! where t,u are active indices. IVEC and JVEC are integer labels,
! denoting sets of coefficients stored on LUSOLV. These are assumed
! to be contravariant representations of the wave operators W1 and W2,
! in the notation of the comments.

use Index_Functions, only: nTri_Elem, nTri3_Elem
use sguga_states, only: SGS
use caspt2_global, only: IDTCEX, LUCIEX
use general_data, only: nAsh, STSym, nLev
use caspt2_module, only: iASym, iSCF, jState, MxCI, nAes, nAshT, nAshT, nConf, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC
real(kind=wp), intent(inout) :: DTU(NASHT,NASHT)

integer(kind=iwp) :: I, ID, ISYM, ISYMT, ITABS, ITLEV, IU, IUABS, IULEV, NOP1, NOP2, NOP3
real(kind=wp) :: OCCNUM, OP0, SCP
integer(kind=iwp), allocatable :: IATOG(:)
real(kind=wp), allocatable :: TRDCI(:), TRDOP1(:), TRDOP2(:), TRDOP3(:), TRDSGM(:), TRDTMP(:)
real(kind=wp), external :: DDOT_
integer(kind=iwp), parameter :: istate=1

! (1): Compute a representation of the operator PCAS*W1T*W2
NOP1 = NASHT**2
NOP2 = nTri_Elem(NOP1)
NOP3 = nTri3_Elem(NOP1)
call MMA_ALLOCATE(TRDOP1,NOP1)
call MMA_ALLOCATE(TRDOP2,NOP2)
call MMA_ALLOCATE(TRDOP3,NOP3)
call MKWWOP(IVEC,JVEC,OP0,TRDOP1,NOP2,TRDOP2,NOP3,TRDOP3)

! (2): Compute the state vector |Temp> = (PCAS*W1T*W2) |0>
! First modify the coefficients, see subroutine MODOP.
call MODOP(TRDOP1,NOP2,TRDOP2,NOP3,TRDOP3)
call MMA_ALLOCATE(TRDTMP,NCONF)
call MMA_ALLOCATE(TRDCI,NCONF)
if (ISCF == 0) then
  ID = IDTCEX(JSTATE)
  call DDAFILE(LUCIEX,2,TRDCI,NCONF,ID)
else
  TRDCI(1) = One
end if
TRDTMP(:) = Zero
call HAM3(OP0,TRDOP1,NOP2,TRDOP2,NOP3,TRDOP3,STSYM,TRDCI,TRDTMP,nCONF)
! No more need for the operators:
call MMA_DEALLOCATE(TRDOP1)
call MMA_DEALLOCATE(TRDOP2)
call MMA_DEALLOCATE(TRDOP3)

! (3) compute <0| E(t,u) W1T W2 |0> as <Sigma_ut| Temp>, and add to DTU

if (ISCF == 0) then
  ! Create reorder table giving the GUGA level, i.e. CI-coupling
  ! ordinal number of each active orbital.
  call mma_allocate(IATOG,NLEV,Label='IATOG')
  ITABS = 0
  do ISYM=1,NSYM
    do I=1,NLEV
      if (SGS(istate)%ISM(I) == ISYM) then
        ITABS = ITABS+1
        IATOG(ITABS) = I
      end if
    end do
  end do
  call MMA_ALLOCATE(TRDSGM,MXCI)
  do ITABS=1,NASHT
    ISYMT = IASYM(ITABS)
    ITLEV = IATOG(ITABS)
    do IU=1,NASH(ISYMT)
      IUABS = NAES(ISYMT)+IU
      IULEV = IATOG(IUABS)
      call GETSGM2(IULEV,ITLEV,STSYM,TRDCI,NCONF,TRDSGM,NCONF)
      SCP = DDOT_(NCONF,TRDSGM,1,TRDTMP,1)
      DTU(ITABS,IUABS) = DTU(ITABS,IUABS)+SCP
    end do
  end do
  call MMA_DEALLOCATE(TRDSGM)
  call MMA_DEALLOCATE(IATOG)
else
  OCCNUM = Two
  if (ISCF == 2) OCCNUM = One
  do ITABS=1,NASHT
    DTU(ITABS,ITABS) = DTU(ITABS,ITABS)+OCCNUM
  end do
end if
! No more need for CI array
call MMA_DEALLOCATE(TRDCI)
! No more need for the TMP state vector
call MMA_DEALLOCATE(TRDTMP)

! (4): Add the correction <0| [W1T,E(tu)] W2 |0>.
call COMMWEW(IVEC,JVEC,DTU)

end subroutine TRDACT
