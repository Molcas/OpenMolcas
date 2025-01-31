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
! Copyright (C) 2011, Jeppe Olsen                                      *
!               2011, Giovanni Li Manni                                *
!***********************************************************************

subroutine ADDDIA_TERM(FACTOR,CVEC,SVEC,IASPGP,IBSPGP,IASM,IBSM)
! Update Sigma vector with diagonal terms for a given block
!     SVEC(IASPGP,IBSPGP) = SVEC(IASPGP,IBSPGP)
!                         + (FACTOR+DIAG(IASPGP,IBSPGP))CVEC(IASPGP,IBSPGP)
!
! Jeppe Olsen and Giovanni Li Manni, September 2011

use stdalloc, only: mma_allocate, mma_deallocate
use strbas, only: NSTSO
use lucia_data, only: ECORE_ORIG, ECORE
use lucia_data, only: IPRDIA
use lucia_data, only: MXNSTR
use lucia_data, only: NTOOB, NACOB, IREOST, IREOTS
use lucia_data, only: NOCTYP
use lucia_data, only: NELEC
use csm_data, only: NSMST

implicit none
real*8 FACTOR
integer IASPGP, IBSPGP, IASM, IBSM
! Input
real*8 CVEC(*)
! Output
real*8 SVEC(*)
integer, allocatable :: LASTR(:), LBSTR(:)
real*8, allocatable :: LSCR(:), LSCR2(:)
real*8, allocatable :: LJ(:), LK(:), LXA(:), LXB(:), LRJKA(:), LH1D(:)
integer, external :: IMNMX
integer NTEST, IATP, IBTP, NAEL, NBEL, NOCTPA, MAXA
real*8 ECOREP, SHIFT, FACTORX

NTEST = 0
NTEST = max(NTEST,IPRDIA)

IATP = 1
IBTP = 2
NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)
NOCTPA = NOCTYP(IATP)

#ifdef _DEBUGPRINT_
if (NTEST >= 10) then
  write(6,*) ' =============================='
  write(6,*) ' ADDDIA_TERM for BK is speaking'
  write(6,*) ' =============================='
  write(6,*) ' NAEL NBEL =',NAEL,NBEL
  write(6,*) ' IASPGP, IBSPGP = ',IASPGP,IBSPGP
end if
#endif
! A bit of scracth
call mma_allocate(LH1D,NTOOB,Label='LH1D')
call mma_allocate(LJ,NTOOB**2,Label='LJ')
call mma_allocate(LK,NTOOB**2,Label='LK')
call mma_allocate(LSCR2,2*NTOOB**2,Label='LSCR2')
call mma_allocate(LXA,NACOB,Label='LXA')
call mma_allocate(LXB,NACOB,Label='LXB')
call mma_allocate(LSCR,2*NACOB,Label='LSCR')
! Space for blocks of strings
call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')

MAXA = IMNMX(NSTSO(IATP)%I,NSMST*NOCTPA,2)
call mma_allocate(LRJKA,MAXA,Label='LRJKA')
! Diagonal of one-body integrals and coulomb and exchange integrals
! Integrals assumed in place so :
call GT1DIA(LH1D)
call GTJK(LJ,LK,NTOOB,LSCR2,IREOTS,IREOST)
! Core energy not included
ECOREP = 0.0d0
call GTJK(LJ,LK,NTOOB,LSCR2,IREOTS,IREOST)

SHIFT = ECORE_ORIG-ECORE
FACTORX = FACTOR+SHIFT

call ADDDIA_TERMS(NAEL,LASTR,NBEL,LBSTR,NACOB,CVEC,SVEC,NSMST,LH1D,LXA,LXB,LSCR,LJ,LK,NSTSO(IATP)%I,NSTSO(IBTP)%I,ECOREP,IPRDIA, &
                  NTOOB,LRJKA,IASPGP,IASM,IBSPGP,IBSM,FACTORX)
! Flush local memory
call mma_deallocate(LH1D)
call mma_deallocate(LJ)
call mma_deallocate(LK)
call mma_deallocate(LSCR2)
call mma_deallocate(LXA)
call mma_deallocate(LXB)
call mma_deallocate(LSCR)
call mma_deallocate(LASTR)
call mma_deallocate(LBSTR)
call mma_deallocate(LRJKA)

end subroutine ADDDIA_TERM
