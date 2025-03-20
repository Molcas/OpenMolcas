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

!#define _DEBUGPRINT_
subroutine ADDDIA_TERM(FACTOR,CVEC,SVEC,IASPGP,IBSPGP,IASM,IBSM)
! Update Sigma vector with diagonal terms for a given block
!     SVEC(IASPGP,IBSPGP) = SVEC(IASPGP,IBSPGP)
!                         + (FACTOR+DIAG(IASPGP,IBSPGP))CVEC(IASPGP,IBSPGP)
!
! Jeppe Olsen and Giovanni Li Manni, September 2011

use lucia_data, only: ECORE, ECORE_ORIG, IREOST, MXNSTR, NACOB, NELEC, NIRREP, NOCTYP, NSTSO, NTOOB
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(in) :: FACTOR, CVEC(*)
real(kind=wp), intent(inout) :: SVEC(*)
integer(kind=iwp), intent(in) :: IASPGP, IBSPGP, IASM, IBSM
integer(kind=iwp) :: IATP, IBTP, MAXA, NAEL, NBEL, NOCTPA
real(kind=wp) :: ECOREP, FACTORX, SHIFT
integer(kind=iwp), allocatable :: LASTR(:), LBSTR(:)
real(kind=wp), allocatable :: LH1D(:), LJ(:), LK(:), LRJKA(:), LXB(:)
integer(kind=iwp), external :: IMNMX

IATP = 1
IBTP = 2
NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)
NOCTPA = NOCTYP(IATP)

#ifdef _DEBUGPRINT_
write(u6,*) ' =============================='
write(u6,*) ' ADDDIA_TERM for BK is speaking'
write(u6,*) ' =============================='
write(u6,*) ' NAEL NBEL =',NAEL,NBEL
write(u6,*) ' IASPGP, IBSPGP = ',IASPGP,IBSPGP
#endif
! A bit of scracth
call mma_allocate(LH1D,NTOOB,Label='LH1D')
call mma_allocate(LJ,NTOOB**2,Label='LJ')
call mma_allocate(LK,NTOOB**2,Label='LK')
call mma_allocate(LXB,NACOB,Label='LXB')
! Space for blocks of strings
call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')

MAXA = IMNMX(NSTSO(IATP)%A,NIRREP*NOCTPA,2)
call mma_allocate(LRJKA,MAXA,Label='LRJKA')
! Diagonal of one-body integrals and coulomb and exchange integrals
! Integrals assumed in place so :
call GT1DIA(LH1D)
call GTJK(LJ,LK,NTOOB,IREOST)
! Core energy not included
ECOREP = Zero
call GTJK(LJ,LK,NTOOB,IREOST)

SHIFT = ECORE_ORIG-ECORE
FACTORX = FACTOR+SHIFT

call ADDDIA_TERMS(NAEL,LASTR,NBEL,LBSTR,NACOB,CVEC,SVEC,NIRREP,LH1D,LXB,LJ,LK,NSTSO(IATP)%A,NSTSO(IBTP)%A,ECOREP,NTOOB,LRJKA, &
                  IASPGP,IASM,IBSPGP,IBSM,FACTORX)
! Flush local memory
call mma_deallocate(LH1D)
call mma_deallocate(LJ)
call mma_deallocate(LK)
call mma_deallocate(LXB)
call mma_deallocate(LASTR)
call mma_deallocate(LBSTR)
call mma_deallocate(LRJKA)

end subroutine ADDDIA_TERM
