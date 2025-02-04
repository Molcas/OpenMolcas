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

subroutine GASDIAT(DIAG,LUDIA,ECORE,ICISTR,I12,IBLTP,NBLOCK,IBLKFO)
! CI diagonal in SD basis for state with symmetry ISM in internal
! space ISPC
!
! GAS version, Winter of 95
!
! Driven by table of TTS blocks, May97

use stdalloc, only: mma_allocate, mma_deallocate
use strbas, only: NSTSO
use lucia_data, only: IPRDIA
use lucia_data, only: PSSIGN
use lucia_data, only: MXNSTR, I_AM_OUT, N_ELIMINATED_BATCHES
use lucia_data, only: IDISK
use lucia_data, only: NTOOB, IREOST, IREOTS, NACOB
use lucia_data, only: NOCTYP
use lucia_data, only: NELEC
use csm_data, only: NSMST
#ifdef _DEBUGPRINT_
use lucia_data, only: IBSPGPFTP
use Definitions, only: u6
#endif

implicit none
! =====
! Input
! =====
integer LUDIA, ICISTR, I12, NBLOCK
real*8 ECORE
integer IBLTP(*)
integer IBLKFO(8,NBLOCK)
! ======
! Output
! ======
real*8 DIAG(*)
integer, allocatable :: LASTR(:), LBSTR(:)
real*8, allocatable :: LSCR2(:)
real*8, allocatable :: LJ(:), LK(:), LXB(:), LH1D(:), LRJKA(:)
integer, external :: IMNMX
integer NTEST, IATP, IBTP, NAEL, NBEL, NOCTPA, MAXA
#ifdef _DEBUGPRINT_
integer NOCTPB, IOCTPA, IOCTPB
#endif

NTEST = 0
NTEST = max(NTEST,IPRDIA)

! Specifications of internal space

IATP = 1
IBTP = 2
NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)
NOCTPA = NOCTYP(IATP)

#ifdef _DEBUGPRINT_
NOCTPB = NOCTYP(IBTP)
! Offsets for alpha and beta supergroups
IOCTPA = IBSPGPFTP(IATP)
IOCTPB = IBSPGPFTP(IBTP)
if (NTEST >= 10) then
  write(u6,*) ' ==============='
  write(u6,*) ' GASDIA speaking'
  write(u6,*) ' ==============='
  write(u6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
  write(u6,*) ' NOCTPA NOCTPB  : ',NOCTPA,NOCTPB
  write(u6,*) ' IOCTPA IOCTPB  : ',IOCTPA,IOCTPB
end if
#endif

! Local memory

call mma_allocate(LJ,NTOOB**2,Label='LJ')
call mma_allocate(LK,NTOOB**2,Label='LK')
call mma_allocate(LSCR2,2*NTOOB**2,Label='LSCR2')
call mma_allocate(LXB,NACOB,Label='LXB')
call mma_allocate(LH1D,NACOB,Label='LH1D')
! Space for blocks of strings
call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
MAXA = IMNMX(NSTSO(IATP)%I,NSMST*NOCTPA,2)
call mma_allocate(LRJKA,MAXA,Label='LRJKA')

! Diagonal of one-body integrals and coulomb and exchange integrals

call GT1DIA(LH1D)
call GTJK(LJ,LK,NTOOB,LSCR2,IREOTS,IREOST)
if (LUDIA > 0) IDISK(LUDIA) = 0
call GASDIAS(NAEL,LASTR,NBEL,LBSTR,NACOB,DIAG,NSMST,LH1D,LXB,LJ,LK,NSTSO(IATP)%I,NSTSO(IBTP)%I,LUDIA,ECORE,PSSIGN,IPRDIA,NTOOB, &
             ICISTR,LRJKA,I12,IBLTP,NBLOCK,IBLKFO,I_AM_OUT,N_ELIMINATED_BATCHES)
! Flush local memory
call mma_deallocate(LJ)
call mma_deallocate(LK)
call mma_deallocate(LSCR2)
call mma_deallocate(LXB)
call mma_deallocate(LH1D)
call mma_deallocate(LASTR)
call mma_deallocate(LBSTR)
call mma_deallocate(LRJKA)

end subroutine GASDIAT
