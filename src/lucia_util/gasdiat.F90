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

!#define _DEBUGPRINT_
subroutine GASDIAT(DIAG,LUDIA,ECORE,ICISTR,I12,IBLTP,NBLOCK,IBLKFO)
! CI diagonal in SD basis for state with symmetry ISM in internal
! space ISPC
!
! GAS version, Winter of 95
!
! Driven by table of TTS blocks, May97

use lucia_data, only: I_AM_OUT, IDISK, IREOST, MXNSTR, N_ELIMINATED_BATCHES, NACOB, NELEC, NIRREP, NOCTYP, NSTSO, NTOOB, PSSIGN
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use lucia_data, only: IBSPGPFTP
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: DIAG(*)
integer(kind=iwp), intent(in) :: LUDIA, ICISTR, I12, IBLTP(*), NBLOCK, IBLKFO(8,NBLOCK)
real(kind=wp), intent(in) :: ECORE
integer(kind=iwp) :: IATP, IBTP, MAXA, NAEL, NBEL, NOCTPA
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IOCTPA, IOCTPB, NOCTPB
#endif
integer(kind=iwp), allocatable :: LASTR(:), LBSTR(:)
real(kind=wp), allocatable :: LH1D(:), LJ(:), LK(:), LRJKA(:), LXB(:)
integer(kind=iwp), external :: IMNMX

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
write(u6,*) ' ==============='
write(u6,*) ' GASDIA speaking'
write(u6,*) ' ==============='
write(u6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
write(u6,*) ' NOCTPA NOCTPB  : ',NOCTPA,NOCTPB
write(u6,*) ' IOCTPA IOCTPB  : ',IOCTPA,IOCTPB
#endif

! Local memory

call mma_allocate(LJ,NTOOB**2,Label='LJ')
call mma_allocate(LK,NTOOB**2,Label='LK')
call mma_allocate(LXB,NACOB,Label='LXB')
call mma_allocate(LH1D,NACOB,Label='LH1D')
! Space for blocks of strings
call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
MAXA = IMNMX(NSTSO(IATP)%A,NIRREP*NOCTPA,2)
call mma_allocate(LRJKA,MAXA,Label='LRJKA')

! Diagonal of one-body integrals and coulomb and exchange integrals

call GT1DIA(LH1D)
call GTJK(LJ,LK,NTOOB,IREOST)
if (LUDIA > 0) IDISK(LUDIA) = 0
call GASDIAS(NAEL,LASTR,NBEL,LBSTR,NACOB,DIAG,NIRREP,LH1D,LXB,LJ,LK,NSTSO(IATP)%A,NSTSO(IBTP)%A,LUDIA,ECORE,PSSIGN,NTOOB,ICISTR, &
             LRJKA,I12,IBLTP,NBLOCK,IBLKFO,I_AM_OUT,N_ELIMINATED_BATCHES)
! Flush local memory
call mma_deallocate(LJ)
call mma_deallocate(LK)
call mma_deallocate(LXB)
call mma_deallocate(LH1D)
call mma_deallocate(LASTR)
call mma_deallocate(LBSTR)
call mma_deallocate(LRJKA)

end subroutine GASDIAT
