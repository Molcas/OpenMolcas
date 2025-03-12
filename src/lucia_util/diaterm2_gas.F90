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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DIATERM2_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,J12,JDC)
! = DIATERM_GAS, just J12 added !
!
! Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
! Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
!
! For the NBLOCKS givem in IBLOCK starting from BLOCK IOFF
!
! Jeppe Olsen, August 1995

use lucia_data, only: ECORE, ECORE_ORIG, IREOST, MXNSTR, NACOB, NELEC, NIRREP, NOCTYP, NSTSO, NTOOB
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use lucia_data, only: IBSPGPFTP, IDC
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: FACTOR
real(kind=wp), intent(_OUT_) :: VEC(*)
integer(kind=iwp), intent(in) :: ITASK, NBLOCK, IBLOCK(8,*), IOFF, J12, JDC
integer(kind=iwp) :: IATP, IBTP, MAXA, NAEL, NBEL, NOCTPA
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IOCTPA, IOCTPB, NOCTPB
#endif
real(kind=wp) :: ECOREP, FACTORX, SHIFT
integer(kind=iwp), allocatable :: LASTR(:), LBSTR(:)
real(kind=wp), allocatable :: LH1D(:), LJ(:), LK(:), LRJKA(:), LXB(:)
integer(kind=iwp), external :: IMNMX

IATP = 1
IBTP = 2
NAEL = NELEC(IATP)
NBEL = NELEC(IBTP)
NOCTPA = NOCTYP(IATP)

!if (JPERT == 0) then
!  ! Use full Hamiltonian
!  I12 = 2
!  IPERTOP = 0
!else
!  ! Use perturbation operator
!  if (IPART == 1) then
!    ! Moller-Plesset partitioning
!    I12 = 1
!    IPERTOP = 1
!  else if (IPART == 2) then
!    ! Epstein-Nesbet Partitioning
!    I12 = 2
!    IPERTOP = 0
!  end if
!end if

#ifdef _DEBUGPRINT_
! Offsets for alpha and beta supergroups
NOCTPB = NOCTYP(IBTP)
IOCTPA = IBSPGPFTP(IATP)
IOCTPB = IBSPGPFTP(IBTP)
write(u6,*) ' ====================='
write(u6,*) ' DIATERM2_GAS speaking'
write(u6,*) ' ====================='
write(u6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
write(u6,*) ' NOCTPA NOCTPB  : ',NOCTPA,NOCTPB
write(u6,*) ' IOCTPA IOCTPB  : ',IOCTPA,IOCTPB
write(u6,*) ' IPART,J12',J12
#endif
! A bit of scracth
call mma_allocate(LJ,NTOOB**2,Label='LJ')
call mma_allocate(LK,NTOOB**2,Label='LK')
call mma_allocate(LXB,NACOB,Label='LX')
call mma_allocate(LH1D,NACOB,Label='LH1D')
! Space for blocks of strings
call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
MAXA = IMNMX(NSTSO(IATP)%A,NIRREP*NOCTPA,2)
call mma_allocate(LRJKA,MAXA,Label='LRJKA')
! Diagonal of one-body integrals and coulomb and exchange integrals
! Integrals assumed in place so :
call GT1DIA(LH1D)
if (J12 == 2) call GTJK(LJ,LK,NTOOB,IREOST)
! Core energy not included
ECOREP = Zero
SHIFT = ECORE_ORIG-ECORE
FACTORX = FACTOR+SHIFT
call DIATERMS_GAS(NAEL,LASTR,NBEL,LBSTR,NACOB,VEC,NIRREP,LH1D,JDC,LXB,LJ,LK,NSTSO(IATP)%A,NSTSO(IBTP)%A,ECOREP,0,0,NTOOB,LRJKA, &
                  J12,IBLOCK(:,IOFF),NBLOCK,ITASK,FACTORX,0,[0])
! Flush local memory
call mma_deallocate(LJ)
call mma_deallocate(LK)
call mma_deallocate(LXB)
call mma_deallocate(LH1D)
call mma_deallocate(LASTR)
call mma_deallocate(LBSTR)
call mma_deallocate(LRJKA)

#ifdef _DEBUGPRINT_
write(u6,*) ' output vector from DIATRM'
call WRTTTS(VEC,IBLOCK(:,IOFF),NBLOCK,NIRREP,NSTSO(IATP)%A,NSTSO(IBTP)%A,IDC)
#endif

end subroutine DIATERM2_GAS
