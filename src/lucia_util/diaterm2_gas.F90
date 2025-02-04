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

subroutine DIATERM2_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,J12,JDC)
! = DIATERM_GAS, just J12 added !
!
! Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
! Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
!
! For the NBLOCKS givem in IBLOCK starting from BLOCK IOFF
!
! Jeppe Olsen, August 1995

use stdalloc, only: mma_allocate, mma_deallocate
use strbas, only: NSTSO
use lucia_data, only: ECORE_ORIG, ECORE
use lucia_data, only: IPRDIA
use lucia_data, only: MXNSTR
use lucia_data, only: NTOOB, IREOST, NACOB
use lucia_data, only: NOCTYP
use lucia_data, only: NELEC
use csm_data, only: NSMST
use Constants, only: Zero
#ifdef _DEBUGPRINT_
use lucia_data, only: IDC, IPERTOP, IBSPGPFTP
use Definitions, only: u6
#endif

implicit none
real*8 FACTOR
integer ITASK, IBLOCK(8,*)
real*8 VEC(*)
integer IOFF, J12, JDC
integer, allocatable :: LASTR(:), LBSTR(:)
real*8, allocatable :: LJ(:), LK(:), LXB(:), LH1D(:), LRJKA(:)
integer, external :: IMNMX
integer NTEST, IATP, IBTP, NAEL, NBEL, NOCTPA, MAXA, NBLOCK
real*8 ECOREP, SHIFT, FACTORX
#ifdef _DEBUGPRINT_
integer NOCTPB, IOCTPA, IOCTPB
#endif

NTEST = 0
NTEST = max(NTEST,IPRDIA)

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
if (NTEST >= 10) then
  write(u6,*) ' ====================='
  write(u6,*) ' DIATERM2_GAS speaking'
  write(u6,*) ' ====================='
  write(u6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
  write(u6,*) ' NOCTPA NOCTPB  : ',NOCTPA,NOCTPB
  write(u6,*) ' IOCTPA IOCTPB  : ',IOCTPA,IOCTPB
  write(u6,*) ' IPART,J12,IPERTOP',J12,IPERTOP
end if
#endif
! A bit of scracth
call mma_allocate(LJ,NTOOB**2,Label='LJ')
call mma_allocate(LK,NTOOB**2,Label='LK')
call mma_allocate(LXB,NACOB,Label='LX')
call mma_allocate(LH1D,NACOB,Label='LH1D')
! Space for blocks of strings
call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
MAXA = IMNMX(NSTSO(IATP)%I,NSMST*NOCTPA,2)
call mma_allocate(LRJKA,MAXA,Label='LRJKA')
! Diagonal of one-body integrals and coulomb and exchange integrals
! Integrals assumed in place so :
call GT1DIA(LH1D)
if (J12 == 2) call GTJK(LJ,LK,NTOOB,IREOST)
! Core energy not included
ECOREP = Zero
SHIFT = ECORE_ORIG-ECORE
FACTORX = FACTOR+SHIFT
call DIATERMS_GAS(NAEL,LASTR,NBEL,LBSTR,NACOB,VEC,NSMST,LH1D,JDC,LXB,LJ,LK,NSTSO(IATP)%I,NSTSO(IBTP)%I,ECOREP,0,0,IPRDIA,NTOOB, &
                  LRJKA,J12,IBLOCK(1,IOFF),NBLOCK,ITASK,FACTORX,0,[0])
!                           IBLOCK,NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)
! Flush local memory
call mma_deallocate(LJ)
call mma_deallocate(LK)
call mma_deallocate(LXB)
call mma_deallocate(LH1D)
call mma_deallocate(LASTR)
call mma_deallocate(LBSTR)
call mma_deallocate(LRJKA)

#ifdef _DEBUGPRINT_
if (NTEST >= 100) then
  write(u6,*) ' output vector from DIATRM'
  call WRTTTS(VEC,IBLOCK(1,IOFF),NBLOCK,NSMST,NSTSO(IATP)%I,NSTSO(IBTP)%I,IDC)
end if
#endif

end subroutine DIATERM2_GAS
