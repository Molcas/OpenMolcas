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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

#include "compiler_features.h"
#if defined(_MOLCAS_MPP_) && defined(_GA_)

subroutine LinDepLag_MPP(lg_BDER,lg_SDER,nAS,nIN,iSym,iCase)
! Parallel LinDepLag
! We always use the canonical orthonormalization.

#ifdef _SCALAPACK_
use scalapack_mod, only: GA_PDSYEVX_
#endif
use caspt2_global, only: idBoriMat, LUSTD
use caspt2_module, only: THRSHS
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: lg_BDER, lg_SDER, nAS, nIN, iSym, iCase
integer(kind=iwp) :: I, IDB, iHi, iLo, J, jHi, jLo, LDB, LDV, lg_B, lg_Lag, lg_S, lg_Vec, mB, mV, myRank
real(kind=wp) :: EVAL, FACT
logical(kind=iwp) :: bStat
real(kind=wp), allocatable :: EIG(:)
#if ! defined (_SCALAPACK_)
real(kind=wp) :: WGRONK(2)
integer(kind=iwp) :: info, NSCRATCH
real(kind=wp), allocatable :: SCRATCH(:), VEC(:)
#endif
#include "global.fh"
#include "mafdecls.fh"

!! Obtain the X matrix
!! First, read S
call PSBMAT_GETMEM('S',lg_S,NAS)
call PSBMAT_READ('S',iCase,iSym,lg_S,NAS)

call mma_allocate(EIG,NAS,Label='EIG')
EIG(:) = Zero

myRank = GA_NodeID()
#ifdef _SCALAPACK_
call PSBMAT_GETMEM('VMAT',lg_Vec,NAS)
call GA_PDSYEVX_(lg_S,lg_Vec,EIG,0)
bSTAT = GA_Destroy(lg_S)
#else
! here for the non-ScaLAPACK version: copy matrix to master process,
! diagonalize using the serial DSYEV routine, and copy the resulting
! eigenvectors back to a global array.  Then distribute the eigenvalues.
if (myRank == 0) then
  call mma_allocate(VEC,NAS*NAS,Label='VEC')
  call GA_Get(lg_S,1,NAS,1,NAS,VEC,NAS)
end if
bSTAT = GA_Destroy(lg_S)
if (myRank == 0) then
  call DSYEV_('V','L',NAS,VEC,NAS,EIG,WGRONK,-1,INFO)
  NSCRATCH = int(WGRONK(1))
  call mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
  call DSYEV_('V','L',NAS,VEC,NAS,EIG,SCRATCH,NSCRATCH,INFO)
  call mma_deallocate(SCRATCH)
end if
call PSBMAT_GETMEM('VMAT',lg_Vec,NAS)
if (myRank == 0) then
  call GA_Put(lg_Vec,1,NAS,1,NAS,VEC,NAS)
  call mma_deallocate(VEC)
end if
call GADGOP(EIG,NAS,'+')
#endif

!! Scale only the independent vectors to avoid
!! any numerically unstable computation
! Form orthonormal transformation vectors by scaling the eigenvectors.
call GA_Sync()
call GA_Distribution(lg_Vec,myRank,iLo,iHi,jLo,jHi)
if (iLo /= 0) then
  call GA_Access(lg_Vec,iLo,iHi,jLo,jHi,mV,LDV)
  if ((jHi-jLo+1) /= NAS) then
    write(u6,*) 'SBDIAG_MPP: error in striping of lg_Vec, ABORT'
    call ABEND()
  end if
  do I=1,jHi-jLo+1 ! NAS
    EVAL = EIG(I)
    if (EVAL < THRSHS) cycle
    FACT = One/sqrt(EVAL)
    call DScal_(iHi-iLo+1,FACT,DBL_MB(mV+LDV*(I-1)),1)
  end do
  call GA_Release_Update(lg_Vec,iLo,iHi,jLo,jHi)
end if
call GA_Sync()

call mma_deallocate(EIG)

! matrix has been prepared

call GA_CREATE_STRIPED('H',NAS,NAS,'Lag',lg_Lag)
call PSBMAT_GETMEM('B',lg_B,NAS)

call GA_Distribution(lg_B,myRank,iLo,iHi,jLo,jHi)
call GA_Access(lg_B,iLo,iHi,jLo,jHi,mB,LDB)
IDB = IDBoriMat(ISYM,ICASE)
call DDAFILE(LUSTD,2,DBL_MB(mB),(iHi-iLo+1)*(jHi-jLo+1),IDB)
call GA_Release_Update(lg_B,iLo,iHi,jLo,jHi)

!! Compute the partial derivative
!! Work(LF)  : B --> lg_B
!! BDER      : D --> lg_BDER
!! Work(LVEC): X^0 and X --> lg_Vec
call GA_DGEMM('N','T',NAS,NAS,NAS,Two,lg_B,lg_BDER,Zero,lg_Lag)
call GA_DGEMM('N','N',NAS,NAS,NAS,One,lg_Lag,lg_VEC,Zero,lg_B)
call GADupl(lg_B,lg_Lag)

call GA_DGEMM('T','N',NAS,NAS,NAS,One,lg_Vec,lg_Lag,Zero,lg_B)
!! At this point,
!! Work(LF) = 2 \mathcal{X}^0 * B * D * \mathcal{X}

!! remove dependent part
!! (linearly indep-indep and dep-dep)
call GA_Distribution(lg_B,myRank,iLo,iHi,jLo,jHi)
if (iLo /= 0) then
  call GA_Access(lg_B,iLo,iHi,jLo,jHi,mV,LDV)
  !if (ilo <= nas-nin) then
  do J=1,nAS-nIN
    do I=1,min(iHi-iLo+1,nAS-nIN-iLo+1)
      DBL_MB(mV+i-1+LDV*(j-1)) = Zero
    end do
  end do
  !end if
  !if (ilo >= nas-nin+1) then
  do J=nAS-nIN+1,nAS
    do I=max(1,nAS-nIN+1-iLo+1),LDV
      DBL_MB(mV+i-1+LDV*(j-1)) = Zero
    end do
  end do
  !end if
  call GA_Release_Update(lg_B,iLo,iHi,jLo,jHi)
end if

!! orthogonal -> non-orthogonal
!! Finalize Eq. (62)
call GA_DGEMM('N','N',NAS,NAS,NAS,One,lg_Vec,lg_B,Zero,lg_Lag)
call GA_DGEMM('N','T',NAS,NAS,NAS,One,lg_Lag,lg_Vec,One,lg_SDER)

call PSBMAT_FREEMEM(lg_Vec)
bSTAT = GA_Destroy(lg_Lag)
call PSBMAT_FREEMEM(lg_B)

#include "macros.fh"
unused_var(bStat)

end subroutine LinDepLag_MPP

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(LinDepLag_MPP)

#endif
