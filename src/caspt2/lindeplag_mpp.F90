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

subroutine LinDepLag_MPP(lg_BDER,lg_SDER,lg_GDER,nAS,nIN,iSym,iCase)
! Parallel LinDepLag

use GA_Wrapper, only: DBL_MB, GA_Destroy, GA_NodeId
#ifdef _SCALAPACK_
use scalapack_mod, only: GA_PDSYEVX_
#endif
use caspt2_global, only: idBoriMat, LUSTD
use caspt2_module, only: IFDORTHO, THRSHN
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: lg_BDER, lg_SDER, lg_GDER, nAS, nIN, iSym, iCase
integer(kind=iwp) :: I, IDB, iHi, iLo, J, jHi, jLo, LDB, LDG, LDS, LDV, lg_B, lg_Lag, lg_S, lg_Vec, M, mB, mG, mS, mV, myRank
real(kind=wp) :: Tmp
logical(kind=iwp) :: bStat
real(kind=wp), allocatable :: EIG(:), SCA(:), SD(:)
#if ! defined (_SCALAPACK_)
real(kind=wp) :: WGRONK(2)
integer(kind=iwp) :: info, NSCRATCH
real(kind=wp), allocatable :: SCRATCH(:), VEC(:)
#endif

!! Obtain the X matrix
!! First, read S
call PSBMAT_GETMEM('S',lg_S,NAS)
call PSBMAT_READ('S',iCase,iSym,lg_S,NAS)

myRank = GA_NodeID()

! Scale factors from the diagonal of S, and scale S in the same way as in SBDIAG
call mma_allocate(SD,NAS,Label='SD')
SD(:) = Zero
call GA_Distribution(lg_S,myRank,iLo,iHi,jLo,jHi)
if (max(iLo,jLo) /= 0) then
  call GA_Access(lg_S,iLo,iHi,jLo,jHi,mS,LDS)
  do I=max(iLo,jLo),min(iHi,jHi)
    SD(I) = DBL_MB(mS+I-iLo+LDS*(I-jLo))
  end do
  call GA_Release(lg_S,iLo,iHi,jLo,jHi)
end if
call GADGOP(SD,NAS,'+')
call mma_allocate(SCA,NAS,Label='SCA')
do I=1,NAS
  if (IFDORTHO) then
    SCA(I) = One
  else if (SD(I) > THRSHN) then
    SCA(I) = (One+real(I,kind=wp)*3.0e-6_wp)/sqrt(SD(I))
  else
    SCA(I) = Zero
  end if
end do
call mma_deallocate(SD)
if (iLo /= 0) then
  call GA_Access(lg_S,iLo,iHi,jLo,jHi,mS,LDS)
  call S_SCALE(NAS,SCA,DBL_MB(mS),iLo,iHi,jLo,jHi,LDS)
  call GA_Release_Update(lg_S,iLo,iHi,jLo,jHi)
end if
call GA_Sync()

call mma_allocate(EIG,NAS,Label='EIG')
EIG(:) = Zero
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

! The first nAS-nIN eigenvectors are the removed ones.
! Divide the rows of the eigenvectors by SCA to get VR, which couples with the densities.
call GA_Sync()
call GA_Distribution(lg_Vec,myRank,iLo,iHi,jLo,jHi)
if (iLo /= 0) then
  call GA_Access(lg_Vec,iLo,iHi,jLo,jHi,mV,LDV)
  if ((jHi-jLo+1) /= NAS) then
    write(u6,*) 'SBDIAG_MPP: error in striping of lg_Vec, ABORT'
    call ABEND()
  end if
  do I=1,iHi-iLo+1
    if (SCA(iLo+I-1) /= Zero) call DScal_(NAS,One/SCA(iLo+I-1),DBL_MB(mV+I-1),LDV)
  end do
  call GA_Release_Update(lg_Vec,iLo,iHi,jLo,jHi)
end if
call GA_Sync()

call GA_CREATE_STRIPED('H',NAS,NAS,'Lag',lg_Lag)
call PSBMAT_GETMEM('B',lg_B,NAS)

call GA_Distribution(lg_B,myRank,iLo,iHi,jLo,jHi)
if ((iLo > 0) .and. (jLo > 0)) then
  call GA_Access(lg_B,iLo,iHi,jLo,jHi,mB,LDB)
  IDB = IDBoriMat(ISYM,ICASE)
  call DDAFILE(LUSTD,2,DBL_MB(mB),(iHi-iLo+1)*(jHi-jLo+1),IDB)
  call GA_Release_Update(lg_B,iLo,iHi,jLo,jHi)
end if

!! Compute the partial derivative
! lg_B   : B
! lg_BDER: D
! lg_Vec : VR and VL
call GA_DGEMM('N','T',NAS,NAS,NAS,Two,lg_B,lg_BDER,Zero,lg_Lag)
call GA_DGEMM('N','N',NAS,NAS,NAS,One,lg_Lag,lg_VEC,Zero,lg_B)

! Add explicit derivative of 2<1|H|0> (lg_GDER) to the S derivative, then replace it with VR^T * G^T * VR
call GA_Distribution(lg_SDER,myRank,iLo,iHi,jLo,jHi)
if (iLo /= 0) then
  call GA_Access(lg_SDER,iLo,iHi,jLo,jHi,mV,LDV)
  call GA_Access(lg_GDER,iLo,iHi,jLo,jHi,mG,LDG)
  do J=1,jHi-jLo+1
    do I=1,iHi-iLo+1
      DBL_MB(mV+I-1+LDV*(J-1)) = DBL_MB(mV+I-1+LDV*(J-1))+DBL_MB(mG+I-1+LDG*(J-1))
    end do
  end do
  call GA_Release(lg_GDER,iLo,iHi,jLo,jHi)
  call GA_Release_Update(lg_SDER,iLo,iHi,jLo,jHi)
end if
call GA_Sync()
call GA_DGEMM('T','N',NAS,NAS,NAS,One,lg_GDER,lg_Vec,Zero,lg_Lag)
call GA_DGEMM('T','N',NAS,NAS,NAS,One,lg_Vec,lg_Lag,Zero,lg_GDER)

! Multiply the rows of lg_Vec by SCA**2 to get VL, which couples with the residual and the derivative of S
call GA_Distribution(lg_Vec,myRank,iLo,iHi,jLo,jHi)
if (iLo /= 0) then
  call GA_Access(lg_Vec,iLo,iHi,jLo,jHi,mV,LDV)
  do I=1,iHi-iLo+1
    call DScal_(NAS,SCA(iLo+I-1)**2,DBL_MB(mV+I-1),LDV)
  end do
  call GA_Release_Update(lg_Vec,iLo,iHi,jLo,jHi)
end if
call GA_Sync()
call GA_DGEMM('T','N',NAS,NAS,NAS,One,lg_Vec,lg_B,Zero,lg_Lag)
!! At this point,
! lg_Lag = 2 VL^T * B * D * VR

! Keep only the rotations between the removed (row) and remaining (column) vectors,
! and add the RHS part with the eigenvalue differences
call GA_Distribution(lg_Lag,myRank,iLo,iHi,jLo,jHi)
if (iLo /= 0) then
  call GA_Access(lg_Lag,iLo,iHi,jLo,jHi,mV,LDV)
  call GA_Access(lg_GDER,iLo,iHi,jLo,jHi,mG,LDG)
  do J=1,nAS
    do I=1,iHi-iLo+1
      if ((iLo+I-1 <= nAS-nIN) .and. (J > nAS-nIN)) then
        DBL_MB(mV+I-1+LDV*(J-1)) = (DBL_MB(mV+I-1+LDV*(J-1))+EIG(iLo+I-1)*DBL_MB(mG+I-1+LDG*(J-1)))/(EIG(J)-EIG(iLo+I-1))
      else
        DBL_MB(mV+I-1+LDV*(J-1)) = Zero
      end if
    end do
  end do
  call GA_Release(lg_GDER,iLo,iHi,jLo,jHi)
  call GA_Release_Update(lg_Lag,iLo,iHi,jLo,jHi)
end if
call GA_Sync()

!! orthogonal -> non-orthogonal
!! Finalize Eq. (62)
call GA_DGEMM('N','N',NAS,NAS,NAS,One,lg_Vec,lg_Lag,Zero,lg_B)
call GA_DGEMM('N','T',NAS,NAS,NAS,One,lg_B,lg_Vec,One,lg_SDER)

! Derivative of the scale factors, which depend on the diagonal elements of S
if (.not. IFDORTHO) then
  call GA_Distribution(lg_SDER,myRank,iLo,iHi,jLo,jHi)
  if (iLo /= 0) then
    call GA_Access(lg_SDER,iLo,iHi,jLo,jHi,mS,LDS)
    call GA_Access(lg_B,iLo,iHi,jLo,jHi,mB,LDB)
    call GA_Access(lg_Vec,iLo,iHi,jLo,jHi,mV,LDV)
    do I=1,iHi-iLo+1
      M = iLo+I-1
      if (SCA(M) == Zero) cycle
      Tmp = Zero
      do J=nAS-nIN+1,nAS
        Tmp = Tmp+DBL_MB(mB+I-1+LDB*(J-1))*EIG(J)*DBL_MB(mV+I-1+LDV*(J-1))
      end do
      DBL_MB(mS+I-1+LDS*(M-1)) = DBL_MB(mS+I-1+LDS*(M-1))-Tmp/(One+real(M,kind=wp)*3.0e-6_wp)**2
    end do
    call GA_Release(lg_Vec,iLo,iHi,jLo,jHi)
    call GA_Release(lg_B,iLo,iHi,jLo,jHi)
    call GA_Release_Update(lg_SDER,iLo,iHi,jLo,jHi)
  end if
  call GA_Sync()
end if

call mma_deallocate(EIG)
call mma_deallocate(SCA)
call PSBMAT_FREEMEM(lg_Vec)
bSTAT = GA_Destroy(lg_Lag)
#include "macros.fh"
unused_var(bStat)
call PSBMAT_FREEMEM(lg_B)

end subroutine LinDepLag_MPP

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(LinDepLag_MPP)

#endif
