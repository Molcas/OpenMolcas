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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

! SVC2010: The following subroutine computes the transformation
! matrices using global arrays rather than replicate arrays.  Currently,
! the maximum memory scales approximately as 3*(NAS**2), which is twice
! what is needed in the replicate routines.  It can be reduced to
! 2*(NAS**2), if the matrix multiplications are rewritten to occur in
! batch mode.  However, unlike in the replicate routine, this amount is
! divided over processors.
#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

subroutine SBDIAG_MPP(ISYM,ICASE,CONDNR,CPU)
! On entry, the DRA metafiles contain the matrices S and B for cases A
! (iCASE=1) en C (iCASE=4).  These symmetric matrices are stored on disk
! in full square format.  The rectangular matrices T are computed, such
! that Sum_(J) [B(I,J)*T(J,MU)] = Sum_(J) [S(I,J)*T(J,MU)*BD(MU)], where
! I,J are in the range (1,NASUP(ISYM,ICASE)) and MU is in the range
! (1,NINDEP(ISYM,ICASE)).  NINDEP is the numerically effective rank of
! S, and the columns of T are orthonormal: Sum_(I,J)
! [S(I,J)*T(J,MU)*T(J,NU)] = Kron(MU,NU).  B is destroyed, and is
! overwritten by BD(MU) and T(I,MU), which is stored in a DRA metafile.

use Index_Functions, only: nTri_Elem
#ifdef _SCALAPACK_
use scalapack_mod, only: GA_PDSYEVX_
#endif
use PrintLevel, only: INSANE
use Para_Info, only: King
use GA_Wrapper, only: DBL_MB, GA_Destroy, GA_NodeId
use EQSOLV, only: IDBMAT, IDSMAT, IDTMAT
use caspt2_global, only: do_grad, do_lindep, idBoriMat, iPrGlb, LUSBT, LUSTD, nStpGrd
use caspt2_module, only: nASup, nISup, Cases, IfDOrtho, ThrShn, ThrShs, nInDep, BMATRIX, BTRANS, BSPECT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iSym, iCase
real(kind=wp), intent(out) :: CondNr, CPU
integer(kind=iwp) :: I, IDB, IDB2, IDS, IDT, IEND, iHi, iLo, IOFF, ISTA, J, jHi, jLo, lDB, LDS, LDV, lg_B, lg_S, lg_ST, lg_T, &
                     lg_V, lg_X, mB, MS, mV, MyRank, NAS, NCOEF, NCOL, NIN, NIS, NTMP
real(kind=wp) :: CPU1, CPU2, CPUE, dTrans, FP, SDiag, SZMAX, SZMIN, TIO, TIOE
logical(kind=iwp) :: bSTAT
character(len=2) :: cCASE, cSYM
real(kind=wp), allocatable :: BD(:), COL(:), COND(:), EIG(:), SCA(:), SD(:), TMP(:), TRANS(:)
#ifndef _SCALAPACK_
integer(kind=iwp) :: Info, NSCRATCH
real(kind=wp) :: WGRONK(2)
real(kind=wp), allocatable :: SCRATCH(:), VEC(:)
#endif
real(kind=wp), external :: PSBMAT_FPRINT, DNRM2_

! Initialize the DRA I/O subsystem with default values.

if ((iCASE /= 1) .and. (iCASE /= 4)) then
  write(u6,*) 'Invalid CASE number used for global SBDIAG, Abort'
  call AbEnd()
end if

write(cCase,'(I2.2)') iCase
write(cSYM,'(I2.2)') iSYM
! Start a long loop over irreps:

CPU = Zero
CONDNR = Zero
NAS = NASUP(ISYM,ICASE)
NIS = NISUP(ISYM,ICASE)
NCOEF = NAS*NIS
if (NCOEF == 0) return

if (IPRGLB >= INSANE) &
  write(u6,'("DEBUG> ",A12,A5,I2,A2,A6,A2,A5,I1)') 'SBDIAG_MPP: ','CASE ',ICASE,' (',CASES(ICASE),') ','SYM ',ISYM

! Allocate memory for the S matrix and read it from disk:
call PSBMAT_GETMEM('S',lg_S,NAS)
call PSBMAT_READ('S',iCase,iSym,lg_S,NAS)
if (IPRGLB >= INSANE) then
  FP = PSBMAT_FPRINT(lg_S,NAS)
  write(u6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM: ',FP
end if

! The S matrices are needed later on by non-global routines.  Take the
! oportunity to save them to LUSBT here.  FIXME: Should be removed once
! full parallelization of use of S matrices is achieved.
if (KING()) then
  NCOL = NAS
  NTMP = nTri_Elem(NAS)
  call mma_allocate(COL,NCOL,Label='COL')
  call mma_allocate(TMP,NTMP,Label='TMP')
  iOFF = 0
  do J=1,NAS
    call GA_Get(lg_S,1,J,J,J,COL,NAS)
    TMP(iOFF+1:iOFF+J) = COL(1:J)
    iOFF = iOFF+J
  end do
  call mma_deallocate(COL)
  IDS = IDSMAT(ISYM,ICASE)
  call DDAFILE(LUSBT,1,TMP,NTMP,IDS)
  call mma_deallocate(TMP)
end if

! Save the diagonal elements from the S matrix for easy access later on.
! FIXME: nicer way to do this?
call mma_allocate(SD,NAS,Label='SD')
SD(:) = Zero
myRank = GA_NodeID()
call GA_Distribution(lg_S,myRank,iLo,iHi,jLo,jHi)
ISTA = max(ILO,JLO)
IEND = min(IHI,JHI)
if (ISTA /= 0) then
  call GA_Access(lg_S,iLo,iHi,jLo,jHi,mS,LDS)
  do I=ISTA,IEND
    SD(I) = DBL_MB(mS+I-ILO+LDS*(I-JLO))
  end do
  call GA_Release(lg_S,iLo,iHi,jLo,jHi)
end if
call GADGOP(SD,NAS,'+')

! Calculate the scaling factors and store them in array SCA.
call mma_allocate(SCA,NAS,Label='SCA')
do I=1,NAS
  SDiag = SD(I)
  if (IFDORTHO) then
    SCA(I) = One
  else
    if (SDiag > THRSHN) then
      SCA(I) = (One+real(I,kind=wp)*3.0e-6_wp)/sqrt(SDiag)
    else
      SCA(I) = Zero
    end if
  end if
end do

! Scale the elements S(I,J) with the factor SCA(I)*SCA(J).
myRank = GA_NodeID()
call GA_Distribution(lg_S,myRank,iLo,iHi,jLo,jHi)
if (iLo /= 0) then
  call GA_Access(lg_S,iLo,iHi,jLo,jHi,mS,LDS)
  call S_SCALE(NAS,SCA,DBL_MB(mS),iLo,iHi,jLo,jHi,LDS)
  call GA_Release_Update(lg_S,iLo,iHi,jLo,jHi)
end if
call GA_Sync()
call TIMING(CPU1,CPUE,TIO,TIOE)

if (IPRGLB >= INSANE) then
  FP = PSBMAT_FPRINT(lg_S,NAS)
  write(u6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM AFTER SCALING: ',FP
end if

call mma_allocate(EIG,NAS,Label='EIG')
EIG(:) = Zero
! Diagonalize the global array S.  Some (old) reports about parallel
! performance recommend PDSYEVX or PDSYEVR as fastest methods if
! eigenvectors are needed (FIXME: should time this).  For the linear
! dependence removal, split eigenvectors in horizontal stripes so that
! each processor has a row window of all column vectors
#ifdef _SCALAPACK_
call PSBMAT_GETMEM('VMAT',lg_V,NAS)
call GA_PDSYEVX_(lg_S,lg_V,EIG,0)
bSTAT = GA_Destroy(lg_S)
#else
! here for the non-ScaLAPACK version: copy matrix to master process,
! diagonalize using the serial DSYEV routine, and copy the resulting
! eigenvectors back to a global array.  Then distribute the eigenvalues.
if (myRank == 0) then
  call mma_allocate(VEC,NAS**2,Label='VEC')
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
call PSBMAT_GETMEM('VMAT',lg_V,NAS)
if (myRank == 0) then
  call GA_Put(lg_V,1,NAS,1,NAS,VEC,NAS)
  call mma_deallocate(VEC)
end if
call GADGOP(EIG,NAS,'+')
#endif

if (IPRGLB >= INSANE) then
  FP = DNRM2_(NAS,EIG,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'SMAT EIGENVALUE NORM: ',FP
end if

NIN = count(EIG(:) >= THRSHS)
NINDEP(ISYM,ICASE) = NIN
if (NIN == 0) then
  call mma_deallocate(SCA)
  call mma_deallocate(EIG)
  call mma_deallocate(SD)
  bSTAT = GA_Destroy(lg_V)
  return
end if

call TIMING(CPU2,CPUE,TIO,TIOE)
CPU = CPU+CPU2-CPU1

call mma_allocate(COND,NIN,Label='COND')
COND(:) = Zero
! Form orthonormal transformation vectors by scaling the eigenvectors.
call GA_Sync()
myRank = GA_NodeID()
call GA_Distribution(lg_V,myRank,iLo,iHi,jLo,jHi)
if (iLo /= 0) then
  call GA_Access(lg_V,iLo,iHi,jLo,jHi,mV,LDV)
  if ((jHi-jLo+1) /= NAS) then
    write(u6,*) 'SBDIAG_MPP: error in striping of lg_V, ABORT'
    call ABEND()
  end if
  call V_SCALE(EIG,SCA(iLo),DBL_MB(mV),iHi-iLo+1,jHi-jLo+1,LDV,NIN,COND)
  call GA_Release_Update(lg_V,iLo,iHi,jLo,jHi)
end if
call GA_Sync()

call mma_deallocate(EIG)
call mma_deallocate(SCA)

! The condition number, after scaling, disregarding linear dep.
! FIXME: adapt to local subroutine for global array lg_V
if (NIN >= 2) then
  call GADGOP(COND,NIN,'+')
  SZMIN = min(1.0e99_wp,minval(COND(:)))
  SZMAX = max(Zero,maxval(COND(:)))
  CONDNR = SZMAX/SZMIN
end if
call mma_deallocate(COND)

! Copy the NIN non-linear dependent eigenvectors to the transformation
! matrix T(NAS,NIN).
call GA_CREATE_STRIPED('H',NAS,NIN,'TMAT',lg_T)
call GA_Copy_Patch('N',lg_V,1,NAS,1,NIN,lg_T,1,NAS,1,NIN)
bStat = GA_Destroy(lg_V)

if (BMATRIX == 'NO') then
  ! In some calculations, we do not use B matrices.
  ! Write the T matrix to disk and exit.  FIXME: This
  ! should be removed when the transformation matrices are stored as disk
  ! resident arrays only.
  if (KING()) then
    call mma_allocate(TRANS,NAS*NIN,Label='TRANS')
    call GA_Get(lg_T,1,NAS,1,NIN,TRANS,NAS)
    if (iPrGlb >= INSANE) then
      dTRANS = dNRM2_(NAS*NIN,TRANS,1)
      write(u6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ',dTRANS
    end if
    IDT = IDTMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
    call mma_deallocate(TRANS)
  end if
  bStat = GA_Destroy(lg_T)
  return
else if (BTRANS /= 'YES') then
  ! In other calculations, the B matrix is used but not transformed.  We
  ! may need the diagonal active energies, i.e. the diagonal values of B
  ! divided by the diagonal values of S. These are placed where the
  ! eigenvalues would go in ordinary CASPT2.
  call PSBMAT_GETMEM('B',lg_B,NAS)
  call PSBMAT_READ('B',iCase,iSym,lg_B,NAS)
  call mma_allocate(BD,NAS,Label='BD')
  BD(:) = Zero
  myRank = GA_NodeID()
  call GA_Distribution(lg_B,myRank,iLo,iHi,jLo,jHi)
  ISTA = max(ILO,JLO)
  IEND = min(IHI,JHI)
  if (ISTA /= 0) then
    call GA_Access(lg_B,iLo,iHi,jLo,jHi,mB,LDB)
    do I=ISTA,IEND
      BD(I) = DBL_MB(mB+I-ILO+LDB*(I-JLO))
    end do
    call GA_Release(lg_B,iLo,iHi,jLo,jHi)
  end if
  call GADGOP(BD,NAS,'+')
  bStat = GA_Destroy(lg_B)
  BD(:) = BD(:)/(SD(:)+1.0e-15_wp)
  IDB = IDBMAT(ISYM,ICASE)
  call DDAFILE(LUSBT,1,BD,NAS,IDB)
  call mma_deallocate(BD)
  ! Write the transformation matrices after the diagonal values of B.
  ! FIXME: This should be removed when the transformation matrices are
  ! stored as disk resident arrays only.
  if (KING()) then
    call mma_allocate(TRANS,NAS*NIN,Label='TRANS')
    call GA_Get(lg_T,1,NAS,1,NIN,TRANS,NAS)
    if (iPrGlb >= INSANE) then
      dTRANS = dNRM2_(NAS*NIN,TRANS,1)
      write(u6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ',dTRANS
    end if
    IDT = IDTMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
    call mma_deallocate(TRANS)
  end if
  bStat = GA_Destroy(lg_T)
  return
end if
call mma_deallocate(SD)

! TRANSFORM B MATRIX TO O-N BASIS. BUT FIRST, SAVE O-N VECTORS.
call PSBMAT_WRITE('T',iCase,iSym,lg_T,NAS*NIN)

if (IPRGLB >= INSANE) then
  FP = PSBMAT_FPRINT(lg_T,NAS)
  write(u6,'(1X,A,ES21.14)') 'EIGENVECTOR NORM BEFORE B TRANS: ',FP
end if

call PSBMAT_GETMEM('B',lg_B,NAS)
call PSBMAT_READ('B',iCase,iSym,lg_B,NAS)

if ((do_grad .or. (nStpGrd == 2)) .and. do_lindep) then
  !! The original B matrix is needed in the LinDepLag subroutine
  IDB2 = idBoriMat(ISYM,ICASE)
  call GA_Distribution(lg_B,myRank,iLo,iHi,jLo,jHi)
  call GA_Access(lg_B,iLo,iHi,jLo,jHi,mB,LDB)
  call DDAFILE(LUSTD,1,DBL_MB(mB),(iHi-iLo+1)*(jHi-jLo+1),IDB2)
  call GA_Release(lg_B,iLo,iHi,jLo,jHi)
end if

if (IPRGLB >= INSANE) write(u6,'(1X,A,ES21.14)') 'BMAT NORM: ',FP

! FIXME: Perform transformation of B using horizontal stripes of B or
! vertical stripes of T to reduce memory usage if necessary as indicated
! by the available memory, which is now scaling as approx. 3*(NAS**2).
call GA_CREATE_STRIPED('H',NAS,NIN,'XMAT',lg_X)
call GA_DGEMM('N','N',NAS,NIN,NAS,One,lg_B,lg_T,Zero,lg_X)
bStat = GA_Destroy(lg_B)

call GA_CREATE_STRIPED('H',NIN,NIN,'BMAT',lg_B)
call GA_DGEMM('T','N',NIN,NIN,NAS,One,lg_T,lg_X,Zero,lg_B)
bStat = GA_Destroy(lg_X)
bStat = GA_Destroy(lg_T)

if (IPRGLB >= INSANE) then
  FP = PSBMAT_FPRINT(lg_B,NIN)
  write(u6,'(1X,A,ES21.14)') 'BMAT NORM AFTER TRANS: ',FP
end if

call TIMING(CPU1,CPUE,TIO,TIOE)

! Diagonalize the transformed B matrix.
call mma_allocate(EIG,NIN,Label='EIG')
EIG(:) = Zero
if (BSPECT /= 'YES') then
  ! Use diagonal approxim., if allowed.
  !call GA_Fill(lg_V, Zero)
  call GA_Zero(lg_V)
  ! FIXME: this original code seemed wrong, using uninitialized SD?
  !IDIAG = 1
  !do I=1,NIN
  !  EIG(I) = B(IDIAG)/SD
  !  IDIAG = IDIAG+1+NIN-I
  !end do
  write(u6,*) 'GLOB_SBDIAG: option not implemented'
  call AbEnd()
else
# ifdef _SCALAPACK_
  call GA_CREATE_STRIPED('H',NIN,NIN,'VMAT',lg_V)
  call GA_PDSYEVX_(lg_B,lg_V,EIG,0)
  bStat = GA_Destroy(lg_B)
# else
  if (myRank == 0) then
    call mma_allocate(VEC,NIN**2,Label='VEC')
    call GA_Get(lg_B,1,NIN,1,NIN,VEC,NIN)
  end if
  bSTAT = GA_Destroy(lg_B)
  if (myRank == 0) then
    call dsyev_('V','L',NIN,VEC,NIN,EIG,WGRONK,-1,INFO)
    NSCRATCH = int(WGRONK(1))
    call mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
    call dsyev_('V','L',NIN,VEC,NIN,EIG,SCRATCH,NSCRATCH,INFO)
    call mma_deallocate(SCRATCH)
  end if
  call GA_Sync()
  call GA_CREATE_STRIPED('H',NIN,NIN,'VMAT',lg_V)
  if (myRank == 0) then
    call GA_Put(lg_V,1,NIN,1,NIN,VEC,NIN)
    call mma_deallocate(VEC)
  end if
  call GADGOP(EIG,NIN,'+')
# endif
end if

if (IPRGLB >= INSANE) then
  FP = DNRM2_(NIN,EIG,1)
  write(u6,'(1X,A,ES21.14)') 'BMAT EIGENVALUE NORM: ',FP
end if

! The eigenvalues are written back at same position as the
! original B matrix, which is destroyed:
IDB = IDBMAT(ISYM,ICASE)
call DDAFILE(LUSBT,1,EIG,NIN,IDB)
call mma_deallocate(EIG)

call TIMING(CPU2,CPUE,TIO,TIOE)
CPU = CPU+CPU2-CPU1

! Finally, we must form the composite transformation matrix: T(NAS,NIN)
! matrix on disk * V(NIN*NIN) matrix in core.  FIXME: for now, asume
! there is enough memory for the full transformation, scaling as
! approx. 3*(NAS**2).  Should be determined by the available memory.
call GA_CREATE_STRIPED('H',NAS,NIN,'XMAT',lg_X)
call PSBMAT_READ('T',iCase,iSym,lg_X,NAS*NIN)
call GA_CREATE_STRIPED('H',NAS,NIN,'TMAT',lg_T)
call GA_DGEMM('N','N',NAS,NIN,NIN,One,lg_X,lg_V,Zero,lg_T)
bStat = GA_Destroy(lg_X)
bStat = GA_Destroy(lg_V)

! Write the composite transformation matrix to disk.
call PSBMAT_WRITE('T',iCase,iSym,lg_T,NAS*NIN)

! Additonally, compute S*T and store in on disk for later use by the RHS
! vector utitlities
call PSBMAT_GETMEM('S',lg_S,NAS)
call PSBMAT_READ('S',iCase,iSym,lg_S,NAS)

call GA_CREATE_STRIPED('H',NAS,NIN,'STMAT',lg_ST)
call GA_DGEMM('N','N',NAS,NIN,NAS,One,lg_S,lg_T,Zero,lg_ST)
bStat = GA_Destroy(lg_S)

call PSBMAT_WRITE('M',iCase,iSym,lg_ST,NAS*NIN)
bStat = GA_Destroy(lg_ST)

! For now, also keep the transformation matrix on disk as a
! replicate array.  FIXME: Should be removed later.
if (KING()) then
  call mma_allocate(TRANS,NAS*NIN,Label='TRANS')
  call GA_Get(lg_T,1,NAS,1,NIN,TRANS,NAS)
  dTRANS = dNRM2_(NAS*NIN,TRANS,1)
  if (iPrGlb >= INSANE) write(u6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ',dTRANS
  IDT = IDTMAT(ISYM,ICASE)
  call DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
  call mma_deallocate(TRANS)
end if

call ga_sync()
bStat = GA_Destroy(lg_T)

#include "macros.fh"
unused_var(bStat)

end subroutine SBDIAG_MPP

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(SBDIAG_MPP)

#endif
