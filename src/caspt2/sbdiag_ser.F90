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

subroutine SBDIAG_SER(ISYM,ICASE,CONDNR,CPU)
! On entry, the file LUSBT contains overlap matrices SMAT at disk
! addresses IDSMAT(ISYM,ICASE), ISYM=1..NSYM, ICASE=1..11, and
! similarly BMAT matrices at IDBMAT(ISYM,ICASE).
! These matrices are stored in a triangular format.
! The rectangular matrices TRANS are computed, such that
!    Sum(J) of BMAT(I,J)*TRANS(J,MU) =
!       Sum(J) of  SMAT(I,J)*TRANS(J,MU)*BDIAG(MU)
! where I,J are in 1..NASUP(ISYM,ICASE)
!        MU is in 1..NINDEP(ISYM,ICASE)
! NINDEP is the numerically effective rank of SMAT, and the
! columns of TRANS are orthonormal:
!    Sum(I,J) of  SMAT(I,J)*TRANS(J,MU)*TRANS(J,NU) = Kron(MU,NU)
! BMAT is destroyed, and is overwritten by BDIAG(MU) and
! TRANS(I,MU), with addresses IDBMAT(ISYM,ICASE) and
! IDTMAT(ISYM,ICASE). Enough file space was thus originally
! reserved on LUSBT to allow this overlay.
! LUSOLV is assumed not to be in use yet, so we use it
! for temporary storage.

use PrintLevel, only: INSANE
use EQSOLV, only: IDBMAT, IDSMAT, IDSTMAT, IDTMAT
use caspt2_global, only: do_grad, do_lindep, idBoriMat, iPrGlb, LUSBT, LUSOLV, LUSTD, nStpGrd
use caspt2_module, only: BMatrix, BSpect, BTrans, Cases, IfDOrtho, nASup, nG3, nInDep, nISup, ThrShn, ThrShs
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, ItoB

implicit none
integer(kind=iwp), intent(in) :: iSym, iCase
real(kind=wp), intent(out) :: CondNr, CPU
integer(kind=iwp) :: I, IDB, IDB2, IDIAG, IDS, IDST, IDT, IDTMP, IDTMP0, IJ, INFO, iPad, J, JOFF, KEND, KSTA, LTRANS1, NAS, NAUX, &
                     NB, NBNEW, NCOEF, NCOL, NIN, NIS, NS, NSCRATCH
real(kind=wp) :: CPU1, CPU2, CPUE, FP, SDIAG, SZ, SZMAX, SZMIN, TIO, TIOE, WGRONK(2)
real(kind=wp), allocatable :: AUX(:), B(:), BD(:), BX(:), EIG(:), S(:), SCA(:), SCRATCH(:), SD(:), ST(:), TRANS(:), VEC(:,:), XBX(:)
real(kind=wp), external :: DNRM2_

SDiag = Zero ! dummy initialize

CPU = Zero
CONDNR = Zero
NAS = NASUP(ISYM,ICASE)
NIS = NISUP(ISYM,ICASE)
NCOEF = NAS*NIS

if (NCOEF == 0) return

if (IPRGLB >= INSANE) &
  write(u6,'("DEBUG> ",A12,A7,I2,A2,A6,A2,A5,I1)') 'SBDIAG_SER: ','CASE ',ICASE,' (',CASES(ICASE),') ','SYM ',ISYM

IDTMP0 = 0
if (do_grad .or. (nStpGrd == 2)) then
  !! correct?
  iPad = ItoB-mod(6*NG3,ItoB)
  IDTMP0 = 6*NG3+iPad
end if

NS = (NAS*(NAS+1))/2
call mma_allocate(S,NS,Label='S')
IDS = IDSMAT(ISYM,ICASE)
call DDAFILE(LUSBT,2,S,NS,IDS)
if (IPRGLB >= INSANE) then
  FP = DNRM2_(NS,S,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM: ',FP
end if

! For some purposes, we need to save the diagonal elements:
if (BMATRIX == 'YES') then
  if (BTRANS /= 'YES') then
    call mma_allocate(SD,NAS,Label='SD')
    IDIAG = 0
    do I=1,NAS
      IDIAG = IDIAG+I
      SD(I) = S(IDIAG)
    end do
  end if
end if

! FIRST, FIND NIN ORTHONORMAL VECTORS BY SCALED SYMMETRIC ON.
! Addition, for the scaled symmetric ON: the S matrix is scaled
! to make the diagonal elements close to 1.
! Extremely small values give scale factor exactly zero.
call mma_allocate(SCA,NAS,Label='SCA')
IDIAG = 0
do I=1,NAS
  IDIAG = IDIAG+I
  SDiag = S(IDIAG)
  if (IFDORTHO) then
    SCA(I) = One
  else
    if (SDiag > THRSHN) then
      ! Small variations of the scale factor were beneficial
      SCA(I) = (One+real(I,kind=wp)*3.0e-6_wp)/sqrt(SDiag)
    else
      SCA(I) = Zero
    end if
  end if
end do
IJ = 0
do J=1,NAS
  do I=1,J
    IJ = IJ+1
    S(IJ) = S(IJ)*SCA(I)*SCA(J)
  end do
end do
! End of addition.
if (IPRGLB >= INSANE) then
  FP = DNRM2_(NS,S,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'SMAT NORM AFTER SCALING: ',FP
end if

! DIAGONALIZE THE SCALED S MATRIX:
call mma_allocate(VEC,NAS,NAS,Label='VEC')
call mma_allocate(EIG,NAS,Label='EIG')

call TIMING(CPU1,CPUE,TIO,TIOE)
IJ = 0
do J=1,NAS
  do I=1,J
    IJ = IJ+1
    VEC(I,J) = S(IJ)
  end do
end do
INFO = 0
call dsyev_('V','L',NAS,VEC,NAS,EIG,WGRONK,-1,INFO)
NSCRATCH = int(WGRONK(1))
call mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
call dsyev_('V','U',NAS,VEC,NAS,EIG,SCRATCH,NSCRATCH,INFO)
call mma_deallocate(SCRATCH)
call mma_deallocate(S)

call TIMING(CPU2,CPUE,TIO,TIOE)
CPU = CPU+CPU2-CPU1
! fingerprint eigenvalues
if (iprglb >= insane) then
  fp = dnrm2_(nas,eig,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'Smat eigval norm: ',fp
end if

! Form orthonormal vectors by scaling eigenvectors
NIN = 0
do I=1,NAS
  if (EIG(I) > THRSHS) then
    NIN = NIN+1
    VEC(:,NIN) = VEC(:,I)/sqrt(EIG(I))
  end if
end do
NINDEP(ISYM,ICASE) = NIN
call mma_deallocate(EIG)
! Addition, for the scaled symmetric ON.
do I=1,NIN
  VEC(:,I) = SCA(:)*VEC(:,I)
end do

call mma_deallocate(SCA)
! The condition number, after scaling, disregarding linear dep.
if (NIN >= 2) then
  SZMIN = 1.0e99_wp
  SZMAX = Zero
  do I=1,NIN
    SZ = DNRM2_(NAS,VEC(:,I),1)
    SZMIN = min(SZMIN,SZ)
    SZMAX = max(SZMAX,SZ)
  end do
  CONDNR = (SZMAX/SZMIN)**2
end if
! End of addition.
if (NIN == 0) then
  call mma_deallocate(VEC)
  return
end if

if (BMATRIX == 'NO') then
  ! In some calculations, we do not use B matrices.
  ! Just write the transformation matrix and branch out:
  IDT = IDTMAT(ISYM,ICASE)
  call DDAFILE(LUSBT,1,VEC,NAS*NIN,IDT)
  call mma_deallocate(VEC)
  if (IPRGLB >= INSANE) write(u6,'("DEBUG> ",A)') 'SBDIAG: skip B matrix'
  return
else if (BTRANS /= 'YES') then
  ! In other calculations, B matrix is used, but not transformed.
  ! We may need the diagonal active energies: the diagonal values of
  ! B divided by the diagonal values of S. These are placed where
  ! the eigenvalues would go in ordinary CASPT2.
  ! NOTE: On LUSBT, the transformation matrices partly overwrite
  ! and destroy the B matrices. The diagonal elements of B must be
  ! extracted before the transformation matrix is written.
  call mma_allocate(BD,NAS,Label='BD')
  NB = (NAS*(NAS+1))/2
  call mma_allocate(B,NB,Label='B')
  IDB = IDBMAT(ISYM,ICASE)
  call DDAFILE(LUSBT,2,B,NB,IDB)
  IDIAG = 0
  do I=1,NAS
    IDIAG = IDIAG+I
    BD(I) = B(IDIAG)
  end do
  call mma_deallocate(B)
  ! Now, the transformation matrix can be written out.
  IDT = IDTMAT(ISYM,ICASE)
  call DDAFILE(LUSBT,1,VEC,NAS*NIN,IDT)
  call mma_deallocate(VEC)
  do I=1,NAS
    SDiag = SD(I)+1.0e-15_wp
    BD(I) = BD(I)/SDiag
  end do
  IDB = IDBMAT(ISYM,ICASE)
  call DDAFILE(LUSBT,1,BD,NAS,IDB)
  call mma_deallocate(SD)
  call mma_deallocate(BD)
  if (IPRGLB >= INSANE) then
    write(u6,'("DEBUG> ",A)') 'SBDIAG: skip B matrix transformation'
    write(u6,'("DEBUG> ",A)') '        but keep B_ii/S_ii values'
  end if
  return
end if

! TRANSFORM B MATRIX TO O-N BASIS. BUT FIRST, SAVE O-N VECTORS.
! USE LUSOLV AS TEMPORARY STORAGE. WE MAY NEED SECTIONING.
! NOTE: SECTIONING MUST BE  PRECISELY THE SAME AS WHEN LATER
! READ BACK (SEE BELOW).
NAUX = min(19,NIN)
IDTMP = IDTMP0
call DDAFILE(LUSOLV,1,VEC,NAS*NAUX,IDTMP)
do KSTA=NAUX+1,NIN,NAUX
  KEND = min(KSTA-1+NAUX,NIN)
  NCOL = 1+KEND-KSTA
  call DDAFILE(LUSOLV,1,VEC(:,KSTA:KEND),NAS*NCOL,IDTMP)
end do
if (IPRGLB >= INSANE) then
  FP = DNRM2_(NAS**2,VEC,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'EIGENVECTOR NORM BEFORE B TRANS: ',FP
end if

! TRANSFORM B. NEW B WILL OVERWRITE AND DESTROY VEC
IDB = IDBMAT(ISYM,ICASE)
NB = NS
call mma_allocate(B,NB,Label='B')
call DDAFILE(LUSBT,2,B,NB,IDB)
if ((do_grad .or. (nStpGrd == 2)) .and. do_lindep) then
  !! The original B matrix is needed in the LinDepLag subroutine
  IDB2 = idBoriMat(ISYM,ICASE)
  call DDAFILE(LUSTD,1,B,NB,IDB2)
end if
if (IPRGLB >= INSANE) then
  FP = DNRM2_(NB,B,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'BMAT NORM: ',FP
end if

call mma_allocate(BX,NAS,Label='BX')
call mma_allocate(XBX,NAS,Label='XBX')
do J=NIN,1,-1
  BX(:) = Zero
# ifdef _CRAY_C90_
  call SSPMV('U',NAS,One,B,VEC(:,J),1,One,BX,1)
# else
  !call DSLMX(NAS,One,B,VEC(:,J),1,BX,1)
  call DSPMV_('U',NAS,One,B,VEC(:,J),1,One,BX,1)
# endif
  ! BX: B * Vector number J.
  XBX(1:J) = Zero
  call DGEMM_('T','N',J,1,NAS,One,VEC,NAS,BX,NAS,Zero,XBX,J)
  ! XBX CONTAINS NOW THE UPPERTRIANGULAR
  ! ELEMENTS OF THE J-th COLUMN OF TRANSFORMED B MATRIX.
  VEC(1:J,J) = XBX(1:J)
end do
call mma_deallocate(BX)
call mma_deallocate(XBX)
call mma_deallocate(B)
! VEC HAS NOW BEEN DESTROYED (OVERWRITTEN BY NEW B).
! COPY TO TRIANGULAR STORAGE.
NBNEW = (NIN*(NIN+1))/2
call mma_allocate(B,NBNEW,Label='B')
do J=1,NIN
  JOFF = (J*(J-1))/2
  B(JOFF+1:JOFF+J) = VEC(1:J,J)
end do
call mma_deallocate(VEC)
if (IPRGLB >= INSANE) then
  FP = DNRM2_(NBNEW,B,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'BMAT NORM AFTER TRANS: ',FP
end if

! DIAGONALIZE THE TRANSFORMED B MATRIX.
call mma_allocate(EIG,NIN,Label='EIG')
call mma_allocate(VEC,NIN,NIN,Label='VEC')
call TIMING(CPU1,CPUE,TIO,TIOE)
! - Alt 0: Use diagonal approxim., if allowed:
if (BSPECT /= 'YES') then
  IDIAG = 0
  write(u6,*) 'Sbdiag: THis code does not make sense!'
  write(u6,*) '        SDiag is not properly defined!'
  call Abend()
  do I=1,NIN
    IDIAG = IDIAG+I
    EIG(I) = B(IDIAG)/SDiag
  end do
else
  IJ = 0
  do J=1,NIN
    do I=1,J
      IJ = IJ+1
      VEC(I,J) = B(IJ)
    end do
  end do
  call DSYEV_('V','U',NIN,VEC,NIN,EIG,WGRONK,-1,INFO)
  NSCRATCH = int(WGRONK(1))
  call mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
  call DSYEV_('V','U',NIN,VEC,NIN,EIG,SCRATCH,NSCRATCH,INFO)
  call mma_deallocate(SCRATCH)
  call mma_deallocate(B)
end if
call TIMING(CPU2,CPUE,TIO,TIOE)
CPU = CPU+CPU2-CPU1
if (IPRGLB >= INSANE) then
  FP = DNRM2_(NIN,EIG,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'BMAT EIGENVALUE NORM: ',FP
end if

! The eigenvalues are written back at same position as the
! original B matrix, which is destroyed:
IDB = IDBMAT(ISYM,ICASE)
call DDAFILE(LUSBT,1,EIG,NIN,IDB)
call mma_deallocate(EIG)

! Finally, we must form the composite transformation matrix,
!  = (NAS*NIN matrix on disk) * (NIN*NIN matrix in core).
! Assume enough space since we got rid of S/B matrices.
! Specifically, assume we have enough space for the two
! full matrices, plus an additional 19 columns of results.
NAUX = min(19,NIN)
call mma_allocate(TRANS,NAS*NIN,Label='TRANS')
call mma_allocate(AUX,NAS*NAUX,Label='AUX')
IDTMP = IDTMP0
call DDAFILE(LUSOLV,2,AUX,NAS*NAUX,IDTMP)
if (BTRANS == 'YES') then
  call DGEMM_('N','N',NAS,NIN,NAUX,One,AUX,NAS,VEC,NIN,Zero,TRANS,NAS)
else
  TRANS(1:NAS*NAUX) = AUX(:)
end if
do KSTA=NAUX+1,NIN,NAUX
  KEND = min(KSTA-1+NAUX,NIN)
  NCOL = 1+KEND-KSTA
  call DDAFILE(LUSOLV,2,AUX,NAS*NCOL,IDTMP)
  if (BTRANS == 'YES') then
    call DGEMM_('N','N',NAS,NIN,NCOL,One,AUX,NAS,VEC(KSTA,1),NIN,One,TRANS,NAS)
  else
    LTRANS1 = NAS*(KSTA-1)
    TRANS(LTRANS1+1:LTRANS1+NAS*NCOL) = AUX(1:NAS*NCOL)
  end if
end do
call mma_deallocate(AUX)
call mma_deallocate(VEC)
IDT = IDTMAT(ISYM,ICASE)
call DDAFILE(LUSBT,1,TRANS,NAS*NIN,IDT)
if (IPRGLB >= INSANE) then
  FP = DNRM2_(NAS*NIN,TRANS,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'TMAT NORM: ',FP
end if

!-SVC: compute S*T and store on disk for later use by RHS vector utilities.
NS = (NAS*(NAS+1))/2
call mma_allocate(S,NS,Label='S')
IDS = IDSMAT(ISYM,ICASE)
call DDAFILE(LUSBT,2,S,NS,IDS)
call mma_allocate(ST,NAS*NIN,Label='ST')
ST(:) = Zero
call TRIMUL(NAS,NIN,One,S,TRANS,NAS,ST,NAS)
call mma_deallocate(S)
call mma_deallocate(TRANS)
IDST = IDSTMAT(ISYM,ICASE)
call DDAFILE(LUSBT,1,ST,NAS*NIN,IDST)
if (IPRGLB >= INSANE) then
  FP = DNRM2_(NAS*NIN,ST,1)
  write(u6,'("DEBUG> ",A,ES21.14)') 'STMAT NORM: ',FP
end if
call mma_deallocate(ST)

end subroutine SBDIAG_SER
