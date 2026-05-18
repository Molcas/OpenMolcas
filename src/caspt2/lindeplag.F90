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

subroutine LinDepLag(BDer,SDer,nAS,nIN,iSym,iCase)
! Compute contributions that arise from the non-invariance effect
! in non-orthogonal -> orthogonal ICB rotations
! See J. Chem. Phys. 2023, 158, 174112. for details, in particular,
! Section II C 3 "Non-invariance with respect to orthogonal..."

use EQSOLV, only: idSMAT
use caspt2_global, only: idBoriMat, LUSBT, LUSTD
use caspt2_module, only: IFDORTHO, THRSHN, THRSHS
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAS, nIN, iSym, iCase
real(kind=wp), intent(inout) :: BDer(nAS,nAS), SDer(nAS,nAS)
integer(kind=iwp) :: I, IDB, IDIAG, idS, IJ, INFO, J, NB, NS, NSCRATCH
real(kind=wp) :: EVAL, FACT, SCAL, SD, WGRONK(2)
real(kind=wp), allocatable :: B(:), EIG(:), F(:,:), LAG(:,:), S(:), SCA(:), SCRATCH(:), SS(:,:), VEC(:,:)

!! Obtain the X matrix
!! First, read S
NS = NAS*(NAS+1)/2
call mma_allocate(S,NS,Label='S')
call mma_allocate(SS,NAS,NAS,Label='SS')
idS = idSMAT(iSym,iCase)
call DDAFILE(LUSBT,2,S,NS,idS)
IJ = 0
do J=1,NAS
  do I=1,J
    IJ = IJ+1
    SS(I,J) = S(IJ)
    SS(J,I) = S(IJ)
  end do
end do

call mma_allocate(VEC,NAS,NAS,Label='VEC')
call mma_allocate(EIG,NAS,Label='EIG')
call mma_allocate(SCA,NAS,Label='SCA')
IDIAG = 0
do I=1,NAS
  IDIAG = IDIAG+I
  SD = S(IDIAG)
  if (IFDORTHO) then
    SCA(I) = One
  else
    if (SD > THRSHN) then
      ! Small variations of the scale factor were beneficial
      SCA(I) = (One+real(I,kind=wp)*3.0e-6_wp)/sqrt(SD)
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

do I=1,NAS
  SCAL = SCA(I)
  call DSCAL_(NAS,SCAL,VEC(I,1),NAS)
end do
call mma_deallocate(SCA)
call mma_deallocate(S)

!! Scale only the independent vectors to avoid
!! any numerically unstable computation
do I=1,NAS
  EVAL = EIG(I)
  if (EVAL < THRSHS) cycle
  FACT = One/sqrt(EVAL)
  call DScal_(nAS,FACT,VEC(1,I),1)
end do

call mma_allocate(LAG,NAS,NAS,Label='LAG')
IDB = IDBoriMat(ISYM,ICASE)
NB = NS
call mma_allocate(B,NB,Label='B')
call DDAFILE(LUSTD,2,B,NB,IDB)
call mma_allocate(F,NAS,NAS,Label='F')
IJ = 0
do J=1,NAS
  do I=1,J
    IJ = IJ+1
    F(I,J) = B(IJ)
    F(J,I) = B(IJ)
  end do
end do

!! Compute the partial derivative
!! F   : B
!! BDER: D
!! VEC : X^0 and X
call DGEMM_('N','T',NAS,NAS,NAS,Two,F,NAS,BDER,NAS,Zero,LAG,NAS)
call DGEMM_('N','N',NAS,NAS,NAS,One,LAG,NAS,VEC,NAS,Zero,F,NAS)
LAG(1:NAS,1:NAS) = F(1:NAS,1:NAS)

call DGEMM_('T','N',NAS,NAS,NAS,One,VEC,NAS,LAG,NAS,Zero,F,NAS)
!! At this point,
!! F = 2 \mathcal{X}^0 * B * D * \mathcal{X}

!! remove dependent part
!! (linearly indep-indep and dep-dep)
F(1:nAS-nIN,1:nAS-nIN) = Zero
F(nAS-nIN+1:nAS,nAS-nIN+1:nAS) = Zero

!! orthogonal -> non-orthogonal
!! Finalize Eq. (62)
call DGEMM_('N','N',NAS,NAS,NAS,One,VEC,NAS,F,NAS,Zero,LAG,NAS)
call DGEMM_('N','T',NAS,NAS,NAS,One,LAG,NAS,VEC,NAS,Zero,F,NAS)

call DaXpY_(nAS*nAS,One,F,1,SDER,1)

call mma_deallocate(LAG)
call mma_deallocate(B)
call mma_deallocate(F)

call mma_deallocate(SS)
call mma_deallocate(EIG)
call mma_deallocate(VEC)

return

end subroutine LinDepLag
