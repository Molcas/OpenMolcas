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

subroutine LinDepLag(BDer,SDer,GDer,nAS,nIN,iSym,iCase)
! Compute contributions that arise from the non-invariance effect
! in non-orthogonal -> orthogonal ICB rotations
! See J. Chem. Phys. 2023, 158, 174112. for details, in particular,
! Section II C 3 "Non-invariance with respect to orthogonal..."
! Both the canonical and the scaled orthonormalization are now supported, as in SBDIAG.

use Index_Functions, only: nTri_Elem
use EQSOLV, only: idSMAT
use caspt2_global, only: idBoriMat, LUSBT, LUSTD
use caspt2_module, only: HZERO, IFDORTHO, THRSHN
use SC_NEVPT2, only: IDBMAT_NEVPT2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAS, nIN, iSym, iCase
real(kind=wp), intent(inout) :: BDer(nAS,nAS), SDer(nAS,nAS)
real(kind=wp), intent(in) :: GDer(nAS,nAS)
integer(kind=iwp) :: I, IDB, idS, IJ, INFO, J, NS, NSCRATCH
real(kind=wp) :: WGRONK(2)
real(kind=wp), allocatable :: B(:), EIG(:), F(:,:), LAG(:,:), S(:), SCA(:), SCRATCH(:), SS(:,:), VEC(:,:)

! Diagonalize S scaled in the same way as in SBDIAG
NS = nTri_Elem(NAS)
call mma_allocate(S,NS,Label='S')
idS = idSMAT(iSym,iCase)
call DDAFILE(LUSBT,2,S,NS,idS)

call mma_allocate(VEC,NAS,NAS,Label='VEC')
call mma_allocate(EIG,NAS,Label='EIG')
call mma_allocate(SCA,NAS,Label='SCA')
do I=1,NAS
  if (IFDORTHO) then
    SCA(I) = One
  else if (S(nTri_Elem(I)) > THRSHN) then
    SCA(I) = (One+real(I,kind=wp)*3.0e-6_wp)/sqrt(S(nTri_Elem(I)))
  else
    SCA(I) = Zero
  end if
end do
IJ = 0
do J=1,NAS
  VEC(1:J,J) = S(IJ+1:IJ+J)*SCA(1:J)*SCA(J)
  IJ = IJ+J
end do
call mma_deallocate(S)
INFO = 0
call dsyev_('V','L',NAS,VEC,NAS,EIG,WGRONK,-1,INFO)
NSCRATCH = int(WGRONK(1))
call mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
call dsyev_('V','U',NAS,VEC,NAS,EIG,SCRATCH,NSCRATCH,INFO)
call mma_deallocate(SCRATCH)

! The first nAS-nIN eigenvectors are the removed (dependent) ones.
! VEC = SCA^{-1} * U couples with the densities.
do I=1,NAS
  if (SCA(I) /= Zero) VEC(I,:) = VEC(I,:)/SCA(I)
end do

call mma_allocate(SS,NAS,NAS,Label='SS')
call mma_allocate(LAG,NAS,NAS,Label='LAG')
call mma_allocate(F,NAS,NAS,Label='F')
! SS = VEC^T * G^T * VEC. Its dep-indep block times the dep eigenvalue is the RHS part of the residual.
call DGEMM_('T','N',NAS,NAS,NAS,One,GDER,NAS,VEC,NAS,Zero,LAG,NAS)
call DGEMM_('T','N',NAS,NAS,NAS,One,VEC,NAS,LAG,NAS,Zero,SS,NAS)

call mma_allocate(B,NS,Label='B')
if ((HZERO == 'DYALL') .and. ((iCase == 1) .or. (iCase == 4))) then
  IDB = IDBMAT_NEVPT2(ISYM,ICASE,1)
  call DDAFILE(LUSBT,2,B,NS,IDB)
else
  IDB = IDBoriMat(ISYM,ICASE)
  call DDAFILE(LUSTD,2,B,NS,IDB)
end if
call SQUARE(B,F,1,NAS,NAS)
call mma_deallocate(B)

! LAG = 2 (SCA*U)^T * B * D * (SCA^{-1}*U), the B part of the residual
call DGEMM_('N','T',NAS,NAS,NAS,Two,F,NAS,BDER,NAS,Zero,LAG,NAS)
call DGEMM_('N','N',NAS,NAS,NAS,One,LAG,NAS,VEC,NAS,Zero,F,NAS)
! VEC = SCA * U couples with the residual and the derivative of S
do I=1,NAS
  VEC(I,:) = VEC(I,:)*SCA(I)**2
end do
call DGEMM_('T','N',NAS,NAS,NAS,One,VEC,NAS,F,NAS,Zero,LAG,NAS)

! Rotations between the dependent (I) and independent (J) vectors
F(:,:) = Zero
do J=nAS-nIN+1,nAS
  do I=1,nAS-nIN
    F(I,J) = (LAG(I,J)+EIG(I)*SS(I,J))/(EIG(J)-EIG(I))
  end do
end do

!! orthogonal -> non-orthogonal
!! Finalize Eq. (62)
call DGEMM_('N','N',NAS,NAS,NAS,One,VEC,NAS,F,NAS,Zero,LAG,NAS)
call DGEMM_('N','T',NAS,NAS,NAS,One,LAG,NAS,VEC,NAS,One,SDER,NAS)

! Derivative of the scale factors, which depend on the diagonal elements of S
if (.not. IFDORTHO) then
  do J=nAS-nIN+1,nAS
    F(1:nAS-nIN,J) = F(1:nAS-nIN,J)*EIG(J)
  end do
  call DGEMM_('N','T',NAS,NAS,NAS,One,F,NAS,VEC,NAS,Zero,LAG,NAS)
  do I=1,NAS
    if (SCA(I) /= Zero) SDER(I,I) = SDER(I,I)-sum(VEC(I,1:nAS-nIN)*LAG(1:nAS-nIN,I))/(One+real(I,kind=wp)*3.0e-6_wp)**2
  end do
end if

call mma_deallocate(SCA)
call mma_deallocate(SS)
call mma_deallocate(LAG)
call mma_deallocate(F)
call mma_deallocate(EIG)
call mma_deallocate(VEC)

end subroutine LinDepLag
