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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine r2elint_sp(rKappa,rMO1,rmo2,FockI,FockA,iDSym,Sgn,Fact,jspin,D,FA)
!***********************************************************************
!
! Constructs the one index transformed Fock-matrixes
! and (pj|kl).
! rKappa : the transformation matrix
! iDSym  : Symmetry of transformation
! Sgn    : 1:real -1:complex
! jspin  : 0:singlet 1:triplet
!
!***********************************************************************

use Symmetry_Info, only: Mul
use MCLR_Data, only: CMO, FIMO, ipCM, ipMat, nA, nCMO, nDens, nMBA, nNA
use input_mclr, only: iMethod, nAsh, nBas, nIsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: rKappa(nDens), Sgn, Fact, D(*), FA(*)
real(kind=wp), intent(inout) :: rMO1(nMba), rmo2(*)
real(kind=wp), intent(out) :: FockI(nDens), FockA(nDens)
integer(kind=iwp), intent(in) :: iDSym, jSpin
integer(kind=iwp) :: i, iA, iB, ip, iS, jA, jB, jS, nDens22
real(kind=wp) :: FacR
logical(kind=iwp) :: lFA, lFI, lMo
real(kind=wp), allocatable :: DA(:), DI(:), FA2(:), FI(:), T1(:), T3(:), T4(:), Tmp2(:)

nDens22 = nDens
do i=1,nSym
  nDens22 = max(nDens22,maxval(nBas(i)*nBas(1:nSym)))
end do
call mma_allocate(T1,nDens22,Label='T1')
call mma_allocate(Tmp2,nDens22,Label='Tmp2')
call mma_allocate(T3,nDens22,Label='T3')
call mma_allocate(T4,nDens22,Label='T4')
call mma_allocate(DI,nCMO,Label='DI')
call mma_allocate(FI,nDens,Label='FI')

FockI(:) = Zero
FockA(:) = Zero
FI(:) = Zero

lFI = .true.
lFa = .false.
lMo = .false.
if (iMethod == 2) then
  call mma_allocate(DA,nCMO,Label='DA')
  call mma_allocate(FA2,nDens,Label='FA2')
  lFa = .true.
  lMo = .true.
else
  call mma_allocate(DA,1,Label='DA')
  call mma_allocate(FA2,1,Label='FA2')
end if

do iS=1,nSym
  do iB=1,nIsh(iS)
    ip = ipCM(iS)+(iB-1)*nBas(iS)+iB-1
    DI(ip) = Two
  end do
end do

if (iMethod == 2) then
  do iS=1,nSym
    do iB=1,nAsh(iS)
      do jB=1,nAsh(iS)
        ip = ipCM(iS)+iB+nIsh(iS)+(jB+nIsh(iS)-1)*nBas(iS)-1
        iA = nA(iS)+iB
        jA = nA(iS)+jB
        DA(ip) = D(iA+(jA-1)*nNA)
      end do
    end do
  end do
end if

! Construct {kappa,(ij|kl)} and the fock matrix contribution
! from one index tranformation of contracted indexes

FacR = Fact
call Read2_2(rmo1,rmo2,FI,FA2,T1,Tmp2,T3,T4,nDens22,DI,DA,rkappa,idsym,Sgn,Facr,jSpin,lFA,lfi,lMo)

! Calculate contribution from uncontracted indexes.

do iS=1,nSym
  jS = Mul(iS,iDSym)
  if (nBas(iS)*nBas(jS) /= 0) then
    call DGEMM_('T','N',nBas(iS),nBas(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),FI(ipMat(iS,jS)),nBas(iS),Zero,FockI(ipMat(iS,jS)), &
                nBas(iS))
    call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sgn*Facr,FIMO(ipCM(iS)),nBas(is),rkappa(ipMat(iS,jS)),nBas(iS),One, &
                FockI(ipMat(iS,jS)),nBas(iS))
    call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,rkappa(ipMat(iS,jS)),nBas(is),FIMO(ipCM(jS)),nBas(jS),One, &
                FockI(ipMat(iS,jS)),nBas(is))
    if (iMethod == 2) then
      call DGEMM_('T','N',nBas(iS),nBas(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),FA2(ipMat(iS,jS)),nBas(iS),Zero, &
                  FockA(ipMat(iS,jS)),nBas(iS))
      call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sgn*Facr,FA(ipCM(iS)),nBas(is),rkappa(ipMat(iS,jS)),nBas(iS),One, &
                  FockA(ipMat(iS,jS)),nBas(iS))
      call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,rkappa(ipMat(iS,jS)),nBas(is),FA(ipCM(jS)),nBas(jS),One, &
                  FockA(ipMat(iS,jS)),nBas(is))
    end if
  end if
end do

call mma_deallocate(FA2)
call mma_deallocate(DA)
call mma_deallocate(FI)
call mma_deallocate(DI)
call mma_deallocate(T4)
call mma_deallocate(T3)
call mma_deallocate(Tmp2)
call mma_deallocate(T1)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine r2elint_sp
