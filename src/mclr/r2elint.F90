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

subroutine r2elint(rKappa,rMO1,rmo2,FockI,FockA,iDSym,Sgn,Fact,jspin)
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

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: CMO, FAMO, FIMO, G1t, ipCM, ipMat, nA, nCMO, nDens, nMBA
use input_mclr, only: CasInt, iMethod, nAsh, nBas, nIsh, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: rKappa(nDens), Sgn, Fact
real(kind=wp), intent(inout) :: rMO1(nMba), rmo2(*)
real(kind=wp), intent(out) :: FockI(nDens), FockA(nDens)
integer(kind=iwp), intent(in) :: iDSym, jSpin
integer(kind=iwp) :: i, iA, iB, ip, ip2, iS, jA, jB, jS, nDens22
real(kind=wp) :: FacR
logical(kind=iwp) :: lFA, lFI, lMo
real(kind=wp), allocatable :: DA(:), DI(:), FA(:), FI(:), T1(:), T3(:), T4(:), Tmp2(:)

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
DI(:) = Zero
lFI = .true.
lFa = .false.
lMo = .false.
if (iMethod == 2) then
  call mma_allocate(DA,nCMO,Label='DA')
  call mma_allocate(FA,nDens,Label='FA')
  lFa = .true.
  lMo = .true.
else
  call mma_allocate(DA,1,Label='DA')
  call mma_allocate(FA,1,Label='FA')
end if
FA(:) = Zero
DA(:) = Zero

do iS=1,nSym
  do iB=1,nIsh(iS)
    ip = ipCM(iS)+(iB-1)*nOrb(iS)+iB-1
    DI(ip) = Two
  end do
end do
if (iMethod == 2) then
  do iS=1,nSym
    do iB=1,nAsh(iS)
      do jB=1,nAsh(iS)
        ip = ipCM(iS)+iB+nIsh(iS)+(jB+nIsh(iS)-1)*nOrb(iS)-1
        iA = nA(iS)+iB
        jA = nA(iS)+jB
        ip2 = iTri(iA,jA)
        DA(ip) = G1t(ip2)
      end do
    end do
  end do
end if

! Construct {kappa,(ij|kl)} and the fock matrix contribution
! from one index tranformation of contracted indexes

FacR = Fact
call Read2_2(rmo1,rmo2,FockI,FockA,T1,Tmp2,T3,T4,nDens22,DI,DA,rkappa,idsym,Sgn,Facr,jSpin,lFA,lfi,lMo)

!call recprt('1 FockI','',FockI,nDens,1)
!call recprt('1 FockA','',FockA,nDens,1)

! Calculate contribution from uncontracted indexes.

do iS=1,nSym
  jS = Mul(iS,iDSym)
  if (nOrb(iS)*nOrb(jS) /= 0) then
    if (.not. CASINT) &
      call DGEMM_('T','N',nOrb(iS),nOrb(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),FI(ipMat(iS,jS)),nBas(iS),Zero, &
                  FockI(ipMat(iS,jS)),nOrb(iS))
    call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Sgn*Facr,FIMO(ipCM(iS)),nOrb(is),rkappa(ipMat(iS,jS)),nOrb(iS),One, &
                FockI(ipMat(iS,jS)),nOrb(iS))
    call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Facr,rkappa(ipMat(iS,jS)),nOrb(is),FIMO(ipCM(jS)),nOrb(jS),One, &
                FockI(ipMat(iS,jS)),nOrb(is))
    if (iMethod == 2) then
      if (.not. CASINT) &
        call DGEMM_('T','N',nOrb(iS),nOrb(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),FA(ipMat(iS,jS)),nBas(iS),Zero, &
                    FockA(ipMat(iS,jS)),nOrb(iS))
      call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Sgn*Facr,FAMO(ipCM(iS)),nOrb(is),rkappa(ipMat(iS,jS)),nOrb(iS),One, &
                  FockA(ipMat(iS,jS)),nOrb(iS))
      call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Facr,rkappa(ipMat(iS,jS)),nOrb(is),FAMO(ipCM(jS)),nOrb(jS),One, &
                  FockA(ipMat(iS,jS)),nOrb(is))
    end if
  end if
end do

!call recprt('2 FockI',' ',FockI,nDens,1)
!call recprt('2 FockA',' ',FockA,nDens,1)

call mma_deallocate(DA)
call mma_deallocate(FA)
call mma_deallocate(FI)
call mma_deallocate(DI)
call mma_deallocate(T4)
call mma_deallocate(T3)
call mma_deallocate(Tmp2)
call mma_deallocate(T1)

end subroutine r2elint
