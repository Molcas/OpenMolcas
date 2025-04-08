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

subroutine r2elint(rKappa,rMO1,rmo2,FockI,FockA,iDSym,sign,Fact,jspin)
!***********************************************************************
!
! Constructs the one index transformed Fock-matrixes
! and (pj|kl).
! rKappa : the transformation matrix
! iDSym  : Symmetry of transformation
! sign   : 1:real -1:complex
! jspin  : 0:singlet 1:triplet
!
!***********************************************************************

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: CMO, G1t, FAMO, FIMO
use MCLR_Data, only: nDens, nMBA, ipCM, ipMat, nA, nCMO
use input_mclr, only: nSym, nAsh, nIsh, nBas, nOrb, iMethod, CasInt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two

implicit none
real*8 rKappa(nDens), rMO1(nMba), rmo2(*), FockI(nDens), FockA(nDens)
integer iDSym, jSpin
real*8 sign, Fact
logical lFI, lFA, lMo
real*8, allocatable :: T1(:), Tmp2(:), T3(:), T4(:), DIL(:), DI(:), DIR(:), FI(:), DAL(:), DAR(:), DA(:), FA(:)
integer nDens22, iAM, iBM, iMem, iS, iB, ip, jB, iA, ip2, jS, jA
real*8 FacR
integer i, j

nDens22 = nDens
iAM = 0
iBM = 0
do i=1,nSym
  do j=1,nSym
    nDens22 = max(nDens22,nBas(i)*nBas(j))
  end do
  iAM = max(iAM,nAsh(i)+nIsh(i))
  iBM = max(iBM,nBAs(i))
end do
imem = nDens22

call mma_allocate(T1,imem,Label='T1')
call mma_allocate(Tmp2,nDens22,Label='Tmp2')
call mma_allocate(T3,nDens22,Label='T3')
call mma_allocate(T4,nDens22,Label='T4')
call mma_allocate(DIL,nDens,Label='DIL')
call mma_allocate(DI,nCMO,Label='DI')
call mma_allocate(DIR,nDens,Label='DIR')
call mma_allocate(FI,nDens,Label='FI')

FockI(:) = Zero
FockA(:) = Zero
FI(:) = Zero
DI(:) = Zero
DIL(:) = Zero
DIR(:) = Zero
lFI = .true.
lFa = .false.
lMo = .false.
if (iMethod == 2) then
  call mma_allocate(DAL,nDens,Label='DAL')
  call mma_allocate(DAR,nDens,Label='DAR')
  call mma_allocate(DA,nCMO,Label='DA')
  call mma_allocate(FA,nDens,Label='FA')
  lFa = .true.
  lMo = .true.
else
  call mma_allocate(DAL,1,Label='DAL')
  call mma_allocate(DAR,1,Label='DAR')
  call mma_allocate(DA,1,Label='DA')
  call mma_allocate(FA,1,Label='FA')
end if
FA(:) = Zero
DA(:) = Zero
DAL(:) = Zero
DAR(:) = Zero

do iS=1,nSym
  do iB=1,nIsh(iS)
    ip = ipCM(iS)+(ib-1)*nOrb(is)+ib-1
    DI(ip) = Two
  end do
end do
if (iMethod == 2) then
  do iS=1,nSym
    do iB=1,nAsh(iS)
      do jB=1,nAsh(iS)
        ip = ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nOrb(is)-1
        iA = nA(is)+ib
        jA = nA(is)+jb
        ip2 = iTri(iA,jA)
        DA(ip) = G1t(ip2)
      end do
    end do
  end do
end if

! Construct {kappa,(ij|kl)} and the fock matrix contribution
! from one index tranformation of contracted indexes

FacR = Fact
call Read2_2(rmo1,rmo2,FockI,FockA,T1,imem,Tmp2,T3,T4,nDens22,DIR,DIL,DI,DAR,DAL,DA,rkappa,idsym,Sign,Facr,jSpin,lFA,lfi,lMo)

!call recprt('1 FockI','',FockI,nDens,1)
!call recprt('1 FockA','',FockA,nDens,1)

! Calculate contribution from uncontracted indexes.

do iS=1,nSym
  jS = Mul(iS,iDSym)
  if (nOrb(iS)*nOrb(jS) /= 0) then
    if (.not. CASINT) &
      call DGEMM_('T','N',nOrb(iS),nOrb(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),FI(ipMat(iS,jS)),nBas(iS),Zero, &
                  FockI(ipMat(iS,jS)),nOrb(iS))
    call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Sign*Facr,FIMO(ipCM(iS)),nOrb(is),rkappa(ipMat(iS,jS)),nOrb(iS),One, &
                FockI(ipMat(iS,jS)),nOrb(iS))
    call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Facr,rkappa(ipMat(iS,jS)),nOrb(is),FIMO(ipCM(jS)),nOrb(jS),One, &
                FockI(ipMat(iS,jS)),nOrb(is))
    if (iMethod == 2) then
      if (.not. CASINT) &
        call DGEMM_('T','N',nOrb(iS),nOrb(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),FA(ipMat(iS,jS)),nBas(iS),Zero, &
                    FockA(ipMat(iS,jS)),nOrb(iS))
      call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Sign*Facr,FAMO(ipCM(iS)),nOrb(is),rkappa(ipMat(iS,jS)),nOrb(iS),One, &
                  FockA(ipMat(iS,jS)),nOrb(iS))
      call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Facr,rkappa(ipMat(iS,jS)),nOrb(is),FAMO(ipCM(jS)),nOrb(jS),One, &
                  FockA(ipMat(iS,jS)),nOrb(is))
    end if
  end if
end do

!call recprt('2 FockI',' ',FockI,nDens,1)
!call recprt('2 FockA',' ',FockA,nDens,1)

call mma_deallocate(DA)
call mma_deallocate(DAR)
call mma_deallocate(DAL)
call mma_deallocate(FA)
call mma_deallocate(FI)
call mma_deallocate(DIR)
call mma_deallocate(DI)
call mma_deallocate(DIL)
call mma_deallocate(T4)
call mma_deallocate(T3)
call mma_deallocate(Tmp2)
call mma_deallocate(T1)

end subroutine r2elint
