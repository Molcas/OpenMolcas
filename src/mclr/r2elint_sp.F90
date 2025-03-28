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

subroutine r2elint_sp(rKappa,rMO1,rmo2,FockI,FockA,iDSym,sign,Fact,jspin,D,FA)
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

use MCLR_Data, only: CMO, FIMO
use MCLR_Data, only: nDens2, nMBA, nNA, ipCM, ipMat, nA, nCMO
use input_mclr, only: nSym, nAsh, nIsh, nBas, iMethod
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two

implicit none
real*8 rKappa(nDens2), rMO1(nMba), rmo2(*), FockI(nDens2), FockA(nDens2)
integer iDSym, jSpin
real*8 sign, Fact
real*8 D(*), FA(*)
logical lFI, lFA, lMo
real*8, allocatable :: T1(:), Tmp2(:), T3(:), T4(:), DIL(:), DI(:), DIR(:), FI(:), FA2(:), DAR(:), DA(:), DAL(:)
integer nDens22, iAM, iBM, iMem, iS, iB, ip, jB, iA, jA, jS
real*8 FacR
! Statement function
integer i, j, irec
irec(i,j) = i+nna*(j-1)

ndens22 = ndens2
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
call mma_allocate(DIL,nDens2,Label='DIL')
call mma_allocate(DI,nCMO,Label='DI')
call mma_allocate(DIR,nDens2,Label='DIR')
call mma_allocate(FI,ndens2,Label='FI')

FockI(:) = Zero
FockA(:) = Zero
FI(:) = Zero
DIL(:) = Zero
DIR(:) = Zero

lFI = .true.
lFa = .false.
lMo = .false.
if (iMethod == 2) then
  call mma_allocate(DAL,nDens2,Label='DAL')
  call mma_allocate(DAR,nDens2,Label='DAR')
  call mma_allocate(DA,nCMO,Label='DA')
  call mma_allocate(FA2,nDens2,Label='FA2')
  lFa = .true.
  lMo = .true.
else
  call mma_allocate(DAL,1,Label='DAL')
  call mma_allocate(DAR,1,Label='DAR')
  call mma_allocate(DA,1,Label='DA')
  call mma_allocate(FA2,1,Label='FA2')
end if

do iS=1,nSym
  do iB=1,nIsh(iS)
    ip = ipCM(iS)+(ib-1)*nBas(is)+ib-1
    DI(ip) = Two
  end do
end do

if (iMethod == 2) then
  do iS=1,nSym
    do iB=1,nAsh(iS)
      do jB=1,nAsh(iS)
        ip = ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nBas(is)-1
        iA = nA(is)+ib
        jA = nA(is)+jb
        DA(ip) = D(irec(iA,jA))
      end do
    end do
  end do
end if

! Construct {kappa,(ij|kl)} and the fock matrix contribution
! from one index tranformation of contracted indexes

FacR = Fact
call Read2_2(rmo1,rmo2,FI,FA2,T1,imem,Tmp2,T3,T4,nDens22,DIR,DIL,DI,DAR,DAL,DA,rkappa,idsym,Sign,Facr,jSpin,lFA,lfi,lMo)

! Calculate contribution from uncontracted indexes.

do iS=1,nSym
  jS = ieor(iS-1,iDSym-1)+1
  if (nBas(iS)*nBas(jS) /= 0) then
    call DGEMM_('T','N',nBas(iS),nBas(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),FI(ipMat(iS,jS)),nBas(iS),Zero,FockI(ipMat(iS,jS)), &
                nBas(iS))
    call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,FIMO(ipCM(iS)),nBas(is),rkappa(ipMat(iS,jS)),nBas(iS),One, &
                FockI(ipMat(iS,jS)),nBas(iS))
    call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,rkappa(ipMat(iS,jS)),nBas(is),FIMO(ipCM(jS)),nBas(jS),One, &
                FockI(ipMat(iS,jS)),nBas(is))
    if (iMethod == 2) then
      call DGEMM_('T','N',nBas(iS),nBas(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(iS),FA2(ipMat(iS,jS)),nBas(iS),Zero, &
                  FockA(ipMat(iS,jS)),nBas(iS))
      call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,FA(ipCM(iS)),nBas(is),rkappa(ipMat(iS,jS)),nBas(iS),One, &
                  FockA(ipMat(iS,jS)),nBas(iS))
      call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,rkappa(ipMat(iS,jS)),nBas(is),FA(ipCM(jS)),nBas(jS),One, &
                  FockA(ipMat(iS,jS)),nBas(is))
    end if
  end if
end do

call mma_deallocate(FA2)
call mma_deallocate(DA)
call mma_deallocate(DAR)
call mma_deallocate(DAL)
call mma_deallocate(FI)
call mma_deallocate(DIR)
call mma_deallocate(DI)
call mma_deallocate(DIL)
call mma_deallocate(T4)
call mma_deallocate(T3)
call mma_deallocate(Tmp2)
call mma_deallocate(T1)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine r2elint_sp
