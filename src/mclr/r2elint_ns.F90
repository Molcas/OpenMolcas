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

subroutine r2elint_ns(rKappa,rMO1,rmo2,FockI,FockA,iDSym,Sgn,Fact,jspin)
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
use MCLR_Data, only: FAMO, FIMO, G1t, ipCM, ipMat, nA, nCMO, nDens, nMBA
use input_mclr, only: iMethod, nAsh, nBas, nIsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: rKappa(nDens), Sgn, Fact
real(kind=wp), intent(inout) :: rMO1(nMBA), rMO2(nMBA)
real(kind=wp), intent(out) :: FockI(nDens), FockA(nDens)
integer(kind=iwp), intent(in) :: iDSym, jSpin
integer(kind=iwp) :: i, iA, iB, ip, ip2, iS, jA, jB, jS, nDens22
real(kind=wp) :: FacR, rdum(1)
logical(kind=iwp) :: lFA, lFI, lMo
real(kind=wp), allocatable :: DA(:), DAL(:), DAR(:), DI(:), DIL(:), DIR(:), FA1(:), FI(:), FI1(:), K1(:), T1(:), T3(:), T4(:), &
                              Tmp2(:)

!                                                                      *
!***********************************************************************
!                                                                      *

nDens22 = nDens
do i=1,nSym
  nDens22 = max(nDens22,maxval(nBas(i)*nBas(1:nSym)))
end do

call mma_allocate(T1,nDens22,Label='T1')
call mma_allocate(Tmp2,nDens22,Label='Tmp2')
call mma_allocate(T3,nDens22,Label='T3')
call mma_allocate(T4,nDens22,Label='T4')
call mma_allocate(DIL,nDens,Label='DIL')
call mma_allocate(DI,nCMO,Label='DI')
call mma_allocate(DIR,nDens,Label='DIR')
call mma_allocate(FI,nDens,Label='FI')
call mma_allocate(FI1,nDens,Label='FI1')
call mma_allocate(K1,nDens,Label='K1')

FockI(:) = Zero
FockA(:) = Zero
FI(:) = Zero
FI1(:) = Zero
K1(:) = Zero
DIR(:) = Zero
DIL(:) = Zero
DI(:) = Zero
lFI = .true.
lFa = .false.
lMo = .false.

if (iMethod == 2) then
  call mma_allocate(DAL,nDens,Label='DAL')
  call mma_allocate(DAR,nDens,Label='DAR')
  call mma_allocate(DA,nCMO,Label='DA')
  call mma_allocate(FA1,nDens,Label='FA1')
  lFa = .true.
  lMo = .true.
else
  call mma_allocate(DAL,1,Label='DAL')
  call mma_allocate(DAR,1,Label='DAR')
  call mma_allocate(DA,1,Label='DA')
  call mma_allocate(FA1,1,Label='FA1')
end if
FA1(:) = Zero
DAL(:) = Zero
DAR(:) = Zero
DA(:) = Zero
! THIS IS THE ENTIRE DENSITY FOR MULTICONF
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
        ip2 = iTri(iA,jA)
        DA(ip) = G1t(ip2)
      end do
    end do
  end do
end if

! Construct {kappa,(ij|kl)} and the fock matrix contribution
! from one index tranformation of contracted indexes

!if (iDsym == 2) then
!  call RecPrt('FI',' ',FI,nDens,1)
!  call RecPrt('DI',' ',DI,nDens,1)
!  stop 10
!end if
FacR = Fact
call Read2_ns(rmo1,rmo2,FockI,FockA,T1,nDens22,Tmp2,T3,T4,DIR,DIL,DI,DAR,DAL,DA,rkappa,idsym,Sgn,Facr,jSpin,lFA,lfi,lMo)
!if (iDsym == 2) then
!  call RecPrt('FI',' ',FI,nDens,1)
!  call RecPrt('DI',' ',DI,nDens,1)
!  stop 10
!end if
do iS=1,nsym
  js = Mul(idsym,is)
  if (nbas(is)*nbas(js) /= 0) call DGETMO(rkappa(ipmat(is,js)),nbas(is),nbas(is),nbas(js),K1(ipmat(js,is)),nbas(js))
end do
K1(:) = -K1(:)
DIR(:) = Zero
DIL(:) = Zero
if (imethod == 2) then
  DAR(:) = Zero
  DAL(:) = Zero
end if
call Read2_ns(rdum,rdum,FI1,FA1,T1,nDens22,Tmp2,T3,T4,DIR,DIL,DI,DAR,DAL,DA,K1,idsym,Sgn,Facr,jSpin,lFA,lfi,.false.)
!if (iDsym == 2) call RecPrt('FI1',' ',FI1,nDens,1)

! Calculate contribution from uncontracted indexes.

do iS=1,nSym
  jS = Mul(iS,iDSym)
  if (nBas(iS)*nBas(jS) /= 0) then
    call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sgn*Facr,FIMO(ipCM(iS)),nBas(is),rkappa(ipMat(iS,jS)),nBas(iS),One, &
                FockI(ipMat(iS,jS)),nBas(iS))
    call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,rkappa(ipMat(iS,jS)),nBas(is),FIMO(ipCM(jS)),nBas(jS),One, &
                FockI(ipMat(iS,jS)),nBas(is))
    if (iMethod == 2) then
      call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sgn*Facr,FAMO(ipCM(iS)),nBas(is),rkappa(ipMat(iS,jS)),nBas(iS),One, &
                  FockA(ipMat(iS,jS)),nBas(iS))
      call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,rkappa(ipMat(iS,jS)),nBas(is),FAMO(ipCM(jS)),nBas(jS),One, &
                  FockA(ipMat(iS,jS)),nBas(is))
    end if
  end if
end do

call mma_deallocate(DA)
call mma_deallocate(DAR)
call mma_deallocate(DAL)
call mma_deallocate(FA1)
call mma_deallocate(FI)
call mma_deallocate(FI1)
call mma_deallocate(K1)
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

end subroutine r2elint_ns
