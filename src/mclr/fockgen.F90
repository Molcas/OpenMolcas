!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine FockGen(d_0,rDens1,rdens2,Fock,FockOut,idSym)
!***********************************************************************
!                                                                      *
!   Constructs active Fock matrix and Q matrix                         *
!                                                                      *
!   Input: rkappa: Rotation matrix                                     *
!          idsym : symmetry of perturbation                            *
!                                                                      *
!                                                                      *
!   Output:MO     : MO integrals                                       *
!          Fock   : Fock matrix (one index transformed integrals)      *
!          MOtilde: MO (one index transformed integrals)               *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use Data_structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use MCLR_Data, only: CMO, FIMO, ipCM, ipMat, nA, nDens, nNA
use input_mclr, only: LuAChoVec, nAsh, nBas, NewCho, nIsh, nOrb, nSym
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: d_0, rDens1(nna,nna), rDens2(*), Fock(nDens), FockOut(*)
integer(kind=iwp) :: idSym
integer(kind=iwp) :: iA, iAA, iAsh, iB, iij, ijS, ikl, iOff, iOff2, iOff3, ip1, ip2, ipF, ipGx, ipM, ipS, iS, iSym, jA, jAA, jAsh, &
                     jS, jSym, kA, kAsh, kS, lAsh, lS, n1, n2, nAG2, nG2
real(kind=wp) :: rd
type(DSBA_type) :: CVa
real(kind=wp), allocatable :: G2x(:), MO(:), Scr(:), Scr1(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
Fock(:) = Zero

n1 = max(0,maxval(nBas(1:nSym)))
n2 = n1**2

if (doDMRG) nash(:) = RGras2(:)  ! yma

!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. NewCho) then  ! Cho-MO

  call mma_allocate(MO,n2,Label='MO')
  call mma_allocate(Scr,n2,Label='Scr')

  do ipS=1,nSym
    do kS=1,nSym
      do iS=1,nSym
        jS = Mul(Mul(ipS,kS),iS)
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Coulomb term: F  = 2(pk|ji)d
        !                kp           ij
        !                                                              *
        !***************************************************************
        !                                                              *
        if ((Mul(ipS,kS) == iDsym) .and. (nBas(ipS)*nIsh(kS) > 0)) then
          do iA=1,nAsh(iS)
            iAA = iA+nIsh(iS)
            do jA=1,nAsh(jS)
              jAA = jA+nIsh(jS)

              call Coul(ipS,kS,iS,jS,iAA,jAA,MO,Scr)

              rD = rDens1(iA+nA(iS),jA+nA(jS))*Two
              Fock(ipMat(ipS,Ks):ipMat(ipS,Ks)+nBas(ipS)*nIsh(kS)-1) = Fock(ipMat(ipS,Ks):ipMat(ipS,Ks)+nBas(ipS)*nIsh(kS)-1)+ &
                rd*MO(1:nBas(ipS)*nIsh(kS))

            end do
          end do
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Exchange term: F = -(pk|ji)d
        !                 pl          kj
        !                                                              *
        !***************************************************************
        !                                                              *
        if ((Mul(ipS,iS) == iDsym) .and. (nBas(ipS) > 0)) then
          do iA=1,nIsh(iS)
            ipF = ipMat(ipS,iS)+nBas(ipS)*(iA-1)
            do jA=1,nAsh(jS)
              jAA = jA+nIsh(jS)

              call Coul(ipS,kS,iS,jS,iA,jAA,MO,Scr)

              ipM = 1+nIsh(kS)*nBas(ipS)
              do kA=1,nAsh(kS)

                rd = rDens1(kA+nA(kS),jA+nA(jS))
                Fock(ipF:ipF+nBas(ipS)-1) = Fock(ipF:ipF+nBas(ipS)-1)-rd*MO(ipM:ipM+nBas(ipS)-1)
                ipM = ipM+nBas(ipS)
              end do

            end do
          end do
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
      end do
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call CreQADD(Fock,rdens2,idsym,MO,Scr,n2)
  call mma_deallocate(Scr)
  call mma_deallocate(MO)

else  ! Cho-Fock

  !*********************************************************************
  !   new Cholesky code                                                *
  !*********************************************************************
  nG2 = 0
  do iSym=1,nSym
    nAG2 = 0
    do jSym=1,nSym
      nAG2 = nAg2+nAsh(jSym)*nAsh(Mul(jSym,iSym))
    end do
    nG2 = nG2+nAG2**2
  end do

  ! Unfold 2-DM

  call mma_allocate(G2x,nG2,Label='G2x')
  ipGx = 0
  do ijS=1,nSym
    do iS=1,nSym
      jS = Mul(is,ijS)
      do kS=1,nSym
        lS = Mul(kS,ijS)
        do kAsh=1,nAsh(ks)
          do lAsh=1,nAsh(ls)
            !ikl = iTri(lAsh+nA(lS),kAsh+nA(kS))
            ikl = nna*(lAsh+nA(lS)-1)+kAsh+nA(kS)
            do iAsh=1,nAsh(is)
              do jAsh=1,nAsh(js)
                !iij = iTri(iAsh+nA(is),jAsh+nA(jS))
                iij = nna*(jAsh+nA(jS)-1)+iAsh+nA(iS)
                ipGx = ipGx+1
                G2x(ipGx) = rdens2(iTri(iij,ikl))
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  ! Get active CMO

  call Allocate_DT(CVa,nAsh,nBas,nSym)
  CVa%A0(:) = Zero

  ioff = 0
  do iS=1,nSym
    ioff2 = ioff+nOrb(iS)*nIsh(iS)
    do iB=1,nAsh(iS)
      ioff3 = ioff2+nOrb(iS)*(iB-1)
      CVa%SB(iS)%A1(iB:iB+nOrb(iS)*nAsh(iS)-1:nAsh(iS)) = CMO(ioff3+1:ioff3+nOrb(iS))
    end do
    ioff = ioff+(nIsh(iS)+nAsh(iS))*nOrb(iS)
  end do

  call mma_allocate(Scr1,n2,2,Label='Scr1')
  Scr1(:,:) = Zero

  call cho_fock_mclr(rdens1,G2x,Scr1(:,1),Scr1(:,2),Fock,CVa,CMO,nIsh,nAsh,LuAChoVec)

  call mma_deallocate(Scr1)
  call Deallocate_DT(CVa)
  call mma_deallocate(G2x)

  call GADSum(Fock,nDens)

end if

!***********************************************************************
!       Common part                                                    *
!***********************************************************************

do iS=1,nSym
  if (nBas(iS) > 0) then
    jS = Mul(is,iDSym)
    do iA=1,nAsh(is)
      do jA=1,nAsh(js)
        rd = rDens1(iA+nA(iS),jA+nA(js))
        ip1 = nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
        ip2 = nBas(iS)*(nIsh(js)+jA-1)+ipmat(is,js)
        Fock(ip2:ip2+nBas(iS)-1) = Fock(ip2:ip2+nBas(iS)-1)+Rd*FIMO(ip1:ip1+nBas(iS)-1)
      end do
    end do
  end if
end do

if (iDsym == 1) then
  do iS=1,nSym
    Fock(ipMat(is,is):ipMat(is,is)+nBas(iS)*nIsh(is)-1) = Fock(ipMat(is,is):ipMat(is,is)+nBas(iS)*nIsh(is)-1)+ &
                                                          Two*d_0*FIMO(ipMat(is,is):ipMat(is,is)+nBas(iS)*nIsh(is)-1)
  end do
end if
do iS=1,nSym
  jS = Mul(iS,idSym)
  if (nBas(is)*nBas(jS) /= 0) &
    call DGeSub(Fock(ipMat(iS,jS)),nBas(iS),'N',Fock(ipMat(jS,iS)),nBas(jS),'T',FockOut(ipMat(iS,jS)),nBas(iS),nBas(iS),nBas(jS))
end do

FockOut(1:nDens) = Two*FockOut(1:nDens)
if (idSym == 1) call AddGrad2(FockOut,d_0)

if (doDMRG) nash(:) = LRras2(:)  ! yma
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine FockGen
