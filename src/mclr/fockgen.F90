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
use Data_structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use MCLR_Data, only: CMO, FIMO
use MCLR_Data, only: nDens2, nNA, ipCM, ipMat, nA
use input_mclr, only: nSym, nAsh, nIsh, nBas, NewCho, LuAChoVec, nOrb
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two

implicit none
real*8 d_0
integer idSym
real*8 Fock(nDens2), FockOut(*), rDens2(*), rDens1(nna,nna)
real*8, allocatable :: MO(:), Scr(:), G2x(:), Scr1(:,:)
type(DSBA_type) :: CVa
integer n1, iS, n2, ipS, kS, jS, iA, iAA, jA, jAA, ipF, ipM, kA, nG2, iSym, nAG2, jSym, kSym, ipGx, ijS, lS, iAsh, jAsh, kAsh, &
        lAsh, iij, iOff, iOff2, iB, iOff3, ip1, ip2, ikl
real*8 rd

!                                                                      *
!***********************************************************************
!                                                                      *
Fock(:) = Zero

n1 = 0
do iS=1,nSym
  n1 = max(n1,nBas(iS))
end do
n2 = n1**2

if (doDMRG) call dmrg_spc_change_mclr(RGras2(1:8),nash)  ! yma

!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. NewCho) then  ! Cho-MO

  call mma_allocate(MO,n2,Label='MO')
  call mma_allocate(Scr,n2,Label='Scr')

  do ipS=1,nSym
    do kS=1,nSym
      do iS=1,nSym
        jS = ieor(ieor(ipS-1,kS-1),iS-1)+1
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Coulomb term: F  = 2(pk|ji)d
        !                kp           ij
        !                                                              *
        !***************************************************************
        !                                                              *
        if ((ieor(ipS-1,kS-1)+1 == iDsym) .and. (nBas(ipS)*nIsh(kS) > 0)) then
          do iA=1,nAsh(iS)
            iAA = iA+nIsh(iS)
            do jA=1,nAsh(jS)
              jAA = jA+nIsh(jS)

              call Coul(ipS,kS,iS,jS,iAA,jAA,MO,Scr)

              rD = rDens1(iA+nA(iS),jA+nA(jS))*Two
              call DaXpY_(nBas(ipS)*nIsh(kS),rd,MO,1,Fock(ipMat(ipS,Ks)),1)

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
        if ((ieor(ipS-1,iS-1)+1 == iDsym) .and. (nBas(ipS) > 0)) then
          do iA=1,nIsh(iS)
            ipF = ipMat(ipS,iS)+nBas(ipS)*(iA-1)
            do jA=1,nAsh(jS)
              jAA = jA+nIsh(jS)

              call Coul(ipS,kS,iS,jS,iA,jAA,MO,Scr)

              ipM = 1+nIsh(kS)*nBas(ipS)
              do kA=1,nAsh(kS)

                rd = rDens1(kA+nA(kS),jA+nA(jS))
                call DaXpY_(nBas(ipS),-rd,MO(ipM),1,Fock(ipF),1)
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
      kSym = ieor(jsym-1,isym-1)+1
      nAG2 = nAg2+nAsh(jSym)*nAsh(kSym)
    end do
    nG2 = nG2+nAG2**2
  end do

  ! Unfold 2-DM

  call mma_allocate(G2x,nG2,Label='G2x')
  ipGx = 0
  do ijS=1,nSym
    do iS=1,nSym
      jS = ieor(is-1,ijS-1)+1
      do kS=1,nSym
        lS = ieor(kS-1,ijS-1)+1
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
      call dcopy_(nOrb(iS),CMO(1+ioff3:),1,CVa%SB(iS)%A1(iB:),nAsh(iS))
    end do
    ioff = ioff+(nIsh(iS)+nAsh(iS))*nOrb(iS)
  end do

  call mma_allocate(Scr1,n2,2,Label='Scr1')
  Scr1(:,:) = Zero

  call cho_fock_mclr(rdens1,G2x,Scr1(:,1),Scr1(:,2),Fock,CVa,CMO,nIsh,nAsh,LuAChoVec)

  call mma_deallocate(Scr1)
  call Deallocate_DT(CVa)
  call mma_deallocate(G2x)

  call GADSum(Fock,nDens2)

end if

!***********************************************************************
!       Common part                                                    *
!***********************************************************************

do iS=1,nSym
  if (nBas(iS) > 0) then
    jS = ieor(is-1,iDSym-1)+1
    do iA=1,nAsh(is)
      do jA=1,nAsh(js)
        rd = rDens1(iA+nA(iS),jA+nA(js))
        ip1 = nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
        ip2 = nBas(iS)*(nIsh(js)+jA-1)+ipmat(is,js)
        call DaXpY_(nBas(iS),Rd,FIMO(ip1),1,Fock(ip2),1)
      end do
    end do
  end if
end do

if (iDsym == 1) then
  do iS=1,nSym
    if (nBas(iS)*nIsh(iS) > 0) call DaXpY_(nBas(iS)*nIsh(is),Two*d_0,FIMO(ipMat(is,is)),1,Fock(ipMat(is,is)),1)
  end do
end if
do iS=1,nSym
  jS = ieor(iS-1,idSym-1)+1
  if (nBas(is)*nBas(jS) /= 0) &
    call DGeSub(Fock(ipMat(iS,jS)),nBas(iS),'N',Fock(ipMat(jS,iS)),nBas(jS),'T',FockOut(ipMat(iS,jS)),nBas(iS),nBas(iS),nBas(jS))
end do

call DScal_(nDens2,Two,FockOut,1)
if (idSym == 1) call AddGrad2(FockOut,d_0)

if (doDMRG) call dmrg_spc_change_mclr(LRras2(1:8),nash)  ! yma
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine FockGen
