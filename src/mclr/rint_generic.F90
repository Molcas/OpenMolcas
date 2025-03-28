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

!#define _DEBUGPRINT_
subroutine RInt_Generic(rkappa,rmos,rmoa,Fock,Q,Focki,Focka,idsym,reco,jspin)
!                          ~
! Constructs  F  = <0|[E  ,H]|0> (+ <0|[[E  , Kappa],H]|0>)
!              pq       pq                pq

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use MCLR_Data, only: W_CMO_Inv => CMO_Inv, W_CMO => CMO, G1t, G2t, FAMO, FIMO
use MCLR_Data, only: nDens2, ipCM, ipMat, ipMatBA, nA, nMBA
#ifdef _DEBUGPRINT_
use Spool, only: LuWr
#endif
use input_mclr, only: NewCho, iMethod, nSym, IsPop, LuAChoVec, LuIChoVec, nAsh, nBas, nIsh, nOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: u6

implicit none
real*8 rkappa(nDens2), rMOs(*), rmoa(*), Fock(nDens2), Q(ndens2), FockI(ndens2), FockA(nDens2)
integer iDSym, jSpin
real*8 reco
logical Fake_CMO2, DoAct
real*8, allocatable :: MT1(:), MT2(:), MT3(:), QTemp(:), Dens2(:), G2x(:)
type(DSBA_Type) CVa(2), DLT(1), DI, DA, Kappa, JI(1), KI, JA, KA, FkI, FkA, QVec, CMO, CMO_Inv
real*8 Fact, Dij
integer iS, iB, jS, nA2, nAct, nG2, iSym, nAG2, jSym, kSym, nAtri, iOff, iOff2, iOff3, iOff4, iOff5, jB, ip2, ipGx, ijS, kS, lS, &
        kAsh, lAsh, ikl, iAsh, jAsh, iij, iRead, ipF, ipFI
#ifdef _DEBUGPRINT_
integer nas
real*8, external :: DDot_
#endif
! Statement function
integer i, j, iTri
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(LuWr,*) 'Focki=',DDot_(nDens2,Focki,1,Focki,1)
write(LuWr,*) 'Focka=',DDot_(nDens2,Focka,1,Focka,1)
#endif

Fact = -One
call dcopy_(ndens2,[Zero],0,Fock,1)
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. NewCho) then  ! Cho-MO

  call mma_allocate(MT1,nmba,Label='MT1')
  call mma_allocate(MT2,nmba,Label='MT2')
  MT1(:) = Zero
  MT2(:) = Zero

  call R2ElInt(rkappa,MT1,MT2,focki,focka,idSym,ReCo,Fact,jspin)
# ifdef _DEBUGPRINT_
  write(LuWr,*) 'MT1=',DDot_(nmba,MT1,1,MT1,1)
  write(LuWr,*) 'MT2=',DDot_(nmba,MT2,1,MT2,1)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Q  = sum(jkl)=(pj|kl)d(ijkl)
  !  pi

  if (iMethod == 2) then
    call mma_allocate(QTemp,ndens2,Label='QTemp')
    call CreQ(Q,MT1,G2t,idsym)
    call CreQ(QTemp,MT2,G2t,idsym)
#   ifdef _DEBUGPRINT_
    write(LuWr,*) 'Q=',DDot_(nDens2,Q,1,Q,1)
    write(LuWr,*) 'QTemp=',DDot_(nDens2,QTemp,1,QTemp,1)
#   endif
    call daxpy_(ndens2,One,QTemp,1,Q,1)
    call mma_deallocate(QTemp)
  end if

else  ! Cho-Fock

  Fake_CMO2 = .false.
  DoAct = .true.

  ! Form inactive density

  call Allocate_DT(DI,nOrb,nOrb,nSym)
  DI%A0(:) = Zero

  do iS=1,nsym
    do iB=1,nIsh(iS)
      DI%SB(iS)%A2(ib,ib) = Two
    end do
  end do

  ! Form AO 1-index transform inactive density

  call mma_allocate(Dens2,nDens2,Label='Dens2')
  Dens2(:) = Zero
  call Allocate_DT(DLT(1),nOrb,nOrb,nSym) ! Note SQ format
  DLT(1)%A0(:) = Zero

  if (idSym /= 1) then
    write(u6,*) 'idSym/=1, idSym=',idsym
    call Abend()
  end if

  do iS=1,nSym
    if (nOrb(iS) /= 0) then
      do jS=1,nSym
        if ((ieor(iS-1,jS-1)+1 == idsym) .and. (nOrb(jS) /= 0)) then
          call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),One,rkappa(ipMat(is,js)),nOrb(iS),DI%SB(js)%A2,nOrb(jS),Zero, &
                      Dens2(ipMat(iS,jS)),nOrb(iS))
          call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),One,Dens2(ipMat(iS,jS)),nOrb(iS),W_CMO(ipCM(is)),nOrb(iS),Zero, &
                      DLT(1)%SB(iS)%A2,nOrb(jS))
          call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),One,DLT(1)%SB(iS)%A2,nOrb(iS),W_CMO(ipCM(js)),nOrb(jS),Zero, &
                      Dens2(ipMat(iS,jS)),nOrb(jS))

          call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),One,DI%SB(js)%A2,nOrb(iS),W_CMO(ipCM(is)),nOrb(iS),Zero,DLT(1)%SB(iS)%A2, &
                      nOrb(jS))
          call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),One,DLT(1)%SB(iS)%A2,nOrb(iS),W_CMO(ipCM(js)),nOrb(jS),Zero,DI%SB(js)%A2, &
                      nOrb(jS))
        end if
      end do
    end if
  end do
  ! Release and allocate again in LT format
  call Deallocate_DT(DLT(1))
  call Allocate_DT(DLT(1),nOrb,nOrb,nSym,aCase='TRI')

  call Fold_Mat(nSym,nOrb,Dens2,DLT(1)%A0)
  DLT(1)%A0(:) = ReCo*DLT(1)%A0(:)

  ! Form active density and MO coefficients

  if (iMethod == 2) then
    na2 = 0
    nAct = 0
    nG2 = 0
    do iSym=1,nSym
      na2 = na2+nAsh(iSym)**2
      nact = nact+nAsh(iSym)
      nAG2 = 0
      do jSym=1,nSym
        kSym = ieor(jsym-1,isym-1)+1
        nAG2 = nAg2+nAsh(jSym)*nAsh(kSym)
      end do
      nG2 = nG2+nAG2**2
    end do
    nAtri = nact*(nact+1)/2
    nAtri = nAtri*(nAtri+1)/2

    call Allocate_DT(CVa(1),nAsh,nOrb,nSym)
    call Allocate_DT(CVa(2),nAsh,nOrb,nSym)
    call Allocate_DT(DA,nAsh,nAsh,nSym)

    ioff = 0
    ioff4 = 1
    do iS=1,nSym
      ioff2 = ioff+nOrb(iS)*nIsh(iS)
      ioff5 = ioff4+nOrb(iS)*nIsh(iS)
      do iB=1,nAsh(iS)
        ioff3 = ioff2+nOrb(iS)*(iB-1)
        CVa(1)%SB(iS)%A2(iB,:) = W_CMO(ioff3+1:ioff3+nOrb(iS))
        do jB=1,nAsh(iS)
          ip2 = itri(nA(is)+ib,nA(is)+jb)
          DA%SB(iS)%A2(iB,jB) = G1t(ip2)
        end do
      end do
      !MGD to check
      call DGEMM_('T','T',nAsh(iS),nOrb(iS),nOrb(iS),One,rkappa(ioff5),nOrb(iS),W_CMO(1+ioff),nOrb(iS),Zero,CVa(2)%SB(iS)%A2, &
                  nAsh(iS))
      ioff = ioff+(nIsh(iS)+nAsh(iS))*nOrb(iS)
      ioff4 = ioff4+nOrb(iS)**2
    end do

    ! Expand 2-body density matrix

    call mma_allocate(G2x,nG2,Label='G2x')
    ipGx = 0
    do ijS=1,nSym
      do iS=1,nSym
        jS = ieor(is-1,ijS-1)+1
        do kS=1,nSym
          lS = ieor(kS-1,ijS-1)+1
          do kAsh=1,nAsh(ks)
            do lAsh=1,nAsh(ls)
              ikl = itri(lAsh+nA(lS),kAsh+nA(kS))
              do iAsh=1,nAsh(is)
                do jAsh=1,nAsh(js)
                  iij = itri(iAsh+nA(is),jAsh+nA(jS))
                  ipGx = ipGx+1
                  G2x(ipGx) = G2t(itri(iij,ikl))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end if

  ! Allocate temp arrays and zero Fock matrices

  call dcopy_(nATri,[Zero],0,rMOs,1)
# ifdef _DEBUGPRINT_
  call RecPrt('DLT',' ',DLT(1)%A0,1,size(DLT(1)%A0))
  call RecPrt('DI ',' ',DI%A0,1,size(DI%A0))
  call RecPrt('DA ',' ',DA%A0,1,size(DA%A0))
  call RecPrt('G2x',' ',G2x,1,size(G2x))
  call RecPrt('rKappa ',' ',rKappa,1,size(rKappa))
# endif

  ! Compute the whole thing

  call Allocate_DT(Kappa,nBas,nBas,nSym,Ref=rKappa)
  call Allocate_DT(JI(1),nBas,nBas,nSym,aCase='TRI')
  JI(1)%A0(:) = Zero
  call Allocate_DT(KI,nBas,nBas,nSym)
  KI%A0(:) = Zero
  call Allocate_DT(JA,nBas,nBas,nSym)
  JA%A0(:) = Zero
  call Allocate_DT(KA,nBas,nBas,nSym)
  KA%A0(:) = Zero
  call Allocate_DT(FkI,nBas,nBas,nSym,Ref=FockI)
  FkI%A0(:) = Zero
  call Allocate_DT(FkA,nBas,nBas,nSym,Ref=FockA)
  FkA%A0(:) = Zero
  call Allocate_DT(QVec,nBas,nAsh,nSym,Ref=Q)
  call Allocate_DT(CMO,nBas,nAsh,nSym,Ref=W_CMO)
  call Allocate_DT(CMO_Inv,nBas,nAsh,nSym,Ref=W_CMO_Inv)
  iread = 2 ! Asks to read the half-transformed Cho vectors

  call CHO_LK_MCLR(DLT,DI,DA,G2x,Kappa,JI,KI,JA,KA,FkI,FkA,rMOs,QVec,CVa,CMO,CMO_inv,nIsh,nAsh,DoAct,Fake_CMO2,LuAChoVec, &
                   LuIChoVec,iread)

  call Deallocate_DT(CMO_Inv)
  call Deallocate_DT(CMO)
  call Deallocate_DT(QVec)
  call Deallocate_DT(FkA)
  call Deallocate_DT(FkI)
  call Deallocate_DT(KA)
  call Deallocate_DT(JA)
  call Deallocate_DT(KI)
  call Deallocate_DT(JI(1))
  call Deallocate_DT(Kappa)
  call GADSum(FockI,nDens2)
  call GADSum(FockA,nDens2)
# ifdef _DEBUGPRINT_
  nas = 0
  do iSym=1,nSym
    write(u6,*) 'iSym=',iSym
    !call RecPrt('FIMO ',' ',FIMO(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    !call RecPrt('FAMO ',' ',FAMO(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('FockI',' ',FockI(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('FockA',' ',FockA(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    !call RecPrt('Q',' ',Q(ipMatba(iSym,iSym)),nOrb(iSym),nAsh(iSym))
    nas = nas+nAsh(iSym)
  end do
  nAtri = nas*(nas+1)/2
  nAtri = nAtri*(nAtri+1)/2
  !call RecPrt('MO1',' ',rMOs,1,nAtri)
# endif

  ! Calculate contribution from uncontracted indexes

  do iS=1,nSym
    jS = ieor(iS-1,iDSym-1)+1
    if (nOrb(iS)*nOrb(jS) /= 0) then
      call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Reco*Fact,FIMO(ipCM(iS)),nOrb(is),rkappa(ipMat(iS,jS)),nOrb(iS),One, &
                  FockI(ipMat(iS,jS)),nOrb(iS))
      call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Fact,rkappa(ipMat(iS,jS)),nOrb(is),FIMO(ipCM(jS)),nOrb(jS),One, &
                  FockI(ipMat(iS,jS)),nOrb(is))
      if (iMethod == 2) then
        call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Reco*Fact,FAMO(ipCM(iS)),nOrb(is),rkappa(ipMat(iS,jS)),nOrb(iS),One, &
                    FockA(ipMat(iS,jS)),nOrb(iS))
        call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Fact,rkappa(ipMat(iS,jS)),nOrb(is),FAMO(ipCM(jS)),nOrb(jS),One, &
                    FockA(ipMat(iS,jS)),nOrb(is))
      end if
    end if
  end do

  ! Deallocate

  call mma_deallocate(Dens2)

  if (iMethod == 2) then
    call mma_deallocate(G2x)
    call Deallocate_DT(CVa(2))
    call Deallocate_DT(CVa(1))
    call deallocate_DT(DA)
  end if

  call deallocate_DT(DLT(1))
  call deallocate_DT(DI)

  call GADSum(Q,nDens2)
  call GADSum(rMOs,nAtri)

# ifdef _DEBUGPRINT_
  nas = 0
  do iSym=1,nSym
    write(u6,*) 'iSym=',iSym
    call RecPrt('FockI',' ',FockI(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('FockA',' ',FockA(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('Q',' ',Q(ipMatba(iSym,iSym)),nOrb(iSym),nAsh(iSym))
    nas = nas+nAsh(iSym)
  end do
  nAtri = nas*(nas+1)/2
  nAtri = nAtri*(nAtri+1)/2
  call RecPrt('MO1',' ',rMOs,1,nAtri)
  call abend()
# endif

end if
!                                                                      *
!***********************************************************************
!                                                                      *

do iS=1,nSym
  jS = ieor(iS-1,idsym-1)+1

  !         I    A
  ! F  = 2 ( F  + F  )
  !  pi       pi   pi

  if (nIsh(iS)*nOrb(jS) > 0) then
    call DaXpY_(nIsh(iS)*nOrb(jS),Two,Focki(ipMat(js,is)),1,Fock(ipMat(js,is)),1)
    if (iMethod == 2) call DaXpY_(nIsh(iS)*nOrb(jS),Two,Focka(ipMat(js,is)),1,Fock(ipMat(js,is)),1)
  end if

  if (nOrb(iS) > 0) then
    do iAsh=1,nAsh(jS)
      do jAsh=1,nAsh(js)
        ipF = ipMat(js,is)+nIsh(js)+jAsh-1
        ipFI = ipMat(is,js)+(nIsh(js)+iAsh-1)*nOrb(is)
        Dij = G1t(itri(iash+nA(js),jAsh+nA(js)))

        !         I
        ! F  = F - F  D
        !  ap   ap  ap ba

        ipF = ipMat(is,js)+(Nish(js)+iAsh-1)*nOrb(is)
        ipFI = ipMat(is,js)+(Nish(js)+jAsh-1)*nOrb(is)

        !         I
        ! F  = F + F  D
        !  pa   pa  pb ab

        call DaXpY_(nOrb(is),Dij,Focki(ipFI),1,Fock(ipF),1)
      end do
    end do
  end if

  ! F  = F  + Q
  !  pa   pa   pa

  if (nAsh(iS)*nOrb(jS) > 0) call DaXpY_(nAsh(is)*nOrb(js),One,Q(ipMatba(js,is)),1,Fock(ipMat(js,is)+nOrb(js)*nIsh(is)),1)

  ! F  = F  - Q
  !  ap   ap   ap

end do
#ifdef _DEBUGPRINT_
write(LuWr,*) 'Fock=',DDot_(nDens2,Fock,1,Fock,1)
#endif

call DYAX(ndens2,Two,Fock,1,Focka,1)
do iS=1,nSym
  js = ieor(is-1,idsym-1)+1
  if (nOrb(is)*nOrb(js) /= 0) &
    call DGESUB(Focka(ipMat(is,js)),nOrb(is),'N',Focka(ipMat(js,is)),nOrb(js),'T',Fock(ipMat(is,js)),nOrb(is),nOrb(is),nOrb(js))
end do
#ifdef _DEBUGPRINT_
write(LuWr,*) 'Fock=',DDot_(nDens2,Fock,1,Fock,1)
#endif

call AddGrad(rKappa,Fock,idsym,Two*fact)
if (.not. newCho) then
  call mma_allocate(MT3,nmba,Label='MT3')
  call DZAXPY(nmba,One,MT1,1,MT2,1,MT3,1)
  call PickMO_MCLR(MT3,rmos,idsym)

  if (ispop /= 0) then
    call DZAXPY(nmba,-One,MT1,1,MT2,1,MT3,1)
    call PickMO_MCLR(MT3,rmoa,idsym)
  end if
  call mma_deallocate(MT3)
  call mma_deallocate(MT2)
  call mma_deallocate(MT1)
end if
#ifdef _DEBUGPRINT_
write(LuWr,*) 'Exit RInt_Generic'
#endif

end subroutine RInt_Generic
