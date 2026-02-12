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

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use MCLR_Data, only: CMO, CMO_Inv, FAMO, FIMO, G1t, G2t, ipCM, ipMat, ipMatBA, nA, nDens, nMBA
use input_mclr, only: iMethod, IsPop, LuAChoVec, LuIChoVec, nAsh, nBas, NewCho, nIsh, nOrb, nSym
#ifdef _DEBUGPRINT_
use Spool, only: LuWr
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: rkappa(nDens), reco
real(kind=wp), intent(_OUT_) :: rMOs(*), rmoa(*)
real(kind=wp), intent(out) :: Fock(nDens)
real(kind=wp), intent(inout) :: Q(nDens), FockI(nDens), FockA(nDens)
integer(kind=iwp), intent(in) :: iDSym, jSpin
integer(kind=iwp) :: iAsh, iB, iij, ijS, ikl, iOff, iOff2, iOff3, iOff4, iOff5, ip2, ipF, ipFI, ipGx, iRead, iS, iSym, jAsh, jB, &
                     jS, jSym, kAsh, kS, lAsh, lS, nAct, nAG2, nAtri, nG2
real(kind=wp) :: Dij, Fact
logical(kind=iwp) :: DoAct, Fake_CMO2
type(DSBA_Type) :: CVa(2), DA, DI, DLT(1), FkA, FkI, JA, JI(1), KA, Kappa, KI, QVec, WCMO, WCMO_Inv
real(kind=wp), allocatable :: Dens2(:), G2x(:), MT1(:), MT2(:), MT3(:), QTemp(:)
#ifdef _DEBUGPRINT_
real(kind=wp), external :: DDot_
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(LuWr,*) 'Focki=',DDot_(nDens,Focki,1,Focki,1)
write(LuWr,*) 'Focka=',DDot_(nDens,Focka,1,Focka,1)
#endif

Fact = -One
Fock(:) = Zero
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
    call mma_allocate(QTemp,nDens,Label='QTemp')
    call CreQ(Q,MT1,G2t,idsym)
    call CreQ(QTemp,MT2,G2t,idsym)
#   ifdef _DEBUGPRINT_
    write(LuWr,*) 'Q=',DDot_(nDens,Q,1,Q,1)
    write(LuWr,*) 'QTemp=',DDot_(nDens,QTemp,1,QTemp,1)
#   endif
    Q(:) = Q(:)+QTemp(:)
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
      DI%SB(iS)%A2(iB,iB) = Two
    end do
  end do

  ! Form AO 1-index transform inactive density

  call mma_allocate(Dens2,nDens,Label='Dens2')
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
        if ((Mul(iS,jS) == idsym) .and. (nOrb(jS) /= 0)) then
          call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),One,rkappa(ipMat(is,js)),nOrb(iS),DI%SB(js)%A2,nOrb(jS),Zero, &
                      Dens2(ipMat(iS,jS)),nOrb(iS))
          call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),One,Dens2(ipMat(iS,jS)),nOrb(iS),CMO(ipCM(is)),nOrb(iS),Zero, &
                      DLT(1)%SB(iS)%A2,nOrb(jS))
          call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),One,DLT(1)%SB(iS)%A2,nOrb(iS),CMO(ipCM(js)),nOrb(jS),Zero, &
                      Dens2(ipMat(iS,jS)),nOrb(jS))

          call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),One,DI%SB(js)%A2,nOrb(iS),CMO(ipCM(is)),nOrb(iS),Zero,DLT(1)%SB(iS)%A2, &
                      nOrb(jS))
          call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),One,DLT(1)%SB(iS)%A2,nOrb(iS),CMO(ipCM(js)),nOrb(jS),Zero,DI%SB(js)%A2, &
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

  nAct = sum(nAsh(1:nSym))
  nAtri = nTri_Elem(nTri_Elem(nAct))
  if (iMethod == 2) then
    nG2 = 0
    do iSym=1,nSym
      nAG2 = 0
      do jSym=1,nSym
        nAG2 = nAG2+nAsh(jSym)*nAsh(Mul(jSym,iSym))
      end do
      nG2 = nG2+nAG2**2
    end do

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
        CVa(1)%SB(iS)%A2(iB,:) = CMO(ioff3+1:ioff3+nOrb(iS))
        do jB=1,nAsh(iS)
          ip2 = iTri(nA(is)+ib,nA(is)+jb)
          DA%SB(iS)%A2(iB,jB) = G1t(ip2)
        end do
      end do
      !MGD to check
      call DGEMM_('T','T',nAsh(iS),nOrb(iS),nOrb(iS),One,rkappa(ioff5),nOrb(iS),CMO(1+ioff),nOrb(iS),Zero,CVa(2)%SB(iS)%A2,nAsh(iS))
      ioff = ioff+(nIsh(iS)+nAsh(iS))*nOrb(iS)
      ioff4 = ioff4+nOrb(iS)**2
    end do

    ! Expand 2-body density matrix

    call mma_allocate(G2x,nG2,Label='G2x')
    ipGx = 0
    do ijS=1,nSym
      do iS=1,nSym
        jS = Mul(is,ijS)
        do kS=1,nSym
          lS = Mul(kS,ijS)
          do kAsh=1,nAsh(ks)
            do lAsh=1,nAsh(ls)
              ikl = iTri(lAsh+nA(lS),kAsh+nA(kS))
              do iAsh=1,nAsh(is)
                do jAsh=1,nAsh(js)
                  iij = iTri(iAsh+nA(is),jAsh+nA(jS))
                  ipGx = ipGx+1
                  G2x(ipGx) = G2t(iTri(iij,ikl))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end if

  ! Allocate temp arrays and zero Fock matrices

  rMOs(1:nATri) = Zero
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
  call Allocate_DT(WCMO,nBas,nAsh,nSym,Ref=CMO)
  call Allocate_DT(WCMO_Inv,nBas,nAsh,nSym,Ref=CMO_Inv)
  iread = 2 ! Asks to read the half-transformed Cho vectors

  call CHO_LK_MCLR(DLT,DI,DA,G2x,Kappa,JI,KI,JA,KA,FkI,FkA,rMOs,QVec,CVa,WCMO,WCMO_inv,nIsh,nAsh,DoAct,Fake_CMO2,LuAChoVec, &
                   LuIChoVec,iread)

  call Deallocate_DT(WCMO_Inv)
  call Deallocate_DT(WCMO)
  call Deallocate_DT(QVec)
  call Deallocate_DT(FkA)
  call Deallocate_DT(FkI)
  call Deallocate_DT(KA)
  call Deallocate_DT(JA)
  call Deallocate_DT(KI)
  call Deallocate_DT(JI(1))
  call Deallocate_DT(Kappa)
  call GADSum(FockI,nDens)
  call GADSum(FockA,nDens)
# ifdef _DEBUGPRINT_
  do iSym=1,nSym
    write(u6,*) 'iSym=',iSym
    !call RecPrt('FIMO ',' ',FIMO(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    !call RecPrt('FAMO ',' ',FAMO(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('FockI',' ',FockI(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('FockA',' ',FockA(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    !call RecPrt('Q',' ',Q(ipMatba(iSym,iSym)),nOrb(iSym),nAsh(iSym))
  end do
  !call RecPrt('MO1',' ',rMOs,1,nAtri)
# endif

  ! Calculate contribution from uncontracted indexes

  do iS=1,nSym
    jS = Mul(iS,iDSym)
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

  call GADSum(Q,nDens)
  call GADSum(rMOs,nAtri)

# ifdef _DEBUGPRINT_
  do iSym=1,nSym
    write(u6,*) 'iSym=',iSym
    call RecPrt('FockI',' ',FockI(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('FockA',' ',FockA(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('Q',' ',Q(ipMatba(iSym,iSym)),nOrb(iSym),nAsh(iSym))
  end do
  call RecPrt('MO1',' ',rMOs,1,nAtri)
  call abend()
# endif

end if
!                                                                      *
!***********************************************************************
!                                                                      *

do iS=1,nSym
  jS = Mul(iS,idsym)

  !         I    A
  ! F  = 2 ( F  + F  )
  !  pi       pi   pi

  if (nIsh(iS)*nOrb(jS) > 0) then
    Fock(ipMat(js,is):ipMat(js,is)+nIsh(iS)*nOrb(jS)-1) = Fock(ipMat(js,is):ipMat(js,is)+nIsh(iS)*nOrb(jS)-1)+ &
                                                          Two*Focki(ipMat(js,is):ipMat(js,is)+nIsh(iS)*nOrb(jS)-1)
    if (iMethod == 2) Fock(ipMat(js,is):ipMat(js,is)+nIsh(iS)*nOrb(jS)-1) = Fock(ipMat(js,is):ipMat(js,is)+nIsh(iS)*nOrb(jS)-1)+ &
                                                                            Two*Focka(ipMat(js,is):ipMat(js,is)+nIsh(iS)*nOrb(jS)-1)
  end if

  if (nOrb(iS) > 0) then
    do iAsh=1,nAsh(jS)
      do jAsh=1,nAsh(js)
        ipF = ipMat(js,is)+nIsh(js)+jAsh-1
        ipFI = ipMat(is,js)+(nIsh(js)+iAsh-1)*nOrb(is)
        Dij = G1t(iTri(iash+nA(js),jAsh+nA(js)))

        !         I
        ! F  = F - F  D
        !  ap   ap  ap ba

        ipF = ipMat(is,js)+(Nish(js)+iAsh-1)*nOrb(is)
        ipFI = ipMat(is,js)+(Nish(js)+jAsh-1)*nOrb(is)

        !         I
        ! F  = F + F  D
        !  pa   pa  pb ab

        Fock(ipF:ipF+nOrb(is)-1) = Fock(ipF:ipF+nOrb(is)-1)+Dij*Focki(ipFI:ipFI+nOrb(is)-1)
      end do
    end do
  end if

  ! F  = F  + Q
  !  pa   pa   pa

  Fock(ipMat(js,is)+nOrb(js)*nIsh(is):ipMat(js,is)+nOrb(js)*(nIsh(is)+nAsh(is))-1) = &
    Fock(ipMat(js,is)+nOrb(js)*nIsh(is):ipMat(js,is)+nOrb(js)*(nIsh(is)+nAsh(is))-1)+ &
    Q(ipMatba(js,is):ipMatba(js,is)+nAsh(is)*nOrb(js)-1)

  ! F  = F  - Q
  !  ap   ap   ap

end do
#ifdef _DEBUGPRINT_
write(LuWr,*) 'Fock=',DDot_(nDens,Fock,1,Fock,1)
#endif

Focka(:) = Two*Fock(:)
do iS=1,nSym
  js = Mul(is,idsym)
  if (nOrb(is)*nOrb(js) /= 0) &
    call DGESUB(Focka(ipMat(is,js)),nOrb(is),'N',Focka(ipMat(js,is)),nOrb(js),'T',Fock(ipMat(is,js)),nOrb(is),nOrb(is),nOrb(js))
end do
#ifdef _DEBUGPRINT_
write(LuWr,*) 'Fock=',DDot_(nDens,Fock,1,Fock,1)
#endif

call AddGrad(rKappa,Fock,idsym,Two*fact)
if (.not. newCho) then
  call mma_allocate(MT3,nmba,Label='MT3')
  MT3(:) = MT1(:)+MT2(:)
  call PickMO_MCLR(MT3,rmos,idsym)

  if (ispop /= 0) then
    MT3(:) = MT2(:)-MT1(:)
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
