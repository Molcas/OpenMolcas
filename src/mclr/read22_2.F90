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

!#define _DEBUGPRINT_
subroutine Read22_2(MO1,Fock,Q,FockI,FockA,Temp2,Scr,Temp3)
!***********************************************************************
!                                                                      *
!   Constructs         everything                                      *
!                                                                      *
!   Output:MO     :MO integrals                                        *
!          Fock   :Fock matrix (one index transformed integrals)       *
!          MOtilde:MO (one index transformed integrals)                *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use MCLR_Data, only: CMO, CMO_Inv, G1t, G2t, Int1, ipCM, ipMat, ipMatBA, LuQDat, nA, nB, nDens
use input_mclr, only: Debug, iAddressQDat, iMethod, LuAChoVec, LuIChoVec, nAsh, nBas, NewCho, nIsh, nOrb, nSym, PotNuc, rIn_Ene, &
                      StepType, TwoStep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: MO1(*), Scr(*)
real(kind=wp), intent(out) :: Fock(nDens), Q(nDens), FockI(nDens), FockA(nDens), Temp2(nDens), Temp3(nDens)
integer(kind=iwp) :: iAsh, iB, iiB, iIJ, ijB, ijB1, ijS, iK, iKK, iKL, iL, iLL, iOff, iOff2, iOff3, ip1, ip2, ipD, ipGx, ipi, ipj, &
                     ipTmp, iS, iStore, iSym, jAsh, jB, jjB, jS, jSym, kAsh, kB, kkB, kS, lAsh, lB, llB, lS, nA2, nAct, nAG2, &
                     nAtri, nG2, nI, nJ, nm, nNB, nNK, nNL
real(kind=wp) :: Fact, rCor, rCora, rCore, rCoreA, rCoreI, rEnergy
logical(kind=iwp) :: DoAct, Fake_CMO2
type(DSBA_Type) :: CVa(2), DA, DI, DLT(1), FkA, FkI, JA, JI(1), KA, Kappa, KI, QVec, WCMO, WCMO_Inv
real(kind=wp), allocatable :: G2x(:)

!                                                                      *
!***********************************************************************
!                                                                      *
FockI(:) = Zero
FockA(:) = Zero
if (TwoStep .and. (StepType == 'RUN2')) then
  nm = sum(nAsh(1:nSym))
  nAtri = nTri_Elem(nTri_Elem(nm))
  Fock(:) = Zero
  Q(:) = Zero
  MO1(1:nAtri) = Zero
  call ddafile(LuQDAT,2,FockA,nDens,iaddressQDAT)
  call ddafile(LuQDAT,2,FockI,nDens,iaddressQDAT)
  call ddafile(LuQDAT,2,Fock,nDens,iaddressQDAT)
  call ddafile(LuQDAT,2,Q,nDens,iaddressQDAT)
  call ddafile(LuQDAT,2,MO1,nAtri,iaddressQDAT)
else

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (.not. NewCho) then  ! Cho-MO algorithm

    do iS=1,nSym
      do jS=1,iS
        ijS = Mul(iS,jS)
        do kS=1,nSym
          do lS=1,ks
            if (nOrb(iS)*nOrb(jS)*nOrb(ks)*nOrb(lS) == 0) cycle
            if (Mul(kS,lS) /= ijS) cycle
            !                                                          *
            !***********************************************************
            !                                                          *
            ijB1 = 0
            do iB=1,nB(iS)
              nnB = nB(jS)
              if (iS == jS) nnB = iB
              do jB=1,nnB

                ! Read a symmetry block

                call COUL(kS,lS,iS,jS,iB,jB,Temp2,Scr)
                ijB1 = ijB1+1
                !                                                      *
                !*******************************************************
                !                                                      *
                ! Add to unrotated inactive fock matrix
                !
                ! Fkl sum(i)     2(ii|kl) -(ki|li)
                !
                ! Coulomb term: F  =2(ii|kl)
                !                kl

                if ((iS == jS) .and. (iB == jB) .and. (iB <= nIsh(iS))) &
                  Focki(ipCM(kS):ipCM(ks)+nOrb(kS)*nOrb(lS)-1) = Focki(ipCM(kS):ipCM(ks)+nOrb(kS)*nOrb(lS)-1)+ &
                                                                 Two*Temp2(1:nOrb(kS)*nOrb(lS))

                !                                                      *
                !*******************************************************
                !                                                      *
                ! Add to unrotated active fock matrix
                !
                ! Coulomb term: F  =2(ij|kl)d   i<j
                !                kl          ij

                if (iMethod == 2) then

                  if (iS == jS) then
                    if (((iB > nIsh(is)) .and. (nAsh(iS) /= 0)) .and. ((jB > nIsh(js)) .and. (nAsh(jS) /= 0))) then
                      ipD = iTri(jB-nIsh(jS)+nA(jS),iB-nIsh(is)+nA(iS))
                      Fact = Two
                      if (iB == jB) Fact = One
                      FockA(ipCM(kS):ipCM(kS)+nOrb(kS)*nOrb(lS)-1) = FockA(ipCM(kS):ipCM(kS)+nOrb(kS)*nOrb(lS)-1)+ &
                                                                     Fact*G1t(ipD)*Temp2(1:nOrb(kS)*nOrb(lS))
                    end if
                  end if

                end if
                !                                                      *
                !*******************************************************
                !                                                      *
              end do  ! jB
            end do    ! iB
            !                                                          *
            !***********************************************************
            !                                                          *
          end do      ! lS
        end do        ! kS
      end do          ! jS
    end do            ! iS
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Construct Q matrix: Q = sum(jkl)(pj|kl)d
    !                      pi                 ijkl

    if (iMethod == 2) then

      call CreQ2(Q,G2t,1,Temp2,Scr,nDens)

      ! Sort out MO (ij|kl)

      do iS=1,nSym
        do jS=1,iS
          do kS=1,iS
            lS = Mul(Mul(iS,jS),kS)
            if ((lS > kS) .or. ((iS == kS) .and. (lS > jS))) cycle

            do iB=1,nAsh(iS)
              iib = ib+nA(iS)
              nnb = nAsh(jS)
              if (iS == jS) nnb = ib
              do jB=1,nnB
                jjb = jb+nA(jS)

                call Coul(kS,lS,iS,jS,iB+nIsh(iS),jB+nIsh(jS),Temp2,Scr)

                nnK = nAsh(kS)
                if (iS == kS) nnK = iB
                do kB=1,nnk
                  kkb = kb+nA(kS)
                  nnL = nAsh(lS)
                  if (kS == lS) nnL = kB
                  if (iib == kkb) nnL = jB
                  do lB=1,nnL
                    llb = lb+nA(lS)

                    ip2 = (lB+nIsh(lS)-1)*nBas(kS)+kB+nIsh(kS)

                    ip1 = iTri(iTri(iib,jjb),iTri(kkb,llb))
                    MO1(ip1) = Temp2(ip2)
                  end do
                end do
              end do
            end do

          end do
        end do
      end do
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    do iS=1,nSym
      kS = iS
      do js=1,nSym
        lS = jS
        if (Mul(is,js) /= Mul(ks,ls)) cycle
        if (nOrb(iS)*nOrb(jS)*nOrb(ks)*nOrb(lS) == 0) cycle
        do LB=1,nB(LS)
          do JB=1,nB(JS)
            !                                                          *
            !***********************************************************
            !                                                          *
            call EXCH(is,js,ks,ls,jb,lb,Temp2,Scr)
            !                                                          *
            !***********************************************************
            !                                                          *
            ! Add to unrotated inactive fock matrix
            !
            ! Exchange term: F  =-(ij|kj)
            !                 ik

            if ((jS == lS) .and. (jB == lB) .and. (jB <= nIsh(jS))) &
              Focki(ipCM(iS):ipCM(iS)+nOrb(iS)*nOrb(kS)-1) = Focki(ipCM(iS):ipCM(iS)+nOrb(iS)*nOrb(kS)-1)-Temp2(1:nOrb(iS)*nOrb(kS))

            !                                                          *
            !***********************************************************
            !                                                          *
            ! Add to unrotated active fock matrix
            !
            ! Exchange term: F  =-1/2(ij|kl)d
            !                 ik             jl

            if (iMethod == 2) then
              if (jS == lS) then
                if (((jB > nIsh(js)) .and. (nAsh(jS) /= 0)) .and. ((lB > nIsh(ls)) .and. (nAsh(lS) /= 0))) then
                  ipD = iTri(lB-nIsh(lS)+nA(lS),jB-nIsh(js)+nA(jS))
                  FockA(ipCM(iS):ipCM(iS)+nOrb(iS)*nOrb(kS)-1) = FockA(ipCM(iS):ipCM(iS)+nOrb(iS)*nOrb(kS)-1)- &
                                                                 Half*G1t(ipD)*Temp2(1:nOrb(iS)*nOrb(kS))
                end if
              end if
            end if
            !                                                          *
            !***********************************************************
            !                                                          *
          end do
        end do
      end do ! jS
    end do ! iS

  else  ! Cho-Fock Algorithm

    Fake_CMO2 = .true.
    DoAct = .true.

    ! Construct inactive density matrix

    Temp2(:) = Zero
    do is=1,nSym
      do iB=1,nIsh(is)
        Temp2(ipCM(iS)+(iB-1)*nOrb(iS)+iB-1) = Two
      end do
    end do

    ! Transform to AO basis

    do iS=1,nSym
      if (nIsh(iS) /= 0) then
        jS = iS
        call DGEMM_('T','T',nIsh(jS),nOrb(iS),nIsh(iS),One,Temp2(ipCM(iS)),nOrb(iS),CMO(ipCM(is)),nOrb(iS),Zero, &
                    Temp3(ipMat(jS,iS)),nOrb(jS))
        call DGEMM_('T','T',nOrb(jS),nOrb(jS),nIsh(iS),One,Temp3(ipMat(jS,iS)),nOrb(iS),CMO(ipCM(js)),nOrb(jS),Zero, &
                    Temp2(ipCM(iS)),nOrb(jS))
      end if
    end do

    call Allocate_DT(DLT(1),nBas,nBas,nSym,aCase='TRI')
    call Fold_Mat(nSym,nOrb,Temp2,DLT(1)%A0)

    ! Form active CMO and density

    nAct = 0
    if (iMethod == 2) then
      na2 = 0
      nG2 = 0
      do iSym=1,nSym
        na2 = na2+nAsh(iSym)**2
        nAct = nAct+nAsh(iSym)
        nAG2 = 0
        do jSym=1,nSym
          nAG2 = nAg2+nAsh(jSym)*nAsh(Mul(jSym,iSym))
        end do
        nG2 = nG2+nAG2**2
      end do
      call Allocate_DT(CVa(1),nAsh,nOrb,nSym)
      CVa(1)%A0(:) = Zero
      call Allocate_DT(CVa(2),nAsh,nOrb,nSym)
      CVa(2)%A0(:) = Zero
      call Allocate_DT(DA,nAsh,nAsh,nSym)

      ioff = 0
      do iSym=1,nSym
        ioff2 = ioff+nOrb(iSym)*nIsh(iSym)
        do ikk=1,nAsh(iSym)
          ioff3 = ioff2+nOrb(iSym)*(ikk-1)
          CVa(1)%SB(iSym)%A2(ikk,:) = CMO(ioff3+1:ioff3+nOrb(iSym))
          ik = ikk+nA(iSym)
          do ill=1,ikk-1
            il = ill+nA(iSym)
            ikl = iTri(ik,il)
            DA%SB(iSym)%A2(ill,ikk) = G1t(ikl)
            DA%SB(iSym)%A2(ikk,ill) = G1t(ikl)
          end do
          ikl = nTri_Elem(ik)
          DA%SB(iSym)%A2(ikk,ikk) = G1t(ikl)
        end do
        ioff = ioff+nOrb(iSym)**2
      end do
      DA%A0(:) = Half*DA%A0(:)

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
    else
      na2 = 1
      nG2 = 1
      call Allocate_DT(CVa(1),[1],[1],1) ! dummy allocation
      call Allocate_DT(CVa(2),[1],[1],1)
      call Allocate_DT(DA,nAsh,nAsh,nSym)
      call mma_allocate(G2x,nG2,Label='G2x')
    end if

    ! Let's go

    call Allocate_DT(JA,nBas,nBas,nSym)
    JA%A0(:) = Zero
    call Allocate_DT(KA,nBas,nBas,nSym)
    KA%A0(:) = Zero

    Q(:) = Zero

    call Allocate_DT(DI,nBas,nBas,nSym,Ref=Temp2)
    call Allocate_DT(JI(1),nBas,nBas,nSym,aCase='TRI',Ref=Temp3)
    JI(1)%A0(:) = Zero
    call Allocate_DT(KI,nBas,nBas,nSym,Ref=Scr)
    KI%A0(:) = Zero
    call Allocate_DT(FkI,nBas,nBas,nSym,Ref=FockI)
    FkI%A0(:) = Zero
    call Allocate_DT(FkA,nBas,nBas,nSym,Ref=FockA)
    FkA%A0(:) = Zero
    call Allocate_DT(QVec,nBas,nAsh,nSym,Ref=Q)
    call Allocate_DT(WCMO,nBas,nBas,nSym,Ref=CMO)
    call Allocate_DT(WCMO_Inv,nBas,nBas,nSym,Ref=CMO_Inv)
    istore = 1 ! Ask to store the half-transformed vectors

    call CHO_LK_MCLR(DLT,DI,DA,G2x,Kappa,JI,KI,JA,KA,FkI,FkA,MO1,QVec,CVa,WCMO,WCMO_inv,nIsh,nAsh,doAct,Fake_CMO2,LuAChoVec, &
                     LuIChoVec,istore)

    nAtri = nTri_Elem(nAct)
    nAtri = nTri_Elem(nAtri)
    MO1(1:nAtri) = Quart*MO1(1:nAtri)
    FkI%A0(:) = -Half*FkI%A0(:)

    call Deallocate_DT(WCMO_Inv)
    call Deallocate_DT(WCMO)
    call Deallocate_DT(QVec)
    call Deallocate_DT(FkA)
    call Deallocate_DT(FkI)
    call Deallocate_DT(KI)
    call Deallocate_DT(JI(1))
    call Deallocate_DT(DI)
    call Deallocate_DT(JA)
    call Deallocate_DT(KA)
    call deallocate_DT(DLT(1))
    call mma_deallocate(G2x)
    call Deallocate_DT(CVa(2))
    call Deallocate_DT(CVa(1))
    call deallocate_DT(DA)

    call GADSum(FockI,nDens)
    call GADSum(FockA,nDens)
    call GADSum(Q,nDens)
    call GADSum(MO1,nAtri)

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
# ifdef _DEBUGPRINT_
  do iSym=1,nSym
    write(u6,*) 'iSym=',iSym
    call RecPrt('FockI',' ',FockI(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('FockA',' ',FockA(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
    call RecPrt('Q',' ',Q(ipMatba(iSym,iSym)),nOrb(iSym),nAsh(iSym))
  end do
  nm = sum(nAsh(1:nSym))
  nAtri = nTri_Elem(nm)
  nAtri = nTri_Elem(nAtri)
  call RecPrt('MO1',' ',MO1,1,nAtri)
# endif
  FockI(:) = FockI(:)+Int1(1:nDens)
  Fock(:) = Zero

  do iS=1,nSym
    if (nOrb(iS) == 0) cycle

    Fock(ipCM(iS):ipCM(iS)+nOrb(iS)*nIsh(is)-1) = Two*FockI(ipCM(iS):ipCM(iS)+nOrb(iS)*nIsh(is)-1)
    if (iMethod == 2) then
      Fock(ipCM(iS):ipCM(iS)+nOrb(iS)*nIsh(is)-1) = Fock(ipCM(iS):ipCM(iS)+nOrb(iS)*nIsh(is)-1)+ &
                                                    Two*FockA(ipCM(iS):ipCM(iS)+nOrb(iS)*nIsh(is)-1)
      Fock(ipCM(iS)+nIsh(is)*nOrb(is):ipCM(iS)+(nIsh(is)+nAsh(is))*nOrb(is)-1) = &
        Q(ipMatba(iS,is):ipMatba(iS,is)+nOrb(iS)*nAsh(is)-1)
      do iAsh=1,nAsh(is)
        ipi = ipCM(iS)+nOrb(is)*(nIsh(is)+iAsh-1)
        do jAsh=1,nAsh(is)
          ipj = ipCM(iS)+nOrb(is)*(nIsh(is)+jAsh-1)
          ni = nA(is)+iAsh
          nj = nA(is)+jAsh
          ipD = iTri(ni,nj)
          Fock(ipj:ipj+nOrb(is)-1) = Fock(ipj:ipj+nOrb(is)-1)+G1t(ipD)*FockI(ipi:ipi+nOrb(is)-1)
        end do
      end do
    end if

  end do

end if
renergy = Zero
rcora = Zero
do iS=1,nSym
  do iB=1,nAsh(is)+nIsh(is)
    rEnergy = rEnergy+Fock(ipCM(is)+nOrb(iS)*(iB-1)+iB-1)
  end do
end do
rcorei = Zero
rcorea = Zero
rcor = Zero
do iS=1,nSym
  iptmp = ipCM(iS)
  do iB=1,nIsh(is)
    rcorei = rcorei+Two*Int1(iptmp)
    rcor = rcor+Two*Focki(iptmp)
    iptmp = iptmp+nOrb(iS)+1
  end do

  do iB=1,nAsh(iS)
    do jB=1,nAsh(iS)
      iiB = nA(iS)+ib
      ijB = nA(iS)+jb
      iij = iTri(iib,ijb)
      iiB = nIsh(iS)+ib
      ijB = nIsh(iS)+jb
      rcorea = rcorea+G1t(iij)*Int1(ipCM(is)-1+nOrb(is)*(iib-1)+ijB)

      rcora = rcora+G1t(iij)*Focki(ipCM(is)+nOrb(is)*(iib-1)+ijB-1)
    end do
  end do
end do
rin_ene = Half*(rcor+rcorei)
rcore = rCorei+rcoreA
if (debug) then
  write(u6,*) 'Checking energy',Half*renergy+potnuc+Half*rcore
  write(u6,*) 'Checking energy',Half*renergy,potnuc,Half*rcore
  write(u6,*)
end if

if (TwoStep .and. (StepType == 'RUN1')) then
  iaddressQDAT = 0
  nm = sum(nAsh(1:nSym))
  nAtri = nTri_Elem(nTri_Elem(nm))
  call ddafile(LuQDAT,1,FockA,nDens,iaddressQDAT)
  call ddafile(LuQDAT,1,FockI,nDens,iaddressQDAT)
  call ddafile(LuQDAT,1,Fock,nDens,iaddressQDAT)
  call ddafile(LuQDAT,1,Q,nDens,iaddressQDAT)
  call ddafile(LuQDAT,1,MO1,nAtri,iaddressQDAT)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Read22_2
