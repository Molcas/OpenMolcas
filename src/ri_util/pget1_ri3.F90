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
! Copyright (C) 1992,2007, Roland Lindh                                *
!***********************************************************************

subroutine PGet1_RI3(PAO,ijkl,nPAO,iCmp,iAO,iAOst,jBas,kBas,lBas,kOp,DSO,DSO_Var,nDSO,ExFac,CoulFac,PMax,V_K,U_K,mV_k,ZpK,nnP1, &
                     nSA,nAct)
!***********************************************************************
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density.                        *
!                                                                      *
!          The indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!                                                                      *
!             Modified for 3-center RI gradients, March 2007           *
!***********************************************************************

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use pso_stuff, only: AOrb, B_PT2, Gamma_On, lPSO, lSA, Thpkl
use Data_Structures, only: V1
use RI_glob, only: BklK, BMP2, CijK, CilK, CMOi, iAdrCVec, iMP2prpt, iUHF, LuBVector, LuCVector, nChOrb, nIJR, nKdens, nYmnij, &
                   tbvec, Yij, Ymnij
#ifdef _DEBUGPRINT_
use pso_stuff, only: D0, iD0Lbl
use RI_glob, only: iOff_Ymnij
#endif
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ijkl, nPAO, iCmp(4), iAO(4), iAOst(4), jBas, kBas, lBas, kOp(4), nDSO, mV_k, nnP1, nSA, nAct(0:7)
real(kind=wp), intent(out) :: PAO(ijkl,nPAO), PMax
real(kind=wp), intent(in) :: DSO(nDSO,nSA), DSO_Var(nDSO), ExFac, CoulFac, V_k(mV_k,nSA), U_k(mV_k), ZpK(nnP1,mV_K,*)
integer(kind=iwp) :: i, i2, i3, i4, iAdr, ij, ijBas, ijk, ik, ik1, ik2, il, il1, il2, ileft, imo, iMO1, iMO2, iMOleft, iMOright, &
                     indexB, Indkl, iOff1, iPAO, irc, iright, iSO, iThpkl, iVec, iVec_, j, jAOj, jC, jik, jmo, jSkip(4), jSO, &
                     jSO_off, jSOj, jSym, k, kAct, kAOk, kmo, kSO, kSOk, kSym, Kth, lAct, lAOl, lBVec, lCVec, lda, lda1, lda2, &
                     lSO, lSOl, lSym, Lth, n2J, nijkl, nik, nj(4), nj2, njk, nk, nKBas, nLBas, nnk, NumOrb(4)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iSym
#endif
real(kind=wp) :: Cpu, Cpu1, Cpu2, ExFac_, Fac, fact, Factor, temp, tmp, Wall, Wall1, Wall2
real(kind=wp), pointer :: Xki(:), Xli(:)
type(V1) :: Xki2(2), Xki3(2), Xli2(2), Xli3(2)
real(kind=wp), external :: Compute_B, dDot_

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call PrMtrx('DSO     ',[iD0Lbl],1,[1],D0)
write(u6,*)
write(u6,*) 'Distribution of Ymnij'
iSym = 1
if (nYmnij(iSym,1) > 0) then
  write(u6,*) 'iSym=',iSym
  do i=iOff_Ymnij(iSym,1)+1,iOff_Ymnij(iSym,1)+nYmnij(iSym,1)
    write(u6,*) 'Ymnij=',Ymnij(1)%A(i)
  end do
end if
write(u6,*) 'jbas,kbas,lbas=',jBas,kBas,lBas
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! DeSymP will compensate for degeneracy due to permutational
! symmetry. We will have to compensate for that here!

call CWTime(Cpu1,Wall1)

iOff1 = nBas(0)
Fac = Quart
PMax = Zero
iPAO = 0

jSym = 1
kSym = 1
lSym = Mul(jSym,kSym)
NumOrb(1) = nChOrb(kSym-1,1)

! Test if we have any exchange contribution of significance

ExFac_ = ExFac
if (ExFac_ /= 0) then

  ! Pick up the number of MOs which passed the threshold test.

  nj2 = 0
  do iSO=1,nKdens
    jSkip(iSO) = 0
    nj(iSO) = nYmnij(jSym,iSO)
    NumOrb(iSO) = nChOrb(kSym-1,iSO)

    ! If all included skip presceening.

    ! trick for skipping unnecessary overhead
    if (-nj(iSO) == NumOrb(iSO)) then
      jSkip(iSO) = 1
      nj(iSO) = NumOrb(iSO)
    end if

    ! If all excluded process only for Coulombic contributions.

    nj2 = nj2+nj(iSO)
  end do
  if ((nj2 == 0) .and. (.not. lPSO)) ExFac_ = Zero
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if ((ExFac_ /= Zero) .and. (NumOrb(1) > 0) .and. (iMP2prpt /= 2) .and. (.not. lPSO) .and. (iUHF == 0)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! HF and Hybrid DFT

  ! number of functions in the kS and lS shell

  nKBas = kBas*iCmp(3)
  nLBas = lBas*iCmp(4)

  kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
  lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)

  ! Pointers to the full list of the X_mu,i elements.

  lda = size(CMOi(1)%SB(1)%A2,1)
  ik = 1+lda*(kSO-1)
  il = 1+lda*(lSO-1)
  Xki(1:) => CMOi(1)%SB(1)%A1(ik:)
  Xli(1:) => CMOi(1)%SB(1)%A1(il:)

  ! Collect the X_mu,i which survived the prescreening.
  ! Replace the pointers above, i.e. Xki, Xli.

  if ((nj(1) <= NumOrb(1)) .and. (jSkip(1) == 0)) then

    ! Note that the X_mu,i are stored as X_i,mu!

    imo = 1
    do k=1,nj(1)
      kmo = Ymnij(1)%A(k) ! CD-MO index

      ! Pick up X_mu,i for all mu's that belong to shell k

      call dcopy_(nKBas,Xki(kmo:),NumOrb(1),Yij(imo,1,1),nj(1))

      ! Pick up X_mu,i for all mu's that belong to shell l

      call dcopy_(nLBas,Xli(kmo:),NumOrb(1),Yij(imo,2,1),nj(1))

      imo = imo+1
    end do
    ! Reset pointers!
    Xki(1:nj(1)*nKBas) => Yij(1:nj(1)*nKBas,1,1)
    Xli(1:nj(1)*nLBas) => Yij(1:nj(1)*nLBas,2,1)
  else if (nj(1) > NumOrb(1)) then
    call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
    call Abend()
  end if

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
    jSO_off = jSO-iOff1

    ! Read a block of C_kl^J

    lCVec = nIJR(kSym,lSym,1)*jBas ! Block size
    iAdr = nIJR(kSym,lSym,1)*(jSO_off-1)+iAdrCVec(jSym,kSym,1)
    call dDaFile(LuCVector(jSym,1),2,CijK,lCVec,iAdr)

    ! Extract only those C_kl^Js for which we deem k and l to
    ! belong to the shell-pair and to be of significance.

    if ((nj(1) <= NumOrb(1)) .and. (jSkip(1) == 0)) then
      ij = 1
      do j=1,nj(1)
        jmo = Ymnij(1)%A(j)
        do i=1,nj(1)
          imo = Ymnij(1)%A(i)
          jC = imo+NumOrb(1)*(jmo-1)
          call dcopy_(jBas,CijK(jC),NumOrb(1)**2,CilK(ij),nj(1)**2)
          ij = ij+1
        end do
      end do
      n2j = nj(1)**2*jBas
      CijK(1:n2j) = CilK(1:n2j)
    end if

    ! Transform according to Eq. 16 (step 4) and generate B_kl^J

    ! E(jK,m) = Sum_i C(i,jK)' * X(i,m)

    call dGEMM_('T','N',nj(1)*jBas,nKBas,nj(1),One,CijK,nj(1),Xki,nj(1),Zero,CilK,nj(1)*jBas)

    ! B(Km,n) = Sum_j E(j,Km)' * X(j,n)

    call dGEMM_('T','N',jBas*nKBas,nLBas,nj(1),One,CilK,nj(1),Xli,nj(1),Zero,BklK,jBas*nKBas)

    do i3=1,iCmp(3)
      kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
      do i4=1,iCmp(4)
        lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

        iPAO = iPAO+1
        nijkl = 0

        do lAOl=0,lBas-1
          lSOl = lSO+lAOl
          do kAOk=0,kBas-1
            kSOk = kSO+kAOk

            indexB = (kAOk+(i3-1)*kBas)*jBas+(lAOl+(i4-1)*lBas)*nKBas*jBas
            Indkl = iTri(kSOk,lSOl)

            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-iOff1
              nijkl = nijkl+1
              indexB = indexB+1

              ! Coulomb contribution: V_k(j)*D(kl)

              temp = CoulFac*V_k(jSOj,1)*DSO(Indkl,1)

              ! Exchange contribution: B(K,m,n)

              temp = temp-ExFac_*Half*BklK(indexB)

              PMax = max(PMax,abs(temp))
              PAO(nijkl,iPAO) = Fac*temp
            end do
          end do
        end do
      end do
    end do
  end do
  nullify(Xki,Xli)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((ExFac_ /= Zero) .and. (NumOrb(1) > 0) .and. (iMP2prpt /= 2) .and. (.not. lPSO) .and. (iUHF == 1)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! UHF and Hybrid UDFT

  ! number of functions in the kS and lS shell

  nKBas = kBas*iCmp(3)
  nLBas = lBas*iCmp(4)

  kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
  lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)

  ! Pointers to the full list of the X_mu,i elements.

  do iSO=1,2
    if (nIJR(kSym,lSym,iSO) /= 0) then

      lda = size(CMOi(iSO)%SB(1)%A2,1)
      ik = 1+lda*(kSO-1)
      il = 1+lda*(lSO-1)
      Xki2(iSO)%A(1:) => CMOi(iSO)%SB(1)%A1(ik:)
      Xli2(iSO)%A(1:) => CMOi(iSO)%SB(1)%A1(il:)

      ! Collect the X_mu,i which survived the prescreening.
      ! Replace the pointers above, i.e. Xki, Xli.

      if ((nj(iSO) <= NumOrb(iSO)) .and. (jSkip(iSO) == 0)) then

        ! Note that the X_mu,i are stored as X_i,mu!

        imo = 1
        do k=1,nj(iSO)
          kmo = Ymnij(iSO)%A(k) ! CD-MO index

          ! Pick up X_mu,i for all mu's that belong to shell k

          call dcopy_(nKBas,Xki2(iSO)%A(kmo:),NumOrb(iSO),Yij(imo,1,iSO),nj(iSO))

          ! Pick up X_mu,i for all mu's that belong to shell l

          call dcopy_(nLBas,Xli2(iSO)%A(kmo:),NumOrb(iSO),Yij(imo,2,iSO),nj(iSO))

          imo = imo+1
        end do
        ! Reset pointers!
        Xki2(iSO)%A(1:nj(iSO)*nKBas) => Yij(1:nj(iSO)*nKBas,1,iSO)
        Xli2(iSO)%A(1:nj(iSO)*nLBas) => Yij(1:nj(iSO)*nLBas,2,iSO)
      else if (nj(iSO) > NumOrb(iSO)) then
        call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
        call Abend()
      end if
    end if
  end do

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
    jSO_off = jSO-iOff1

    Factor = Zero

    do iSO=1,2
      if ((nIJR(kSym,lSym,iSO) /= 0) .and. (nj(iSO) /= 0)) then

        ! Read a block of C_kl^J

        lCVec = nIJR(kSym,lSym,iSO)*jBas ! Block size
        iAdr = nIJR(kSym,lSym,iSO)*(jSO_off-1)+iAdrCVec(jSym,kSym,iSO)
        call dDaFile(LuCVector(jSym,iSO),2,CijK,lCVec,iAdr)

        ! Extract only those C_kl^Js for which we deem k and l to
        ! belong to the shell-pair and to be of significance.

        if ((nj(iSO) <= NumOrb(iSO)) .and. (jSkip(iSO) == 0)) then
          ij = 1
          do j=1,nj(iSO)
            jmo = Ymnij(iSO)%A(j)
            do i=1,nj(iSO)
              imo = Ymnij(iSO)%A(i)
              jC = imo+NumOrb(iSO)*(jmo-1)
              call dcopy_(jBas,CijK(jC),NumOrb(iSO)**2,Cilk(ij),nj(iSO)**2)
              ij = ij+1
            end do
          end do
          n2j = nj(iSO)**2*jBas
          CijK(1:n2j) = CilK(1:n2j)
        end if

        ! Transform according to Eq. 16 (step 4) and generate B_kl^J

        ! E(jK,m) = Sum_i C(i,jK)' * X(i,m)

        call dGEMM_('T','N',nj(iSO)*jBas,nKBas,nj(iSO),One,CijK,nj(iSO),Xki2(iSO)%A,nj(iSO),Zero,CilK,nj(iSO)*jBas)

        ! B(Km,n) = Sum_j E(j,Km)' * X(j,n)

        call dGEMM_('T','N',jBas*nKBas,nLBas,nj(iSO),One,CilK,nj(iSO),Xli2(iSO)%A,nj(iSO),Factor,BklK,jBas*nKBas)
        Factor = One
      end if
    end do

    do i3=1,iCmp(3)
      kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
      do i4=1,iCmp(4)
        lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

        iPAO = iPAO+1
        nijkl = 0

        do lAOl=0,lBas-1
          lSOl = lSO+lAOl
          do kAOk=0,kBas-1
            kSOk = kSO+kAOk

            indexB = (kAOk+(i3-1)*kBas)*jBas+(lAOl+(i4-1)*lBas)*nKBas*jBas
            Indkl = iTri(kSOk,lSOl)

            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-iOff1
              nijkl = nijkl+1
              indexB = indexB+1

              ! Coulomb contribution: V_k(j)*D(kl)

              temp = CoulFac*V_k(jSOj,1)*DSO(Indkl,1)

              ! Exchange contribution: B(K,m,n)

              temp = temp-ExFac_*BklK(indexB)

              PMax = max(PMax,abs(temp))
              PAO(nijkl,iPAO) = Fac*temp
            end do
          end do
        end do
      end do
    end do
  end do
  do iSO=1,2
    nullify(Xki2(iSO)%A,Xli2(iSO)%A)
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((ExFac_ /= Zero) .and. (NumOrb(1) > 0) .and. (iMP2prpt /= 2) .and. lPSO .and. (.not. LSA)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! CASSCF case

  ! number of functions in the kS and lS shell

  nKBas = kBas*iCmp(3)
  nLBas = lBas*iCmp(4)

  kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
  lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)

  ! Pointers to the full list of the X_mu,i elements.

  lda = size(CMOi(1)%SB(1)%A2,1)
  ik = 1+lda*(kSO-1)
  il = 1+lda*(lSO-1)
  Xki(1:) => CMOi(1)%SB(1)%A1(ik:)
  Xli(1:) => CMOi(1)%SB(1)%A1(il:)

  ! Collect the X_mu,i which survived the prescreening.
  ! Replace the pointers above, i.e. Xki, Xli.

  if ((nj(1) <= NumOrb(1)) .and. (jSkip(1) == 0) .and. (nj(1) /= 0)) then

    ! Note that the X_mu,i are stored as X_i,mu!

    imo = 1
    do k=1,nj(1)
      kmo = Ymnij(1)%A(k)

      ! Pick up X_mu,i for all mu's that belong to shell k

      call dcopy_(nKBas,Xki(kmo:),NumOrb(1),Yij(imo,1,1),nj(1))

      ! Pick up X_mu,i for all mu's that belong to shell l

      call dcopy_(nLBas,Xli(kmo:),NumOrb(1),Yij(imo,2,1),nj(1))

      imo = imo+1
    end do
    ! Reset pointers!
    Xki(1:nj(1)*nKBas) => Yij(1:nj(1)*nKBas,1,1)
    Xli(1:nj(1)*nLBas) => Yij(1:nj(1)*nLBas,2,1)
  else if (nj(1) > NumOrb(1)) then
    call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
    call Abend()
  end if

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
    jSO_off = jSO-iOff1

    ! Read a block of C_kl^J

    lCVec = nIJR(kSym,lSym,1)*jBas ! Block size
    iAdr = nIJR(kSym,lSym,1)*(jSO_off-1)+iAdrCVec(jSym,kSym,1)
    call dDaFile(LuCVector(jSym,1),2,CijK,lCVec,iAdr)

    ! Extract only those C_kl^Js for which we deem k and l to
    ! belong to the shell-pair and to be of significance.

    if (nj(1) /= 0) then
      if ((nj(1) <= NumOrb(1)) .and. (jSkip(1) == 0)) then
        ij = 1
        do j=1,nj(1)
          jmo = Ymnij(1)%A(j)
          do i=1,nj(1)
            imo = Ymnij(1)%A(i)
            jC = imo+NumOrb(1)*(jmo-1)
            call dcopy_(jBas,CijK(jC),NumOrb(1)**2,CilK(ij),nj(1)**2)
            ij = ij+1
          end do
        end do
        n2j = nj(1)**2*jBas
        CijK(1:n2j) = CilK(1:n2j)
      end if

      ! Transform according to Eq. 16 (step 4) and generate B_kl^J

      ! E(jK,m) = Sum_i C(i,jK)' * X(i,m)

      call dGEMM_('T','N',nj(1)*jBas,nKBas,nj(1),One,CijK,nj(1),Xki,nj(1),Zero,CilK,nj(1)*jBas)

      ! B(Km,n) = Sum_j E(j,Km)' * X(j,n)

      call dGEMM_('T','N',jBas*nKBas,nLBas,nj(1),One,CilK,nj(1),Xli,nj(1),Zero,BklK,jBas*nKBas)

    else
      BklK(1:jBas*nKBas*nLBas) = Zero
    end if

    ! Active term

    do jAOj=0,jBas-1
      jSOj = jSO_off+jAOj
      do i4=1,iCmp(4)
        lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
        do lAOl=0,lBas-1
          lSOl = lSO+lAOl
          do kAct=1,nAct(kSym-1)
            tmp = ddot_(kact,Zpk(iTri(kAct,1),jSOj,1),1,AOrb(1)%SB(1)%A2(:,lSOl),1)

            do lAct=kAct+1,nAct(lSym-1)
              tmp = tmp+Zpk(iTri(lAct,kAct),jSOj,1)*AOrb(1)%SB(1)%A2(lAct,lSOl)
            end do
            CilK(kAct) = tmp
          end do

          do i3=1,iCmp(3)
            kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
            iThpkl = jAOj+(i3-1)*kBas*jBas+(lAOl+(i4-1)*lBas)*nKBas*jBas+1
            lda = size(AOrb(1)%SB(1)%A2,1)
            ik = 1+lda*(kSO-1)
            call dGeMV_('T',nAct(kSym-1),kBas,One,AOrb(1)%SB(1)%A1(ik:),nAct(kSym-1),Cilk,1,Zero,Thpkl(iThpkl),jBas)
          end do
        end do
      end do
    end do

    do i3=1,iCmp(3)
      kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
      do i4=1,iCmp(4)
        lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

        iPAO = iPAO+1
        nijkl = 0

        do lAOl=0,lBas-1
          lSOl = lSO+lAOl
          do kAOk=0,kBas-1
            kSOk = kSO+kAOk

            iThpkl = (kAOk+(i3-1)*kBas)*jBas+(lAOl+(i4-1)*lBas)*nKBas*jBas
            indexB = iThpkl
            Indkl = iTri(kSOk,lSOl)

            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-iOff1
              nijkl = nijkl+1
              indexB = indexB+1
              iThpkl = iThpkl+1

              ! Coulomb contribution: V_k(j)*D(kl)

              temp = CoulFac*V_k(jSOj,1)*DSO(Indkl,1)

              ! Exchange contribution: B(K,m,n)

              temp = temp-ExFac_*Half*Bklk(indexB)

              ! Active space contribution: Sum_p Z(p,K)*Th(p,m,n)

              temp = temp+Thpkl(iThpkl)

              PMax = max(PMax,abs(temp))
              PAO(nijkl,iPAO) = Fac*temp
            end do
          end do
        end do
      end do
    end do
  end do
  nullify(Xki,Xli)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((ExFac_ /= Zero) .and. (iMP2prpt /= 2) .and. lPSO .and. lSA) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! SA-CASSCF case

  ! number of functions in the kS and lS shell

  nKBas = kBas*iCmp(3)
  nLBas = lBas*iCmp(4)

  kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
  lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)

  ! Pointers to the full list of the X_mu,i elements.

  do iSO=1,2
    if (nIJR(kSym,lSym,iSO) /= 0) then
      iMOleft = iSO
      iMOright = iSO+2

      lda1 = size(CMOi(iSO)%SB(1)%A2,1)
      lda2 = size(CMOi(iSO+2)%SB(1)%A2,1)
      ik1 = 1+lda1*(kSO-1)
      ik2 = 1+lda2*(kSO-1)
      il1 = 1+lda1*(lSO-1)
      il2 = 1+lda2*(lSO-1)

      Xki2(iSO)%A(1:) => CMOi(iSO+2)%SB(1)%A1(ik2:)
      Xki3(iSO)%A(1:) => CMOi(iSO)%SB(1)%A1(ik1:)
      Xli2(iSO)%A(1:) => CMOi(iSO)%SB(1)%A1(il1:)
      Xli3(iSO)%A(1:) => CMOi(iSO+2)%SB(1)%A1(il2:)

      ! Collect the X_mu,i which survived the prescreening.
      ! Replace the pointers above, i.e. Xki, Xli.

      if ((nj(iMOright) <= NumOrb(iMOright)) .and. (jSkip(iMOright) == 0)) then

        ! Note that the X_mu,i are stored as X_i,mu!

        imo = 1
        do k=1,nj(iMOright)
          kmo = Ymnij(iMOright)%A(k) ! CD-MO index

          ! Pick up X_mu,i for all mu's that belong to shell k

          call dcopy_(nKBas,Xki2(iSO)%A(kmo:),NumOrb(iMOright),Yij(imo,1,iMOright),nj(iMOright))

          call dcopy_(nLBas,Xli3(iSO)%A(kmo:),NumOrb(iMOright),Yij(imo,2,iMOright),nj(iMOright))

          imo = imo+1
        end do
        ! Reset pointers!
        nk = nj(iMOright)
        Xki2(iSO)%A(1:nk*nKBas) => Yij(1:nk*nKBas,1,iMOright)
        Xli3(iSO)%A(1:nk*nLBas) => Yij(1:nk*nLBas,2,iMOright)
      else if (nj(iMOright) > NumOrb(iMOright)) then
        call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
        call Abend()
      end if

      if ((nj(iMOleft) <= NumOrb(iMOleft)) .and. (jSkip(iMOleft) == 0)) then

        imo = 1
        do k=1,nj(iMOleft)
          kmo = Ymnij(iMOleft)%A(k) ! CD-MO index

          ! Pick up X_mu,i for all mu's that belong to shell l

          call dcopy_(nLBas,Xli2(iSO)%A(kmo:),NumOrb(iMOleft),Yij(imo,2,iMOleft),nj(iMOleft))

          call dcopy_(nKBas,Xki3(iSO)%A(kmo:),NumOrb(iMOleft),Yij(imo,1,iMOleft),nj(iMOleft))
          imo = imo+1
        end do
        nk = nj(iMOleft)
        Xli2(iSO)%A(1:nk*nLBas) => Yij(1:nk*nLBas,2,iMOleft)
        Xki3(iSO)%A(1:nk*nKBas) => Yij(1:nk*nKBas,1,iMOleft)
      else if (nj(iMOleft) > NumOrb(iMOleft)) then
        call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
        call Abend()
      end if
    end if
  end do

  do i2=1,iCmp(2)

    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
    jSO_off = jSO-iOff1

    Factor = Zero

    do iSO=1,2
      iMOleft = iSO
      iMOright = iSO+2
      if ((nIJR(kSym,lSym,iSO) /= 0) .and. (nj(iMOright) /= 0) .and. (nj(iMOleft) /= 0)) then

        ! Read a block of C_kl^J

        lCVec = nIJR(kSym,lSym,iSO)*jBas ! Block size
        iAdr = nIJR(kSym,lSym,iSO)*(jSO_off-1)+iAdrCVec(jSym,kSym,iSO)
        call dDaFile(LuCVector(jSym,iSO),2,CijK,lCVec,iAdr)

        ! Extract only those C_kl^Js for which we deem k and l to
        ! belong to the shell-pair and to be of significance.

        !MGD skipped jSkip() since not used and complicated in this case
        if ((nj(iMOright) <= NumOrb(iMOright)) .or. (nj(iMOleft) <= NumOrb(iMOleft))) then
          ij = 1
          do j=1,nj(iMOleft)
            jmo = Ymnij(iMOleft)%A(j)
            do i=1,nj(iMOright)
              imo = Ymnij(iMOright)%A(i)
              jC = imo+NumOrb(iMOright)*(jmo-1)
              call dcopy_(jBas,CijK(jC),NumOrb(iMOright)*NumOrb(iMOleft),CilK(ij),nj(iMOright)*nj(iMOleft))
              ij = ij+1
            end do
          end do
          n2j = nj(iMOright)*nj(iMOleft)*jBas
          CijK(1:n2j) = CilK(1:n2j)
        end if

        ! Transform according to Eq. 16 (step 4) and generate B_kl^J

        ! E(jK,m) = Sum_i C(i,jK)' * X(i,m)

        call dGEMM_('T','N',nj(iMOleft)*jBas,nKBas,nj(iMOright),One,CijK,nj(iMOright),Xki2(iSO)%A,nj(iMOright),Zero,CilK, &
                    nj(iMOleft)*jBas)

        ! B(Km,n) = Sum_j E(j,Km)' * X(j,n)

        call dGEMM_('T','N',jBas*nKBas,nLBas,nj(iMOleft),One,CilK,nj(iMOleft),Xli2(iSO)%A,nj(iMOleft),Factor,BklK,jBas*nKBas)
        Factor = One

        ! Add transpose

        ! Transpose Cijk->Cjik
        do ijBas=1,jBas
          nnk = nj(iMOleft)*nj(iMOright)*(ijBas-1)

          do ileft=1,nj(iMOleft)
            njk = nj(iMOright)*(ileft-1)+nnk

            do iright=1,nj(iMOright)
              nik = nj(iMOleft)*(iright-1)+nnk

              ijk = iright+njk
              jik = ileft+nik

              CilK(jik) = CijK(ijk)
            end do
          end do
        end do

        ! E(iK,m) = Sum_j C(j,iK)' * X(j,m)

        call dGEMM_('T','N',nj(iMOright)*jBas,nKBas,nj(iMOleft),One,CilK,nj(iMOleft),Xki3(iSO)%A,nj(iMOleft),Zero,CijK, &
                    nj(iMOright)*jBas)

        ! B(Km,n) = Sum_j E(i,Km)' * X(i,n)

        call dGEMM_('T','N',jBas*nKBas,nLBas,nj(iMOright),One,CijK,nj(iMOright),Xli3(iSO)%A,nj(iMOright),Factor,BklK,jBas*nKBas)
      end if
    end do

    ! Active term

    Thpkl(1:jBas*nKBas*nLBas) = Zero
    do iVec=1,4
      iMO1 = 1
      iMO2 = 1
      iVec_ = iVec
      fact = One
      if (iVec == 2) iMO2 = 2
      if (iVec == 3) fact = Two
      if (iVec == 4) then
        iMO1 = 2
        iVec_ = 2
      end if

      do jAOj=0,jBas-1
        jSOj = jSO_off+jAOj
        do i4=1,iCmp(4)
          lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
          do lAOl=0,lBas-1
            !lSOl = lSO+lAOl-1
            lSOl = lSO+lAOl
            do kAct=1,nAct(kSym-1)
              tmp = ddot_(kact,Zpk(iTri(kAct,1),jSOj,iVec_),1,AOrb(iMO1)%SB(1)%A2(:,lSOl),1)
              do lAct=kAct+1,nAct(lSym-1)
                tmp = tmp+Zpk(iTri(lAct,kAct),jSOj,iVec_)*AOrb(iMO1)%SB(1)%A2(lAct,lSOl)
              end do
              CilK(kAct) = tmp
            end do

            do i3=1,iCmp(3)
              kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
              iThpkl = jAOj+(i3-1)*kBas*jBas+(lAOl+(i4-1)*lBas)*nKBas*jBas+1
              lda = size(AOrb(iMO2)%SB(1)%A2,1)
              ik = 1+lda*(kSO-1)
              call dGeMV_('T',nAct(kSym-1),kBas,fact,AOrb(iMO2)%SB(1)%A1(ik:),nAct(kSym-1),Cilk,1,One,Thpkl(iThpkl),jBas)
            end do
          end do
        end do
      end do
    end do

    do i3=1,iCmp(3)
      kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
      do i4=1,iCmp(4)
        lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

        iPAO = iPAO+1
        nijkl = 0

        do lAOl=0,lBas-1
          lSOl = lSO+lAOl
          Lth = lAOl+(i4-1)*lBas+1
          do kAOk=0,kBas-1
            kSOk = kSO+kAOk
            Kth = kAOk+(i3-1)*kBas+1

            iThpkl = (kAOk+(i3-1)*kBas)*jBas+(lAOl+(i4-1)*lBas)*nKBas*jBas
            indexB = iThpkl
            Indkl = iTri(kSOk,lSOl)

            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-iOff1
              nijkl = nijkl+1
              indexB = indexB+1
              iThpkl = iThpkl+1

              ! SA-CASSCF Coulomb contribution

              temp = CoulFac*(V_k(jSOj,1)*DSO(Indkl,2)+V_k(jSOj,2)*DSO(Indkl,1)+V_k(jSOj,3)*DSO(Indkl,4)+V_k(jSOj,4)*DSO(Indkl,3)+ &
                              V_k(jSOj,1)*DSO(Indkl,5)+V_k(jSOj,5)*DSO(Indkl,1))

              ! Exchange contribution: B(K,m,n)

              temp = temp-Factor*ExFac_*Half*BklK(indexB)

              ! Active space contribution: Sum_p Z(p,K)*Th(p,m,n)

              temp = temp+Thpkl(iThpkl)

              if (Gamma_On) temp = temp+B_PT2(jSOj,Kth,Lth) ! For CASPT2

              PMax = max(PMax,abs(temp))
              PAO(nijkl,iPAO) = Fac*temp
            end do
          end do
        end do
      end do
    end do
  end do
  do iSO=1,2
    nullify(Xki2(iSO)%A,Xki3(iSO)%A,Xli2(iSO)%A,Xli3(iSO)%A)
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((ExFac_ /= Zero) .and. (NumOrb(1) > 0) .and. (iMP2prpt == 2)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! MP2 case

  nKBas = kBas*iCmp(3)
  nLBas = lBas*iCmp(4)

  kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
  lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)

  lda = size(CMOi(1)%SB(1)%A2,1)
  ik = 1+lda*(kSO-1)
  il = 1+lda*(lSO-1)
  Xki(1:) => CMOi(1)%SB(1)%A1(ik:)
  Xli(1:) => CMOi(1)%SB(1)%A1(il:)

  if ((nj(1) <= NumOrb(1)) .and. (jSkip(1) == 0)) then
    imo = 1
    do k=1,nj(1)
      kmo = Ymnij(1)%A(k)
      call dcopy_(nKBas,Xki(kmo:),NumOrb(1),Yij(imo,1,1),nj(1))
      call dcopy_(nLBas,Xli(kmo:),NumOrb(1),Yij(imo,2,1),nj(1))
      imo = imo+1
    end do
    Xki(1:nj(1)*nKBas) => Yij(1:nj(1)*nKBas,1,1)
    Xli(1:nj(1)*nLBas) => Yij(1:nj(1)*nLBas,2,1)
  else if (nj(1) > NumOrb(1)) then
    call WarningMessage(2,'Pget1_RI3: nj > NumOrb.')
    call Abend()
  end if

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
    jSO_off = jSO-iOff1

    lCVec = nIJR(kSym,lSym,1)*jBas
    iAdr = nIJR(kSym,lSym,1)*(jSO_off-1)+iAdrCVec(jSym,kSym,1)
    call dDaFile(LuCVector(jSym,1),2,CijK,lCVec,iAdr)

    if ((nj(1) <= NumOrb(1)) .and. (jSkip(1) == 0)) then
      ij = 1
      do j=1,nj(1)
        jmo = Ymnij(1)%A(j)
        do i=1,nj(1)
          imo = Ymnij(1)%A(i)
          jC = imo+NumOrb(1)*(jmo-1)
          call dcopy_(jBas,CijK(jC),NumOrb(1)**2,CilK(ij),nj(1)**2)
          ij = ij+1
        end do
      end do
      n2j = nj(1)**2*jBas
      CijK(1:n2j) = CilK(1:n2j)
    end if

    ! C(jK,m) = sum_i C(i,jK)' * X(i,m)

    call dGEMM_('T','N',nj(1)*jBas,nKBas,nj(1),One,CijK,nj(1),Xki,nj(1),Zero,CilK,nj(1)*jBas)

    ! B(Km,n) = sum_j C(j,Km)' * X(j,n)

    call dGEMM_('T','N',jBas*nKBas,nLBas,nj(1),One,CilK,nj(1),Xli,nj(1),Zero,BklK,jBas*nKBas)

    lBVec = nBas(0)*nBas(0)*jBas
    do i=1,2
      iAdr = 1+nBas(0)*nBas(0)*(jSO_off-1)
      call dDaFile(LuBVector(i),2,Bmp2(:,i),lBVec,iAdr)
    end do

    do i3=1,iCmp(3)
      kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
      do i4=1,iCmp(4)
        lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

        iPAO = iPAO+1
        nijkl = 0

        do lAOl=0,lBas-1
          lSOl = lSO+lAOl
          do kAOk=0,kBas-1
            kSOk = kSO+kAOk

            indexB = (kAOk+(i3-1)*kBas)*jBas+(lAOl+(i4-1)*lBas)*nKBas*jBas
            Indkl = iTri(kSOk,lSOl)

            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-iOff1
              nijkl = nijkl+1
              indexB = indexB+1

              ! Coulomb contribution: V_k(j)*D(kl)

              temp = CoulFac*(V_k(jSOj,1)*DSO(Indkl,1)+U_k(jSOj)*DSO(Indkl,1)+V_k(jSOj,1)*(DSO_Var(Indkl)-DSO(indkl,1))+ &
                              Compute_B(irc,kSOk,lSOl,jAOj,iOff1,2))

              ! Exchange contribution: B(K,m,n)

              temp = temp-ExFac_*Half*(BklK(indexB)+Compute_B(irc,kSOk,lSOl,jAOj,iOff1,1))
              !temp = temp-ExFac_*Half*Compute_B(irc,kSOk,lSOl,jAOj,iOff1,1)

              PMax = max(PMax,abs(temp))
              PAO(nijkl,iPAO) = Fac*temp
            end do
          end do
        end do
      end do
    end do
  end do
  nullify(Xki,Xli)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Pure DFT or case when no exhange

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
    do i3=1,iCmp(3)
      kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
      do i4=1,iCmp(4)
        lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

        iPAO = iPAO+1
        nijkl = 0

        do lAOl=0,lBas-1
          lSOl = lSO+lAOl
          do kAOk=0,kBas-1
            kSOk = kSO+kAOk

            Indkl = iTri(kSOk,lSOl)

            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-iOff1
              nijkl = nijkl+1

              ! Coulomb contribution: V_k(j)*D(kl)

              temp = CoulFac*V_k(jSOj,1)*DSO(Indkl,1)
              !temp = Zero

              PMax = max(PMax,abs(temp))
              PAO(nijkl,iPAO) = Fac*temp
            end do
          end do
        end do
      end do
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *

if (iPAO /= nPAO) then
  write(u6,*) ' Error in PGet1_RI3!'
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt(' In PGet1_RI3:PAO ',' ',PAO,ijkl,nPAO)
do i=1,ijkl
  write(u6,*) DDot_(nPAO,PAO(i,1),ijkl,PAO(i,1),ijkl)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(Cpu2,Wall2)
Cpu = Cpu2-Cpu1
Wall = Wall2-Wall1
tbvec(1) = tbvec(1)+Cpu
tbvec(2) = tbvec(2)+Wall

return

end subroutine PGet1_RI3
