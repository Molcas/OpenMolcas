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

subroutine Read2_ns(rMO1,rMO2,FockI,FockA,Temp1,nDens22,Temp2,Temp3,Temp4,DI13,DI24,DI,DA13,DA24,DA,rkappa,idsym,Signa,Fact,jSpin, &
                    lfat,lfit,lMOt)
!***********************************************************************
!                                        ~     ~                       *
!   Monster routine for construction of Fock, MO                       *
!   Handles non (anti) symmetric orbital rotations                     *
!                                                                      *
!   Input: rkappa : Rotation matrix                                    *
!          idsym  : symmetry of perturbation                           *
!          DIR    : Inactive One electron density                      *
!          DAR    : Inactive One electron density                      *
!                                                                      *
!   Scrtch:Temp1, Temp2,Temp3,Temp4,DR,DL                              *
!                                                                      *
!   Output:                                                            *
!         FockI(A):Fock matrix (one index transformed integrals)       *
!                  ~~                                                  *
!         rMO1    (pj|kl)\     Added together this gives oneindex      *
!                     ~~  -->  transformed integrals. They are         *
!         rMO2    (pj|kl)/     separated to make it easy to go         *
!                              to spindependent perturbations          *
!                                                                      *
!   Remember Coulomb type integrals are used to construct              *
!   exchange part Fock matrix and exchange integrals to construct      *
!   Coulomb part.                                                      *
!                                                                      *
!   Sgn  =  1                                                          *
!                                                                      *
!   Sgn  = -1  {I,K}=KI+signIK                                         *
!                                                                      *
!   jspin =  0  Fock matrixes and MO's needed for singlet              *
!               perturbations                                          *
!                                                                      *
!   jspin =  1  Fock matrixes and MO's needed for triplet              *
!               perturbations                                          *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipCM, ipMat, ipMO, nB, nCMO, nDens, nMBA
use input_mclr, only: iMethod, nAsh, nBas, nIsh, nSym
use Constants, only: Zero, Half, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: rmo1(nMBA), rmo2(nMBA), FockI(nDens), FockA(nDens)
integer(kind=iwp), intent(in) :: nDens22, iDSym, jSpin
real(kind=wp), intent(out) :: Temp1(nDens22), Temp2(nDens), Temp3(nDens), Temp4(nDens), DI13(nDens), DI24(nDens), DA13(nDens), &
                              DA24(nDens)
real(kind=wp), intent(in) :: DI(nCMO), DA(nCMO), rkappa(nDens), Signa, Fact
logical(kind=iwp), intent(in) :: lFAt, lFIT, lmot
integer(kind=iwp) :: iB, iiB, ijA, ijS, ilA, ip1, ip2, ip3, ip4, ipA, ipD, ipF, ipS, iS, jB, jjB, jS, kS, lB, lS, nNB
real(kind=wp) :: Sgn
logical(kind=iwp) :: singlet

!                                                                      *
!***********************************************************************
!                                                                      *
!  (mn|pq)=sum(o) T  (on|pq) + sgn*T  (mo|pq)+T (mn|oq) +sgn*T  (mn|po)
!                  mo               no         po             qo
!
!   DL = sum(po) D  T   C     (13)
!     bj          ij pi  bp
!
!                     *
!   DR  = sum(qo) D  T  C     (24)
!     ib           ij jp bp
!                                                                      *
!***********************************************************************
!                                                                      *
if (jspin == 1) then
  Singlet = .false.
else if (jspin == 0) then
  Singlet = .true.
else
  Singlet = .false.
  write(u6,*) 'Error jspin=/=1,0'
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!                    t
! Construct (C* kappa * D) & (C * kappa) and store it with the
! general index as the first index and the occupied as the second.
! The general index is transformed to AO index (contravariant).
!
!                                                                      *
!***********************************************************************
!                                                                      *
do iS=1,nSym
  if (nBas(iS) /= 0) then
    do jS=1,nSym
      if ((Mul(iS,jS) == idsym) .and. (nB(jS) /= 0)) then

        call DGEMM_('N','N',nBas(iS),nB(jS),nB(jS),signa,rkappa(ipMat(is,js)),nBas(iS),DI(ipCM(js)),nBas(jS),Zero, &
                    DI24(ipMat(iS,jS)),nBas(iS))

        call DGEMM_('T','N',nBas(iS),nB(jS),nB(jS),One,rkappa(ipMat(js,is)),nBas(jS),DI(ipCM(js)),nBas(jS),Zero, &
                    DI13(ipMat(iS,jS)),nBas(iS))

        if (iMethod == 2) then

          call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),signa,rkappa(ipMat(is,js)),nBas(iS),DA(ipCM(js)),nBas(jS),Zero, &
                      DA24(ipMat(iS,jS)),nBas(iS))

          call DGEMM_('T','N',nBas(iS),nBas(jS),nB(jS),One,rKappa(ipMat(js,iS)),nBas(jS),DA(ipCM(js)),nBas(js),Zero, &
                      DA13(ipMat(iS,jS)),nBas(iS))

        end if

      end if
    end do
  end if
end do

Sgn = One
!                                                                      *
!***********************************************************************
!                                                                      *
!   Read in integrals                                                  *
!                                                                      *
!***********************************************************************
!                                                                      *
do iS=1,nSym
  do jS=1,iS
    ijS = Mul(iS,jS)
    do kS=1,nSym
      do lS=1,ks
        if (nBas(iS)*nBas(js)*nBas(kS)*nBAs(ls) /= 0) then
          if (Mul(kS,lS) == ijS) then
            do iB=1,nB(iS)
              nnB = nB(jS)
              if (iS == jS) nnB = iB
              do jB=1,nnB

                call COUL(lS,kS,iS,jS,iB,jB,Temp2,Temp3)
                !*******************************************************
                !                                                      *
                !  A C T I V E    F O C K M A T R I X                  *
                !                                                      *
                !*******************************************************
                !                       ~                                                 IJ  I
                ! Fpj = sum(i,l) Dil (pl|ij) = sum(ilr) Dil Klr (pr|ij) = Fqj = sum(rI)  I  DL
                !                                                                         qr  r

                if (lFAT) then

                  if (Mul(ls,is) == idsym) then
                    ipD = ipMat(lS,iS)+nBas(lS)*(ib-1)
                    ipF = ipMat(kS,jS)+nBas(kS)*(jB-1)
                    call dGeMV_('T',nBas(lS),nBas(kS),-Fact*Sgn*Half,Temp2,nBas(lS),DA24(ipD),1,One,FockA(ipF),1)
                  end if

                  !               IJ  I
                  ! FqJ=sum(rI)  I  DL
                  !               rq  r

                  if ((kS /= ls) .and. (Mul(kS,is) == iDSym)) then
                    ipD = ipMat(kS,iS)+nBas(kS)*(ib-1)
                    ipF = ipMat(lS,jS)+nBas(lS)*(jB-1)
                    call dGeMV_('N',nBas(lS),nBas(kS),-Fact*Sgn*Half,Temp2,nBas(lS),DA24(ipD),1,One,FockA(ipF),1)
                  end if

                  !               JI  J
                  ! FqI=sum(rJ)  I  DL
                  !               qr  r

                  if ((iS /= jS) .or. (iB /= jb)) then
                    if (Mul(lS,js) == iDSym) then
                      ipD = ipMat(lS,jS)+nBas(lS)*(jb-1)
                      ipF = ipMat(kS,iS)+nBas(kS)*(iB-1)
                      call dGeMV_('T',nBas(lS),nBas(kS),-Fact*Sgn*Half,Temp2,nBas(lS),DA24(ipD),1,One,FockA(ipF),1)
                    end if

                    !                JI  J
                    ! Fqr=sum(rJ) = I  DL
                    !                rq  r

                    if ((kS /= ls) .and. (Mul(kS,js) == iDSym)) then
                      ipD = ipMat(kS,jS)+nBas(kS)*(jb-1)
                      ipF = ipMat(lS,iS)+nBas(lS)*(iB-1)
                      call dGeMV_('N',nBas(lS),nBas(kS),-Fact*Sgn*Half,Temp2,nBas(lS),DA24(ipD),1,One,FockA(ipF),1)
                    end if

                  end if
                end if
                !*******************************************************
                !                                                      *
                !  I N A C T I V E    F O C K M A T R I X              *
                !                                                      *
                !*******************************************************
                if (lFIT) then
                  !               IJ  I
                  ! Fqj=sum(rI)  I  DL
                  !               qr  r

                  if (Mul(ls,is) == idsym) then
                    ipD = ipMat(lS,iS)+nBas(lS)*(ib-1)
                    ipF = ipMat(kS,jS)+nBas(kS)*(jB-1)
                    call dGeMV_('T',nBas(lS),nBas(kS),-Fact*Sgn*Half,Temp2,nBas(lS),DI24(ipD),1,One,FockI(ipF),1)
                  end if

                  !               IJ  I
                  ! FqJ=sum(rI)  I  DL
                  !               rq  r

                  if ((kS /= ls) .and. (Mul(kS,is) == iDSym)) then
                    ipD = ipMat(kS,iS)+nBas(kS)*(ib-1)
                    ipF = ipMat(lS,jS)+nBas(lS)*(jB-1)
                    call dGeMV_('N',nBas(lS),nBas(kS),-Fact*Sgn*Half,Temp2,nBas(lS),DI24(ipD),1,One,FockI(ipF),1)
                  end if

                  !               JI  J
                  ! FqI=sum(rJ)  I  DL
                  !               qr  r

                  if ((iS /= jS) .or. (iB /= jb)) then
                    if (Mul(lS,js) == iDSym) then
                      ipD = ipMat(lS,jS)+nBas(lS)*(jb-1)
                      ipF = ipMat(kS,iS)+nBas(kS)*(iB-1)
                      call dGeMV_('T',nBas(lS),nBas(kS),-Fact*Sgn*Half,Temp2,nBas(lS),DI24(ipD),1,One,FockI(ipF),1)
                    end if

                    !                JI  J
                    ! Fqr=sum(rJ) = I  DL
                    !                rq  r

                    if ((kS /= ls) .and. (Mul(kS,js) == iDSym)) then
                      ipD = ipMat(kS,jS)+nBas(kS)*(jb-1)
                      ipF = ipMat(lS,iS)+nBas(lS)*(iB-1)
                      call dGeMV_('N',nBas(lS),nBas(kS),-Fact*Sgn*Half,Temp2,nBas(lS),DI24(ipD),1,One,FockI(ipF),1)
                    end if

                  end if
                end if
                !                                                      *
                !*******************************************************
                !                                                      *
                ! M O !!!
                !  ~           ~
                ! (pj|kl) &  (pj|kl)
                !                                                      *
                !*******************************************************
                !                                                      *
                if (LMOT .and. (jB > nIsh(js)) .and. (ib > nish(is)) .and. (nBas(kS)*nBas(lS) /= 0)) then
                  iib = ib-nIsh(is)
                  jjb = jb-nIsh(js)
                  call DGETMO(Temp2,nbas(ls),nbas(ls),nBas(kS),Temp4,nbas(ks))

                  ipS = Mul(kS,idsym)
                  ip1 = ipMO(ls,js,is)+nBas(ips)*nAsh(ls)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip2 = ipMO(ips,js,is)+nBas(ls)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip3 = ipMO(ls,is,js)+nBas(ips)*nAsh(ls)*((jjB-1)*nAsh(iS)+iib-1)
                  ip4 = ipMO(ips,is,js)+nBas(ls)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                  if (nBas(ips)*nBas(lS) > 0) &
                    call DGEMM_('N','N',nBas(ips),nBas(ls),nBas(ks),One,rKappa(ipMat(ips,ks)),nBas(ipS),Temp4,nBas(ks),Zero,Temp3, &
                                nBas(ips))
                  !  ~
                  ! (pl|ij)

                  if (nBas(ipS)*nAsh(ls) > 0) &
                    call DGEACC(Fact,Temp3(nbas(ips)*nish(ls)+1),nbas(ips),'N',rmo1(ip1),nbas(ips),nbas(ips),nash(ls))
                  if (nBas(lS)*nAsh(ipS) > 0) &
                    call DGEACC(Fact,Temp3(nish(ips)+1),nbas(ips),'T',rmo2(ip2),nbas(ls),nbas(ls),nash(ips))
                  !  ~
                  ! (pl|ji)

                  if ((is /= js) .or. (ib /= jb)) then
                    if (nBas(ipS)*nAsh(lS) > 0) &
                      call Dgeacc(Fact,Temp3(nbas(ips)*nish(ls)+1),nbas(ips),'N',rmo1(ip3),nbas(ips),nbas(ips),nash(ls))
                    if (nBas(ls)*nAsh(ipS) > 0) &
                      call Dgeacc(Fact,Temp3(nish(ips)+1),nbas(ips),'T',rmo2(ip4),nbas(ls),nbas(ls),nash(ips))
                  end if
                  if (ks /= ls) then
                    ipS = Mul(lS,idsym)
                    ip1 = ipMO(ks,js,is)+nBas(ips)*nAsh(ks)*((iiB-1)*nAsh(jS)+jjb-1)
                    ip2 = ipMO(ips,js,is)+nBas(ks)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                    ip3 = ipMO(ks,is,js)+nBas(ips)*nAsh(ks)*((jjB-1)*nAsh(iS)+iib-1)
                    ip4 = ipMO(ips,is,js)+nBas(ks)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                    if (nBas(ipS)*nBas(kS) > 0) &
                      call DGEMM_('N','T',nBas(ips),nBas(ks),nBas(ls),One,rKappa(ipmat(ips,ls)),nBas(ips),Temp4,nBas(ks),Zero, &
                                  Temp3,nBas(ips))
                    !  ~
                    ! (pk|ij)

                    if (nBas(ipS)*nAsh(ks) > 0) &
                      call Dgeacc(Fact,Temp3(nbas(ips)*nish(ks)+1),nbas(ips),'N',rmo1(ip1),nbas(ips),nbas(ips),nash(ks))
                    if (nBas(ks)*nAsh(ips) > 0) &
                      call Dgeacc(Fact,Temp3(nish(ips)+1),nbas(ips),'T',rmo2(ip2),nbas(ks),nbas(ks),nash(ips))

                    !  ~
                    ! (pk|ji)

                    if ((is /= js) .or. (ib /= jb)) then
                      if (nBas(ipS)*nAsh(ks) > 0) &
                        call DGEACC(Fact,Temp3(nbas(ips)*nish(ks)+1),nbas(ips),'N',rmo1(ip3),nbas(ips),nbas(ips),nash(ks))
                      if (nBas(ks)*nAsh(ips) > 0) &
                        call DGEACC(Fact,Temp3(nish(ips)+1),nbas(ips),'T',rmo2(ip4),nbas(ks),nbas(ks),nash(ips))
                    end if
                  end if ! kl

                  ipS = Mul(lS,idsym)
                  ip1 = ipMO(ips,js,is)+nBas(kS)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip2 = ipMO(ks,js,is)+nBas(ipS)*nAsh(ks)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip3 = ipMO(ips,is,js)+nBas(kS)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                  ip4 = ipMO(ks,is,js)+nBas(ipS)*nAsh(ks)*((jjB-1)*nAsh(iS)+iib-1)
                  if (nBas(ks)*nBas(ips) > 0) &
                    call DGEMM_('N','N',nBas(ks),nbas(ips),nBas(ls),One,Temp4,nBas(kS),rkappa(ipMat(ls,ips)),nBas(ls),Zero,Temp3, &
                                nBas(ks))
                  !   ~
                  ! (pl|ji)

                  if (nBas(ks)*nAsh(ipS) > 0) &
                    call DGEACC(Fact*Signa,Temp3(1+nbas(ks)*nish(ips)),nbas(ks),'N',rmo1(ip1),nbas(ks),nbas(ks),nash(ips))
                  if (nBas(ips)*nAsh(ks) > 0) &
                    call DGEACC(Fact*Signa,Temp3(nish(ks)+1),nbas(ks),'T',rmo2(ip2),nbas(ips),nbas(ips),nash(ks))
                  !   ~
                  ! (pl|ij)

                  if ((is /= js) .or. (ib /= jb)) then
                    if (nBas(ks)*nAsh(ips) > 0) &
                      call DGEACC(Fact*Signa,Temp3(1+nbas(ks)*nish(ips)),nbas(ks),'N',rmo1(ip3),nbas(ks),nbas(ks),nash(ips))
                    if (nBas(ips)*nAsh(ks) > 0) &
                      call DGEACC(Fact*Signa,Temp3(nish(ks)+1),nbas(ks),'T',rmo2(ip4),nbas(ips),nbas(ips),nash(ks))
                  end if

                  if (ks /= ls) then
                    ipS = Mul(kS,idsym)
                    ip1 = ipMO(ips,js,is)+nBas(ls)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                    ip2 = ipMO(ls,js,is)+nBas(ips)*nAsh(ls)*((iiB-1)*nAsh(jS)+jjb-1)
                    ip3 = ipMO(ips,is,js)+nBas(ls)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                    ip4 = ipMO(ls,is,js)+nBas(ips)*nAsh(ls)*((jjB-1)*nAsh(iS)+iib-1)
                    if (nBas(ls)*nBas(ips) > 0) &
                      call DGEMM_('T','N',nBas(ls),nBas(ips),nBas(ks),One,Temp4,nBas(ks),rKappa(ipMat(ks,ips)),nBas(ks),Zero, &
                                  Temp3,nBas(ls))
                    !   ~
                    ! (pk|ij)

                    if (nBas(ls)*nAsh(ips) > 0) &
                      call Dgeacc(Fact*Signa,Temp3(1+nbas(ls)*nish(ips)),nbas(ls),'N',rmo1(ip1),nbas(ls),nbas(ls),nash(ips))
                    if (nBas(ips)*nAsh(ls) > 0) &
                      call Dgeacc(Fact*Signa,Temp3(1+nish(ls)),nbas(ls),'T',rmo2(ip2),nbas(ips),nbas(ips),nash(ls))
                    !   ~
                    ! (pk|ji)

                    if ((is /= js) .or. (ib /= jb)) then
                      if (nBas(ls)*nAsh(ips) > 0) &
                        call DGeAcc(Fact*Signa,Temp3(1+nbas(ls)*nish(ips)),nbas(ls),'N',rmo1(ip3),nbas(ls),nbas(ls),nash(ips))
                      if (nBas(ips)*nAsh(ls) > 0) &
                        call DGeAcc(Fact*Signa,Temp3(1+nish(ls)),nbas(ls),'T',rmo2(ip4),nbas(ips),nbas(ips),nash(ls))
                    end if
                  end if ! (kl)
                  ! ^
                  ! |
                  ! NOT DEBUGGED?
                end if ! lmo

              end do
            end do
          end if
        end if
      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do iS=1,nSym
  do jS=1,nSym
    ijS = Mul(iS,jS)
    do kS=1,nSym
      do lS=1,nSym
        if (nBas(is)*nBAs(js)*nBas(ks)*nBas(ls) /= 0) then
          if (Mul(kS,lS) == ijS) then
            do lB=1,nB(lS)
              do jB=1,nB(jS)
                call EXCH(is,js,ks,ls,jb,lb,Temp1,Temp2)

                !*******************************************************
                !                                                      *
                !  I N A C T I V E    F O C K M A T R I X              *
                !                                                      *
                !*******************************************************
                if (lFIT) then
                  !           ~
                  ! Fpi=Dkl(pi|kl)=Dkl Kkr (pi|rl)=
                  !
                  !        il   l
                  ! F =   I   DR
                  !  pi    pr   r

                  if (singlet) then
                    if (Mul(ks,ls) == idsym) then
                      ipD = ipMat(kS,lS)+nBas(kS)*(lb-1)
                      ipF = ipMat(iS,jS)+nBas(iS)*(jB-1)
                      call dGeMV_('N',nBas(iS),nBas(kS),Fact,Temp1,nBas(iS),DI13(ipD),1,One,FockI(ipF),1)

                      !           ~
                      ! Fpi=Dkl(pi|kl)=Dkl Klr (pi|lr)=(Dkl Klr (pi|rl)
                      !
                      !        il   l
                      ! F =   I   DL
                      !  pi    pr   r

                      ipD = ipMat(kS,lS)+nBas(kS)*(lb-1)
                      ipF = ipMat(iS,jS)+nBas(iS)*(jB-1)
                      call dGeMV_('N',nBas(iS),nBas(kS),Fact*Sgn,Temp1,nBas(iS),DI24(ipD),1,One,FockI(ipF),1)
                    end if
                  end if
                  !           ~
                  ! Fpl=Dkj(pj|kl)=Dkj Kkr (pj|rl)
                  !
                  !             jl   j
                  !            I   DR
                  !             pr   r
                  if (Mul(ks,js) == idsym) then
                    ipD = ipMat(kS,jS)+nBas(kS)*(jb-1)
                    ipF = ipMat(iS,lS)+nBas(iS)*(lB-1)
                    call dGeMV_('N',nBas(iS),nBas(kS),-Fact*Half,Temp1,nBas(iS),DI13(ipD),1,One,FockI(ipF),1)

                  end if
                end if
                !*******************************************************
                !                                                      *
                !  A C T I V E    F O C K M A T R I X                  *
                !                                                      *
                !*******************************************************
                if (lFAT) then
                  if (singlet) then
                    !           ~
                    ! Fpi=Dkl(pi|kl)=Dkl Kkr (pi|rl)=
                    !
                    !        il   l
                    ! F =   I   DR
                    !  pi    pr   r

                    if (Mul(ks,ls) == idsym) then
                      ipD = ipMat(kS,lS)+nBas(kS)*(lb-1)
                      ipF = ipMat(iS,jS)+nBas(iS)*(jB-1)
                      call dGeMV_('N',nBas(iS),nBas(kS),Fact,Temp1,nBas(iS),DA13(ipD),1,One,FockA(ipF),1)

                      !           ~
                      ! Fpi=Dkl(pi|kl)=Dkl Krl (pi|rl)=
                      !
                      !        il   l
                      ! F =   I   DL
                      !  pi    pr   r

                      ipD = ipMat(kS,lS)+nBas(kS)*(lb-1)
                      ipF = ipMat(iS,jS)+nBas(iS)*(jB-1)
                      call dGeMV_('N',nBas(iS),nBas(kS),Fact*Sgn,Temp1,nBas(iS),DA24(ipD),1,One,FockA(ipF),1)
                    end if
                  end if
                  !           ~
                  ! Fpl=Dkj(pj|kl)=Dkj Kkr (pj|rl)
                  !
                  !             jl   j
                  !            I   DR
                  !             pr   r
                  if (Mul(kS,js) == idsym) then
                    ipD = ipMat(kS,jS)+nBas(kS)*(jb-1)
                    ipF = ipMat(iS,lS)+nBas(iS)*(lB-1)
                    call dGeMV_('N',nBas(iS),nBas(kS),-Fact*Half,Temp1,nBas(iS),DA13(ipD),1,One,FockA(ipF),1)

                  end if
                end if
                !                                                      *
                !*******************************************************
                !                                                      *
                ! M O !!!
                !     ~           ~
                ! (pj|kl) &  (pj|kl)
                !                                                      *
                !*******************************************************
                !                                                      *
                if (lmot .and. (nAsh(js)*nAsh(ls) /= 0) .and. ((jb > nish(js)) .and. (lB > nish(ls)))) then
                  Temp3(1:nBas(iS)*nBas(kS)) = Temp1(1:nBas(iS)*nBas(kS))
                  !     ~            ~
                  ! (pj|kl) &   (pj|lk)
                  !
                  !  JL                      JL
                  ! I   k    (iJ|kL)      & I   k    (iJLk)
                  !  ir  kr                  ir  kr
                  ips = Mul(kS,iDsym)
                  if (nBas(is)*nAsh(ipS) /= 0) &
                    call DGEMM_('N','T',nBas(iS),nAsh(ipS),nBas(kS),One,Temp3,nBas(iS),rKappa(ipMat(ips,ks)+nish(ips)),nBas(ips), &
                                Zero,Temp4,nBas(iS))
                  ija = jB-nIsh(jS)
                  ila = lB-nIsh(lS)
                  !     ~
                  ! (pj|kl)

                  ip2 = ipMO(js,ips,ls)+nBas(iS)*(ija-1)+nBas(is)*nAsh(js)*nAsh(ips)*(ilA-1)
                  ip3 = ipMO(js,ls,ips)+nBas(iS)*(ija-1)+nBas(is)*nAsh(js)*(ilA-1)

                  ip1 = 1
                  do ipa=1,nAsh(ips)
                    rMO1(ip2:ip2+nBas(iS)-1) = rMO1(ip2:ip2+nBas(iS)-1)+Fact*Temp4(ip1:ip1+nBas(iS)-1)
                    rMO2(ip3:ip3+nBas(iS)-1) = rMO2(ip3:ip3+nBas(iS)-1)+Fact*Temp4(ip1:ip1+nBas(iS)-1)
                    ip2 = ip2+nBas(is)*nAsh(js)
                    ip3 = ip3+nBas(is)*nAsh(js)*nash(ls)
                    ip1 = ip1+nBas(is)
                  end do
                  !      ~
                  ! (pj|lk)

                  if (nBas(iS)*nAsh(ipS) /= 0) &
                    call DGEMM_('N','N',nBas(iS),nAsh(ipS),nBas(kS),One,Temp3,nBas(iS),rKappa(ipMat(ks,ips)+nbas(ks)*nish(ips)), &
                                nBas(ks),Zero,Temp4,nBas(iS))
                  ip2 = ipMO(js,ls,ips)+nBas(iS)*(ija-1)+nBas(is)*nAsh(js)*(ilA-1)
                  ip3 = ipMO(js,ips,ls)+nBas(iS)*(ija-1)+nBas(is)*nAsh(js)*nash(ips)*(ilA-1)
                  ip1 = 1
                  do ipa=1,nAsh(ips)
                    rMO1(ip2:ip2+nBas(iS)-1) = rMO1(ip2:ip2+nBas(iS)-1)+Fact*signa*Temp4(ip1:ip1+nBas(iS)-1)
                    rMO2(ip3:ip3+nBas(iS)-1) = rMO2(ip3:ip3+nBas(iS)-1)+Fact*signa*Temp4(ip1:ip1+nBas(iS)-1)
                    ip2 = ip2+nBas(is)*nAsh(js)*nAsh(lS)
                    ip3 = ip3+nBas(is)*nAsh(js)
                    ip1 = ip1+nBas(is)
                  end do
                end if

              end do
            end do
          end if
        end if
      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Read2_ns
