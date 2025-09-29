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

subroutine Read2_2(rMO1,rMO2,FockI,FockA,Temp1,Temp2,Temp3,Temp4,nDens22,DI,DA,rkappa,idsym,Signa,Fact,jSpin,lfat,lfit,lMOt)
!*******************************************************************
!                                        ~     ~                   *
!   Monster routine for construction of Fock, MO                   *
!                                                                  *
!   Input: rkappa : Rotation matrix                                *
!          idsym  : symmetry of perturbation                       *
!          DIR    : Inactive One electron density                  *
!          DAR    : Inactive One electron density                  *
!                                                                  *
!   Scrtch:Temp1, Temp2,Temp3,Temp4,DR,DL                          *
!                                                                  *
!   Output:                                                        *
!         FockI(A):Fock matrix (one index transformed integrals)   *
!                  ~~                                              *
!         rMO1    (pj|kl)\     Added together this gives oneindex  *
!                     ~~  -->  transformed integrals. They are     *
!         rMO2    (pj|kl)/     separated to make it easy to go     *
!                              to spindependent perturbations      *
!                                                                  *
!   Remember Coulomb type integrals are used to construct          *
!   exchange part Fock matrix and exchange integrals to construct  *
!   Coulomb part.                                                  *
!                                                                  *
!   Sgn  =  1                                                      *
!                                                                  *
!   Sgn  = -1  {I,K}=KI+signIK                                     *
!                                                                  *
!   jspin =  0  Fock matrixes and MO's needed for singlet          *
!               perturbations                                      *
!                                                                  *
!   jspin =  1  Fock matrixes and MO's needed for triplet          *
!               perturbations                                      *
!                                                                  *
!*******************************************************************

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipCM, ipMat, ipMO, nB, nCMO, nDens, nMBA
use input_mclr, only: iMethod, nAsh, nIsh, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: rmo1(nMba), rmo2(nmba), FockI(nDens), FockA(nDens)
integer(kind=iwp), intent(in) :: nDens22, iDSym, jSpin
real(kind=wp), intent(out) :: Temp1(nDens22), Temp2(nDens22), Temp3(nDens22), Temp4(nDens22)
real(kind=wp), intent(in) :: DI(nCMO), DA(nCMO), rkappa(nDens), Signa, Fact
logical(kind=iwp), intent(in) :: lFAt, lFIT, lmot
integer(kind=iwp) :: iB, iiB, ijA, ijS, ilA, ip1, ip2, ip3, ip4, ipA, ipD, ipF, ipi, ipS, iS, jB, jjB, jS, kS, lB, lS, n, nNB
real(kind=wp) :: Sgn
logical(kind=iwp) :: singlet
real(kind=wp), allocatable :: DA13(:), DA24(:), DI13(:), DI24(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! (mn|pq)=sum(o) T  (on|pq) + sgn*T  (mo|pq)+T (mn|oq) +sgn*T  (mn|po)
!                 mo               no         po             qo
!
!   DL = sum(po) D  T   C     (13)
!     bj          ij pi  bp
!
!                     *
!   DR  = sum(qo) D  T  C     (24)
!     ib           ij jp bp

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
! Construct (C* kappa * D) & (C * kappa) and store it with
! the general index as the 1st index and the occupied as the 2nd.
! The general index is transformed to AO index (contravariant).
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(DI13,nDens,Label='DI13')
call mma_allocate(DI24,nDens,Label='DI24')
if (iMethod == 2) then
  call mma_allocate(DA13,nDens,Label='DA13')
  call mma_allocate(DA24,nDens,Label='DA24')
end if
do iS=1,nSym
  if (nOrb(iS) /= 0) then
    do jS=1,nSym
      if ((Mul(iS,jS) == idsym) .and. (nB(jS) /= 0)) then

        call DGEMM_('N','N',nOrb(iS),nB(jS),nOrb(jS),signa,rkappa(ipMat(is,js)),nOrb(iS),DI(ipCM(js)),nOrb(jS),Zero, &
                    DI24(ipMat(iS,jS)),nOrb(iS))

        call DGEMM_('T','N',nOrb(iS),nB(jS),nOrb(jS),One,rkappa(ipMat(js,is)),nOrb(jS),DI(ipCM(js)),nOrb(jS),Zero, &
                    DI13(ipMat(iS,jS)),nOrb(iS))

        if (iMethod == 2) then

          call DGEMM_('N','N',nOrb(iS),nB(jS),nOrb(jS),signa,rkappa(ipMat(is,js)),nOrb(iS),DA(ipCM(js)),nOrb(jS),Zero, &
                      DA24(ipMat(iS,jS)),nOrb(iS))

          call DGEMM_('T','N',nOrb(iS),nB(jS),nOrb(jS),One,rKappa(ipMat(js,iS)),nOrb(jS),DA(ipCM(js)),nOrb(js),Zero, &
                      DA13(ipMat(iS,jS)),nOrb(iS))

        end if
      end if
    end do
  end if
end do

Sgn = One
!                                                                      *
!***********************************************************************
!                                                                      *
!     Read in integrals                                                *
!                                                                      *
!***********************************************************************
!                                                                      *
do iS=1,nSym
  do jS=1,iS

    if (nOrb(iS)*nOrb(jS) /= 0) then

      ijS = Mul(iS,jS)
      do kS=1,nSym
        do lS=1,ks
          if (nOrb(kS)*nOrb(lS) /= 0) then

            if (Mul(kS,lS) == ijS) then
              do iB=1,nB(iS)
                nnB = nB(jS)
                if (iS == jS) nnB = iB
                do jB=1,nnB

                  call COUL(lS,kS,IS,jS,IB,JB,Temp2,Temp4)

                  !*****************************************************
                  !                                                    *
                  !  A C T I V E    F O C K M A T R I X                *
                  !                                                    *
                  !*****************************************************
                  !                       ~                                                 IJ  I
                  ! Fpj = sum(i,l) Dil (pl|ij) = sum(ilr) Dil Klr (pr|ij) = Fqj = sum(rI)  I  DL
                  !                                                                         qr  r

                  if (lFAT) then

                    if (Mul(ls,is) == idsym) then
                      ipD = ipMat(lS,iS)+nOrb(lS)*(ib-1)
                      ipF = ipMat(kS,jS)+nOrb(kS)*(jB-1)
                      call dGeMV_('T',nOrb(lS),nOrb(kS),-Fact*Sgn*Half,Temp2,nOrb(lS),DA24(ipD),1,One,FockA(ipF),1)
                    end if

                    !               IJ  I
                    ! FqJ=sum(rI)  I  DL
                    !               rq  r

                    if ((kS /= ls) .and. (Mul(kS,is) == iDSym)) then
                      ipD = ipMat(kS,iS)+nOrb(kS)*(ib-1)
                      ipF = ipMat(lS,jS)+nOrb(lS)*(jB-1)
                      call dGeMV_('N',nOrb(lS),nOrb(kS),-Fact*Sgn*Half,Temp2,nOrb(lS),DA24(ipD),1,One,FockA(ipF),1)
                    end if

                    !               JI  J
                    ! FqI=sum(rJ)  I  DL
                    !               qr  r

                    if ((iS /= jS) .or. (iB /= jb)) then
                      if (Mul(lS,js) == iDSym) then
                        ipD = ipMat(lS,jS)+nOrb(lS)*(jb-1)
                        ipF = ipMat(kS,iS)+nOrb(kS)*(iB-1)
                        call dGeMV_('T',nOrb(lS),nOrb(kS),-Fact*Sgn*Half,Temp2,nOrb(lS),DA24(ipD),1,One,FockA(ipF),1)
                      end if

                      !                JI  J
                      ! Fqr=sum(rJ) = I  DL
                      !                rq  r

                      if ((kS /= ls) .and. (Mul(kS,js) == iDSym)) then
                        ipD = ipMat(kS,jS)+nOrb(kS)*(jb-1)
                        ipF = ipMat(lS,iS)+nOrb(lS)*(iB-1)
                        call dGeMV_('N',nOrb(lS),nOrb(kS),-Fact*Sgn*Half,Temp2,nOrb(lS),DA24(ipD),1,One,FockA(ipF),1)
                      end if

                    end if
                  end if
                  !*****************************************************
                  !                                                    *
                  !  I N A C T I V E    F O C K M A T R I X            *
                  !                                                    *
                  !*****************************************************
                  if (lFIT) then
                    !               IJ  I
                    ! Fqj=sum(rI)  I  DL
                    !               qr  r

                    if (Mul(ls,is) == idsym) then
                      ipD = ipMat(lS,iS)+nOrb(lS)*(ib-1)
                      ipF = ipMat(kS,jS)+norb(kS)*(jB-1)
                      call dGeMV_('T',nOrb(lS),nOrb(kS),-Fact*Sgn*Half,Temp2,nOrb(lS),DI24(ipD),1,One,FockI(ipF),1)
                    end if

                    !               IJ  I
                    ! FqJ=sum(rI)  I  DL
                    !               rq  r

                    if ((kS /= ls) .and. (Mul(kS,is) == iDSym)) then
                      ipD = ipMat(kS,iS)+nOrb(kS)*(ib-1)
                      ipF = ipMat(lS,jS)+norb(lS)*(jB-1)
                      call dGeMV_('N',nOrb(lS),nOrb(kS),-Fact*Sgn*Half,Temp2,nOrb(lS),DI24(ipD),1,One,FockI(ipF),1)
                    end if

                    !               JI  J
                    ! FqI=sum(rJ)  I  DL
                    !               qr  r

                    if ((iS /= jS) .or. (iB /= jb)) then
                      if (Mul(lS,js) == iDSym) then
                        ipD = ipMat(lS,jS)+nOrb(lS)*(jb-1)
                        ipF = ipMat(kS,iS)+norb(kS)*(iB-1)
                        call dGeMV_('T',nOrb(lS),nOrb(kS),-Fact*Sgn*Half,Temp2,nOrb(lS),DI24(ipD),1,One,FockI(ipF),1)
                      end if

                      !                JI  J
                      ! Fqr=sum(rJ) = I  DL
                      !                rq  r

                      if ((kS /= ls) .and. (Mul(kS,js) == iDSym)) then
                        ipD = ipMat(kS,jS)+nOrb(kS)*(jb-1)
                        ipF = ipMat(lS,iS)+norb(lS)*(iB-1)
                        call dGeMV_('N',nOrb(lS),nOrb(kS),-Fact*Sgn*Half,Temp2,nOrb(lS),DI24(ipD),1,One,FockI(ipF),1)
                      end if

                    end if
                  end if

                  !*****************************************************
                  !
                  ! The rest of the Fock matrix is taken from {F,K}
                  ! and exchange type integrals
                  !
                  !*****************************************************
                  !*****************************************************
                  !
                  ! M O !!!
                  !  ~           ~
                  ! (pj|kl) &  (pj|kl)
                  !
                  !*****************************************************
                  ! NOT DEBUGGED?
                  ! |
                  ! V
                  if (LMOT .and. (jB > nIsh(js)) .and. (ib > nish(is)) .and. (nOrb(ls)*nOrb(ks) > 0)) then
                    iib = ib-nIsh(is)
                    jjb = jb-nIsh(js)

                    call DGETMO(Temp2,nOrb(ls),nOrb(ls),nOrb(kS),Temp4,nOrb(ks))

                    ipS = Mul(kS,idsym)
                    n = nOrb(ips)*nAsh(lS)
                    if (n > 0) then
                      ipI = norb(ks)*nIsh(ls)+1
                      ip1 = ipMO(ls,js,is)+n*((iiB-1)*nAsh(jS)+jjb-1)
                      ip4 = ipMO(ls,is,js)+n*((jjB-1)*nAsh(iS)+iib-1)
                      call DGEMM_('N','N',nOrb(ips),nAsh(ls),nOrb(ks),One,rKappa(ipMat(ips,ks)),nOrb(ipS),Temp4(ipI),nOrb(ks), &
                                  Zero,Temp3,nOrb(ips))
                      !  ~
                      ! (pl|ij)

                      rmo1(ip1:ip1+n-1) = rmo1(ip1:ip1+n-1)+Fact*Temp3(1:n)

                      !  ~
                      ! (pl|ji)

                      if ((is /= js) .or. (ib /= jb)) rmo1(ip4:ip4+n-1) = rmo1(ip4:ip4+n-1)+Fact*temp3(1:n)
                    end if ! nOrb(ips)*nAsh(lS)

                    if (ks /= ls) then
                      ipS = Mul(lS,idsym)
                      n = nOrb(ips)*nAsh(ks)
                      if (n > 0) then
                        ip1 = ipMO(ks,js,is)+n*((iiB-1)*nAsh(jS)+jjb-1)
                        ip4 = ipMO(ks,is,js)+n*((jjB-1)*nAsh(iS)+iib-1)
                        call DGEMM_('N','T',nOrb(ips),nAsh(ks),nOrb(ls),One,rKappa(ipmat(ips,ls)),nOrb(ips),Temp4(nIsh(ks)+1), &
                                    nOrb(ks),Zero,Temp3,nOrb(ips))
                        !  ~
                        ! (pk|ij)

                        rmo1(ip1:ip1+n-1) = rmo1(ip1:ip1+n-1)+Fact*temp3(1:n)

                        !  ~
                        ! (pk|ji)

                        if ((is /= js) .or. (ib /= jb)) rmo1(ip4:ip4+n-1) = rmo1(ip4:ip4+n-1)+Fact*temp3(1:n)
                      end if !(nOrb(ips)
                    end if ! kl

                    ipS = Mul(lS,idsym)
                    n = nOrb(ks)*nAsh(ips)
                    if (n > 0) then
                      ip2 = ipMO(ips,js,is)+n*((iiB-1)*nAsh(jS)+jjb-1)
                      ip3 = ipMO(ips,is,js)+n*((jjB-1)*nAsh(iS)+iib-1)
                      call DGEMM_('N','N',nOrb(ks),nAsh(ips),nOrb(ls),One,Temp4,nOrb(kS),rkappa(ipMat(ls,ips)+nOrb(ls)*nIsh(ips)), &
                                  nOrb(ls),Zero,Temp3,nOrb(ks))
                      !   ~
                      ! (pl|ji)

                      rmo1(ip2:ip2+n-1) = rmo1(ip2:ip2+n-1)+Fact*Signa*temp3(1:n)

                      !   ~
                      ! (pl|ij)

                      if ((is /= js) .or. (ib /= jb)) rmo1(ip3:ip3+n-1) = rmo1(ip3:ip3+n-1)+Fact*Signa*temp3(1:n)
                    end if ! nOrb(ips)

                    if (ks /= ls) then
                      ipS = Mul(kS,idsym)
                      n = nOrb(ls)*nAsh(ips)
                      if (n > 0) then
                        ip1 = ipMO(ips,js,is)+n*((iiB-1)*nAsh(jS)+jjb-1)
                        ip4 = ipMO(ips,is,js)+n*((jjB-1)*nAsh(iS)+iib-1)
                        call DGEMM_('T','N',nOrb(ls),nAsh(ips),nOrb(ks),One,Temp4,nOrb(ks), &
                                    rKappa(ipMat(ks,ips)+nOrb(ks)*nIsh(ips)),nOrb(ks),Zero,Temp3,nOrb(ls))
                        !   ~
                        ! (pk|ij)

                        rmo1(ip1:ip1+n-1) = rmo1(ip1:ip1+n-1)+Fact*Signa*temp3(1:n)

                        !   ~
                        ! (pk|ji)

                        if ((is /= js) .or. (ib /= jb)) rmo1(ip4:ip4+n-1) = rmo1(ip4:ip4+n-1)+Fact*Signa*temp3(1:n)
                      end if !(nOrb)
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

    end if

  end do
end do

!***********************************************************************

do iS=1,nSym
  do jS=1,nSym

    if (nOrb(iS)*nOrb(jS) /= 0) then

      ijS = Mul(iS,jS)
      do kS=1,nSym
        do lS=1,nSym

          if (Mul(kS,lS) == ijS) then
            if (nOrb(kS)*nOrb(lS) /= 0) then

              do lB=1,nB(lS)
                do jB=1,nB(jS)

                  call EXCH(is,js,ks,ls,jb,lb,Temp1,Temp2)

                  !*****************************************************
                  !                                                    *
                  !  I N A C T I V E    F O C K M A T R I X            *
                  !                                                    *
                  !*****************************************************
                  if (lFIT) then
                    !           ~
                    ! Fpi=Dkl(pi|kl)=Dkl Kkr (pi|rl)=
                    !
                    !        il   l
                    ! F =   I   DR
                    !  pi    pr   r

                    if (singlet) then
                      if (Mul(ks,ls) == idsym) then
                        ipD = ipMat(kS,lS)+nOrb(kS)*(lb-1)
                        ipF = ipMat(iS,jS)+norb(iS)*(jB-1)
                        call dGeMV_('N',nOrb(iS),nOrb(kS),Fact,Temp1,nOrb(iS),DI13(ipD),1,One,FockI(ipF),1)

                        !           ~
                        ! Fpi=Dkl(pi|kl)=Dkl Klr (pi|lr)=(Dkl Klr (pi|rl)
                        !
                        !        il   l
                        ! F =   I   DL
                        !  pi    pr   r

                        ipD = ipMat(kS,lS)+nOrb(kS)*(lb-1)
                        ipF = ipMat(iS,jS)+nOrb(iS)*(jB-1)
                        call dGeMV_('N',nOrb(iS),nOrb(kS),Fact*Sgn,Temp1,nOrb(iS),DI24(ipD),1,One,FockI(ipF),1)
                      end if
                    end if
                    !   ~
                    ! Fpl=Dkj(pj|kl)=Dkj Kkr (pj|rl)
                    !
                    !             jl   j
                    !            I   DR
                    !             pr   r
                    if (Mul(ks,js) == idsym) then
                      ipD = ipMat(kS,jS)+nOrb(kS)*(jb-1)
                      ipF = ipMat(iS,lS)+nOrb(iS)*(lB-1)
                      call dGeMV_('N',nOrb(iS),nOrb(kS),-Fact*Half,Temp1,nOrb(iS),DI13(ipD),1,One,FockI(ipF),1)

                    end if
                  end if
                  !*****************************************************
                  !                                                    *
                  !  A C T I V E    F O C K M A T R I X                *
                  !                                                    *
                  !*****************************************************
                  if (lFAT) then
                    if (singlet) then
                      !           ~
                      ! Fpi=Dkl(pi|kl)=Dkl Kkr (pi|rl)=
                      !
                      !        il   l
                      ! F =   I   DR
                      !  pi    pr   r

                      if (Mul(ks,ls) == idsym) then
                        ipD = ipMat(kS,lS)+nOrb(kS)*(lb-1)
                        ipF = ipMat(iS,jS)+nOrb(iS)*(jB-1)
                        call dGeMV_('N',nOrb(iS),nOrb(kS),Fact,Temp1,nOrb(iS),DA13(ipD),1,One,FockA(ipF),1)

                        !           ~
                        ! Fpi=Dkl(pi|kl)=Dkl Krl (pi|rl)=
                        !
                        !        il   l
                        ! F =   I   DL
                        !  pi    pr   r

                        ipD = ipMat(kS,lS)+nOrb(kS)*(lb-1)
                        ipF = ipMat(iS,jS)+nOrb(iS)*(jB-1)
                        call dGeMV_('N',nOrb(iS),nOrb(kS),Fact*Sgn,Temp1,nOrb(iS),DA24(ipD),1,One,FockA(ipF),1)
                      end if
                    end if
                    !           ~
                    ! Fpl=Dkj(pj|kl)=Dkj Kkr (pj|rl)
                    !
                    !             jl   j
                    !            I   DR
                    !             pr   r
                    if (Mul(kS,js) == idsym) then
                      ipD = ipMat(kS,jS)+nOrb(kS)*(jb-1)
                      ipF = ipMat(iS,lS)+nOrb(iS)*(lB-1)
                      call dGeMV_('N',nOrb(iS),nOrb(kS),-Fact*Half,Temp1,nOrb(iS),DA13(ipD),1,One,FockA(ipF),1)

                    end if
                  end if

                  !*****************************************************
                  !
                  ! M O !!!
                  !     ~           ~
                  ! (pj|kl) &  (pj|kl)
                  !
                  !*****************************************************

                  if (lmot .and. (nAsh(js)*nAsh(ls) /= 0) .and. ((jb > nish(js)) .and. (lB > nish(ls)))) then

                    Temp3(1:nOrb(iS)*nOrb(kS)) = Temp1(1:nOrb(iS)*nOrb(kS))

                    !     ~            ~
                    ! (pj|kl) &   (pj|lk)
                    !
                    !  JL                      JL
                    ! I   k    (iJ|kL)      & I   k    (iJLk)
                    !  ir  kr                  ir  kr
                    ips = Mul(kS,iDsym)
                    if (nOrb(iS)*nAsh(ipS) /= 0) &
                      call DGEMM_('N','T',nOrb(iS),nAsh(ipS),nOrb(kS),One,Temp3,nOrb(iS),rKappa(ipMat(ips,ks)+nIsh(ips)), &
                                  nOrb(ips),Zero,Temp4,nOrb(iS))
                    ija = jB-nIsh(jS)
                    ila = lB-nIsh(lS)
                    !     ~
                    ! (pj|kl)

                    ip2 = ipMO(js,ips,ls)+nOrb(iS)*(ija-1)+nOrb(is)*nAsh(js)*nAsh(ips)*(ilA-1)
                    ip1 = 1
                    do ipa=1,nAsh(ips)
                      rMO2(ip2:ip2+nOrb(iS)-1) = rMO2(ip2:ip2+nOrb(iS)-1)+Fact*Temp4(ip1:ip1+nOrb(iS)-1)
                      ip2 = ip2+nOrb(is)*nAsh(js)
                      ip1 = ip1+nOrb(is)
                    end do
                    !      ~
                    ! (pj|lk)

                    if (nOrb(is)*nAsh(ipS) /= 0) &
                      call DGEMM_('N','N',nOrb(iS),nAsh(ipS),nOrb(kS),One,Temp3,nOrb(iS),rKappa(ipMat(ks,ips)+nOrb(ks)*nIsh(ipS)), &
                                  nOrb(ks),Zero,Temp4,nOrb(iS))
                    ip3 = ipMO(js,ls,ips)+nOrb(iS)*(ija-1)+nOrb(is)*nAsh(js)*(ilA-1)
                    ip1 = 1
                    do ipa=1,nAsh(ips)
                      rMO2(ip3:ip3+nOrb(iS)-1) = rMO2(ip3:ip3+nOrb(iS)-1)+Fact*signa*Temp4(ip1:ip1+nOrb(iS)-1)
                      ip3 = ip3+nOrb(is)*nAsh(js)*nAsh(lS)
                      ip1 = ip1+nOrb(is)
                    end do
                  end if

                end do
              end do
            end if
          end if
        end do
      end do
    end if
  end do
end do

call mma_deallocate(DI13)
call mma_deallocate(DI24)
if (iMethod == 2) then
  call mma_deallocate(DA13)
  call mma_deallocate(DA24)
end if

end subroutine Read2_2
