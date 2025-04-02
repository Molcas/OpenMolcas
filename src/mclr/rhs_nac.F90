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

subroutine RHS_NAC(Fock,SLag)

use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: ipin, ipnout, opout, W
use MCLR_Data, only: ipCI, nConf1, ipMat, n1Dens, n2Dens, nDens2
use MCLR_Data, only: NSSA
use MCLR_Data, only: XISPSM
use CandS, only: ICSM, ISSM
use input_mclr, only: ntAsh, PT2, State_Sym, nSym, nRoots, nConf, nBas, nCSF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp

implicit none
real*8 Fock(*)
real*8 :: SLag(*)
integer ng1, ng2, i, j, k, l, ij, kl, ijkl, ij2, kl2, ijkl2
integer nConfL, nConfR, iRC, LuDens
real*8 factor
integer iSLag !,jR,kR
real*8, allocatable :: G1q(:), G1m(:), G1r(:), G2q(:), G2r(:), CIL(:), CIR(:), T(:), F(:)

!                                                                      *
!***********************************************************************
!                                                                      *
ng1 = nTri_Elem(ntAsh)
ng2 = nTri_Elem(ng1)
if (PT2) then
  call mma_allocate(G1q,n1dens,Label='G1q')
  call mma_allocate(G2q,n2dens,Label='G2q')
else
  call mma_allocate(G1q,nG1,Label='G1q')
  call mma_allocate(G2q,nG2,Label='G2q')
end if
call mma_allocate(G1m,nG1,Label='G1m')
call mma_allocate(G1r,n1dens,Label='G1r')
G1r(:) = Zero
call mma_allocate(G2r,n2dens,Label='G2r')
G2r(:) = Zero

! Calculate one- and two-particle transition matrices
! from the CI vectors of the two NAC states
! (code copied from CIdens_SA, same symmetry)

nConfL = max(nconf1,nint(xispsm(State_sym,1)))
nConfR = max(nconf1,nint(xispsm(State_sym,1)))
call mma_allocate(CIL,nConfL,Label='CIL')
call mma_allocate(CIR,nConfR,Label='CIR')
if (PT2) then
  call PT2_SLag()
else
  call ipIn(ipCI)
  call CSF2SD(W(ipCI)%A(1+(NSSA(2)-1)*nconf1),CIL,State_sym)
  call opout(ipCI)
  call CSF2SD(W(ipCI)%A(1+(NSSA(1)-1)*nconf1),CIR,State_sym)
  call opout(ipCI)
  call ipnout(-1)
  icsm = 1
  issm = 1
  call Densi2_mclr(2,G1r,G2r,CIL,CIR,0,0,0,n1dens,n2dens)
end if
call mma_deallocate(CIL)
call mma_deallocate(CIR)

! Symmetrize densities
! For the one-particle density, save the antisymmetric part too

ij = 0
do i=0,ntAsh-1
  do j=0,i-1
    ij = ij+1
    G1q(ij) = (G1r(1+i*ntAsh+j)+G1r(1+j*ntAsh+i))*Half
    ! Note that the order of subtraction depends on how the matrix
    ! will be used when contracting with derivative integrals
    ! This is found to give the correct results:
    G1m(ij) = (G1r(1+j*ntAsh+i)-G1r(1+i*ntAsh+j))*Half
  end do
  ij = ij+1
  G1q(ij) = G1r(1+i*ntAsh+i)
  G1m(ij) = Zero
end do

! The anti-symmetric RDM is contructed somewhere in the CASPT2
! module. It will be read from disk in out_pt2.f.
if (PT2) call DCopy_(ng1,[Zero],0,G1m,1)

do i=1,ntAsh**2
  j = nTri_Elem(i)
  G2r(j) = Half*G2r(j)
end do
do i=0,ntAsh-1
  do j=0,i-1
    ij = nTri_Elem(i)+j
    do k=0,ntAsh-1
      do l=0,k
        kl = nTri_Elem(k)+l
        if (ij >= kl) then
          factor = Quart
          if (ij == kl) factor = Half
          ijkl = nTri_Elem(ij)+kl
          ij2 = i*ntAsh+j
          kl2 = k*ntAsh+l
          G2q(1+ijkl) = factor*G2r(1+nTri_Elem(ij2)+kl2)
          ij2 = max(j*ntAsh+i,l*ntAsh+k)
          kl2 = min(j*ntAsh+i,l*ntAsh+k)
          G2q(1+ijkl) = G2q(1+ijkl)+factor*G2r(1+nTri_Elem(ij2)+kl2)
          if (k /= l) then
            ij2 = i*ntAsh+j
            kl2 = l*ntAsh+k
            G2q(1+ijkl) = G2q(1+ijkl)+factor*G2r(1+nTri_Elem(ij2)+kl2)
            if (ij /= kl) then
              ij2 = max(j*ntAsh+i,k*ntAsh+l)
              kl2 = min(j*ntAsh+i,k*ntAsh+l)
              G2q(1+ijkl) = G2q(1+ijkl)+factor*G2r(1+nTri_Elem(ij2)+kl2)
            end if
          end if
        end if
      end do
    end do
  end do
  ij = nTri_Elem(i+1)-1
  do k=0,ntAsh-1
    do l=0,k
      kl = nTri_Elem(k)+l
      if (ij >= kl) then
        factor = Half
        if (ij == kl) factor = One
        ijkl = nTri_Elem(ij)+kl
        ij2 = i*ntAsh+i
        kl2 = k*ntAsh+l
        G2q(1+ijkl) = factor*G2r(1+nTri_Elem(ij2)+kl2)
        if (k /= l) then
          kl2 = l*ntAsh+k
          G2q(1+ijkl) = G2q(1+ijkl)+factor*G2r(1+nTri_Elem(ij2)+kl2)
        end if
      end if
    end do
  end do
end do

! Write the symmetric densities in the runfile and
! the antisymmetric one-particle in MCLRDENS

iRC = 0
LuDens = 20
call DaName(LuDens,'MCLRDENS')
call dDaFile(LuDens,1,G1m,ng1,iRC)
call DaClos(LuDens)
call Put_dArray('D1mo',G1q,ng1)
call Put_dArray('P2mo',G2q,ng2)

! Store transition Fock matrix

do i=1,ntAsh
  do j=1,ntAsh
    G1r(ntAsh*(j-1)+i) = G1q(iTri(i,j))
  end do
end do
do i=1,ntAsh
  do j=1,ntAsh
    ij = iTri(i,j)
    ij2 = ntAsh*(i-1)+j
    do k=1,ntAsh
      do l=1,ntAsh
        kl = iTri(k,l)
        kl2 = ntAsh*(k-1)+l
        factor = One
        if ((ij >= kl) .and. (k == l)) factor = Two
        if ((ij < kl) .and. (i == j)) factor = Two
        ijkl = iTri(ij,kl)
        ijkl2 = iTri(ij2,kl2)
        G2r(ijkl2) = factor*G2q(ijkl)
      end do
    end do
  end do
end do

! Note: 1st arg = zero for no inactive density (TDM)
call mma_allocate(T,nDens2,Label='T')
call mma_allocate(F,nDens2,Label='F')
call FockGen(Zero,G1r,G2r,T,Fock,1)
call TCMO(T,1,-2)
ij = 0
do k=1,nSym
  do i=0,nBas(k)-1
    do j=0,i-1
      ij = ij+1
      F(ij) = T(ipMat(k,k)+nBas(k)*j+i)+T(ipMat(k,k)+nBas(k)*i+j)
    end do
    ij = ij+1
    F(ij) = T(ipMat(k,k)+nBas(k)*i+i)
  end do
end do
call Put_dArray('FockOcc',F,nDens2)

call mma_deallocate(T)
call mma_deallocate(F)
call mma_deallocate(G1r)
call mma_deallocate(G2r)
call mma_deallocate(G2q)
call mma_deallocate(G1m)
call mma_deallocate(G1q)

contains

subroutine PT2_SLag()

  ! Almost the same to the subroutine in rhs_sa,
  ! but slightly modified

  integer jR, kR
  real*8 vSLag

  !iR = iRLXRoot
  do jR=1,nRoots
    do kR=1,jR
      vSLag = Zero
      !write(u6,*) 'jr,kr= ',jr,kr
      !write(u6,*) vslag
      iSLag = jR+nRoots*(kR-1)
      vSLag = SLag(iSLag)
      !write(u6,*) vslag

      call CSF2SD(W(ipCI)%A(1+(jR-1)*nconf1),CIL,1)
      !call opout(ipCI)
      call CSF2SD(W(ipCI)%A(1+(kR-1)*nconf1),CIR,1)
      !call opout(ipCI)
      !call ipnout(-1)
      !icsm = 1
      !issm = 1

      if (abs(vSLag) > 1.0e-10_wp) then
        call Densi2_mclr(2,G1q,G2q,CIL,CIR,0,0,0,n1dens,n2dens)
        call DaXpY_(n1dens,vSLag,G1q,1,G1r,1)
        call DaXpY_(n2dens,vSLag,G2q,1,G2r,1)
      end if

      if (kR /= jR) then
        iSLag = kR+nRoots*(jR-1)
        vSLag = SLag(iSLag)
        if (abs(vSLag) > 1.0e-10_wp) then
          call Densi2_mclr(2,G1q,G2q,CIR,CIL,0,0,0,n1dens,n2dens)
          call DaXpY_(n1dens,vSLag,G1q,1,G1r,1)
          call DaXpY_(n2dens,vSLag,G2q,1,G2r,1)
        end if
      end if
    end do
  end do
  nConf = ncsf(1) !! nconf is overwritten somewhere in densi2_mclr

end subroutine PT2_SLag

end subroutine RHS_NAC
