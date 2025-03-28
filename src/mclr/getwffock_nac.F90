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
! Copyright (C) 2021, Paul B Calio                                     *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Based on cmsbk.f from Jie J. Bao                               *
! Additional work from  rhs_nac                                  *
! ****************************************************************

subroutine GetWFFock_NAC(FOccMO,bk,R,nTri,P2MOt,NG2)
! Partially readapted from rhs_sa

use ipPage, only: W
use MCLR_Data, only: nDens2, nConf1, ipCI, nNA
use MCLR_Data, only: NACSTATES
use MCLR_Data, only: LuJob
use MCLR_Data, only: XISPSM
use input_mclr, only: nRoots, ntAsh, State_Sym, iTOC, nCSF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Quart

implicit none
! Input
real*8, dimension(nRoots**2) :: R
integer nTri, NG2
! Output
real*8, dimension(nDens2) :: FOccMO
real*8, dimension(nDens2) :: bk
real*8, dimension(nG2) :: P2MOt
! Auxiliaries
real*8, dimension(:), allocatable :: FinCI
! FinCI: CI Vectors in final CMS state basis
real*8, dimension(1) :: rdum
real*8, dimension(:), allocatable :: Fock, T, G1r, G2r, G2rt, CIL, CIR, G1q, G2q, G1qs, G2qs, G1m
real*8, dimension(:), allocatable :: DMatAO, D5, D6
integer I, J, iTri, K, NCSFs
real*8 Fact
integer iB, jB, kB, lB, iDkl, iRijkl
integer IJ, KL, IJKL, IJ2, KL2
real*8 factor
integer nG1, nConfL, nConfR, ijkl2, iRC, LuDens, iDij, iRij, iRkl, iIJKL, jDisk
! Statement function
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *
ng1 = itri(ntash,ntash)
ng2 = itri(ng1,ng1)

call mma_allocate(FinCI,nconf1*nroots,Label='FinCI')
call mma_allocate(Fock,ndens2,Label='Fock')
call mma_allocate(T,ndens2,Label='T')
call mma_allocate(G1q,ng1,Label='G1q')
call mma_allocate(G1m,ng1,Label='G1m')
call mma_allocate(G2q,ng2,Label='G2q')
call mma_allocate(G1r,ntash**2,Label='G1r')
call mma_allocate(G2r,itri(ntash**2,ntash**2),Label='G2r')
call mma_allocate(G2rt,itri(ntash**2,ntash**2),Label='G2rt')
! Rotate CI vectors back to those for reference states
NCSFs = NCSF(state_sym)
call DGEMM_('n','n',NCSFS,nRoots,nRoots,One,W(ipCI)%Vec,NCSFs,R,nRoots,Zero,FinCI,nCSFs)
nConfL = max(ncsf(state_sym),nint(xispsm(state_sym,1)))
nConfR = max(ncsf(state_sym),nint(xispsm(state_sym,1)))

call mma_allocate(CIL,nConfL)
call mma_allocate(CIR,nConfR)

I = NACstates(1)
J = NACstates(2)
call CSF2SD(FinCI(1+(J-1)*NCSFs),CIL,state_sym)
call CSF2SD(FinCI(1+(I-1)*NCSFs),CIR,state_sym)
call Densi2_mclr(2,G1r,G2rt,CIL,CIR,0,0,0,ntash**2,itri(ntash**2,ntash**2))

! Copied from rhs_nac
ij = 0
do iB=0,nnA-1
  do jB=0,iB-1
    ij = ij+1
    G1q(ij) = (G1r(1+iB*ntAsh+jB)+G1r(1+jB*ntAsh+iB))*Half
    ! Note that the order of subtraction depends on how the matrix
    ! will be used when contracting with derivative integrals
    ! This is found to give the correct results:
    G1m(ij) = (G1r(1+jB*ntAsh+iB)-G1r(1+iB*ntAsh+jB))*Half
  end do
  ij = ij+1
  G1q(ij) = G1r(1+iB*ntAsh+iB)
  G1m(ij) = Zero
end do

do iB=1,ntash
  do jB=1,ntash
    G1r(ib+(jb-1)*ntash) = G1q(itri(ib,jb))
  end do
end do

do iB=1,ntAsh**2
  jB = itri(iB,iB)
  G2rt(jB) = Half*G2rt(jB)
end do
do iB=0,ntAsh-1
  do jB=0,iB-1
    ij = iB*(iB+1)/2+jB
    do kB=0,ntAsh-1
      do lB=0,kB
        kl = kB*(kB+1)/2+lB
        if (ij >= kl) then
          factor = Quart
          if (ij == kl) factor = Half
          ijkl = ij*(ij+1)/2+kl
          ij2 = iB*ntAsh+jB
          kl2 = kB*ntAsh+lB
          G2q(1+ijkl) = factor*G2rt(1+ij2*(ij2+1)/2+kl2)
          ij2 = max(jB*ntAsh+iB,lB*ntAsh+kB)
          kl2 = min(jB*ntAsh+iB,lB*ntAsh+kB)
          G2q(1+ijkl) = G2q(1+ijkl)+factor*G2rt(1+ij2*(ij2+1)/2+kl2)
          if (kB /= lB) then
            ij2 = iB*ntAsh+jB
            kl2 = lB*ntAsh+kB
            G2q(1+ijkl) = G2q(1+ijkl)+factor*G2rt(1+ij2*(ij2+1)/2+kl2)
            if (ij /= kl) then
              ij2 = max(jB*ntAsh+iB,kB*ntAsh+lB)
              kl2 = min(jB*ntAsh+iB,kB*ntAsh+lB)
              G2q(1+ijkl) = G2q(1+ijkl)+factor*G2rt(1+ij2*(ij2+1)/2+kl2)
            end if
          end if
        end if
      end do
    end do
  end do
  ij = iB*(iB+1)/2+iB
  do kB=0,ntAsh-1
    do lB=0,kB
      kl = kB*(kB+1)/2+lB
      if (ij >= kl) then
        factor = Half
        if (ij == kl) factor = One
        ijkl = ij*(ij+1)/2+kl
        ij2 = iB*ntAsh+iB
        kl2 = kB*ntAsh+lB
        G2q(1+ijkl) = factor*G2rt(1+ij2*(ij2+1)/2+kl2)
        if (kB /= lB) then
          kl2 = lB*ntAsh+kB
          G2q(1+ijkl) = G2q(1+ijkl)+factor*G2rt(1+ij2*(ij2+1)/2+kl2)
        end if
      end if
    end do
  end do
end do

do iB=1,ntAsh
  do jB=1,ntAsh
    ij = iTri(iB,jB)
    ij2 = ntAsh*(iB-1)+jB
    do kB=1,ntAsh
      do lB=1,ntAsh
        kl = iTri(kB,lB)
        kl2 = ntAsh*(kB-1)+lB
        factor = One
        if ((ij >= kl) .and. (kB == lB)) factor = Two
        if ((ij < kl) .and. (iB == jB)) factor = Two
        ijkl = iTri(ij,kl)
        ijkl2 = iTri(ij2,kl2)
        G2r(ijkl2) = factor*G2q(ijkl)
      end do
    end do
  end do
end do

call FockGen(Zero,G1r,G2r,FOccMO,bk,1)

! D1MOt: CMS-PDFT 1RDM for computing 1-electron gradient
call Put_DArray('D1MOt',G1q,ng1)

iRC = 0
LuDens = 20
call DaName(LuDens,'MCLRDENS')
call dDaFile(LuDens,1,G1m,ng1,iRC)
call DaClos(LuDens)

do iB=1,ntash
  do jB=1,ntash
    iDij = iTri(ib,jB)
    iRij = jb+(ib-1)*ntash
    do kB=1,ntash
      do lB=1,ntash
        iDkl = iTri(kB,lB)
        iRkl = lb+(kb-1)*ntash
        fact = One
        if ((iDij >= iDkl) .and. (kB == lB)) fact = Half
        if ((iDij < iDkl) .and. (iB == jB)) fact = Half
        iijkl = itri(iDij,iDkl)
        iRijkl = itri(iRij,iRkl)
        G2q(iijkl) = Fact*G2r(iRijkl)
      end do
    end do
  end do
end do

call Get_dArray_chk('P2MOt',P2MOt,ng2)
call DaXpY_(ng2,One,G2q,1,P2MOt,1)

! Done with the info from CMS final state

! Doing some computation for computing non-active-active 2RDM in
! integral_util/prepp.f
call mma_allocate(D5,nTri)
call mma_allocate(D6,nTri)
! D5: Used in ptrans_sa when isym == jsym (PDFT parts cancel WF
!     parts for intermediate states)
! D6: Used in ptrans_sa when isym /= jsym (sum of inactive parts of
! intermediate-state 1RDMs cancels that of the final state)
call mma_allocate(DMatAO,nTri)
call Get_DArray('MSPDFTD6',D6,nTri)
call GetDMatAO(G1q,DMatAO,ng1,nTri)
call DaXpY_(nTri,One,DMatAO,1,D6,1)
call DCopy_(nTri,DMatAO,1,D5,1)
call Put_DArray('MSPDFTD5',D5,nTri)
call Put_DArray('MSPDFTD6',D6,nTri)
call mma_deallocate(D5)
call mma_deallocate(D6)
call mma_deallocate(DMatAO)
! Beginning of the info for CMS intermediate states
jdisk = itoc(3)
call mma_allocate(G1qs,ng1*nRoots)
call mma_allocate(G2qs,ng2*nRoots)
do K=1,nRoots
  call dDaFile(LUJOB,2,G1q,ng1,jDisk)
  call dDaFile(LUJOB,0,rdum,ng1,jDisk)
  call dDaFile(LUJOB,2,G2q,Ng2,jDisk)
  call dDaFile(LUJOB,0,rdum,Ng2,jDisk)
  call dcopy_(ng1,G1q,1,G1qs((K-1)*ng1+1),1)
  call dcopy_(ng2,G2q,1,G2qs((K-1)*ng2+1),1)
  call mma_allocate(DMatAO,ntri)
  call GetDMatAO(G1q,DMatAO,ng1,nTri)
  call mma_deallocate(DMatAO)

  do iB=1,ntash
    do jB=1,ntash
      G1r(ib+(jb-1)*ntash) = G1q(itri(ib,jb))
    end do
  end do

  do iB=1,ntash
    do jB=1,ntash
      iDij = iTri(ib,jB)
      iRij = jb+(ib-1)*ntash
      do kB=1,ntash
        do lB=1,ntash
          iDkl = iTri(kB,lB)
          iRkl = lb+(kb-1)*ntash
          fact = One
          if ((iDij >= iDkl) .and. (kB == lB)) fact = Two
          if ((iDij < iDkl) .and. (iB == jB)) fact = Two
          iijkl = itri(iDij,iDkl)
          iRijkl = itri(iRij,iRkl)
          G2r(iRijkl) = Fact*G2q(iijkl)
        end do
      end do
    end do
  end do

  call FockGen(Zero,G1r,G2r,T,Fock,1)
  call Daxpy_(nDens2,-R((I-1)*nRoots+K)*R((J-1)*nRoots+K),Fock,1,bk,1)
  call Daxpy_(nDens2,-R((I-1)*nRoots+K)*R((J-1)*nRoots+K),T,1,FOccMO,1)
  call DaXpY_(ng2,-R((I-1)*nRoots+K)*R((J-1)*nRoots+K),G2q,1,P2MOt,1)
end do
call Put_DArray('D1INTER',G1qs,ng1*nRoots)
call Put_DArray('P2INTER',G2qs,ng2*nRoots)
call mma_deallocate(G1qs)
call mma_deallocate(G2qs)
call mma_deallocate(Fock)
call mma_deallocate(T)
call mma_deallocate(G1r)
call mma_deallocate(G1m)
call mma_deallocate(G2r)
call mma_deallocate(G2rt)
call mma_deallocate(G1q)
call mma_deallocate(G2q)
call mma_deallocate(CIL)
call mma_deallocate(CIR)
call mma_deallocate(FinCI)

end subroutine GetWFFock_NAC
