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

subroutine GetWFFock_NAC(FOccMO,bk,R,nTri,P2MOt,nG2)
! Partially readapted from rhs_sa

use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: W
use MCLR_Data, only: ipCI, LuJob, NACSTATES, nConf1, nDens, nNA, XISPSM
use input_mclr, only: iTOC, nCSF, nRoots, ntAsh, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTri, nG2
real(kind=wp), intent(out) :: FOccMO(nDens), bk(nDens), P2MOt(nG2)
real(kind=wp), intent(in) :: R(nRoots,nRoots)
integer(kind=iwp) :: I, iB, iDij, iDkl, iIJKL, IJ, IJ2, IJKL, ijkl2, iRC, iRij, iRijkl, iRkl, J, jB, jDisk, K, kB, KL, KL2, lB, &
                     LuDens, nConfL, nConfR, NCSFs, nG1
real(kind=wp) :: Fact, factor, rdum(1)
real(kind=wp), allocatable :: CIL(:), CIR(:), D5(:), D6(:), FinCI(:), Fock(:), G1m(:), G1q(:), G1qs(:,:), G1r(:), G2q(:), &
                              G2qs(:,:), G2r(:), G2rt(:), T(:)

! FinCI: CI Vectors in final CMS state basis

!                                                                      *
!***********************************************************************
!                                                                      *
ng1 = nTri_Elem(ntash)
if (ng2 /= nTri_Elem(ng1)) call Abend()

call mma_allocate(FinCI,nconf1*nroots,Label='FinCI')
call mma_allocate(Fock,nDens,Label='Fock')
call mma_allocate(T,nDens,Label='T')
call mma_allocate(G1q,ng1,Label='G1q')
call mma_allocate(G1m,ng1,Label='G1m')
call mma_allocate(G2q,ng2,Label='G2q')
call mma_allocate(G1r,ntash**2,Label='G1r')
call mma_allocate(G2r,nTri_Elem(ntash**2),Label='G2r')
call mma_allocate(G2rt,nTri_Elem(ntash**2),Label='G2rt')
! Rotate CI vectors back to those for reference states
NCSFs = NCSF(state_sym)
call DGEMM_('n','n',NCSFS,nRoots,nRoots,One,W(ipCI)%A,NCSFs,R,nRoots,Zero,FinCI,nCSFs)
nConfL = max(ncsf(state_sym),nint(xispsm(state_sym,1)))
nConfR = max(ncsf(state_sym),nint(xispsm(state_sym,1)))

call mma_allocate(CIL,nConfL)
call mma_allocate(CIR,nConfR)

I = NACstates(1)
J = NACstates(2)
call CSF2SD(FinCI(1+(J-1)*NCSFs),CIL,state_sym)
call CSF2SD(FinCI(1+(I-1)*NCSFs),CIR,state_sym)
call Densi2_mclr(2,G1r,G2rt,CIL,CIR,0,0,0,ntash**2,nTri_Elem(ntash**2))

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
    G1r(ib+(jb-1)*ntash) = G1q(iTri(ib,jb))
  end do
end do

do iB=1,ntAsh**2
  jB = nTri_Elem(iB)
  G2rt(jB) = Half*G2rt(jB)
end do
do iB=0,ntAsh-1
  do jB=0,iB-1
    ij = nTri_Elem(iB)+jB
    do kB=0,ntAsh-1
      do lB=0,kB
        kl = nTri_Elem(kB)+lB
        if (ij >= kl) then
          factor = Quart
          if (ij == kl) factor = Half
          ijkl = nTri_Elem(ij)+kl
          ij2 = iB*ntAsh+jB
          kl2 = kB*ntAsh+lB
          G2q(1+ijkl) = factor*G2rt(1+nTri_Elem(ij2)+kl2)
          ij2 = max(jB*ntAsh+iB,lB*ntAsh+kB)
          kl2 = min(jB*ntAsh+iB,lB*ntAsh+kB)
          G2q(1+ijkl) = G2q(1+ijkl)+factor*G2rt(1+nTri_Elem(ij2)+kl2)
          if (kB /= lB) then
            ij2 = iB*ntAsh+jB
            kl2 = lB*ntAsh+kB
            G2q(1+ijkl) = G2q(1+ijkl)+factor*G2rt(1+nTri_Elem(ij2)+kl2)
            if (ij /= kl) then
              ij2 = max(jB*ntAsh+iB,kB*ntAsh+lB)
              kl2 = min(jB*ntAsh+iB,kB*ntAsh+lB)
              G2q(1+ijkl) = G2q(1+ijkl)+factor*G2rt(1+nTri_Elem(ij2)+kl2)
            end if
          end if
        end if
      end do
    end do
  end do
  ij = nTri_Elem(iB)+iB
  do kB=0,ntAsh-1
    do lB=0,kB
      kl = nTri_Elem(kB)+lB
      if (ij >= kl) then
        factor = Half
        if (ij == kl) factor = One
        ijkl = nTri_Elem(ij)+kl
        ij2 = iB*ntAsh+iB
        kl2 = kB*ntAsh+lB
        G2q(1+ijkl) = factor*G2rt(1+nTri_Elem(ij2)+kl2)
        if (kB /= lB) then
          kl2 = lB*ntAsh+kB
          G2q(1+ijkl) = G2q(1+ijkl)+factor*G2rt(1+nTri_Elem(ij2)+kl2)
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
        iijkl = iTri(iDij,iDkl)
        iRijkl = iTri(iRij,iRkl)
        G2q(iijkl) = Fact*G2r(iRijkl)
      end do
    end do
  end do
end do

call Get_dArray_chk('P2MOt',P2MOt,ng2)
P2MOt(:) = P2MOt(:)+G2q(:)

! Done with the info from CMS final state

! Doing some computation for computing non-active-active 2RDM in
! integral_util/prepp.f
call mma_allocate(D5,nTri)
call mma_allocate(D6,nTri)
! D5: Used in ptrans_sa when isym == jsym (PDFT parts cancel WF
!     parts for intermediate states)
! D6: Used in ptrans_sa when isym /= jsym (sum of inactive parts of
! intermediate-state 1RDMs cancels that of the final state)
call Get_DArray('MSPDFTD6',D6,nTri)
call GetDMatAO(G1q,D5,ng1,nTri)
D6(:) = D6(:)+D5(:)
call Put_DArray('MSPDFTD5',D5,nTri)
call Put_DArray('MSPDFTD6',D6,nTri)
call mma_deallocate(D5)
call mma_deallocate(D6)
! Beginning of the info for CMS intermediate states
jdisk = itoc(3)
call mma_allocate(G1qs,ng1,nRoots)
call mma_allocate(G2qs,ng2,nRoots)
do K=1,nRoots
  call dDaFile(LUJOB,2,G1q,ng1,jDisk)
  call dDaFile(LUJOB,0,rdum,ng1,jDisk)
  call dDaFile(LUJOB,2,G2q,Ng2,jDisk)
  call dDaFile(LUJOB,0,rdum,Ng2,jDisk)
  G1qs(:,K) = G1q(:)
  G2qs(:,K) = G2q(:)

  do iB=1,ntash
    do jB=1,ntash
      G1r(ib+(jb-1)*ntash) = G1q(iTri(ib,jb))
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
          iijkl = iTri(iDij,iDkl)
          iRijkl = iTri(iRij,iRkl)
          G2r(iRijkl) = Fact*G2q(iijkl)
        end do
      end do
    end do
  end do

  call FockGen(Zero,G1r,G2r,T,Fock,1)
  bk(:) = bk(:)-R(K,I)*R(K,J)*Fock(:)
  FOccMO(:) = FOccMO(:)-R(K,I)*R(K,J)*T(:)
  P2MOt(:) = P2MOt(:)-R(K,I)*R(K,J)*G2q(:)
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
