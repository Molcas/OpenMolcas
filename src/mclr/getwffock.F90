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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine GetWFFock(FOccMO,bk,R,nTri,P2MOt,NG2)
! Partially readapted from rhs_sa

use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: W
use MCLR_Data, only: ipCI, IRLXROOT, LuJob, nConf1, nDens, nNA, XISPSM
use input_mclr, only: iTOC, nCSF, nRoots, ntAsh, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTri, nG2
real(kind=wp), intent(out) :: FOccMO(nDens), bk(nDens), P2MOt(nG2)
real(kind=wp), intent(in) :: R(nRoots,nRoots)
integer(kind=iwp) :: I, iA, iB, iDij, iDkl, iIJKL, ij1, iRij, iRijkl, iRkl, jA, jB, JDisk, K, kA, kB, kl1, kl2, lA, lB, nConfL, &
                     NCSFs, nG1
real(kind=wp) :: Fact, rdum(1)
real(kind=wp), allocatable :: CIL(:), CIR(:), D5(:), D6(:), DIAO(:), DMatAO(:), FinCI(:), Fock(:), G1q(:), G1qs(:), G1r(:), &
                              G2q(:), G2qs(:), G2r(:), G2rt(:), T(:)

!                                                                      *
!***********************************************************************
!                                                                      *
ng1 = nTri_Elem(ntash)
if (ng2 /= nTri_Elem(ng1)) call Abend()

call mma_allocate(FinCI,nconf1*nroots,Label='FinCI')
call mma_allocate(Fock,nDens,Label='Fock')
call mma_allocate(T,nDens,Label='T')
call mma_allocate(G1q,ng1,Label='G1q')
call mma_allocate(G2q,ng2,Label='G2q')
call mma_allocate(G1r,ntash**2,Label='G1r')
call mma_allocate(G2r,nTri_Elem(ntash**2),Label='G2r')
call mma_allocate(G2rt,nTri_Elem(ntash**2),Label='G2rt')
! Rotate CI vectors back to those for reference states
NCSFs = NCSF(state_sym)
call DGEMM_('n','n',NCSFS,nRoots,nRoots,One,W(ipCI)%A,NCSFs,R,nRoots,Zero,FinCI,nCSFs)
nConfL = max(ncsf(state_sym),nint(xispsm(state_sym,1)))

call mma_allocate(CIL,nConfL)
call mma_allocate(CIR,nConfL)

I = IRlxRoot
call CSF2SD(FinCI(1+(I-1)*NCSFs),CIL,state_sym)
CIR(:) = CIL(:)
call Densi2_mclr(2,G1r,G2rt,CIL,CIR,0,0,0,ntash**2,nTri_Elem(ntash**2))
do iA=1,nnA
  do jA=1,nnA
    do kA=1,nnA
      do la=1,nnA
        ij1 = nnA*(iA-1)+ja
        !ij2 = nna*(ja-1)+ia
        kl1 = nnA*(ka-1)+la
        kl2 = nna*(la-1)+ka
        if ((iA == jA) .or. (kA == la)) then
          G2r(iTri(ij1,kl1)) = G2rt(iTri(ij1,kl1))
        else
          G2r(iTri(ij1,kl1)) = (G2rt(iTri(ij1,kl1))+G2rt(iTri(ij1,kl2)))*Half
        end if
      end do
    end do
  end do
end do
call FockGen(One,G1r,G2r,FOccMO,bk,1)

do iB=1,ntash
  do jB=1,iB
    G1q(iTri(ib,jb)) = G1r(ib+(jb-1)*ntash)
  end do
end do
! D1MOt: CMS-PDFT 1RDM for computing 1-electron gradient
call Put_DArray('D1MOt',G1q,ng1)
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
call mma_allocate(DMatAO,nTri)
call mma_allocate(DIAO,nTri)
call Get_DArray('MSPDFTD5',DIAO,nTri)
call Get_DArray('MSPDFTD6',D6,nTri)
call GetDMatAO(G1q,DMatAO,ng1,nTri)
D6(:) = D6(:)+DMatAO(:)
D5(:) = DMatAO(:)+Half*DIAO(:)
call Put_DArray('MSPDFTD5',D5,nTri)
call Put_DArray('MSPDFTD6',D6,nTri)
call mma_deallocate(D5)
call mma_deallocate(D6)
call mma_deallocate(DMatAO)
call mma_deallocate(DIAO)
! Beginning of the info for CMS intermediate states

jdisk = itoc(3)
call mma_allocate(G1qs,ng1*nRoots)
call mma_allocate(G2qs,ng2*nRoots)
do K=1,nRoots
  call dDaFile(LUJOB,2,G1q,ng1,jDisk)
  call dDaFile(LUJOB,0,rdum,ng1,jDisk)
  call dDaFile(LUJOB,2,G2q,Ng2,jDisk)
  call dDaFile(LUJOB,0,rdum,Ng2,jDisk)
  G1qs((K-1)*ng1+1:K*ng1) = G1q(:)
  G2qs((K-1)*ng2+1:K*ng2) = G2q(:)
  call mma_allocate(DMatAO,ntri)
  call GetDMatAO(G1q,DMatAO,ng1,nTri)
  call mma_deallocate(DMatAO)
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

  call FockGen(One,G1r,G2r,T,Fock,1)
  bk(:) = bk(:)-R(K,I)**2*Fock(:)
  FOccMO(:) = FOccMO(:)-R(K,I)**2*T(:)
  P2MOt(:) = P2MOt(:)-R(K,I)**2*G2q(:)
end do
call Put_DArray('D1INTER',G1qs,ng1*nRoots)
call Put_DArray('P2INTER',G2qs,ng2*nRoots)
call mma_deallocate(G1qs)
call mma_deallocate(G2qs)
call mma_deallocate(Fock)
call mma_deallocate(T)
call mma_deallocate(G1r)
call mma_deallocate(G2r)
call mma_deallocate(G2rt)
call mma_deallocate(G1q)
call mma_deallocate(G2q)
call mma_deallocate(CIL)
call mma_deallocate(CIR)
call mma_deallocate(FinCI)

end subroutine GetWFFock
