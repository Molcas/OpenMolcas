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

subroutine rhs_sa(Fock,SLag,ipS2)

use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: W
use MCLR_Data, only: Int1, ipCI, ipCM, ipMat, ISTATE, LuJob, n1Dens, n2Dens, nA, nConf1, nDens, nNA, XISPSM
use input_mclr, only: Debug, iRoot, iTOC, nAsh, nBas, nConf, nCSF, nIsh, nOrb, nRoots, nSym, ntAsh, PT2
use PCM_grad, only: do_RF, DSSAO, DSSMO, PCMSSAO, PCMSSMO, PCM_grad_CLag, PrepPCM2
use ISRotation, only: InvEne, InvSCF
use dmrginfo, only: DoDMRG, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Fock(nDens)
real(kind=wp), intent(in) :: SLag(nRoots,nRoots)
integer(kind=iwp), intent(in) :: ipS2
integer(kind=iwp) :: i, iB, iDij, iDkl, ii, iiB, iIJ, iIJKL, ij, ij2, ijB, ijkl, iR, iRij, iRijkl, iRkl, iS, j, jB, jDisk, jj, jR, &
                     k, kB, kl, kl2, kR, l, lB, nConfL, nConfR, nG1, nG2
real(kind=wp) :: Fact, Factor, rCoreA, rCoreI, rdum(1), rEnergy, vSLag
real(kind=wp), allocatable :: CIL(:), CIR(:), F(:), G1q(:), G1r(:), G2q(:), G2r(:), T(:)

!                                                                      *
!***********************************************************************
!                                                                      *

if (doDMRG) call dmrg_dim_change_mclr(RGras2,ntash,0)  ! yma

ng1 = nTri_Elem(ntash)
ng2 = nTri_Elem(ng1)

call mma_allocate(T,nDens,Label='T')
call mma_allocate(F,nDens,Label='F')
call mma_allocate(G1q,ng1,Label='G1q')
call mma_allocate(G2q,ng2,Label='G2q')
call mma_allocate(G1r,ntash**2,Label='G1r')
call mma_allocate(G2r,nTri_Elem(ntash**2),Label='G2r')

! Pick up densities from JobIph file

iR = iroot(istate)
jdisk = itoc(3)
do ii=1,iR-1
  call dDaFile(LUJOB,0,rdum,ng1,jDisk)
  call dDaFile(LUJOB,0,rdum,ng1,jDisk)
  call dDaFile(LUJOB,0,rdum,Ng2,jDisk)
  call dDaFile(LUJOB,0,rdum,Ng2,jDisk)
end do
call dDaFile(LUJOB,2,G1q,ng1,jDisk)
call dDaFile(LUJOB,0,rdum,ng1,jDisk)
call dDaFile(LUJOB,2,G2q,Ng2,jDisk)
call dDaFile(LUJOB,0,rdum,Ng2,jDisk)

! Add SLag (rotations of states) contributions from the partial
! derivative of the CASPT2 energy. G1q and G2q are modified.
! The modified density will be used in out_pt2.f and ptrans_sa.f
! If the target energy is not invariant, rotations that are in
! the SLag array is evaluated here (even if everything zero).
if (PT2 .or. (InvSCF .and. (.not. InvEne))) then
  do iB=1,ntash
    do jB=1,ntash
      G1r(ib+(jb-1)*ntash) = G1q(iTri(ib,jb))
    end do
  end do

  nConfL = max(nconf1,nint(xispsm(1,1)))
  nConfR = max(nconf1,nint(xispsm(1,1)))
  call mma_allocate(CIL,nConfL,Label='CIL')
  call mma_allocate(CIR,nConfR,Label='CIR')
  !iR = iRLXRoot
  do jR=1,nRoots
    call CSF2SD(W(ipCI)%A(1+(jR-1)*nconf1),CIL,1)
    do kR=1,jR ! jR-1
      vSLag = SLag(jR,kR)
      if (abs(vSLag) <= 1.0e-10_wp) cycle

      call CSF2SD(W(ipCI)%A(1+(jR-1)*nconf1),CIL,1)
      !call opout(ipCI)
      call CSF2SD(W(ipCI)%A(1+(kR-1)*nconf1),CIR,1)
      !call opout(ipCI)
      !call ipnout(-1)
      !icsm = 1
      !issm = 1
      call Densi2_mclr(2,G1r,G2r,CIL,CIR,0,0,0,n1Dens,n2Dens)
      ! For RDM1
      ij = 0
      do i=0,ntAsh-1
        do j=0,i-1
          ij = ij+1
          G1q(ij) = G1q(ij)+(G1r(1+i*ntAsh+j)+G1r(1+j*ntAsh+i))*Half*vSLag
        end do
        ij = ij+1
        G1q(ij) = G1q(ij)+G1r(1+i*ntAsh+i)*vSLag
      end do
      ! For RDM2
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
                factor = Quart*vSLag
                if (ij == kl) factor = Half*vSLag
                ijkl = nTri_Elem(ij)+kl
                ij2 = i*ntAsh+j
                kl2 = k*ntAsh+l
                G2q(1+ijkl) = G2q(1+ijkl)+factor*G2r(1+nTri_Elem(ij2)+kl2)
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
              factor = Half*vSLag
              if (ij == kl) factor = One*vSLag
              ijkl = nTri_Elem(ij)+kl
              ij2 = i*ntAsh+i
              kl2 = k*ntAsh+l
              G2q(1+ijkl) = G2q(1+ijkl)+factor*G2r(1+nTri_Elem(ij2)+kl2)
              if (k /= l) then
                kl2 = l*ntAsh+k
                G2q(1+ijkl) = G2q(1+ijkl)+factor*G2r(1+nTri_Elem(ij2)+kl2)
              end if
            end if
          end do
        end do
      end do
    end do
  end do
  call mma_deallocate(CIL)
  call mma_deallocate(CIR)
  nConf = ncsf(1)
end if

call Put_dArray('P2mo',G2q,ng2)
call Put_dArray('D1mo',G1q,ng1)

do iB=1,ntash
  do jB=1,ntash
    G1r(iB+(jB-1)*ntash) = G1q(iTri(iB,jB))
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

if (do_RF) then
  ! Update state-specific quantities using rotated G1r
  DSSMO(1:ntAsh,1:ntAsh) = reshape(G1r(1:ntAsh**2),[ntAsh,ntAsh])
  call PrepPCM2(1,DSSMO,DSSAO,PCMSSAO,PCMSSMO)
  ! compute CLag again with the rotated SS density
  ! There is no need to construct the derivative twice, but
  ! we anyway have to evaluate the derivatives of energy and
  ! Lagrangian separately, so there is no much benefit to change
  call PCM_grad_CLag(2,ipCI,ipS2)
end if

if (doDMRG) call dmrg_dim_change_mclr(RGras2,nna,0)  ! yma

call FockGen(One,G1r,G2r,T,Fock,1)
!do iS=1,nsym
!  call RecPrt(' ',' ',fock(ipMat(is,is)),nbas(is),nbas(is))
!end do

if (.not. debug) then !yma debug ??
  renergy = Zero
  do ii=1,nsym
    do jj=1,nbas(ii)
      renergy = renergy+T(ipmat(ii,ii)+jj-1+nbas(ii)*(jj-1))
    end do
  end do

  rcorei = Zero
  rcorea = Zero
  do iS=1,nSym
    do iB=1,nIsh(is)
      rcorei = rcorei+Two*Int1(ipCM(iS)+nOrb(iS)*(iB-1)+iB-1)
    end do

    do iB=1,nAsh(iS)
      do jB=1,nAsh(iS)
        iiB = nA(iS)+ib
        ijB = nA(iS)+jb
        iij = iTri(iib,ijb)
        iiB = nIsh(iS)+ib
        ijB = nIsh(iS)+jb
        rcorea = rcorea+G1q(iij)*Int1(ipCM(is)+nOrb(is)*(iib-1)+ijB-1)
      end do
    end do
  end do
  !rcore = rCorei+rcoreA
  !write(u6,*) 'In rhs_sa'
  !write(u6,*) 'Checking energy',Half*renergy+potnuc+Half*rcore !yma
  !write(u6,*) 'Checking energy',Half*renergy,potnuc,rcore      !yma
  !write(u6,*)
end if
!do iS=1,nsym
!  call RecPrt(' ',' ',fock(ipMat(is,is)),nbas(is),nbas(is))
!end do

call mma_deallocate(G1q)
call mma_deallocate(G2q)

call TCMO(T,1,-2)
ijb = 0
do is=1,nsym
  do ib=1,nbas(is)
    do jb=1,ib-1
      ijb = ijb+1
      F(ijb) = T(ipmat(is,is)+nbas(is)*(JB-1)+IB-1)+T(ipmat(is,is)+nbas(is)*(IB-1)+JB-1)
    end do
    ijb = ijb+1
    F(ijb) = T(ipmat(is,is)+nbas(is)*(iB-1)+IB-1)
  end do
end do
F(ijb+1:) = Zero
call Put_dArray('FockOcc',F,nDens)

!call recprt('RHS',' ',fock,nDens,1)

call mma_deallocate(G1r)
call mma_deallocate(G2r)
call mma_deallocate(T)
call mma_deallocate(F)

end subroutine rhs_sa
