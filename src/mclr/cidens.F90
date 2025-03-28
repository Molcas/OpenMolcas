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

subroutine CIDens(response,iLS,iRS,iL,iR,rP,rD)

use ipPage, only: W
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use MCLR_Data, only: nConf1, n1Dens, n2Dens, nNA
use MCLR_Data, only: XISPSM
use MCLR_Data, only: NOCSF
use CandS, only: ICSM, ISSM
use input_mclr, only: TimeDep, nCSF
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2

implicit none
logical Response
integer iLS, iRS, iL, iR
real*8 rP(*), rD(*)
real*8, allocatable :: De(:), Pe(:), CIL(:), CIR(:)
integer i, j, itri
integer nDim, nConfL, nConfR
integer IA, JA, KA, LA, ij1, ij2, kl1, kl2
! Statement function
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

! LS = CI
!
! Ok, we once more want to hide Jeppe's routines from
! the eyes of the world, so everyone believes that I have done
! all the work.
! If we have spin dependent Hamiltonian we will work with
! SD in all parts of the program, no CSF is necessary,
! otherwise we will do the optimization in CSF's to increase
! convergence.
!
! Input:
!
! response: true if the density should be used for response calculations
! iLS: CI Coeff for left state
! iRS: CI Coeff for right state
! iL : Symmetry of left state
! iR : Symmetry of right state
!
!           +       +
! iS=1 E  =a  a  + a  a
!       pq  ap aq   Bp Bq
!
!            +       +
! iS=-1 T  =a  a  - a a
!        pq  ap aq   Bp Bq
!
! Output:
!
!  rP : Two Electron Density
!  rD : One Electron Density

if (nconf1 == 0) return

if (doDMRG) then  !yma
  call dmrg_dim_change_mclr(LRras2(1:8),ndim,0)
  call dmrg_dim_change_mclr(LRras2(1:8),nna,0)
  n1dens = ndim**2
  n2dens = n1dens*(n1dens+1)/2
end if

call mma_allocate(De,n1dens,Label='De')
call mma_allocate(Pe,n2dens,Label='Pe')
De(:) = Zero
Pe(:) = Zero

if (nocsf == 0) then
  nConfL = max(ncsf(il),nint(xispsm(il,1)))
  nConfR = max(ncsf(iR),nint(xispsm(iR,1)))
  call mma_allocate(CIL,nConfL,Label='CIL')
  call ipin1(iLS,nconfL)
  call CSF2SD(W(iLS)%Vec,CIL,iL)
  call opout(iLs)
  call mma_allocate(CIR,nConfR,Label='CIR')
  call ipin1(iRS,nconfR)
  call CSF2SD(W(iRS)%Vec,CIR,iR)
  call opout(irs)
  call ipnout(-1)
  icsm = iR
  issm = iL
  call Densi2_mclr(2,De,Pe,CIL,CIR,0,0,0,n1dens,n2dens)

  if (.not. timedep) then
    if (response) then
      do iA=1,nnA
        do jA=1,nnA
          do kA=1,nnA
            do la=1,nnA
              ij1 = nnA*(iA-1)+ja
              ij2 = nna*(ja-1)+ia
              kl1 = nnA*(ka-1)+la
              kl2 = nna*(la-1)+ka
              rp(itri(ij1,kl1)) = Pe(itri(ij1,kl1))+Pe(itri(ij2,kl2))
            end do
          end do
        end do
      end do
      do iA=1,nnA
        do jA=1,nnA
          ij1 = nnA*(iA-1)+ja
          ij2 = nna*(ja-1)+ia
          rD(ij1) = De(ij1)+De(ij2)
        end do
      end do

    else
      call dcopy_(n2dens,Pe,1,rp,1)
      call dcopy_(n1dens,De,1,rD,1)
    end if
  else
    call dcopy_(n2dens,Pe,1,rp,1)
    call dcopy_(n1dens,De,1,rD,1)
    iCSM = iL
    iSSM = iR
    call Densi2_mclr(2,De,Pe,CIR,CIL,0,0,0,n1dens,n2dens)
    call daxpy_(n2Dens,-One,Pe,1,rp,1)
    call daxpy_(n1Dens,-One,De,1,rD,1)
  end if
  call mma_deallocate(CIL)
  call mma_deallocate(CIR)
else
  issm = iL
  icsm = iR
  call ipin(iLS)
  call ipin(iRS)
  call Densi2_mclr(2,De,Pe,W(iLS)%Vec,W(iRS)%Vec,0,0,0,n1dens,n2dens)
  if (.not. timedep) then
    if (response) then
      do iA=1,nnA
        do jA=1,nnA
          do kA=1,nnA
            do la=1,nnA
              ij1 = nnA*(iA-1)+ja
              ij2 = nna*(ja-1)+ia
              kl1 = nnA*(ka-1)+la
              kl2 = nna*(la-1)+ka
              rp(itri(ij1,kl1)) = Pe(itri(ij1,kl1))+Pe(itri(ij2,kl2))
            end do
          end do
        end do
      end do
      do iA=1,nnA
        do jA=1,nnA
          ij1 = nnA*(iA-1)+ja
          ij2 = nna*(ja-1)+ia
          rD(ij1) = De(ij1)+De(ij2)
        end do
      end do
    else
      call dcopy_(n2dens,Pe,1,rp,1)
      call dcopy_(n1dens,De,1,rD,1)
    end if
  else
    call dcopy_(n2dens,Pe,1,rp,1)
    call dcopy_(n1dens,De,1,rD,1)
    iCSM = iL
    iSSM = iR
    call ipin(iRS)
    call ipin(iLS)
    call Densi2_mclr(2,De,Pe,W(iRS)%Vec,W(iLS)%Vec,0,0,0,n1dens,n2dens)
    call daxpy_(n2Dens,-One,Pe,1,rp,1)
    call daxpy_(n1Dens,-One,De,1,rD,1)
  end if
end if
call mma_deallocate(Pe)
call mma_deallocate(De)

if (doDMRG) then  ! yma
  call dmrg_dim_change_mclr(RGras2(1:8),ndim,0)
  call dmrg_dim_change_mclr(RGras2(1:8),nna,0)
  n1dens = ndim**2
  n2dens = n1dens*(n1dens+1)/2
end if

end subroutine CIDens
