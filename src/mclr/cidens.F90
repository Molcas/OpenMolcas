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

use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: ipin, W
use MCLR_Data, only: n1Dens, n2Dens, nConf1, nNA
use CI_Interfaces_MCLR, only: CI_MCLR_TDens
use input_mclr, only: nCSF, TimeDep
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
logical(kind=iwp), intent(in) :: Response
integer(kind=iwp), intent(in) :: iLS, iRS, iL, iR
real(kind=wp), intent(_OUT_) :: rP(*), rD(*)
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp) :: IA, ij1, ij2, JA, KA, kl1, kl2, LA, nDim
real(kind=wp), allocatable :: De(:), Pe(:)

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
  call dmrg_dim_change_mclr(LRras2,ndim,0)
  call dmrg_dim_change_mclr(LRras2,nna,0)
  n1Dens = ndim**2
  n2Dens = nTri_Elem(n1Dens)
end if

call mma_allocate(De,n1Dens,Label='De')
call mma_allocate(Pe,n2Dens,Label='Pe')
De(:) = Zero
Pe(:) = Zero

! (through the CI interface; L of symmetry iL, R of symmetry iR)
call ipin(iLS)
call ipin(iRS)
call CI_MCLR_TDens(iState,2,iL,iR,nCSF(iL),W(iLS)%A(1:nCSF(iL)),nCSF(iR),W(iRS)%A(1:nCSF(iR)),n1Dens,De,n2Dens,Pe)
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
            rp(iTri(ij1,kl1)) = Pe(iTri(ij1,kl1))+Pe(iTri(ij2,kl2))
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
    rp(1:n2Dens) = Pe(:)
    rD(1:n1Dens) = De(:)
  end if
else
  rp(1:n2Dens) = Pe(:)
  rD(1:n1Dens) = De(:)
  call ipin(iRS)
  call ipin(iLS)
  call CI_MCLR_TDens(iState,2,iR,iL,nCSF(iR),W(iRS)%A(1:nCSF(iR)),nCSF(iL),W(iLS)%A(1:nCSF(iL)),n1Dens,De,n2Dens,Pe)
  rp(1:n2Dens) = rp(1:n2Dens)-Pe(:)
  rD(1:n1Dens) = rD(1:n1Dens)-De(:)
end if
call mma_deallocate(Pe)
call mma_deallocate(De)

if (doDMRG) then  ! yma
  call dmrg_dim_change_mclr(RGras2,ndim,0)
  call dmrg_dim_change_mclr(RGras2,nna,0)
  n1Dens = ndim**2
  n2Dens = nTri_Elem(n1dens)
end if

end subroutine CIDens
