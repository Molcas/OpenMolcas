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

!#define _DEBUGPRINT_
subroutine CIDens_TD(iCI,iS,rP,rD)

use CI_Info, only: CI_St
use ipPage, only: ipin, W
use MCLR_Data, only: ipCI, n1Dens, n2Dens, nConf1
use CI_Operations, only: CI_Lib_Dens
use input_mclr, only: nCSF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Index_Functions, only: iTri
use input_mclr, only: ntAsh
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: iCI, iS
real(kind=wp), intent(_OUT_) :: rP(*), rD(*)
integer(kind=iwp) :: iSy
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, ijkl, j, k, l
#endif
real(kind=wp), allocatable :: De(:), Pe(:)

! LS = CI
!
! Ok we once more want to hide Jeppe's routines from
! the eyes of the world, so everyone belives that I have done
! all the work.
! If we have spin dependent Hamiltonian we will work with
! SD in all parts of the program, no CSF is necessary,
! otherwise we will do the optimazation in CSF's to increase
! convergence.
!
! Input:
!
! response: true if the density should be used for response
!           calculations
! LS : CI Coeff for left state
! RS : CI Coeff for right state
! iL : Symmetry of left state
! iR : Symmetry of right state
!
!           +       +
! iS=1 E  =a  a  + a  a     ! Singlet operator
!       pq  ap aq   Bp Bq
!
!            +       +
! iS=-1 T  =a  a  - a a     ! Triplet operator
!        pq  ap aq   Bp Bq
!
! Output:
!
!  rP : Two Electron Density
!  rD : One Electron Density

!call ipin(iCI)
!write(u6,*) 'iCI*iCI',ddot_(2*nConf1,W(iCI)%A,1,W(iCI)%A,1)

if (nconf1 == 0) return

call mma_allocate(De,n1Dens,Label='De')
call mma_allocate(Pe,n2Dens,Label='Pe')
rD(1:n1Dens) = Zero
rP(1:n2Dens) = Zero

! (through the CI interface) iCI is as long as ipcid for the timedep case!
iSy = CI_St(iState)%STSym
! <0|E_pq|P> & <0|e_pqrs|P>: the bra is the reference
call ipin(iCI)
call ipin(ipCI)
call CI_Lib_Dens(iState,2,0,iSy,iS,ncsf(iSy),W(ipCI)%A(1:ncsf(iSy)),ncsf(iS),W(iCI)%A(1:ncsf(iS)),n1Dens,De,n2Dens,Pe)
rp(1:n2Dens) = Pe(:)
rD(1:n1Dens) = De(:)
! minus <P'|E_pq|0> & <P'|e_pqrs|0>, P' the second half of iCI
call ipin(iCI)
call ipin(ipCI)
call CI_Lib_Dens(iState,2,0,iS,iSy,ncsf(iS),W(iCI)%A(nconf1+1:nconf1+ncsf(iS)),ncsf(iSy),W(ipCI)%A(1:ncsf(iSy)),n1Dens,De, &
                   n2Dens,Pe)
rp(1:n2Dens) = rp(1:n2Dens)-Pe(:)
rD(1:n1Dens) = rD(1:n1Dens)-De(:)
# ifdef _DEBUGPRINT_
do i=1,ntash
  do j=1,ntash
    do k=1,ntash
      do l=1,ntash
        ijkl = iTri(ntash*(j-1)+i,k+(l-1)*ntash)
        write(u6,'(I1,I1,I1,I1,F12.6)') i,j,k,l,rp(ijkl)
      end do
    end do
  end do
end do
# endif

call mma_deallocate(Pe)
call mma_deallocate(De)

end subroutine CIDens_TD
