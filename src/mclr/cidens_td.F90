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

use ipPage, only: ipin, ipnout, W
use MCLR_Data, only: ipCI, n1Dens, n2Dens, nConf1, NOCSF, XISPSM
use CandS, only: ICSM, ISSM
use input_mclr, only: nCSF, State_Sym
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
integer(kind=iwp), intent(in) :: iCI, iS
real(kind=wp), intent(_OUT_) :: rP(*), rD(*)
integer(kind=iwp) :: nC, nConfL, nConfR
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, ijkl, j, k, l
#endif
real(kind=wp), allocatable :: CIL(:), CIR(:), De(:), Pe(:)

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

if (nocsf == 0) then
  nConfL = max(ncsf(iS),nint(xispsm(iS,1)))
  nConfR = max(ncsf(State_SYM),nint(xispsm(STATE_SYM,1)))
  nC = max(nconfL,nconfR)
  call mma_allocate(CIL,nC,Label='CIL')
  call mma_allocate(CIR,nC,Label='CIR')

  ! iCI is as long as ipcid for the timedep case!

  call ipin(iCI)
  call ipin(ipCI)
  call CSF2SD(W(iCI)%A,CIR,iS)
  call CSF2SD(W(ipCI)%A,CIL,State_SYM)

  !write(u6,*) 'ipL*ipL',ddot_(nConfL,CIL,1,CIL,1)
  !write(u6,*) 'ipR*ipR',ddot_(nConfR,CIR,1,CIR,1)
  !call RecPrt('CIL',' ',CIL,nConfL,1)
  !call RecPrt('CIR',' ',CIR,nConfR,1)

  call ipnout(-1)
  icsm = iS
  issm = STATE_SYM

  ! <P|E_pq|0> & <P|e_pqrs|0> -> ipDe & ipP
  ! ipL is the bra side vector
  De(:) = Zero
  Pe(:) = Zero
  call Densi2_mclr(2,De,Pe,CIL,CIR,0,0,0,n1Dens,n2Dens)

  !write(u6,*) 'De*De',ddot_(n1Dens,De,1,De,1)
  !call RecPrt('De',' ',De,n1Dens,1)

  rp(1:n2Dens) = Pe(:)
  rD(1:n1Dens) = De(:)

  !write(u6,*) 'rD*rD',ddot_(n1Dens,rD,1,rD,1)
  !call RecPrt('rD',' ',rD,n1Dens,1)

  call ipin(iCI)
  call ipin(ipCI)
  call CSF2SD(W(iCI)%A(1+nconf1),CIL,iS)
  call CSF2SD(W(ipci)%A,CIR,State_SYM)

  !write(u6,*) 'CIL*CIL',ddot_(nConfL,CIL,1,CIL,1)
  !write(u6,*) 'CIR*CIR',ddot_(nConfR,CIR,1,CIR,1)
  !call RecPrt('CIL',' ',CIL,nConfL,1)
  !call RecPrt('CIR',' ',CIR,nConfR,1)

  call ipnout(-1)
  issm = iS
  icsm = STATE_SYM
  De(:) = Zero
  Pe(:) = Zero
  call Densi2_mclr(2,De,Pe,CIL,CIR,0,0,0,n1Dens,n2Dens)

  !write(u6,*) 'De*De',ddot_(n1Dens,De,1,De,1)
  !Call RecPrt('De',' ',De,n1Dens,1)

  rp(1:n2Dens) = rp(1:n2Dens)-Pe(:)
  rD(1:n1Dens) = rD(1:n1Dens)-De(:)

  !rP(1:n2Dens) = -rP(1:n2Dens)
  !rD(1:n1Dens) = -rD(1:n1Dens)

  !write(u6,*) 'rD*rD',ddot_(n1Dens,rD,1,rD,1)
  !call RecPrt('rD',' ',rD,n1Dens,1)

  call mma_deallocate(CIL)
  call mma_deallocate(CIR)

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
end if

call mma_deallocate(Pe)
call mma_deallocate(De)

end subroutine CIDens_TD
