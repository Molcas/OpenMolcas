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
! Copyright (C) 1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine H0(rdia,MP1,MP2,MQ,isym,nprciv,TimeDep)
! The explicit part of the Hamiltonian for the CI preconditioner: the P
! space of at most MP2 CSFs and its block of H from the CI library
! (CI_PSpace), shifted by the energy, squared and factorized (H0S, H0F;
! SBIDT: the CSFs of the P space). Only a P2 space (MP1 = MQ = 0).

use Index_Functions, only: nTri_Elem
use CI_Operations, only: CI_PSpace
use MCLR_Data, only: FIMO, H0F, H0S, Int2, nGP, SBIDT
use input_mclr, only: nCSF
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: rdia(*)
integer(kind=iwp), intent(in) :: MP1, MP2, MQ, iSym
integer(kind=iwp), intent(out) :: nprciv
logical(kind=iwp), intent(in) :: TimeDep
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp) :: i, iRC, MXP
real(kind=wp) :: ENA
real(kind=wp), allocatable :: H0T(:)
real(kind=wp), external :: E2, E2_TD

if ((MP1 /= 0) .or. (MQ /= 0)) then
  write(u6,*) 'H0: only a P2 space is implemented (MP1 = MQ = 0)'
  call Abend()
end if
MXP = MP2
if (TimeDep) then
  EnA = E2_td(FIMO,Int2,0,-1)
else
  EnA = E2(FIMO,Int2,0,-1)
end if

call mma_allocate(SBIDT,MXP,Label='SBIDT')
call mma_allocate(H0S,MXP**2,Label='H0S')
call mma_allocate(H0F,MXP,Label='H0F')
call mma_allocate(H0T,nTri_Elem(MXP),Label='H0T')
call CI_PSpace(iState,iSym,nCSF(iSym),rdia,MXP,nprciv,SBIDT,size(H0T),H0T)
do i=1,nprciv
  H0T(nTri_Elem(i)) = H0T(nTri_Elem(i))-ENA
end do
if (NGP) call mkp1(nprciv,SBIDT,H0T,rdia)
call square(H0T,H0S,1,NPRCIV,NPRCIV)
call mma_deallocate(H0T)
irc = 0
call dgetrf_(NPRCIV,NPRCIV,H0S,NPRCIV,H0F,irc)
if (irc /= 0) then
  write(u6,*) 'Sorry but you have an singular ci matrix'
  write(u6,*) 'Set ExpDimension and restart mclr'
  call Abend()
end if

end subroutine H0
