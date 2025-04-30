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

subroutine DMInvKap(rMFact,iMFact,rIn,rOut,rtemp,isym,iter)
!***********************************************************************
!                                                                      *
!     _____     -1                                                     *
!     Kappa  = M  Kappa                                                *
!          ip   pq     iq                                              *
!                                                                      *
!                                                                      *
!     In: rMFact        Factorized preconditioner (diagonal part       *
!                       of the electronic hessian that couple          *
!                       rotations with one common index)               *
!     Out: rOut         Orbital rotaotion                              *
!                                                                      *
!     iSym              Symmetry of rotation                           *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: Mul
use Spool, only: LuWr
use MCLR_Data, only: ipMat, nDens, nDensC, SA
use input_mclr, only: nAsh, nIsh, nOrb, nRs1, nRs2, nRs3, nSym, PT2
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: rMFact(*), rIn(nDensC)
integer(kind=iwp), intent(in) :: iMFact(*), iter
real(kind=wp), intent(out) :: rOut(nDensC), rtemp(nDens)
integer(kind=iwp), intent(inout) :: iSym
integer(kind=iwp) :: ii, ip1, ip2, ipi, iRC, iS, jS, nd
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
ip1 = 1
ipi = 1

if (doDMRG) then  ! yma
  nash(:) = RGras2(:)
  nrs2(:) = RGras2(:)
end if

call Uncompress2(rIn,rtemp,isym)
!                                                                      *
!***********************************************************************
!                                                                      *
do jS=1,nSym
  iS = Mul(js,iSym)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !    kappa
  !         ip
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iI=1,nIsh(js)
    nD = nOrb(is)-nIsh(is)
    if (nd /= 0) then
      ip2 = ipMat(is,js)+nOrb(is)*(iI-1)
      irc = 0
      call dgetrs_('N',ND,1,rMFact(ip1),nd,iMFact(ipi),rtemp(ip2),nd,irc)
      if (irc /= 0) then
        write(u6,*) 'Error in DGETRS called from dminvkap'
        call Abend()
      end if
      ip1 = ip1+nD**2
      ipi = ipi+nD
    end if
  end do

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !    kappa
  !         ap
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iI=1,nAsh(js)
    nD = nOrb(is)-nAsh(is)
    if (SA .or. PT2) then
      if (iI <= nRs1(jS)) then
        nD = nOrb(is)-nRs1(js)
      else if (iI <= nRs1(jS)+nRs2(jS)) then
        nD = nOrb(is)-nRs2(js)
      else if (iI <= nRs1(jS)+nRs2(jS)+nRs3(jS)) then
        nD = nOrb(is)-nRs3(js)
      end if
    end if
    if (nd /= 0) then
      ip2 = ipMat(is,js)+nOrb(is)*(iI-1+nIsh(js))
      irc = 0
      call dgetrs_('N',ND,1,rMFact(ip1),nd,iMFact(ipi),rtemp(ip2),nd,irc)
      if (irc /= 0) then
        write(u6,*) 'Error in DGETRS called from dminvkap'
        call Abend()
      end if
      ip1 = ip1+nD**2
      ipi = ipi+nD
    end if
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call Compress2(rtemp,rOut,isym)

if (doDMRG) nash(:) = LRras2(:)
!                                                                      *
!***********************************************************************
!                                                                      *
! Warn if the trial vector becomes large

if ((ddot_(nDensC,rOut,1,rOut,1) > 100.0_wp) .and. (iter == 1)) then
  write(LuWr,*) '****************************************'
  write(LuWr,*) '*                                      *'
  write(LuWr,*) '*           WARNING!!                  *'
  write(LuWr,*) '* Elements in the E^[2] matrix small!! *'
  write(LuWr,*) '* The calculation might diverge.       *'
  write(LuWr,*) '*                                      *'
  write(LuWr,*) '* Check your active space!!!!          *'
  write(LuWr,*) '*                                      *'
  write(LuWr,*) '* Make sure degenerate orbitals do not *'
  write(LuWr,*) '* belong to different spaces.          *'
  write(LuWr,*) '* Note that no LR code can handle      *'
  write(LuWr,*) '* 2.0 occupancy in active orbitals!!   *'
  write(LuWr,*) '****************************************'
end if
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine DMInvKap
