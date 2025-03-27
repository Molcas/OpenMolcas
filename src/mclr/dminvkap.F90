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

subroutine DMInvKap(rMFact,rIn,nrIn,rOut,nrOut,rtemp,nrtemp,isym,iter)
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
!     In,Out rOut       Orbital rotaotion                              *
!                                                                      *
!     iSym              Symmetry of rotation                           *
!                                                                      *
!***********************************************************************

use Spool, only: LuWr
use MCLR_Data, only: SA
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use Definitions, only: wp, u6

implicit none
integer nrIn, nrOut, nrTemp, iSym, iter
real*8 rOut(nrOut), rMFact(*), rIn(nrIn), rtemp(nrTemp)

!                                                                      *
!***********************************************************************
!                                                                      *
call DMInvKap_Internal(rMFact)

! This is to allow type punning without an explicit interface
contains

subroutine DMInvKap_Internal(rMFact)

  use iso_c_binding
  use MCLR_Data, only: ipMat, nDensC
  use input_mclr, only: nSym, PT2, nAsh, nIsh, nOrb, nRs1, nRs2, nRs3

  implicit none
  real*8, target :: rMFact(*)
  integer, pointer :: iMFact(:)
  integer ip1, iS, jS, ii, nd, ip2, iRC
  real*8, external :: DDot_
  ip1 = 1

  if (doDMRG) then  ! yma
    call dmrg_spc_change_mclr(RGras2(1:8),nash)
    call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
  end if

  call DCopy_(nDensC,rIn,1,rOut,1)
  call Uncompress2(rIn,rtemp,isym)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do jS=1,nSym
    iS = ieor(js-1,iSym-1)+1
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !    kappa
    !         ip
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    do iI=1,nIsh(js)
      nD = nOrb(is)-nIsh(is)
      if (nd /= 0) then
        ip2 = ipMat(is,js)+nOrb(is)*(iI-1)
        irc = 0
        call c_f_pointer(c_loc(rMFact(ip1+nd**2)),iMFact,[ND])
        call dgetrs_('N',ND,1,rMFact(ip1),nd,iMFact,rtemp(ip2),nd,irc)
        nullify(iMFact)
        if (irc /= 0) then
          write(u6,*) 'Error in DGETRS called from dminvkap'
          call Abend()
        end if
        ip1 = ip1+nD*(nD+1)
      end if
    end do

    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !    kappa
    !         ap
    !                                                                  *
    !*******************************************************************
    !                                                                  *
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
        call c_f_pointer(c_loc(rMFact(ip1+nd**2)),iMFact,[ND])
        call dgetrs_('N',ND,1,rMFact(ip1),nd,iMFact,rtemp(ip2),nd,irc)
        nullify(iMFact)
        if (irc /= 0) then
          write(u6,*) 'Error in DGETRS called from dminvkap'
          call Abend()
        end if
        ip1 = ip1+nD*(nD+1)
      end if
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call Compress2(rtemp,nrtemp,rOut,nrOut,isym)

  if (doDMRG) call dmrg_spc_change_mclr(LRras2(1:8),nash)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Warn if the trial vector becomes large

  if ((ddot_(ndensc,rout,1,rout,1) > 100.0_wp) .and. (iter == 1)) then
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
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine DMInvKap_Internal

end subroutine DMInvKap
