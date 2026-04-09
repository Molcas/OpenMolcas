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
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine OutRAS(iKapDisp,iCiDisp)
!*******************************************************************
!                                                                  *
! Contracts the response coefficient to the hessian                *
!                                                                  *
! Input                                                            *
!       iKapDisp : Disk locations of solutions to respons equation *
!       iCIDisp  : Disk locations of CI Soulutions to response     *
!                                                                  *
! Author: Anders Bernhardsson, 1996                                *
!         Theoretical Chemistry, University of Lund                *
!*******************************************************************

use Symmetry_Info, only: Mul
use MckDat, only: sLength
use gugx, only: SGS, CIS, EXS
use MCLR_Data, only: DspVec, lDisp, LuTEMP, nConf1, nDens, nDensC
use input_mclr, only: iMethod, iSpin, kPrint, nActEl, nConf, nCSF, nDisp, nElec3, nHole1, nRS1, nRS2, nRS3, nSym, nTPert, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iKapDisp(nDisp), iCiDisp(nDisp)
integer(kind=iwp) :: iDIs, iDisk, iDisp, iLen, iOpt, iPert, iRC, iSym, iSymL, jDisp, kDisp, Length, nConfm, Pstate_sym
logical(kind=iwp) :: CI
character(len=8) :: Label
real(kind=wp), allocatable :: CIp1(:), Kap1(:), Kap2(:), Kap3(:)

!----------------------------------------------------------------------*
!
! Ok construct hessian
!
!----------------------------------------------------------------------*

write(u6,*)
write(u6,*) '      Writing response to disk in Split guga GUGA format'
write(u6,*)

idisp = 0
do iSym=1,nSym
  call Setup_MCLR(iSym)
  PState_SYM = Mul(State_Sym,iSym)
  nconfM = ncsf(PState_Sym)
  nconf1 = ncsf(PState_Sym)
  CI = .false.
  if ((iMethod == 2) .and. (nconf1 > 0)) CI = .true.
  if (CI .and. (nconf1 == 1) .and. (isym == 1)) CI = .false.

  ! Allocate areas for scratch and state variables

  call mma_allocate(Kap1,nDens,Label='Kap1')
  call mma_allocate(Kap2,nDens,Label='Kap2')
  call mma_allocate(Kap3,nDens,Label='Kap3')
  if (CI) then
    call mma_allocate(CIp1,nconfM,Label='CIp1')
    call InCSFSD(Pstate_sym,State_sym)
  end if
  do jDisp=1,lDisp(iSym)
    iDisp = iDisp+1
    if (btest(ntpert(idisp),4)) then
      kdisp = DspVec(idisp)

      iDisk = iKapDisp(iDisp)
      if (iDisk /= -1) then
        Length = nDensC
        call dDaFile(LuTemp,2,Kap1,Length,iDisk)
        call Uncompress(Kap1,Kap3,isym)
        if (CI) then
          ilen = nconfM
          idis = iCIDisp(iDisp)
          call dDaFile(LuTemp,2,CIp1,iLen,iDis)
        end if
        call GASync()
      else
        call GASync()
        Length = nDensC
        Kap1(1:Length) = Zero
        call GADSum(Kap1,Length)
        if (CI) then
          Length = nconfM
          CIp1(1:Length) = Zero
          call GADSum(CIp1,Length)
        end if
      end if
      call GASync()
      call TCMO(Kap3,isym,-1)
      irc = nDens
      Label = 'KAPPA'
      iopt = ibset(0,sLength)
      isyml = ibset(0,isym-1)
      ipert = kdisp
      write(u6,'(A,I5," jDisp: ",I5," and iSym:",I5)') 'Writing KAPPA and CI in mclr for iDisp:',iDisp,jDisp,iSym
      call dWrMCk(iRC,iOpt,Label,ipert,Kap3,isyml)
      if (irc /= 0) call SysAbendMsg('outras','Error in wrmck','label=KAPPA')
      irc = nconfM
      iopt = ibset(0,sLength)
      Label = 'CI'
      isyml = ibset(0,isym-1)
      ipert = kdisp

      if (btest(kprint,3)) write(u6,*) 'Perturbation ',ipert
      if (CI) then
        call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,CIp1,0,pstate_sym,State_Sym)
        NCSF(1:nSym) = CIS%NCSF(1:nSym)
        NCONF = CIS%NCSF(pstate_sym)
        call mkGuga_Free(SGS,CIS,EXS)
      end if

      if ((imethod == 2) .and. (.not. CI) .and. (nconfM == 1)) CIp1(1) = Zero
      call dWrMCk(iRC,iOpt,Label,ipert,CIp1,isyml)
      if (irc /= 0) call SysAbendMsg('outras','Error in wrmck',' ')
    end if
    !*******************************************************************

  end do

  ! Free areas for scratch and state variables

  if (CI) call mma_deallocate(CIp1)
  call mma_deallocate(Kap3)
  call mma_deallocate(Kap2)
  call mma_deallocate(Kap1)
end do

end subroutine OutRAS
