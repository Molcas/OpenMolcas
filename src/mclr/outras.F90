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

use MckDat, only: sLength
use gugx, only: SGS, CIS, EXS
use stdalloc, only: mma_allocate, mma_deallocate
use MCLR_Data, only: nConf1, nDensC, nDens2
use MCLR_Data, only: DspVec, lDisp
use MCLR_Data, only: LuTEMP
use input_mclr, only: nDisp, nSym, State_Sym, iMethod, nCSF, nConf, iMethod, iSpin, kPrint, nActEl, nElec3, nHole1, nRS1, nRS2, &
                      nRS3, nTPert

implicit none
integer iKapDisp(nDisp), iCiDisp(nDisp)
character(len=8) Label
integer Pstate_sym
logical CI
real*8, allocatable :: Kap1(:), Kap2(:), Kap3(:), CIp1(:)
integer iDisp, iSym, nConfm, jDisp, kDisp, iDisk, Len, iLen, iDIs, iRC, iOpt, iSymL, iPert

!----------------------------------------------------------------------*
!
! Ok construct hessian
!
!----------------------------------------------------------------------*

write(6,*)
write(6,*) '      Writing response to disk in Split guga GUGA format'
write(6,*)

idisp = 0
do iSym=1,nSym
  call Setup_MCLR(iSym)
  PState_SYM = ieor(State_Sym-1,iSym-1)+1
  nconfM = ncsf(PState_Sym)
  nconf1 = ncsf(PState_Sym)
  CI = .false.
  if ((iMethod == 2) .and. (nconf1 > 0)) CI = .true.
  if (CI .and. (nconf1 == 1) .and. (isym == 1)) CI = .false.

  ! Allocate areas for scratch and state variables

  call mma_allocate(Kap1,nDens2,Label='Kap1')
  call mma_allocate(Kap2,nDens2,Label='Kap2')
  call mma_allocate(Kap3,nDens2,Label='Kap3')
  if (CI) then
    call mma_allocate(CIp1,nconfM,Label='CIp1')
    call InCSFSD(Pstate_sym,State_sym,.true.)
  end if
  do jDisp=1,lDisp(iSym)
    iDisp = iDisp+1
    if (iand(ntpert(idisp),2**4) == 16) then
      kdisp = DspVec(idisp)

      iDisk = iKapDisp(iDisp)
      if (iDisk /= -1) then
        Len = nDensC
        call dDaFile(LuTemp,2,Kap1,Len,iDisk)
        call Uncompress(Kap1,Kap3,isym)
        if (CI) then
          ilen = nconfM
          idis = iCIDisp(iDisp)
          call dDaFile(LuTemp,2,CIp1,iLen,iDis)
        end if
        call GASync()
      else
        call GASync()
        Len = nDensC
        call FZero(Kap1,Len)
        call GADSum(Kap1,Len)
        if (CI) then
          len = nconfM
          call FZero(CIp1,Len)
          call GADSum(CIp1,Len)
        end if
      end if
      call GASync()
      call TCMO(Kap3,isym,-1)
      irc = ndens2
      Label = 'KAPPA'
      iopt = ibset(0,sLength)
      isyml = 2**(isym-1)
      ipert = kdisp
      write(6,'(A,I5," jDisp: ",I5," and iSym:",I5)') 'Writing KAPPA and CI in mclr for iDisp:',iDisp,jDisp,iSym
      call dWrMCk(iRC,iOpt,Label,ipert,Kap3,isyml)
      if (irc /= 0) call SysAbendMsg('outras','Error in wrmck','label=KAPPA')
      irc = nconfM
      iopt = ibset(0,sLength)
      Label = 'CI'
      isyml = 2**(isym-1)
      ipert = kdisp

      if (iand(kprint,8) == 8) write(6,*) 'Perturbation ',ipert
      if (CI) then
        call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,CIp1,0,pstate_sym,State_Sym)
        NCSF(1:nSym) = CIS%NCSF(1:nSym)
        NCONF = CIS%NCSF(pstate_sym)
        call mkGuga_Free(SGS,CIS,EXS)
      end if

      if ((imethod == 2) .and. (.not. CI) .and. (nconfM == 1)) CIp1(1) = 0.0d0
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
