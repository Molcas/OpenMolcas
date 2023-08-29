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

subroutine Cho_P_OpenVR(iOpt)
!
! Purpose: open (iOpt=1) or close (iOpt=2) local and global storage
!          files.

use Para_Info, only: Is_Real_Par, nProcs
use Cholesky, only: Cho_AdrVec, Cho_Fake_Par, Cho_Real_Par, LuCho, LuCho_G, LuPri, LuRed_G, LuRst_G, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iOpt
integer(kind=iwp) :: ID, iSym
character(len=6) :: FNRst, FNVec(8)
character(len=5) :: FNRed
character(len=*), parameter :: SecNam = 'Cho_P_OpenVR'

! Local files.
! ------------

if (Cho_Real_Par) then
  ID = 1
else
  ID = 2
end if
call Cho_OpenVR(iOpt,ID)

! Global files for restart info and reduced set indices.
! ------------------------------------------------------

if (Cho_Real_Par) then
  if (iOpt == 1) then
    LuRed_G = 7
    FNRed = 'CHRED'
    call DAName_MF_WA(LuRed_G,FNRed)
    LuRst_G = 7
    FNRst = 'CHORST'
    call DAName_MF_WA(LuRst_G,FNRst)
    do iSym=1,nSym
      LuCho_G(iSym) = 7
      write(FNVec(iSym),'(A5,I1)') 'CHVEC',iSym
      call DaName_MF_WA(LuCho_G(iSym),FNVec(iSym))
    end do
  else if (iOpt == 2) then
    if (LuRed_G > 0) then
      call DAClos(LuRed_G)
      LuRed_G = 0
    end if
    if (LuRst_G > 0) then
      call DAClos(LuRst_G)
      LuRst_G = 0
    end if
    do iSym=1,nSym
      if (LuCho_G(iSym) > 0) then
        call DaClos(LuCho_G(iSym))
        LuCho_G(iSym) = 0
      end if
    end do
  else
    write(Lupri,*) SecNam,': iOpt out of bounds: ',iOpt
    call Cho_Quit('Error in '//SecNam,104)
  end if
else
  if (CHO_FAKE_PAR .and. (nProcs > 1) .and. Is_Real_Par()) then
    if (iOpt == 1) then
      if (CHO_ADRVEC == 1) then
        do iSym=1,nSym
          LuCho_G(iSym) = 7
          write(FNVec(iSym),'(A5,I1)') 'CHVCL',iSym
          call DaName_MF_WA(LuCho_G(iSym),FNVec(iSym))
        end do
      else if (CHO_ADRVEC == 2) then
        do iSym=1,nSym
          LuCho_G(iSym) = 7
          write(FNVec(iSym),'(A5,I1)') 'CHVCL',iSym
          call DaName_MF(LuCho_G(iSym),FNVec(iSym))
        end do
      else
        call Cho_Quit('CHO_ADRVEC out of bounds in '//SecNam,102)
        LuCho_G(1:nSym) = 0
      end if
      ! Swap units so that
      !    LuCho_G points to 'CHVEC'
      !    LuCho   points to 'CHVCL'
      call iSwap(nSym,LuCho,1,LuCho_G,1)
    else if (iOpt == 2) then
      do iSym=1,nSym
        if (LuCho_G(iSym) > 0) then
          call DaClos(LuCho_G(iSym))
          LuCho_G(iSym) = 0
        end if
      end do
    else
      write(Lupri,*) SecNam,': iOpt out of bounds: ',iOpt
      call Cho_Quit('Error in '//SecNam,104)
    end if
  end if
end if

end subroutine Cho_P_OpenVR
