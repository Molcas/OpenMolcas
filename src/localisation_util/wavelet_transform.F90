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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine Wavelet_Transform(irc,CMO,nSym,nBas,nFro,nOrb2Loc,inv,Silent,xNrm)
! Author: F. Aquilante
!
! Purpose: wavelet transform of the MO basis (inv=0)
!          "       backtransform (inv=1)

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nOrb2Loc(nSym), inv
real(kind=wp), intent(inout) :: CMO(*)
real(kind=wp), intent(out) :: xNrm
logical(kind=iwp), intent(in) :: Silent
integer(kind=iwp) :: iSym, kOff1, kOff2, l_Scr, njOrb
type(DSBA_Type) :: C
real(kind=wp), allocatable :: Scr(:)
character(len=*), parameter :: SecNam = 'Wavelet_Transform'
integer(kind=iwp), external :: Log2
real(kind=wp), external :: ddot_

irc = 0
xNrm = Zero
if (.not. Silent) then
  if (inv == 0) write(u6,'(/,1X,A)') 'Wavelet transform of the MOs'
  if (inv == 1) write(u6,'(/,1X,A)') 'Inverse wavelet transform of the MOs'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',nFro(:)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to transform:',nOrb2Loc(:)
end if

call Allocate_DT(C,nBas,nBas,nSym,label='C',Ref=CMO)

if (inv == 1) then
  ! Inverse wavelet transform
  njOrb = Log2(nOrb2Loc(1))
  l_Scr = nBas(1)*(2**njOrb)
  do iSym=2,nSym
    njOrb = Log2(nOrb2Loc(iSym))
    l_Scr = max(l_Scr,nBas(iSym)*(2**njOrb))
  end do
  call mma_allocate(Scr,l_Scr,label='Scratch')
  do iSym=1,nSym
    if (nOrb2Loc(iSym) > 0) then
      kOff1 = nFro(iSym)+1
      kOff2 = kOff1
      njOrb = Log2(nOrb2Loc(iSym))
      do while (njOrb >= 1)
        call Inv_FWT_Haar(nBas(iSym),njOrb,Scr,C%SB(iSym)%A2(:,kOff2:))
        njOrb = 2**njOrb
        kOff2 = kOff2+njOrb
        njOrb = Log2(nOrb2Loc(iSym)-njOrb)
      end do
      xNrm = xNrm+dDot_(nBas(iSym)*nOrb2Loc(iSym),C%SB(iSym)%A2(:,kOff1:),1,C%SB(iSym)%A2(:,kOff1:),1)
      if (irc /= 0) then
        irc = 1
        xNrm = -huge(xNrm)
        call FreeMem()
        return
      end if
    end if
  end do
  xNrm = sqrt(xNrm)
  call mma_deallocate(Scr)
else
  njOrb = Log2(nOrb2Loc(1))
  l_Scr = nBas(1)*(2**njOrb-1)
  do iSym=2,nSym
    njOrb = Log2(nOrb2Loc(iSym))
    l_Scr = max(l_Scr,nBas(iSym)*(2**njOrb-1))
  end do
  call mma_allocate(Scr,l_Scr,label='Scratch')
  do iSym=1,nSym
    if (nOrb2Loc(iSym) > 0) then
      kOff1 = nFro(iSym)+1
      kOff2 = kOff1
      njOrb = Log2(nOrb2Loc(iSym))
      do while (njOrb >= 1)
        call FWT_Haar(nBas(iSym),njOrb,Scr,C%SB(iSym)%A2(:,kOff2:))
        njOrb = 2**njOrb
        kOff2 = kOff2+njOrb
        njOrb = Log2(nOrb2Loc(iSym)-njOrb)
      end do
      xNrm = xNrm+dDot_(nBas(iSym)*nOrb2Loc(iSym),C%SB(iSym)%A2(:,kOff1:),1,C%SB(iSym)%A2(:,kOff1:),1)
      if (irc /= 0) then
        irc = 1
        xNrm = -huge(xNrm)
        call FreeMem()
        return
      end if
    end if
  end do
  xNrm = sqrt(xNrm)
  call mma_deallocate(Scr)
end if

call FreeMem()

return

contains

subroutine FreeMem()

  call Deallocate_DT(C)
  if (allocated(Scr)) call mma_deallocate(Scr)

end subroutine FreeMem

end subroutine Wavelet_Transform
