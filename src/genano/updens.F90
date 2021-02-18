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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine UpDens()

use Genano_globals, only: nSym, nBas, kSet, tDsym, pDsym, wSet, Cmo, Occ, Eps, thr, wc0, wc1, BasName
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IndCmo, IndNme, IndOcc, iOrb, iSym
real(kind=wp) :: e, o, w

pDsym(:) = Zero
IndOcc = 1
IndCmo = 1
IndNme = 1
do iSym=1,nSym
  !write(u6,'(a,i5)') ' iSym:   ',iSym
  do iOrb=1,nBas(iSym)
    !write(u6,'(a,i5)') '   iOrb:   ',iOrb
    !write(u6,'(a,i5)') '   IndOcc: ',IndOcc
    !write(u6,'(a,i5)') '   IndCmo: ',IndCmo
    o = Occ(IndOcc)
    e = Eps(IndOcc)
    w = wc0-wc1*e
    w = max(w,One)
    !write(u6,'(a,2i3,3f12.6)') 'iSym.iOrb,o,e,w',iSym,iOrb,o,e,w
    if (abs(o) > thr) then
      call UpOrb(nBas(iSym),o,w,Cmo(IndCmo),BasName(IndNme))
    end if
    IndCmo = IndCmo+nBas(iSym)
    IndOcc = IndOcc+1
  end do
  IndNme = IndNme+nBas(iSym)
end do
!write(u6,*) '*** Partial density matrix in UpDens ***'
!iBlk = 0
!do iLqn=0,MxLqn
!  do iShell=-iLqn,iLqn
!    iBlk = iBlk+1
!    if (nPrim(iLqn) > 0) then
!      write(u6,'(a,2i5)') ' Block',iLqn,iShell
!      call Triprt(' ','(6f12.6)',pDsym(iSymBk(iBlk)),nPrim(iLqn))
!    end if
!  end do
!end do
tDsym(:) = tDsym(:)+wSet(kSet)*pDsym(:)
!write(u6,*) '*** Total density matrix in UpDens ***'
!iBlk = 0
!do iLqn=0,MxLqn
!  do iShell=-iLqn,iLqn
!    iBlk = iBlk+1
!    if (nPrim(iLqn) > 0) then
!      write(u6,'(a,2i5)') ' Block',iLqn,iShell
!      call Triprt(' ','(6f12.6)',tDsym(iSymBk(iBlk)),nPrim(iLqn))
!    end if
!  end do
!end do

return

end subroutine UpDens
