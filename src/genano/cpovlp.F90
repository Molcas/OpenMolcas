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

subroutine CpOvlp(from,to)

use Genano_globals, only: MxLqn, nSym, nBas, iSymBk, Center, BasName, symlab
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: from(*)
real(kind=wp), intent(out) :: to(*)
integer(kind=iwp) :: iBas, iBasX, iBlk, ijBas, ind, indx, iSym, jBas, jBasX, jndx
logical(kind=iwp) :: oki, okj

do iBlk=1,(MxLqn+1)**2
  indx = 0
  ijBas = 0
  iBasX = 0
  do iSym=1,nSym
    do iBas=1,nBas(iSym)
      iBasX = iBasX+1
      oki = BasName(iBasX)(1:len(Center)) == Center
      oki = oki .and. (BasName(iBasX)(len(Center)+1:) == symlab(iBlk))
      if (oki) indx = indx+1
      jndx = 0
      jBasX = iBasX-iBas
      do jBas=1,iBas
        jBasX = jBasX+1
        ijBas = ijBas+1
        okj = BasName(jBasX)(1:len(Center)) == Center
        okj = okj .and. (BasName(jBasX)(len(Center)+1:) == symlab(iBlk))
        if (okj) jndx = jndx+1
        if (oki .and. okj) then
          ind = jndx+indx*(indx-1)/2+iSymBk(iBlk)-1
          to(ind) = from(ijBas)
        end if
      end do
    end do
  end do
end do
!write(u6,*) '*** Overlap matrix in CpOvlp ***'
!iBlk = 0
!do iLqn=0,MxLqn
!  do iShell=-iLqn,iLqn
!    iBlk = iBlk+1
!    if (nPrim(iLqn) > 0) then
!      write(u6,'(a,2i5)') ' Block',iLqn,iShell
!      call Triprt(' ','(6F12.6)',to(iSymBk(iBlk)),nPrim(iLqn))
!    end if
!  end do
!end do

return

end subroutine CpOvlp
