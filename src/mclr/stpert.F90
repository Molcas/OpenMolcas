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

subroutine StPert()

use Index_Functions, only: nTri_Elem
use MckDat, only: sNew
use MCLR_Data, only: FnMck, Hss, lDisp, LuMck, SwLbl
use Orbital_Spaces, only: nAsh, nIsh, nSym
use input_mclr, only: McKinley, nBas, nDisp, nTPert, PT2
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: idum(1), iDummer, iOpt, iRC, iS, nHss
character(len=288) :: Header
character(len=16) :: Label
character(len=8) :: MckLbl

iDummer = 0
nHss = 0
do iS=1,nSym
  nHss = nHss+nTri_Elem(lDisp(is))
end do
call mma_allocate(Hss,nHss,Label='Hss')
Hss(:) = Zero

if (.not. Mckinley) then
  irc = -1
  iopt = ibset(0,sNew)
  call OPNMCK(irc,iopt,FNMCK,LUMCK)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error opening MCKINT'
    call Abend()
  end if
  irc = -1
  iopt = 0
  LABEL = 'SEWARD'
  if (PT2) LABEL = 'PT2LAG'
  MckLbl = 'PERT'
  call cWrMck(iRC,iOpt,MckLbl,1,LABEL,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'NDISP'
  idum(1) = ndisp
  call WrMck(iRC,iOpt,MckLbl,1,idum,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'TDISP'
  call WrMck(iRC,iOpt,MckLbl,1,ntpert,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'Title'
  Header = ' '
  call cWrMck(iRC,iOpt,MckLbl,1,Header,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'nSym'
  idum(1) = nSym
  call WrMck(iRC,iOpt,MckLbl,1,idum,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'nBas'
  call WrMck(iRC,iOpt,MckLbl,1,nBas,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'ldisp'
  call WrMck(iRC,iOpt,MckLbl,1,ldisp,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'chdisp'
  call cWrMck(iRC,iOpt,MckLbl,1,swlbl(1),iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'NISH'
  call WrMck(iRC,iOpt,MckLbl,1,nish,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
  irc = -1
  iopt = 0
  MckLbl = 'NASH'
  call WrMck(iRC,iOpt,MckLbl,1,nash,iDummer)
  if (irc /= 0) then
    write(u6,*) 'StPert: Error writing to MCKINT'
    write(u6,'(A,A)') 'MckLbl=',MckLbl
    call Abend()
  end if
end if

end subroutine StPert
