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

subroutine InpCtl_GenANO()

use Genano_globals, only: nSets, iProj, kRfSet, nPrim, nCore, kSet, isUHF, thr, wSet, wc0, wc1, rowise, lftdeg, rydgen, Center, &
                          Title
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: LuIn, i, err
logical(kind=iwp) :: RdHead, done
character(len=72) :: Key, KeyIn
integer(kind=iwp), external :: IsFreeUnit

Title = 'Atom'
RdHead = .false.
nSets = 1
Center = 'ANO '
thr = 1.0e-8_wp
wc0 = Zero
wc1 = Zero
kRfSet = 1
rowise = .false.
lftdeg = .false.
rydgen = .false.
iProj = 0
nPrim(:) = 0
nCore(:) = 0
call mma_allocate(wSet,1,label='wSet')
wSet(1) = One
kSet = 0
isUHF = 0
!-----------------------------------------------------------------------
LuIn = IsFreeUnit(11)
call SpoolInp(LuIn)
call RdNlst(LuIn,'GENANO')
!-----------------------------------------------------------------------
done = .false.
do
  read(LuIn,'(a)',iostat=err) KeyIn
  if (err /= 0) exit
  Key = KeyIn
  !write(u6,*) 'echo> ',KeyIn
  call zlcase(Key)
  if (Key(1:1) == '*') then
    continue
  else if (Key == ' ') then
    continue
  else if (Key == 'title') then
    RdHead = .true.
  else if (Key == 'sets') then
    RdHead = .false.
    read(LuIn,*,iostat=err) nSets
    if (err /= 0) then
      write(u6,*) 'Error while reading input, keyword: ',trim(Key)
      call Quit_OnUserError()
    end if
    if (allocated(wSet)) call mma_deallocate(wSet)
    call mma_allocate(wSet,nSets,label='wSet')
    wSet(:) = One/nSets
  else if (Key == 'center') then
    RdHead = .false.
    read(LuIn,'(a)',iostat=err) Center
    if (err /= 0) then
      write(u6,*) 'Error while reading input, keyword: ',trim(Key)
      call Quit_OnUserError()
    end if
    call Upcase(Center)
  else if ((Key == 'no threshold') .or. (Key == 'nothreshold')) then
    read(LuIn,*,iostat=err) thr
    if (err /= 0) then
      write(u6,*) 'Error while reading input, keyword: ',trim(Key)
      call Quit_OnUserError()
    end if
  else if ((Key == 'row wise') .or. (Key == 'rowwise')) then
    RdHead = .false.
    rowise = .true.
  else if ((Key == 'lift degeneracy') .or. (Key == 'liftdegeneracy')) then
    RdHead = .false.
    lftdeg = .true.
  else if (Key == 'rydberg') then
    RdHead = .false.
    rydgen = .true.
  else if (Key == 'weights') then
    RdHead = .false.
    read(LuIn,*) (wSet(i),i=1,nSets)
  else if ((Key == 'orbitalweights') .or. (Key == 'orbital weights')) then
    read(LuIn,*) wc0,wc1
    !write(u6,'(a,f12.6)') 'inpctl: wc0',wc0
    !write(u6,'(a,f12.6)') 'inpctl: wc1',wc1
  else if (Key == 'project') then
    RdHead = .false.
    iProj = 1
  else if (Key == 'project 1') then
    RdHead = .false.
    iProj = 1
  else if (Key == 'project 2') then
    RdHead = .false.
    iProj = 2
  else if (Key == 'end of input') then
    RdHead = .false.
    done = .true.
  else if (RdHead) then
    if (Title == 'Atom') Title = KeyIn
    write(u6,*) trim(KeyIn)
  else
    write(u6,*) '*** Illegal keyword: ',KeyIn
    call Quit_OnUserError()
  end if
  if (done) exit
end do

return

end subroutine InpCtl_GenANO
