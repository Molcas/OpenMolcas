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
! Copyright (C) 1991, Roland Lindh                                     *
!               2008, Giovanni Ghigo                                   *
!***********************************************************************

subroutine Def_CtoF(lNew)
!***********************************************************************
!                                                                      *
!     Author: Giovanni Ghigo, Dep. of General and Organic Chemistry    *
!             University of Torino, ITALY                              *
!             July 2008                                                *
!     Adapted from  DefInt by                                          *
!             Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: lNew
integer(kind=iwp) :: i1, i2, i3, iEnd, iFrst, jEnd, Lu_UDIC, mCntr, msAtom, nCntr, neq, nGo, nTemp
real(kind=wp) :: Tmp
character(len=120) :: Line, Temp
character(len=16) :: filnam
character(len=8) :: Frmt, Labels
character(len=6) :: Typ
integer(kind=iwp), allocatable :: Ind(:,:)
real(kind=wp), allocatable :: Mass(:,:), Temp2(:,:), xyz(:,:)

nTemp = len(Temp)
write(Frmt,'(A,I3.3,A)') '(F',nTemp,'.0)'

Lu_UDIC = 91
filnam = 'UDIC'
call molcas_open(Lu_UDIC,filnam)
!open(Lu_UDIC,File=filnam,Form='Formatted',Status='OLD')
rewind(Lu_UDIC)
write(u6,*)
write(u6,*) '****************************************************************'
if (lNew) then
  write(u6,*) '* New value of the internal coordinate to follow               *'
else
  write(u6,*) '* Original value of the internal coordinate to follow          *'
end if
write(u6,*) '****************************************************************'

! Step 1. BSet up the b vectors from which we will define the
! internal coordinates.

!iBVct = 0
read(Lu_UDIC,'(A)') Line
Temp = Line
call UpCase(Temp)

! Move the label of the internal coordinate

neq = index(Line,'=')
if (neq == 0) then
  call WarningMessage(2,'Error in Def_CTOF')
  write(u6,'(A)') '***********************************'
  write(u6,'(A)') ' Syntax error in line :            '
  write(u6,'(A)') Line(1:33),'...'
  write(u6,'(A)') '***********************************'
  call Quit_OnUserError()
else
  iFrst = 1
  call NxtWrd(Line,iFrst,iEnd)
  jEnd = iEnd
  if (Line(iEnd:iEnd) == '=') jEnd = jEnd-1
  if (jEnd-iFrst+1 > 8) then
    call WarningMessage(2,'Error in Def_CTOF')
    write(u6,'(A)') '***********************************'
    write(u6,'(A)') ' Syntax error in line :            '
    write(u6,'(A)') Line(1:33),'...'
    write(u6,'(A,A)') Line(iFrst:jEnd),' has more than 8 character'
    write(u6,'(A)') '***********************************'
    call Quit_OnUserError()
  end if
  Labels = Line(iFrst:jEnd)
end if

! Construct the corresponding transformation vector

mCntr = 0
if (index(Temp,'CART') /= 0) then
  nCntr = 1
  nGo = index(Temp,'CART')
  nGo = nGo-1+index(Temp(nGo:nTemp),' ')
  if (index(Temp(nGo:nTemp),'X') /= 0) then
    nGo = nGo-1+index(Temp(nGo:nTemp),'X')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    Typ = 'X     '
  else if (index(Temp(nGo:nTemp),'Y') /= 0) then
    nGo = nGo-1+index(Temp(nGo:nTemp),'Y')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    Typ = 'Y     '
  else if (index(Temp(nGo:nTemp),'Z') /= 0) then
    nGo = nGo-1+index(Temp(nGo:nTemp),'Z')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    Typ = 'Z     '
  else
    nGo = -1
    call WarningMessage(2,'Error in Def_CTOF')
    write(u6,*) 'DefInt: wrong cartesian type'
    write(u6,'(A,A)') 'Temp=',Temp
    call Quit_OnUserError()
  end if
else if (index(Temp,'BOND') /= 0) then
  nCntr = 2
  nGo = index(Temp,'BOND')
  nGo = nGo-1+index(Temp(nGo:nTemp),' ')
  Typ = 'STRTCH'
else if (index(Temp,'LANGLE(2)') /= 0) then
  nCntr = 3
  nGo = index(Temp,'LANGLE(2)')
  nGo = nGo-1+index(Temp(nGo:nTemp),' ')
  Typ = 'LBEND2'
else if (index(Temp,'LANGLE(1)') /= 0) then
  nCntr = 3
  nGo = index(Temp,'LANGLE(1)')
  nGo = nGo-1+index(Temp(nGo:nTemp),' ')
  Typ = 'LBEND1'
else if (index(Temp,'ANGL') /= 0) then
  nCntr = 3
  nGo = index(Temp,'ANGL')
  nGo = nGo-1+index(Temp(nGo:nTemp),' ')
  Typ = 'BEND  '
else if (index(Temp,'DIHE') /= 0) then
  nCntr = 4
  nGo = index(Temp,'DIHE')
  nGo = nGo-1+index(Temp(nGo:nTemp),' ')
  Typ = 'TRSN  '
else if (index(Temp,'OUTO') /= 0) then
  nCntr = 4
  nGo = index(Temp,'OUTO')
  nGo = nGo-1+index(Temp(nGo:nTemp),' ')
  Typ = 'OUTOFP'
else if (index(Temp,'DISS') /= 0) then
  i1 = index(Line,'(')
  i2 = index(Line,'+')
  i3 = index(Line,')')
  if ((i1 >= i2) .or. (i2 >= i3)) then
    call WarningMessage(2,'Error in Def_CTOF')
    write(u6,*) ' Line contains syntax error!'
    write(u6,'(A)') Line
    write(u6,*) i1,i2,i3
    call Quit_OnUserError()
  end if
  nGo = i3+1
  Temp = Line(i1+1:i2-1)
  read(Temp,Frmt) Tmp
  nCntr = nint(Tmp)
  Temp = Line(i2+1:i3-1)
  read(Temp,Frmt) Tmp
  mCntr = nint(Tmp)
  Typ = 'DISSOC'
else
  nGo = -1
  call WarningMessage(2,'Error in Def_CTOF')
  write(u6,*) ' Line contains syntax error!'
  write(u6,'(A)') Line
  call Quit_OnUserError()
end if

msAtom = nCntr+mCntr
call mma_allocate(xyz,3,msAtom,Label='xyz')
call mma_allocate(Temp2,3,msAtom,Label='Temp2')
call mma_allocate(Ind,2,msAtom,Label='Ind')
call mma_allocate(Mass,2,msAtom,Label='Mass')

call CllCtoF(Line(nGo:nTemp),nCntr,mCntr,xyz,Temp2,Ind,Typ,Mass,Labels)

call mma_deallocate(Mass)
call mma_deallocate(Ind)
call mma_deallocate(Temp2)
call mma_deallocate(xyz)

close(Lu_UDIC)

return

end subroutine Def_CtoF
