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

subroutine PrintTitles(LuVal,nShowMOs,isDensity,nMOs,GRef,isEner,Occ,iType,Crypt,NZ,E,ifpartial,isLine,isSphere,isColor,isLuscus, &
                       nCoor,nBlocks,nInc)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LuVal, nShowMOs, nMOs, GRef(*), iType(*), NZ(*), nCoor, nBlocks, nInc
logical(kind=iwp), intent(in) :: isDensity, isEner, ifpartial, isLine, isSphere, isColor, isLuscus
real(kind=wp), intent(in) :: Occ(*), E(*)
character(len=7), intent(in) :: Crypt
integer(kind=iwp) :: i, ib, j, Sizeof8 !, iActOrb
character(len=128) :: Line
character(len=10) :: LineT
character :: bb
#include "macros.fh"
unused_var(nBlocks)
unused_var(nInc)

Sizeof8 = 8
!iActOrb = 0
LineT = 'GridName= '
if (isLine) LineT = '#GridName='
!write(u6,*) 'here',nInc, nBlocks, nCoor
if (isLuscus) then
  !NFIRST = nInc

  !if (nBlocks == 1) NFIRST=nCoor

  !  ii = 0
  !  do i=1,nShowMOs
  !    write(Line,'(A17,1000i12)') ' File_pointers = ',((nFirst*(j-1)+ii)*Sizeof8,j=1,nBlocks)
  !    call PRINTLINE(LUVAL,LINE,12*nBlocks+17,.false.)
  !    ii = ii+nFirst*nBlocks
  !  end do
  write(Line,'(A,i22)') ' ORBOFF = ',nCoor*nShowMOs*Sizeof8
  !     & + 28*nShowMOs-20
  ! ?? -20?
  call PRINTLINE(LUVAL,LINE,32,.false.)
end if
do i=1,nShowMOs-merge(1,0,isDensity)-merge(1,0,isSphere)-merge(1,0,isColor)
  j = GRef(i)
  if (isEner) then
    !if ((Occ(j) > Zero) .and. (Occ(j) < Two)) then
    !  iActOrb = iActOrb+1
    !  if (isLuscus) then
    !    write(LINE,1010) 'VB orbital',iActOrb,' (',VBocc,')'
    !    call PRINTLINE(LUVAL,LINE,45,.false.)
    !  else
    !    write(line,'(2a,i4,5x,a,f4.2,a)') LineT,'VB orbital',iActOrb,' (',VBocc,')'
    !    call PrintLine(LuVal,line,38,.false.)
    !  end if
    !else
    ib = iType(j)
    bb = ' '
    if ((ib > 0) .and. (ib < 8)) bb = Crypt(ib:ib)
    if (isLuscus) then
      write(LINE,1000) NZ(j),NZ(j+nMOs),E(j),Occ(j),bb
      call PRINTLINE(LUVAL,LINE,72,.false.)
    else
      write(line,'(a,i2,i5,f12.4," (",f4.2,")",1x,a)') LineT,NZ(j),NZ(j+nMOs),E(j),Occ(j),bb
      call PrintLine(LuVal,line,38,.false.)
    end if
    !end if
  else
    !if ((Occ(j) > Zero) .and. (Occ(j) < Two)) then
    !  iActOrb = iActOrb+1
    !  if (isLuscus) then
    !    write(LINE,1010) iActOrb,VBocc
    !    call PRINTLINE(LUVAL,LINE,45,.false.)
    !  else
    !    write(line,'(2a,i4,5x,a,f4.2,a)') LineT,'VB orbital',iActOrb,' (',VBocc,')'
    !    call PrintLine(LuVal,line,30,.false.)
    !  end if
    !else
    ib = iType(j)
    bb = ' '
    if ((ib > 0) .and. (ib < 8)) bb = Crypt(ib:ib)
    if (isLuscus) then
      write(LINE,1020) NZ(j),NZ(j+nMOs),Occ(j),bb
      call PRINTLINE(LUVAL,LINE,53,.false.)
    else
      write(line,'(a,i2,i5," (",f8.6,")",1x,a)') LineT,NZ(j),NZ(j+nMOs),Occ(j),bb
      call PrintLine(LuVal,line,30,.false.)
    end if
    !end if
  end if
end do
if (isSphere) then
  write(line,'(a,a)') LineT,'  Sphere '
  call PrintLine(LuVal,line,19,.false.)
end if
if (isColor) then
  write(line,'(a,a)') LineT,'0 Color  '
  call PrintLine(LuVal,line,19,.false.)
end if
if (isDensity) then
  if (ifpartial) then
    if (isLuscus) then
      LINE = ' GridName= Density (partial)'
      call PRINTLINE(LUVAL,LINE,28,.false.)
    else
      write(line,'(a,a)') LineT,'  Density (partial)'
      call PrintLine(LuVal,line,29,.false.)
    end if
  else
    if (isLuscus) then
      LINE = ' GridName= Density'
      call PRINTLINE(LUVAL,LINE,18,.false.)
    else
      write(line,'(a,a)') LineT,'  Density'
      call PrintLine(LuVal,line,19,.false.)
    end if
  end if
end if
if (isLuscus) then
  LINE = ' <DENSITY>'
  call PRINTLINE(LUVAL,LINE,10,.false.)
end if

return

1000 format(1X,'GridName= Orbital sym=',i2,' index=',i5,' Energ=',F12.4,' occ=',F4.2,' type=',a1)
!1010 format(1X,'GridName= VB_orbital iActOrb= ',I4,' occ= ',F4.2)
1020 format(1X,'GridName= Orbital sym=',I2,' index=',I5,' occ=',F4.2,' type=',A1)

end subroutine PrintTitles
