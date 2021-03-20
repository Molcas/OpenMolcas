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

subroutine PrintTitles(LuVal,nShowMOs,isDensity,nMOs,iWipGRef,isEner,WipOcc,iWipType,Crypt,iWipNZ,WipE,VBocc,ifpartial,isLine, &
                       isSphere,isColor,ISLUSCUS,nCoor,nBlocks,nInc)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LuVal, nShowMOs, isDensity, nMOs, iWipGRef(*), isEner, iWipType(*), iWipNZ(*), ifpartial, isLine, &
                                 isSphere, isColor, ISLUSCUS, nCoor, nBlocks, nInc
real(kind=iwp), intent(in) :: WipOcc(*), WipE(*), VBocc
character(len=7), intent(in) :: Crypt
integer(kind=iwp) :: i, iActOrb, ib, j, Sizeof8
character(len=12000) :: Line
character(len=10) :: LineT
character :: bb

Sizeof8 = 8
iActOrb = 0
LineT = 'GridName= '
if (isLine == 1) LineT = '#GridName='
!write(u6,*) 'here',nInc, nBlocks, nCoor
if (isLUSCUS == 1) then
  !NFIRST = nInc

  !if (nBlocks == 1) NFIRST=nCoor

  !  ii = 0
  !  do i=1,nShowMOs
  !    write(Line,'(A17,1000i12)') ' File_pointers = ',((nFirst*(j-1)+ii)*Sizeof8,j=1,nBlocks)
  !    call PRINTLINE(LUVAL,LINE,12*nBlocks+17,0)
  !    ii = ii+nFirst*nBlocks
  !  end do
  write(Line,'(A,i22)') ' ORBOFF = ',nCoor*nShowMOs*Sizeof8
  !     & + 28*nShowMOs-20
  ! ?? -20?
  call PRINTLINE(LUVAL,LINE,32,0)
end if
do i=1,nShowMOs-isDensity-isSphere-isColor
  j = iWipGRef(i)
  if (isEner == 1) then
    if (.not.(0 == 1 .and. WipOcc(j) > Zero .and. WipOcc(j) < Two)) then
      ib = iWipType(j)
      bb = ' '
      if (ib > 0 .and. ib < 8) bb = Crypt(ib:ib)
      if (ISLUSCUS == 1) then
        write(LINE,1000) iWipNZ(j),iWipNZ(j+nMOs),WipE(j),WipOcc(j),bb
        call PRINTLINE(LUVAL,LINE,72,0)
      else
        write(line,'(a,i2,i5,f12.4," (",f4.2,")",1x,a)') LineT,iWipNZ(j),iWipNZ(j+nMOs),WipE(j),WipOcc(j),bb
        call PrintLine(LuVal,line,38,0)
      end if
1000  format(1X,'GridName= Orbital sym=',i2,' index=',i5,' Energ=',F12.4,' occ=',F4.2,' type=',a1)
    else
      iActOrb = iActOrb+1
      if (ISLUSCUS == 1) then
        write(LINE,1010) 'VB orbital',iActOrb,' (',VBocc,')'
1010    format(1X,'GridName= VB_orbital iActOrb= ',I4,' occ= ',F4.2)
        call PRINTLINE(LUVAL,LINE,45,0)
      else
        write(line,'(2a,i4,5x,a,f4.2,a)') LineT,'VB orbital',iActOrb,' (',VBocc,')'
        call PrintLine(LuVal,line,38,0)
      end if
    end if
  else
    if (.not.(0 == 1 .and. WipOcc(j) > Zero .and. WipOcc(j) < Two)) then
      ib = iWipType(j)
      bb = ' '
      if (ib > 0 .and. ib < 8) bb = Crypt(ib:ib)
      if (ISLUSCUS == 1) then
        write(LINE,1020) iWipNZ(j),iWipNZ(j+nMOs),WipOcc(j),bb
1020    format(1X,'GridName= Orbital sym=',I2,' index=',I5,' occ=',F4.2,' type=',A1)
        call PRINTLINE(LUVAL,LINE,53,0)
      else
        write(line,'(a,i2,i5," (",f8.6,")",1x,a)') LineT,iWipNZ(j),iWipNZ(j+nMOs),WipOcc(j),bb
        call PrintLine(LuVal,line,30,0)
      end if
    else
      iActOrb = iActOrb+1
      if (ISLUSCUS == 1) then
        write(LINE,1010) iActOrb,VBocc
        call PRINTLINE(LUVAL,LINE,45,0)
      else
        write(line,'(2a,i4,5x,a,f4.2,a)') LineT,'VB orbital',iActOrb,' (',VBocc,')'
        call PrintLine(LuVal,line,30,0)
      end if
    end if
  end if
end do
if (isSphere == 1) then
  write(line,'(a,a)') LineT,'  Sphere '
  call PrintLine(LuVal,line,19,0)
end if
if (isSphere == 1) then
  write(line,'(a,a)') LineT,'0 Color  '
  call PrintLine(LuVal,line,19,0)
end if
if (isDensity == 1) then
  if (ifpartial == 0) then
    if (ISLUSCUS == 1) then
      LINE = ' GridName= Density'
      call PRINTLINE(LUVAL,LINE,18,0)
    else
      write(line,'(a,a)') LineT,'  Density'
      call PrintLine(LuVal,line,19,0)
    end if
  else
    if (ISLUSCUS == 1) then
      LINE = ' GridName= Density (partial)'
      call PRINTLINE(LUVAL,LINE,28,0)
    else
      write(line,'(a,a)') LineT,'  Density (partial)'
      call PrintLine(LuVal,line,29,0)
    end if
  end if
end if
if (ISLUSCUS == 1) then
  LINE = ' <DENSITY>'
  call PRINTLINE(LUVAL,LINE,10,0)
end if

return
! Avoid unused argumet warnings
if (.false.) then
  call Unused_integer(nBlocks)
  call Unused_integer(nInc)
end if

end subroutine PrintTitles
