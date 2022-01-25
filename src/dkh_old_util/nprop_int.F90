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

integer function nProp_Int(Do_Index,Index,nIndex)

implicit real*8(a-h,o-z)
#include "real.fh"
character*8 Label
logical Do_Index
integer index(4,nIndex)
dimension nint(1)

!                                                                      *
!***********************************************************************
!                                                                      *
nProp_Int = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan the ONEINT file for multipole moment operators

!write(6,*) ' Starting scan of ONEINT for multipole moments'
do iMltpl=1,99
  nComp = (iMltpl+1)*(iMltpl+2)/2

  write(Label,'(a,i2)') 'MLTPL ',iMltpl
  irc = -1
  iopt = 1
  iComp = 1
  call iRdOne(irc,iopt,Label,iComp,nInt,iSmLbl)
  if (irc /= 0) Go To 110
  if (Do_Index) then
    do iComp=1,nComp
      nProp_Int = nProp_Int+1
      index(1,nProp_Int) = 1
      index(2,nProp_Int) = iMltpl
      index(3,nProp_Int) = iComp
      index(4,nProp_Int) = 0
    end do
  else
    nProp_Int = nProp_Int+nComp
  end if
end do
110 continue
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan 'ONEINT' for electric field integrals

!write(6,*) ' Starting scan of ONEINT for various elec. field integrals'
do iEF=0,2
  nComp = (iEF+1)*(iEF+2)/2

  maxCen = 9999
  do iCent=1,maxCen
    write(Label,'(a,i1,i5)') 'EF',iEF,iCent
    irc = -1
    iopt = 1
    iComp = 1
    call iRdOne(irc,iopt,Label,iComp,nInt,iSmLbl)
    if (irc /= 0) Go To 201
    if (Do_Index) then
      do iComp=1,nComp
        nProp_Int = nProp_Int+1
        index(1,nProp_Int) = 2
        index(2,nProp_Int) = iEF
        index(3,nProp_Int) = iComp
        index(4,nProp_Int) = iCent
      end do
    else
      nProp_Int = nProp_Int+nComp
    end if
  end do
201 continue
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan 'ONEINT' for contact term integrals

!write(6,*) ' Starting scan of ONEINT for various contact term integrals'
nComp = 1

maxCen = 9999
do iCent=1,maxCen
  write(Label,'(a,i5)') 'Cnt',iCent
  irc = -1
  iopt = 1
  iComp = 1
  call iRdOne(irc,iopt,Label,iComp,nInt,iSmLbl)
  if (irc /= 0) Go To 301
  if (Do_Index) then
    do iComp=1,nComp
      nProp_Int = nProp_Int+1
      index(1,nProp_Int) = 3
      index(2,nProp_Int) = 1
      index(3,nProp_Int) = 1
      index(4,nProp_Int) = iCent
    end do
  else
    nProp_Int = nProp_Int+nComp
  end if
end do
301 continue
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan 'ONEINT' for magnetic integrals (copied from ONEREL)

nComp = 9
iSmLbl_ = 255
maxCen = 9999
do iCent=1,maxCen
  write(Label,'(A,I3)') 'MAGXP',iCent
  irc = -1
  iopt = 1
  iComp = 1
  call iRdOne(irc,iopt,Label,iComp,nInt,iSmLbl_)
  if (irc /= 0) Go To 401
  if (Do_Index) then
    do iComp=1,nComp
      nProp_Int = nProp_Int+1
      index(1,nProp_Int) = 4
      index(2,nProp_Int) = 0
      index(3,nProp_Int) = iComp
      index(4,nProp_Int) = iCent
    end do
  else
    nProp_Int = nProp_Int+nComp
  end if
end do
401 continue
!                                                                      *
!***********************************************************************
!                                                                      *
!if (Do_Index) call iVcPrt('Index',' ',Index,4*nProp_Int)

return

end function nProp_Int
