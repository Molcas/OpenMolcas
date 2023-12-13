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

function nProp_Int(Do_Index,Idx,nIdx)

use OneDat, only: sOpSiz
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nProp_Int
logical(kind=iwp), intent(in) :: Do_Index
integer(kind=iwp), intent(in) :: nIdx
integer(kind=iwp), intent(inout) :: Idx(4,nIdx)
integer(kind=iwp) :: iCent, iComp, iEF, iMltpl, iopt, irc, iSmLbl, iSmLbl_, maxCen, nComp, n_Int(1)
character(len=8) :: Label

!                                                                      *
!***********************************************************************
!                                                                      *
nProp_Int = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan the ONEINT file for multipole moment operators

!write(u6,*) ' Starting scan of ONEINT for multipole moments'
do iMltpl=1,99
  nComp = (iMltpl+1)*(iMltpl+2)/2

  write(Label,'(a,i2)') 'MLTPL ',iMltpl
  irc = -1
  iopt = ibset(0,sOpSiz)
  iComp = 1
  call iRdOne(irc,iopt,Label,iComp,n_Int,iSmLbl)
  if (irc /= 0) exit
  if (Do_Index) then
    do iComp=1,nComp
      nProp_Int = nProp_Int+1
      Idx(1,nProp_Int) = 1
      Idx(2,nProp_Int) = iMltpl
      Idx(3,nProp_Int) = iComp
      Idx(4,nProp_Int) = 0
    end do
  else
    nProp_Int = nProp_Int+nComp
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan 'ONEINT' for electric field integrals

!write(u6,*) ' Starting scan of ONEINT for various elec. field integrals'
do iEF=0,2
  nComp = (iEF+1)*(iEF+2)/2

  maxCen = 9999
  do iCent=1,maxCen
    write(Label,'(a,i1,i5)') 'EF',iEF,iCent
    irc = -1
    iopt = ibset(0,sOpSiz)
    iComp = 1
    call iRdOne(irc,iopt,Label,iComp,n_Int,iSmLbl)
    if (irc /= 0) exit
    if (Do_Index) then
      do iComp=1,nComp
        nProp_Int = nProp_Int+1
        Idx(1,nProp_Int) = 2
        Idx(2,nProp_Int) = iEF
        Idx(3,nProp_Int) = iComp
        Idx(4,nProp_Int) = iCent
      end do
    else
      nProp_Int = nProp_Int+nComp
    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan 'ONEINT' for contact term integrals

!write(u6,*) ' Starting scan of ONEINT for various contact term integrals'
nComp = 1

maxCen = 9999
do iCent=1,maxCen
  write(Label,'(a,i5)') 'Cnt',iCent
  irc = -1
  iopt = ibset(0,sOpSiz)
  iComp = 1
  call iRdOne(irc,iopt,Label,iComp,n_Int,iSmLbl)
  if (irc /= 0) exit
  if (Do_Index) then
    do iComp=1,nComp
      nProp_Int = nProp_Int+1
      Idx(1,nProp_Int) = 3
      Idx(2,nProp_Int) = 1
      Idx(3,nProp_Int) = 1
      Idx(4,nProp_Int) = iCent
    end do
  else
    nProp_Int = nProp_Int+nComp
  end if
end do
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
  iopt = ibset(0,sOpSiz)
  iComp = 1
  call iRdOne(irc,iopt,Label,iComp,n_Int,iSmLbl_)
  if (irc /= 0) exit
  if (Do_Index) then
    do iComp=1,nComp
      nProp_Int = nProp_Int+1
      Idx(1,nProp_Int) = 4
      Idx(2,nProp_Int) = 0
      Idx(3,nProp_Int) = iComp
      Idx(4,nProp_Int) = iCent
    end do
  else
    nProp_Int = nProp_Int+nComp
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!if (Do_Index) call iVcPrt('Idx',' ',Idx,4*nProp_Int)

return

end function nProp_Int
