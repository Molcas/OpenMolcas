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

subroutine Read_h0(nSize,nBas,ip_h0,Restart)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: nSize, nBas, ip_h0
logical(kind=iwp) :: Restart
#include "WrkSpc.fh"
character(len=8) :: Label
integer(kind=iwp) :: iComp, iOpt0, iOpt1, iRc, iSyLbl, nInts(1)

!                                                                      *
!***********************************************************************
!                                                                      *
iOpt0 = 0
iOpt1 = 1

call Allocate_Work(ip_h0,nsize)

iComp = 1
iSyLbl = 1
Label = 'OneHam  '
iRc = -1
if (Restart) then
  call Get_dArray('LoProp H0',Work(ip_h0),nSize)
else
  call iRdOne(iRc,iOpt1,Label,iComp,nInts,iSyLbl)
  if (iRc /= 0) then
    write(u6,*) 'Read_h0: Error reading ONEINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  if (nInts(1)+4 /= nSize) then
    write(u6,*) 'Local_Polar: nInts+4.ne.nSize',nInts(1)+4,nSize
    call Abend()
  end if
  iRc = -1
  call RdOne(iRc,iOpt0,Label,iComp,Work(ip_h0),iSyLbl)
  call Put_dArray('LoProp H0',Work(ip_h0),nSize)
end if
!call TriPrt('H0 ',' ',Work(ip_h0),nBas)
!                                                                      *
!***********************************************************************
!                                                                      *
return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nBas)

end subroutine Read_h0
