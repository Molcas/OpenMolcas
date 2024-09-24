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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine Start3(CMO,TrM,mBB,nD,OneHam,Ovrlp,mBT)
!***********************************************************************
!                                                                      *
!     purpose: Get starting orbitals from density matrix read as input.*
!                                                                      *
!***********************************************************************

use InfSCF, only: nBB, nBO, nBT, nSym, nBas
use Constants, only: Half

implicit none
integer mBB, nD, mBT
real*8 CMO(mBB,nD), TrM(mBB,nD), OneHam(mBT), Ovrlp(mBT), Dens(mBT,nD)
! Define local variables
integer nBasX(8), i, iD, iSym, nSymX
character(len=8) Location
real*8 ra, rb
#include "SysDef.fh"

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

Location = 'Start3'

! Compute transformation matrix
do iD=1,nD
  call TrGen(TrM(1,iD),nBB,Ovrlp,OneHam,nBT)
  call DCopy_(nBO,TrM(1,iD),1,CMO(1,iD),1)
end do

! read old system definitions
! check for compatibility of the calculations

!call get_iScalar('nSym',nSymX)
call Peek_iScalar('nSym',nSymX)
if (nSymX /= nSym) then
  call SysWarnMsg(Location,'Error inconsistent number of Irreps',' ')
  call SysCondMsg('nSymX=nSym',nSymX,'/=',nSym)
end if
call Get_iArray('nBas',nBasX,nSymX)
do iSym=1,nSym
  if (nBasX(iSym) /= nBas(iSym)) then
    call SysWarnMsg(Location,'Error inconsistent nBas',' ')
    call SysCondMsg('nBasX(iSym)=nBas (iSym)',nBasX(iSym),'/=',nBas(iSym))
  end if
end do

! read old density matrix
call Get_dArray_chk('D1AO',Dens(1,1),nBT)
if (nD == 2) then
  call Get_dArray_chk('D1sao',Dens(1,2),nBT)
  ! now we need to fix interface - actually we read a+b,a-b
  do i=1,nBT
    ra = Half*(Dens(i,1)+Dens(i,2))
    rb = Half*(Dens(i,1)-Dens(i,2))
    Dens(i,1) = ra
    Dens(i,2) = rb
  end do

end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine Start3
