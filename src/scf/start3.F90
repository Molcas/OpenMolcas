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

use InfSCF, only: nBas, nBB, nBO, nBT, nSym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mBB, nD, mBT
real(kind=wp), intent(out) :: CMO(mBB,nD), TrM(mBB,nD)
real(kind=wp), intent(in) :: OneHam(mBT), Ovrlp(mBT)
integer(kind=iwp) :: iD, iSym, nBasX(8), nSymX
character(len=8) :: Location

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

Location = 'Start3'

! Compute transformation matrix
do iD=1,nD
  call TrGen(TrM(:,iD),nBB,Ovrlp,OneHam,nBT)
  CMO(1:nBO,iD) = TrM(1:nBO,iD)
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
!call mma_allocate(Dens,nBT,nD,Label='Dens')
!call Get_dArray_chk('D1AO',Dens(:,1),nBT)
!if (nD == 2) then
!  call Get_dArray_chk('D1sao',Dens(:,2),nBT)
!  ! now we need to fix interface - actually we read a+b,a-b
!  call mma_allocate(Tmp,nBT,Label='Tmp')
!  Tmp(:) = Half*(Dens(:,1)-Dens(:,2))
!  Dens(:,1) = Half*(Dens(:,1)+Dens(:,2))
!  Dens(:,2) = Tmp(:)
!  call mma_deallocate(Tmp)
!end if
!call mma_deallocate(Dens)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine Start3
