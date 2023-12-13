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
! Copyright (C) 2010,2012,2017, Francesco Aquilante                    *
!               2015,2017, Alexander Zech                              *
!***********************************************************************

subroutine cWrap_DrvNQ(KSDFT,nFckDim,Func,D_DS,nh1,nD_DS,Do_Grad,Grad,nGrad,DFTFOCK,F_corr)

use OFembed, only: Do_Core
use nq_Info, only: Dens_I, Grad_I, mBas, mIrrep, nAsh, nFro, nIsh, Tau_I
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: KSDFT
integer(kind=iwp), intent(in) :: nFckDim, nh1, nD_DS, nGrad
real(kind=wp), intent(out) :: Func
real(kind=wp), intent(in) :: D_DS(nh1,nD_DS)
logical(kind=iwp), intent(in) :: Do_Grad
real(kind=wp), intent(inout) :: Grad(nGrad), F_corr(nh1,nFckDim)
character(len=4), intent(in) :: DFTFOCK
integer(kind=iwp) :: nOrbA
logical(kind=iwp) :: Do_MO, Do_TwoEl, F_nAsh

!                                                                      *
!***********************************************************************
!                                                                      *
Func = Zero
Dens_I = Zero
Grad_I = Zero
Tau_I = Zero
Do_MO = .false.
Do_TwoEl = .false.

call Get_iScalar('nSym',mIrrep)
call Get_iArray('nBas',mBas(0),mIrrep)
call Get_iArray('nFro',nFro(0),mIrrep)
call Get_iArray('nIsh',nIsh(0),mIrrep)
call qpg_dArray('nAsh',F_nAsh,nOrbA)
if ((.not. F_nAsh) .or. (nOrbA == 0)) then
  nAsh(0:mIrrep-1) = 0
else
  call Get_iArray('nAsh',nAsh(0),mIrrep)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
Do_Core = .true.
call Driver(KSDFT,Do_Grad,Func,Grad,nGrad,Do_MO,Do_TwoEl,D_DS,F_corr,nh1,nFckDim,DFTFOCK)
Do_Core = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine cWrap_DrvNQ
