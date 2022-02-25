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

subroutine Wrap_DrvNQ(KSDFT,F_DFT,nFckDim,Func,D_DS,nh1,nD_DS,Do_Grad,Grad,nGrad,DFTFOCK)

use nq_Info

implicit real*8(a-h,o-z)
character*(*) KSDFT
integer nh1, nFckDim, nD_DS
real*8 F_DFT(nh1,nFckDim), D_DS(nh1,nD_DS), Func
logical Do_Grad
real*8 Grad(nGrad)
character*4 DFTFOCK
#include "real.fh"
#include "debug.fh"
logical Do_MO, Do_TwoEl, F_nAsh

!                                                                      *
!***********************************************************************
!                                                                      *
! DFT functionals, compute integrals over the potential

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
call qpg_iArray('nAsh',F_nAsh,nOrbA)
if ((.not. F_nAsh) .or. (nOrbA == 0)) then
  call Izero(nAsh(0),mIrrep)
else
  call Get_iArray('nAsh',nAsh(0),mIrrep)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Driver(KSDFT,Do_Grad,Func,Grad,nGrad,Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nFckDim,DFTFOCK)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Wrap_DrvNQ
