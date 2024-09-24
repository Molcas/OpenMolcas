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

!#define _DEBUGPRINT_
subroutine Get_Ecorr_dft(nh1,Grad,nGrad,DFTFOCK,F_DFT,D_DS,nBT,nD,KSDFT,Ec_AB)

use OFembed, only: dFMD, Do_Core
use nq_Info, only: Dens_I, Grad_I, Tau_I
use Constants, only: Zero, One

implicit none
integer nh1, nGrad, nBT, nD
#ifdef _DEBUGPRINT_
integer i
#endif
real*8 Grad(nGrad)
character(LEN=4) DFTFOCK
character(LEN=80) KSDFT
real*8 :: F_DFT(nBT,nD), D_DS(nBT,nD)
real*8 Ec_AB
logical Do_MO, Do_TwoEl, Do_Grad
real*8 dFMD_, Func
integer nFckDim

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
Do_Grad = .false.

nFckDim = 2
dFMD_ = dFMD
dFMD = One
!                                                                      *
!***********************************************************************
!                                                                      *
Do_Core = .true.
call Driver(KSDFT,Do_Grad,Func,Grad,nGrad,Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nFckDim,DFTFOCK)
Do_Core = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
Ec_AB = Func

#ifdef _DEBUGPRINT_
write(6,*) ' Correlation energy: ',Ec_AB
write(6,*)
write(6,*) ' Correlation potentials: (itri,F_alpha,F_beta)'
write(6,*)
do i=1,nh1
  write(6,'(i4,3f22.16)') i,F_DFT(i,1),F_DFT(i,2)
end do
#endif

dFMD = dFMD_

return

end subroutine Get_Ecorr_dft
