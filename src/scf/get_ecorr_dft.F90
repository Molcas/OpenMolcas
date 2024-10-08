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
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nh1, nGrad, nBT, nD
real(kind=wp), intent(inout) :: Grad(nGrad), F_DFT(nBT,nD)
character(len=4), intent(in) :: DFTFOCK
real(kind=wp), intent(in) :: D_DS(nBT,nD)
character(len=80), intent(in) :: KSDFT
real(kind=wp), intent(out) :: Ec_AB
integer(kind=iwp) :: nFckDim
real(kind=wp) :: dFMD_
logical(kind=iwp) :: Do_Grad, Do_MO, Do_TwoEl
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
! DFT functionals, compute integrals over the potential

Ec_AB = Zero
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
call Driver(KSDFT,Do_Grad,Ec_AB,Grad,nGrad,Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nFckDim,DFTFOCK)
Do_Core = .false.
!                                                                      *
!***********************************************************************
!                                                                      *

#ifdef _DEBUGPRINT_
write(u6,*) ' Correlation energy: ',Ec_AB
write(u6,*)
write(u6,*) ' Correlation potentials: (itri,F_alpha,F_beta)'
write(u6,*)
do i=1,nh1
  write(u6,'(i4,3f22.16)') i,F_DFT(i,1),F_DFT(i,2)
end do
#endif

dFMD = dFMD_

return

end subroutine Get_Ecorr_dft
