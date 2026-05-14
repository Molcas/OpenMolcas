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
! Copyright (C) 1995, Martin Schuetz                                   *
!               2022, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Mk_FockAO(nIter_)
!***********************************************************************
!                                                                      *
!     purpose: Update Fock matrix from actual OneHam & TwoHam matrices *
!                                                                      *
!***********************************************************************

use InfSCF, only: FockAO, OneHam, TwoHam, Vxc
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nIter_
integer(kind=iwp) :: i2Hm, iD, nD, NumDT
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nDT
#endif

NumDT = size(TwoHam,3)
nD = size(FockAO,2)

if (nIter_ == 1) then
  i2Hm = 1
else
  i2Hm = NumDT
end if
#ifdef _DEBUGPRINT_
nDT = size(FockAO,1)
write(u6,*) 'i2Hm=',i2Hm
call NrmClc(OneHam,nDT,'UpdFck','OneHam')
call NrmClc(TwoHam(1,1,i2Hm),nDT*nD,'UpdFck','T in i2Hm')
call NrmClc(FockAO,(nDT)*nD,'UpdFck','FockAO')
call NrmClc(Vxc(1,1,i2Hm),nDT*nD,'UpdFck','V in i2Hm  ')
#endif
do iD=1,nD

  ! F = h + G(D), where D is the inter-/extrapolated
  ! densitity matrix. This part is exact.
  !
  ! G(D) = Sum_i c_i G(D_i)
  !
  ! Here we add the contribution to the Fock matrix from
  ! the external field (the none linear or bi-linear parts).
  ! Note that Vxc(D_a+D_b)=/= Vxc(D_a) + Vxc(D_b). However,
  ! using a first order Taylor expansion we can demonstate that
  ! this is correct to first order. Hence, we pick up the
  ! inter-/extra-polated external potential here as generated
  ! by OptClc.f.
  !
  ! Proof:
  !
  ! D = Sum_i c_i * D_i; Sum_i c_i = 1
  !
  ! First order Taylor expansion of Vxc around D_0 (arbitrary)
  !
  ! Vxc(D) = Vxc(D_0) + dVxc(D_0)/dD (D-D_0)
  !
  ! Vxc(D) = Vxc(D_0) + Sum_i C_i dVxc(D_0)/dD (D_i_D_0)
  !
  ! Vxc(D) = Vxc(D_0) + Sum_i C_i (Vxc(D_i) - Vxc(D_0)
  !
  ! Vxc(D) = Sum_i c_i Vxc(D_i), which is produced by OptClc.f

  FockAO(:,iD) = OneHam(:)+TwoHam(:,iD,i2Hm)+Vxc(:,iD,i2Hm)

# ifdef _DEBUGPRINT_
  write(u6,*) 'Fock'
  !write(u6,'(5f12.6)') (FockAO(ivv,iD),ivv=1,nDT)
  call NrmClc(FockAO(1,iD),nDT,'Fock  ','UpdFck ')
  call NrmClc(TwoHam(1,iD,i2Hm),nDT,'TwoHam','UpdFck ')
  call NrmClc(Vxc(1,iD,i2Hm),nDT,'Vxc','UpdFck ')
# endif

end do

end subroutine Mk_FockAO
