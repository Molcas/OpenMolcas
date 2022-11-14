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
! Copyright (C) 2011, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LDF_Fock_CoulombError(PrintNorm,ComputeF,Mode,PackedD,PackedF,nD,FactC,ip_D,F)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: compute Coulomb error
!
!      F(uv)-Ftilde(uv)= FactC * sum_kl { (uv|kl)-[uv|kl] } * D(kl)
!
!      where [uv|kl] are the LDF integrals.
!
! If ComputF: the LDF Fock matrix is computed here
! Else: on input, F should be the Fock matrix computed from
! LDF integrals (replaced with the error on exit)!

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: PrintNorm, ComputeF, PackedD, PackedF
integer(kind=iwp), intent(in) :: Mode, nD, ip_D(nD)
real(kind=wp), intent(inout) :: FactC(nD), F(*)
logical(kind=iwp) :: Timing, Add
integer(kind=iwp) :: IntegralOption, ipF, lF, iD
real(kind=wp), allocatable :: myF(:)
real(kind=wp) :: ThrPS(2)
real(kind=wp), external :: ddot_
#include "localdf_bas.fh"

if (ComputeF) then
  IntegralOption = 0
  Timing = .false.
  ThrPS(1) = Zero
  ThrPS(2) = Zero
  Add = .false.
  call LDF_Fock_CoulombOnly(IntegralOption,Timing,Mode,ThrPS,Add,PackedD,PackedF,nD,FactC,ip_D,F)
end if

if (PackedF) then
  lF = nBas_Valence*(nBas_Valence+1)/2
else
  lF = nBas_Valence**2
end if
call mma_allocate(myF,nD*lF,label='myF')
IntegralOption = 222 ! use conventional integrals
Timing = .false.
ThrPS(1) = Zero
ThrPS(2) = Zero
Add = .false.
call LDF_Fock_CoulombOnly(IntegralOption,Timing,Mode,ThrPS,Add,PackedD,PackedF,nD,FactC,ip_D,myF)
do iD=1,nD
  ipF = (iD-1)*lF+1
  call dAXPY_(lF,-One,myF(ipF),1,F(ipF),1)
  call dScal_(lF,-One,F(ipF),1)
end do
call mma_deallocate(myF)

if (PrintNorm) then
  do iD=1,nD
    ipF = (iD-1)*lF+1
    write(u6,'(A,I10,A,1P,D20.10)') 'Norm of Coulomb error for density',iD,':',sqrt(dDot_(lF,F(ipF),1,F(ipF),1))
  end do
  call xFlush(u6)
end if

end subroutine LDF_Fock_CoulombError
