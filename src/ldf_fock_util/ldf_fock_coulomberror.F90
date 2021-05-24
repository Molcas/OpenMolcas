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

subroutine LDF_Fock_CoulombError(PrintNorm,ComputeF,Mode,PackedD,PackedF,nD,FactC,ip_D,ip_F)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: compute Coulomb error
!
!      F(uv)-Ftilde(uv)= FactC * sum_kl { (uv|kl)-[uv|kl] } * D(kl)
!
!      where [uv|kl] are the LDF integrals.
!
! If ComputF: the LDF Fock matrix is computed here
! Else: on input, ip_F should point to the Fock matrix computed from
! LDF integrals (replaced with the error on exit)!

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
logical(kind=iwp), intent(in) :: PrintNorm, ComputeF, PackedD, PackedF
integer(kind=iwp), intent(in) :: Mode, nD, ip_D(nD), ip_F(nD)
real(kind=wp), intent(in) :: FactC(nD)
logical(kind=iwp) :: Timing, Add
integer(kind=iwp) :: IntegralOption, ip_myF, l_myF, ipF, lF, iD
real(kind=wp) :: ThrPS(2)
real(kind=r8), external :: ddot_
#include "WrkSpc.fh"
#include "localdf_bas.fh"

if (ComputeF) then
  IntegralOption = 0
  Timing = .false.
  ThrPS(1) = Zero
  ThrPS(2) = Zero
  Add = .false.
  call LDF_Fock_CoulombOnly(IntegralOption,Timing,Mode,ThrPS,Add,PackedD,PackedF,nD,FactC,ip_D,ip_F)
end if

if (PackedF) then
  lF = nBas_Valence*(nBas_Valence+1)/2
else
  lF = nBas_Valence**2
end if
l_myF = nD
call GetMem('myFPtr','Allo','Inte',ip_myF,l_myF)
do iD=1,nD
  call GetMem('myF','Allo','Real',ipF,lF)
  iWork(ip_myF-1+iD) = ipF
end do
IntegralOption = 222 ! use conventional integrals
Timing = .false.
ThrPS(1) = Zero
ThrPS(2) = Zero
Add = .false.
call LDF_Fock_CoulombOnly(IntegralOption,Timing,Mode,ThrPS,Add,PackedD,PackedF,nD,FactC,ip_D,iWork(ip_myF))
do iD=1,nD
  ipF = iWork(ip_myF-1+iD)
  call dAXPY_(lF,-One,Work(ipF),1,Work(ip_F(iD)),1)
  call dScal_(lF,-One,Work(ip_F(iD)),1)
end do
do iD=1,nD
  ipF = iWork(ip_myF-1+iD)
  call GetMem('myF','Free','Real',ipF,lF)
end do
call GetMem('myFPtr','Free','Inte',ip_myF,l_myF)

if (PrintNorm) then
  do iD=1,nD
    write(u6,'(A,I10,A,1P,D20.10)') 'Norm of Coulomb error for density',iD,':',sqrt(dDot_(lF,Work(ip_F(iD)),1,Work(ip_F(iD)),1))
  end do
  call xFlush(u6)
end if

end subroutine LDF_Fock_CoulombError
