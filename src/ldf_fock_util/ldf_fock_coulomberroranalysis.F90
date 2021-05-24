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

subroutine LDF_Fock_CoulombErrorAnalysis(ComputeF,Mode,PackedD,PackedF,nD,FactC,ip_D,ip_F)
! Thomas Bondo Pedersen, January 2011.
!
! Purpose: analyze Coulomb error
!
!      F(uv)-Ftilde(uv)= FactC * sum_kl { (uv|kl)-[uv|kl] } * D(kl)
!
!      where [uv|kl] are the LDF integrals, and compare to the upper bound.
!
! If ComputF: the LDF Fock matrix is computed here
! Else: on input, ip_F should point to the Fock matrix computed from
! LDF integrals (replaced with the error on exit)!

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
logical(kind=iwp), intent(in) :: ComputeF, PackedD, PackedF
integer(kind=iwp), intent(in) :: Mode, nD, ip_D(nD), ip_F(nD)
real(kind=wp), intent(in) :: FactC(nD)
real(kind=wp) :: Stat(7,3), RMS1, RMS2, RMS3
logical(kind=iwp) :: Add, Packed_myF
integer(kind=iwp) :: ip_myF, l_myF, ipF, lF, iD, i
character(len=29), parameter :: SecNam = 'LDF_Fock_CoulombErrorAnalysis'
logical(kind=iwp), parameter :: PrintNorm = .false.
real(kind=r8), external :: ddot_
#include "WrkSpc.fh"
#include "localdf_bas.fh"

if (nD < 1) return

! Compute error
call LDF_Fock_CoulombError(PrintNorm,ComputeF,Mode,PackedD,PackedF,nD,FactC,ip_D,ip_F)

! Compute upper bound
Add = .false.
Packed_myF = PackedF
l_myF = nD
call GetMem('CEAmyFP','Allo','Inte',ip_myF,l_myF)
if (Packed_myF) then
  lF = nBas_Valence*(nBas_Valence+1)/2
else
  lF = nBas_Valence**2
end if
do iD=1,nD
  call GetMem('CEAmyF','Allo','Real',ipF,lF)
  iWork(ip_myF-1+iD) = ipF
end do
call LDF_Fock_CoulombUpperBound_Full(PrintNorm,Add,PackedD,Packed_myF,nD,FactC,ip_D,iWork(ip_myF))

! Analysis
call Cho_Head('Coulomb Error','-',80,u6)
do iD=1,nD
  call Statistics(Work(iWork(ip_myF-1+iD)),lF,Stat(1,1),1,2,3,4,5,6,7)
  RMS1 = dDot_(lF,Work(iWork(ip_myF-1+iD)),1,Work(iWork(ip_myF-1+iD)),1)
  call Statistics(Work(ip_F(iD)),lF,Stat(1,2),1,2,3,4,5,6,7)
  RMS2 = dDot_(lF,Work(ip_F(iD)),1,Work(ip_F(iD)),1)
  do i=1,lF
    Work(iWork(ip_myF-1+iD)-1+i) = Work(iWork(ip_myF-1+iD)-1+i)-abs(Work(ip_F(iD)-1+i))
  end do
  call Statistics(Work(iWork(ip_myF-1+iD)),lF,Stat(1,3),1,2,3,4,5,6,7)
  RMS3 = dDot_(lF,Work(iWork(ip_myF-1+iD)),1,Work(iWork(ip_myF-1+iD)),1)
  write(u6,'(/,2X,A,I10,A)') 'Coulomb error for density',iD,' (Upper bound,Actual,Diff):'
  write(u6,'(2X,A,1P,3D20.10)') 'Average error......',Stat(1,1),Stat(1,2),Stat(1,3)
  write(u6,'(2X,A,1P,3D20.10)') 'Abs average error..',Stat(2,1),Stat(2,2),Stat(2,3)
  write(u6,'(2X,A,1P,3D20.10)') 'Min error..........',Stat(3,1),Stat(3,2),Stat(3,3)
  write(u6,'(2X,A,1P,3D20.10)') 'Max error..........',Stat(4,1),Stat(4,2),Stat(4,3)
  write(u6,'(2X,A,1P,3D20.10)') 'Max abs error......',Stat(5,1),Stat(5,2),Stat(5,3)
  write(u6,'(2X,A,1P,3D20.10)') 'Variance...........',Stat(6,1),Stat(6,2),Stat(6,3)
  write(u6,'(2X,A,1P,3D20.10)') 'Norm...............',sqrt(RMS1),sqrt(RMS2),sqrt(RMS3)
  if (lF > 0) then
    RMS1 = sqrt(RMS1/real(lF,kind=wp))
    RMS2 = sqrt(RMS2/real(lF,kind=wp))
    RMS3 = sqrt(RMS3/real(lF,kind=wp))
  else
    RMS1 = Zero
    RMS2 = Zero
    RMS3 = Zero
  end if
  write(u6,'(2X,A,1P,3D20.10)') 'RMS error..........',RMS1,RMS2,RMS3
  call xFlush(u6)
  if ((Stat(5,1)-Stat(5,2)) < Zero) then
    if (abs(Stat(5,1)-Stat(5,2)) > 1.0e-6_wp) then
      call WarningMessage(2,SecNam//': max abs error is greater than upper bound!')
      call LDF_Quit(1)
    end if
  end if
end do

! Deallocations
do iD=1,nD
  ipF = iWork(ip_myF-1+iD)
  call GetMem('CEAmyF','Free','Real',ipF,lF)
end do
call GetMem('CEAmyFP','Free','Inte',ip_myF,l_myF)

end subroutine LDF_Fock_CoulombErrorAnalysis
