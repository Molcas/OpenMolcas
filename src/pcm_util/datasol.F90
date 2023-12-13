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

subroutine DataSol(IDSolv)
! Database of optical and physical data for various solvent.

use Solvent_Data, only: Init_Solvent_Data, SolvData
use rctfld_module, only: DerEPS, Eps, EPS_USER, EpsInf, EpsInf_USER, kT, MXA, nTT, rDiff, rSolv, rWT, TCE, vMol
use Constants, only: Zero, One
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IDSolv
integer(kind=iwp) :: i

call Init_Solvent_Data()

!Tabs = SolvData(IDSolv)%Tabs
Eps = SolvData(IDSolv)%Eps
EpsInf = SolvData(IDSolv)%EpsInf
DerEps = SolvData(IDSolv)%DerEps
RSolv = SolvData(IDSolv)%RSolv
VMol = SolvData(IDSolv)%VMol
TCE = SolvData(IDSolv)%TCE
!STen = SolvData(IDSolv)%STen
!DSTen = SolvData(IDSolv)%DSTen
!CMF = SolvData(IDSolv)%CMF
! Atomic parameters for dispersion and repulsion
!Rho = SolvData(IDSolv)%Rho
if (size(SolvData(IDSolv)%Atoms) > MxA) then
  call WarningMessage(2,'DataSol: num. solv. atoms > MxA')
  call Abend()
end if
do i=1,size(SolvData(IDSolv)%Atoms)
  if (SolvData(IDSolv)%Atoms(i)%NTT == 0) then
    !NATyp = i-1
    exit
  end if
  NTT(i) = SolvData(IDSolv)%Atoms(i)%NTT
  RDiff(i) = SolvData(IDSolv)%Atoms(i)%RDiff
  KT(i) = SolvData(IDSolv)%Atoms(i)%KT
  RWT(i) = SolvData(IDSolv)%Atoms(i)%RWT
end do

! Use user specified value of the dielectric constant

if (Eps_User /= -One) Eps = Eps_User
if (EpsInf_User /= Zero) EpsInf = EpsInf_User

return

end subroutine DataSol
