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
! Copyright (C) Luca De Vico                                           *
!***********************************************************************

subroutine Get_Nuc_Charge_All_(Coord_Unique,Charges_Unique,nUnique_Atoms,Charges_All,nAll_Atoms)

use Symmetry_Info, only: iOper, nIrrep, Symmetry_Info_Get
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nUnique_Atoms, nAll_Atoms
real(kind=wp) :: Coord_Unique(3,nUnique_Atoms), Charges_Unique(nUnique_Atoms), Charges_All(nAll_Atoms)
integer(kind=iwp) :: Active = 0, iAll_Atom, iChAtom, iChxyz, iCo, iCoSet(0:7,0:7), iGen(3), iStab(0:7), iUnique_Atom, MaxDCR, &
                     nCoSet, nGen, nStab
real(kind=wp) :: Charge_Old

!write(u6,*) 'Enter Get_Nuc_Charge_All_'
!                                                                      *
!***********************************************************************
!                                                                      *
if (Active == 0) then
  call Symmetry_Info_Get()
  Active = 1
end if
!write(u6,*) 'Get_Nuc_Charge_All_: nIrrep=',nIrrep
!write(u6,*) 'Get_Nuc_Charge_All_: iOper=',(iOper(i),i=0,nIrrep-1)
!                                                                      *
!***********************************************************************
!                                                                      *
nGen = 0
if (nIrrep == 2) nGen = 1
if (nIrrep == 4) nGen = 2
if (nIrrep == 8) nGen = 3
if (nGen >= 1) iGen(1) = iOper(1)
if (nGen >= 2) iGen(2) = iOper(2)
if (nGen == 3) iGen(3) = iOper(4)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of all nuclear charges

iAll_Atom = 0
MaxDCR = 0
do iUnique_Atom=1,nUnique_Atoms
  iChAtom = iChxyz(Coord_Unique(1,iUnique_Atom),iGen,nGen)
  call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
  nCoSet = nIrrep/nStab

  Charge_Old = Charges_Unique(iUnique_Atom)

  do iCo=0,nCoSet-1
    iAll_Atom = iAll_Atom+1
    Charges_All(iAll_Atom) = Charge_Old
  end do

end do

!call RecPrt('Charges_Unique',' ',Charges_Unique,1,nUnique_Atoms)
!call RecPrt('Charges_All',' ',Charges_All,1,nAll_Atoms)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Get_Nuc_Charge_All_
