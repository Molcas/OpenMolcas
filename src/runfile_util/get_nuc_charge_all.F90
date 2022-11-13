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
!  Get_Nuc_Charge_All
!
!> @brief
!>   Get nuclear charges from RUNFILE
!> @author L. De Vico
!>
!> @details
!> Place nuclear charges (in a.u.) into array \p Charges_All(*).
!> Based on ::Get_Coord_All
!>
!> @param[out] Charges_All Array of charges
!> @param[in]  nAtoms_All  Number of atoms
!***********************************************************************

subroutine Get_Nuc_Charge_All(Charges_All,nAtoms_All)

use Symmetry_Info, only: iOper, nIrrep, Symmetry_Info_Get
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms_All
real(kind=wp), intent(out) :: Charges_All(nAtoms_All)
integer(kind=iwp) :: Active = 0, iAll_Atom, iChAtom, iCo, iCoSet(0:7,0:7), iGen(3), iStab(0:7), iUnique_Atom, MaxDCR, nAtoms, &
                     nAtoms_Allx, nCoSet, nGen, nStab
real(kind=wp) :: Charge_Old
real(kind=wp), allocatable :: CMu(:), CU(:,:)
integer(kind=iwp), external :: iChxyz

call Get_nAtoms_All(nAtoms_Allx)
if (nAtoms_All /= nAtoms_Allx) then
  write(u6,*) 'Get_Nuc_Charge_All: nAtoms_All /= nAtoms_Allx'
  write(u6,*) 'nAtoms_All=',nAtoms_All
  write(u6,*) 'nAtoms_Allx=',nAtoms_Allx
  call Abend()
end if

call Get_iScalar('Unique atoms',nAtoms)

call mma_allocate(CU,3,nAtoms,label='CU')
call Get_dArray('Unique Coordinates',CU,3*nAtoms)

call mma_allocate(CMu,nAtoms,label='CMu')
call Get_dArray('Nuclear charge',CMu,nAtoms)

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
do iUnique_Atom=1,nAtoms
  iChAtom = iChxyz(CU(:,iUnique_Atom),iGen,nGen)
  call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
  nCoSet = nIrrep/nStab

  Charge_Old = CMu(iUnique_Atom)

  do iCo=0,nCoSet-1
    iAll_Atom = iAll_Atom+1
    Charges_All(iAll_Atom) = Charge_Old
  end do

end do

!call RecPrt('CMu',' ',CMu,1,nAtoms)
!call RecPrt('Charges_All',' ',Charges_All,1,nAtoms_All)
!                                                                      *
!***********************************************************************
!                                                                      *

call mma_deallocate(CU)
call mma_deallocate(CMu)

return

end subroutine Get_Nuc_Charge_All
