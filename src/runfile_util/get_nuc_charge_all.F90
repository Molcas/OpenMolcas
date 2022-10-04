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

implicit real*8(a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
real*8 Charges_All(nAtoms_All)

call Get_nAtoms_All(nAtoms_Allx)
if (nAtoms_All /= nAtoms_Allx) then
  write(6,*) 'Get_Nuc_Charge_All: nAtoms_All /= nAtoms_Allx'
  write(6,*) 'nAtoms_All=',nAtoms_All
  write(6,*) 'nAtoms_Allx=',nAtoms_Allx
  call Abend()
end if

call Get_iScalar('Unique atoms',nAtoms)

call Allocate_Work(ipCU,3*nAtoms)
call Get_dArray('Unique Coordinates',Work(ipCU),3*nAtoms)

call Allocate_Work(ipCMu,nAtoms)
call Get_dArray('Nuclear charge',Work(ipCMu),nAtoms)

call Get_Nuc_Charge_All_(Work(ipCU),Work(ipCMu),nAtoms,Charges_All,nAtoms_All)

call Free_Work(ipCMu)
call Free_Work(ipCU)

return

end subroutine Get_Nuc_Charge_All
