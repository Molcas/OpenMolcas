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

subroutine Rd1Int_MCLR()
!***********************************************************************
!                                                                      *
!     Read header and matrices from the one-electron integral file     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use input_mclr, only: nSym, nBas, iMethod, nAtoms, ntBas, ntBTri, ntBSqr, AtLbl, ChIrr, Header1I, PotNuc, Coor

implicit none
#include "LenIn.fh"
character(len=8) Method
integer iSym

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!debug = .true.

call Get_cArray('Relax Method',Method,8)
if (Method == 'RHF-SCF') then
  iMethod = 1
else if ((Method == 'RASSCF') .or. (Method == 'CASSCF')) then
  iMethod = 2
else if ((Method == 'RASSCFSA') .or. (Method == 'CASSCFSA')) then
  iMethod = 2
else if (Method == 'CASPT2') then
  !iMethod = 3
  iMethod = 2 ! tentative
else if (Method == 'MBPT2') then
  iMethod = 4
else if (Method == 'MCPDFT') then
  iMethod = 2
else if (Method == 'MSPDFT') then
  iMethod = 2
end if
!---  read file header  -----------------------------------------------*
call Get_cArray('Seward Title',Header1I,144)
!---  read number of symm. species ------------------------------------*
call get_iScalar('nSym',nSym)
!---  read number of basis functions per symmetry species -------------*
call Get_iArray('nBas',nBas,nSym)
!---  read nuclear potential ------------------------------------------*
!call Get_PotNuc(PotNuc)
call Get_dScalar('PotNuc',Potnuc)
!---  read number of atoms --------------------------------------------*
call Get_iScalar('Unique atoms',nAtoms)
!---  read atom labels ------------------------------------------------*
call Get_cArray('Unique Atom Names',AtLbl,LenIn*nAtoms)
!---  read atom coordinates -------------------------------------------*
call Get_dArray('Unique Coordinates',Coor,3*nAtoms)
!---  read labels of the irreps ---------------------------------------*
call Get_cArray('Irreps',ChIrr,24)
!----------------------------------------------------------------------*
!     Precompute the total sum of variables and size of matrices       *
!----------------------------------------------------------------------*
ntBas = 0
ntBtri = 0
ntBsqr = 0
do iSym=1,nSym
  ntBas = ntBas+nBas(iSym)
  ntBtri = ntBtri+nTri_Elem(nBas(iSym))
  ntBsqr = ntBsqr+nBas(iSym)*nBas(iSym)
end do
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine Rd1Int_MCLR
