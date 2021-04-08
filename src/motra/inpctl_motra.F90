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
! Copyright (C) 1991, Markus P. Fuelscher                              *
!***********************************************************************

subroutine InpCtl_Motra()
!***********************************************************************
!                                                                      *
! Purpose:                                                             *
! Read all information required                                        *
!                                                                      *
!**** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************

use motra_global, only: CMO, iAutoCut, iPrint, iRFpert, nTot2, Ovlp
use stdalloc, only: mma_allocate

implicit none

!----------------------------------------------------------------------*
! Read the content of the one electron integral file                   *
!----------------------------------------------------------------------*
call Rd1Int_Motra()
!----------------------------------------------------------------------*
! Read auxiliary  input                                                *
!----------------------------------------------------------------------*
call RdInp_Motra()
!----------------------------------------------------------------------*
! Read Reaction field and add to one-electron integrals                *
!----------------------------------------------------------------------*
if (iRFpert == 1) call RdRfld()
!----------------------------------------------------------------------*
! Read the MO coefficients and occupations                             *
!----------------------------------------------------------------------*
call mma_allocate(CMO,nTot2,label='CMO')
call RdCmo_motra(CMO,Ovlp)
!----------------------------------------------------------------------*
! Delete orbitals with occupations samller than a given value          *
!----------------------------------------------------------------------*
if (iAutoCut == 1) call AutoCut()
!----------------------------------------------------------------------*
! Print the input and orbital definitions                              *
!----------------------------------------------------------------------*
if (iPrint >= 0) call PrInp(CMO)
!----------------------------------------------------------------------*
! Normal termination                                                   *
!----------------------------------------------------------------------*

return

end subroutine InpCtl_Motra
