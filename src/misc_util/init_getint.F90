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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************
!*************************************************************
! Initializes info needed by the Cholesky integral generator
!
! F. Aquilante
!*************************************************************

subroutine INIT_GETINT(RC)

use GetInt_mod, only: LuCVec, nBas, NumCho, pq1, nRS, nPQ
use Definitions, only: iwp
use Index_Functions, only: nTri_Elem
use RICD_Info, only: Do_DCCD

implicit none
integer(kind=iwp), intent(out) :: RC
integer(kind=iwp) :: nSym

rc = 0
call get_iscalar('nSym',nSym)
call get_iarray('nBas',nBas,nSym)
call INIT_NumCV(NumCho,nSym)

If (Do_DCCD) Then
   If (NumCho(1) < 1) Then
      Write (6,*) 'Init_GetInt: NumCho(1) < 1'
      Call Abend()
   End If

   nPQ=nTri_Elem(nBas(1))
   nRS=nPQ
End If

LuCVec(1) = -1
LuCVec(2) = -1

pq1 = 0

return

end subroutine INIT_GETINT
