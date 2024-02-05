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
Module RASSCF_LUCIA
Private
      INTEGER, Public:: C_POINTER, Lucia_base, Lucia_length,            &
     &                  kvec3_length, iSigma_on_disk,ini_h0,            &
     &                  Memory_Needed_Lucia
      INTEGER, Public:: LW6,LW7,LW8,LW9,LW_RF1,LW_RF2
End Module RASSCF_LUCIA
