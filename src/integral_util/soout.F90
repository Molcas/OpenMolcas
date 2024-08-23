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
      Subroutine SOOUT(label,cnt_ico,phase_ico)
      Implicit None
#include "Molcas.fh"
      Integer cnt_ico(0:7,*),phase_ico(0:7,*)
      Character(LEN=LENIN8) Label(MaxBfn+MaxBfn_Aux)
!
      Call SOCtl_mod(Label,Maxbfn+MaxBfn_Aux,Cnt_ico,Phase_ico)
!
      End Subroutine SOOUT
