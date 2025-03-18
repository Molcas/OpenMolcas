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
      Subroutine SA_PREC(S,rdia)
      use ipPage, only: W
      use MCLR_Data, only: ipCI
      use input_mclr, only: nRoots,ERASSCF
      Implicit None
      Real*8 S(nroots**2,nroots),rdia(*)

      Integer irc,i
      Integer, External:: ipIN

      irc=ipin(ipci)
      Do i=1,nroots
         Call SA_PREC2(rdia,S(1,i),W(ipci)%Vec,ERASSCF(i))
      End Do
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(irc)
#endif
      End Subroutine SA_PREC
