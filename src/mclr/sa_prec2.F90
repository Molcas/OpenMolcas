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
      Subroutine SA_PREC2(rdia,S,CI,ENE)
      use Constants, only: Zero
      use input_mclr, only: nRoots,nCSF,State_Sym
      Implicit None
      Real*8 rdia(*),S(nroots,nroots),CI(*)
      Real*8 ENE

      Integer i,j,k
      Real*8 dnum

      Do i=0,nroots-1
       Do j=0,nroots-1
         S(i+1,j+1)=Zero
         Do k=1,ncsf(State_Sym)
          dnum=rdia(k)-Ene
          dnum=Sign(Max(Abs(dnum),1.0d-16),dnum)
          S(i+1,j+1)=S(i+1,j+1)+                                        &
     &    CI(i*ncsf(State_Sym)+k)*CI(j*ncsf(State_Sym)+k)/dnum
         End Do
       End Do
      End Do
      Call MatInvert(S,nroots)
      end Subroutine SA_PREC2
