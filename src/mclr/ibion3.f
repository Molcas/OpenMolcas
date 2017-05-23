************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Integer Function iBion(M,N)
C
C     Compute binomial coefficients recursively
C
      Real*8 Temp1,Temp2,Temp3,Temp4
C
      If ( N.lt.0 .or. M.lt.N ) then
         write(6,*) 'Wrong params is iBion',M,N
         Call Abend
      End If
      If ( N.eq.0 .or. M.eq.0 ) then
         iBion=1
         Return
      Else If ( N.eq.1 ) then
         iBion=M
         Return
      Else If ( N.eq.2 ) then
         iBion=M*(M-1)/2
         Return
      End If
      K1=Max((M-N),N)
      K2=Min((M-N),N)
      Temp1=1.0D0
      Do K=1,K2
         Temp2=DBLE(K1+K)
         Temp3=DBLE(K)
         Temp4=Temp2/Temp3
         Temp1=Temp1*Temp4
      End Do
      iBion=nInt(Temp1)
*
      Return
      End
