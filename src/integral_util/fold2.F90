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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************
      Subroutine Fold2(nSym,nBas,A,B)

      Implicit None

      Integer nSym,nBas(*)
      Real*8 A(*) , B(*)

      Integer iOff1, iOff2, iSym, mBas, iBas, jBas

      iOff1 = 0
      iOff2 = 0
      Do iSym = 1, nSym
        mBas = nBas(iSym)
        Do iBas= 1, mBas
          Do jBas = 1 , iBas-1
            B(iOff2+jBas) =   A(iOff1+jBas)
          End Do
          B(iOff2+iBas) =  A(iOff1+iBas)
          iOff1 = iOff1 + mBas
          iOff2 = iOff2 + iBas
        End Do
      End Do

      Return
      end Subroutine Fold2
