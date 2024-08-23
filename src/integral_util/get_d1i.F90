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
      Subroutine Get_D1I(CMO,D1It,D1I,nish,nbas,nsym)
      use Constants, only: Zero, Two
      Implicit None
      Integer nSym
      Real*8 CMO(*), D1I(*),D1It(*)
      Integer nbas(nsym),nish(nsym)

      Integer iOff1, iSym, iBas, iOrb, iOff2, i, j, k
      Real*8 Sum

      iOff1 = 0
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nIsh(iSym)
        If ( iBas.ne.0 ) then
          iOff2 = iOff1
          Do i = 1,iBas
            Do j = 1,iBas
              Sum = Zero
              Do k = 0,iOrb-1
                Sum = Sum + Two * CMO(iOff1+k*iBas+i)                   &
     &                          * CMO(iOff1+k*iBas+j)
              End Do
              D1I(iOff2 + j) = Sum
            End Do
            iOff2 = iOff2 + iBas
          End Do
          iOff1 = iOff1 + iBas*iBas
        End If
      End Do
      Call Fold2(nsym,nBas,D1I,D1It)
      Return
      End Subroutine Get_D1I
