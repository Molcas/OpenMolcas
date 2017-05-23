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
      Function arCos(Arg)
      Implicit Real*8 (a-h,o-z)
      Real*8 ArCos
      Character*72 Warning
#include "real.fh"
      A=Arg
      Delta=1.0D-12
      IF (ABS(A).GT.One) Then
         Write(Warning,3) A
3        FORMAT(1X,'Warning argument of aCos= ',1F21.18)
         If (ABS(A).lt.One+Delta) Then
            Call WarningMessage(1,Warning)
            A=Sign(One,A)
         Else
            Call WarningMessage(2,Warning)
            Call Abend()
         End If
      End If
*
      ArCos=ACos(A)
      Return
      End
