************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      Logical Function RinT_(iOpT,nOpT,iOpR)
************************************************************************
*                                                                      *
* Object: to return .true. if R is in {T}.                             *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             June '90                                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Integer iOpT(nOpT)
      RinT_ = .False.
      Do 10 i = 1, nOpT
         If (iOpT(i).eq.iOpR) Then
            RinT_ = .True.
            Return
         End If
 10   Continue
      Return
      End
