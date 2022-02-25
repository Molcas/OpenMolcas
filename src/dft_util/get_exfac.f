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
      Function Get_ExFac(KSDFT)
************************************************************************
*     Return the factor which determines how much "exact exchange" that*
*     should be included.                                              *
************************************************************************
      Use Functionals, Only: Get_Func_ExFac
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Get_ExFac
      Character*(*) KSDFT
      Character(LEN=80)  cTmp
*                                                                      *
************************************************************************
*                                                                      *
*     Write functional to run file.
*
      If (KSDFT.ne.'Overlap') Then
         cTmp=KSDFT
         Call Put_cArray('DFT functional',cTmp,LEN(cTmp))
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (KSDFT(1:2).eq.'T:' .or. KSDFT(1:3).eq.'FT:') Then
         Get_ExFac=Zero
         Return
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     We bring in only cases where it is different from zero.
      Select Case(KSDFT)
*                                                                      *
************************************************************************
*                                                                      *
*     CASDFT                                                           *
*                                                                      *
      Case ('CASDFT')
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     SCF                                                              *
*                                                                      *
      Case ('SCF')
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     CS                                                               *
*                                                                      *
      Case ('CS')
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
      Case Default
         Get_ExFac = Get_Func_ExFac(KSDFT)
*                                                                      *
************************************************************************
*                                                                      *
      End Select
*                                                                      *
************************************************************************
*                                                                      *
      End
