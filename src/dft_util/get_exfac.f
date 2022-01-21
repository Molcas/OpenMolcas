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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "hflda.fh"
      Real*8 Get_ExFac
      Character*(*) KSDFT
      Character*16  cTmp
      logical l_casdft
*                                                                      *
************************************************************************
*                                                                      *
      Get_ExFac=One
c      Get_ExFac=HFLDA
*
*     Write functional to run file.
*
      If (KSDFT.ne.'Overlap') Then
         cTmp=KSDFT
         Call Put_cArray('DFT functional',cTmp,16)
      End If
*                                                                      *
************************************************************************
* Global variable for MCPDFT                                           *

       l_casdft = KSDFT.eq.'TLSDA'   .or.
     &            KSDFT.eq.'TLSDA5'  .or.
     &            KSDFT.eq.'TBLYP'   .or.
     &            KSDFT.eq.'TSSBSW'  .or.
     &            KSDFT.eq.'TSSBD'   .or.
     &            KSDFT.eq.'TS12G'   .or.
     &            KSDFT.eq.'TPBE'    .or.
     &            KSDFT.eq.'FTPBE'   .or.
     &            KSDFT.eq.'TOPBE'   .or.
     &            KSDFT.eq.'FTOPBE'  .or.
     &            KSDFT.eq.'TREVPBE' .or.
     &            KSDFT.eq.'FTREVPBE'.or.
     &            KSDFT.eq.'FTLSDA'  .or.
     &            KSDFT.eq.'FTBLYP'
      If (l_casdft) Then
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
*     BR89G1                                                           *
*                                                                      *
      Case ('BR89G1')
         Get_ExFac=0.22D0
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP                                                            *
*                                                                      *
      Case ('B3LYP ')
         Get_ExFac=One-0.80D0
*                                                                      *
************************************************************************
*                                                                      *
*     O3LYP                                                            *
*                                                                      *
      Case ('O3LYP ')
         Get_ExFac = 0.1161D0
*                                                                      *
************************************************************************
*                                                                      *
*     B2PLYP                                                           *
*                                                                      *
      Case ('B2PLYP')
         Get_ExFac=0.530D0
*                                                                      *
************************************************************************
*                                                                      *
*     O2PLYP                                                           *
*                                                                      *
      Case ('O2PLYP')
         Get_ExFac=0.50D0
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP5                                                           *
*                                                                      *
      Case ('B3LYP5')
         Get_ExFac=One-0.80D0
*                                                                      *
************************************************************************
*                                                                      *
*     CASDFT, SCF, CS, M06-HF, TLYP                                    *
*                                                                      *
      Case ('CASDFT','SCF','CS','M06HF','TLYP')
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     S12H                                                             *
*                                                                      *
      Case ('S12H')
         Get_ExFac=0.25d0
*                                                                      *
************************************************************************
*                                                                      *
*     M06                                                              *
*                                                                      *
      Case ('M06 ')
         Get_ExFac=0.27D0
*                                                                      *
************************************************************************
*                                                                      *
*     M06-2X                                                           *
*                                                                      *
      Case ('M062X')
         Get_ExFac=0.54D0
*
************************************************************************
*                                                                      *
*     PBE0                                                             *
*                                                                      *
      Case ('PBE0')
         Get_ExFac=0.25D0
*                                                                      *
************************************************************************
*                                                                      *
      Case Default
         Get_ExFac=0.0D0
*                                                                      *
************************************************************************
*                                                                      *
      End Select
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
