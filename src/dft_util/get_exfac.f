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
      use libxc_parameters
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Get_ExFac
      Character*(*) KSDFT
      Character*16  cTmp
      logical l_casdft
*                                                                      *
************************************************************************
*                                                                      *
      Get_ExFac=One
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
*     M06-HF                                                           *
*                                                                      *
      Case ('M06HF')
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     TLYP                                                             *
*                                                                      *
      Case ('TLYP')
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP                                                            *
*                                                                      *
      Case ('B3LYP ')
         Get_ExFac=LibXC_ExFac(XC_HYB_GGA_XC_B3LYP)
*                                                                      *
************************************************************************
*                                                                      *
*     O3LYP                                                            *
*                                                                      *
      Case ('O3LYP ')
         Get_ExFac=LibXC_ExFac(XC_HYB_GGA_XC_O3LYP)
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
         Get_ExFac=LibXC_ExFac(XC_HYB_GGA_XC_B3LYP5)
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
*
************************************************************************
*                                                                      *
*    BR3P86                                                            *
*                                                                      *
      Case ('BR3P86')
         Get_ExFac=0.22D0
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
      Contains
      Function LibXC_ExFac(funcid)
      Real*8 :: LibXC_ExFac
      Integer(kind=LibxcInt) :: funcid
      call xc_f03_func_init(xc_func(1),funcid,int(1, kind=LibxcInt))
      LibXC_ExFac=xc_f03_hyb_exx_coef(xc_func(1))
      call xc_f03_func_end(xc_func(1))
      Return
      End Function LibXC_ExFac
      End
