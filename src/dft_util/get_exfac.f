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
*     Return the factor which determine how much "exact exchange" that *
*     should be included.                                              *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "hflda.fh"
      Real*8 Get_ExFac
      Character*(*) KSDFT
      Character*16  cTmp
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
*                                                                      *
*      LSDA LDA SVWN                                                   *
*                                                                      *
       If (KSDFT.eq.'LSDA ' .or.
     &     KSDFT.eq.'LDA '  .or.
     &     KSDFT.eq.'SVWN ') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*      MC-PDFT                                                         *
*                                                                      *
       Else If (KSDFT(1:5).eq.'TLSDA'.or. !GLM
     &          KSDFT(1:5).eq.'TBLYP'.or.
     &          KSDFT(1:4).eq.'TSSB'.or.
     &          KSDFT(1:6).eq.'FTLSDA'.or.
     &          KSDFT(1:6).eq.'FTBLYP'.or.
     &          KSDFT(1:5).eq.'FTPBE'.or.
     &          KSDFT(1:8).eq.'FTREVPBE'.or.
     &          KSDFT(1:7).eq.'TREVPBE'.or.
     &          KSDFT(1:4).eq.'TPBE') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*      TST                                                             *
*                                                                      *
       Else If (KSDFT.eq.'TST' ) Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*      LSDA5 LDA5 SVWN5                                                *
*                                                                      *
       Else If (KSDFT.eq.'LSDA5' .or.
     &          KSDFT.eq.'LDA5'  .or.
     &          KSDFT.eq.'SVWN5') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     HFB                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFB') Then
         Get_ExFac=Zero
************************************************************************
*                                                                      *
*     HFO                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFO') Then
         Get_ExFac=Zero
************************************************************************
*                                                                      *
*     HFG                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFG') Then
         Get_ExFac=Zero
************************************************************************
*                                                                      *
*     HFB86                                                            *
*                                                                      *
       Else If (KSDFT.eq.'HFB86') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*      HFS                                                             *
*                                                                      *
       Else If (KSDFT.eq.'HFS') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*      XALPHA                                                          *
*                                                                      *
       Else If (KSDFT.eq.'XALPHA') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     Overlap                                                          *
*                                                                      *
      Else If (KSDFT.eq.'Overlap') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     BWIG                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BWIG') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     BLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BLYP') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     OLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'OLYP') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     KT3                                                              *
*                                                                      *
      Else If (KSDFT.eq.'KT3') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     KT2                                                              *
*                                                                      *
      Else If (KSDFT.eq.'KT2') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     BPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BPBE') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     GLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'GLYP') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     GPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'GPBE') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     B86LYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B86LYP') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     B86PBE                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B86PBE') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     OPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'OPBE') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     TLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'TLYP') Then
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     NLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'NLYP') Then
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     NEWF0                                                            *
*                                                                      *
      Else If (KSDFT.eq.'NEWF0') Then
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     NEWF1                                                            *
*                                                                      *
      Else If (KSDFT.eq.'NEWF1') Then
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP                                                            *
*                                                                      *
      Else If (KSDFT.eq.'B3LYP ') Then
         Get_ExFac=One-0.80D0
*                                                                      *
************************************************************************
*                                                                      *
*     O3LYP                                                            *
*                                                                      *
      Else If (KSDFT.eq.'O3LYP ') Then
         Get_ExFac = 0.1161D0
*                                                                      *
************************************************************************
*                                                                      *
*     B2PLYP                                                           *
*                                                                      *
      Else If (KSDFT(1:6).eq.'B2PLYP') Then
         Get_ExFac=0.530D0
*                                                                      *
************************************************************************
*                                                                      *
*     O2PLYP                                                           *
*                                                                      *
      Else If (KSDFT(1:6).eq.'O2PLYP') Then
         Get_ExFac=0.50D0
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP5                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B3LYP5') Then
         Get_ExFac=One-0.80D0
*                                                                      *
************************************************************************
*                                                                      *
*     CASDFT                                                           *
*                                                                      *
      Else If (KSDFT.eq.'CASDFT') Then
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     SCF                                                              *
*                                                                      *
      Else If (KSDFT.eq.'SCF') Then
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     PAM                                                              *
*                                                                      *
      Else If (KSDFT(1:3).eq.'PAM') Then
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     CS                                                               *
*                                                                      *
      Else If (KSDFT(1:2).eq.'CS') Then
         Get_ExFac=One
*                                                                      *
************************************************************************
*                                                                      *
*     PBE                                                              *
*                                                                      *
      Else If (KSDFT(1:4).eq.'PBE') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     REVPBE                                                           *
*                                                                      *
      Else If (KSDFT(1:6).eq.'REVPBE') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     SSB                                                              *
*                                                                      *
      Else If (KSDFT(1:4).eq.'SSB') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     PBEsol                                                           *
*                                                                      *
      Else If (KSDFT(1:6).eq.'PBESOL') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     RGE2                                                             *
*                                                                      *
      Else If (KSDFT(1:6).eq.'RGE2') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     TCA                                                              *
*                                                                      *
      Else If (KSDFT(1:6).eq.'PTCA') Then
         Get_ExFac=Zero

*                                                                      *
************************************************************************
*                                                                      *
*     M06-L                                                            *
*                                                                      *
      Else If (KSDFT(1:4).eq.'M06L') Then
         Get_ExFac=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     M06                                                              *
*                                                                      *
      Else If (KSDFT(1:4).eq.'M06 ') Then
         Get_ExFac=0.27D0
*                                                                      *
************************************************************************
*                                                                      *
*     M06-2X                                                           *
*                                                                      *
      Else If (KSDFT(1:5).eq.'M062X') Then
         Get_ExFac=0.54D0
*                                                                      *
************************************************************************
*                                                                      *
*     M06-HF                                                           *
*                                                                      *
      Else If (KSDFT(1:5).eq.'M06HF ') Then
         Get_ExFac=1.0D0
*
************************************************************************
*                                                                      *
*     PBE0                                                             *
*                                                                      *
      Else If (KSDFT(1:4).eq.'PBE0') Then
         Get_ExFac=0.25D0
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Call WarningMessage(2,
     &               'Unknown DFT functional;')
         Write (6,*) 'KSDFT=',KSDFT
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
