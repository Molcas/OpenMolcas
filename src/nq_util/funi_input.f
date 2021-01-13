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
      Subroutine Funi_Input(LuRd)
      Implicit Real*8 (a-h,o-z)
*
#include "nq_info.fh"
#include "real.fh"
      Character*180 Get_Ln,Key,KWord
      External Get_Ln
      Logical Check
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      Check(i,j)=iAnd(i,2**(j-1)).ne.0
*                                                                      *
************************************************************************
*                                                                      *
*
      mask_111110=62
      mask_111101=61
      mask_111011=59
      mask_111010=58
*
*     KeyWord directed input
*
 999  Continue
      Key = Get_Ln(LuRd)
*     Write (*,*) ' Processing:',Key
      KWord = Key
      Call UpCase(KWord)
      If (KWord(1:4).eq.'RTHR') Go To 100
      If (KWord(1:4).eq.'GRID') Go To 101
      If (KWord(1:4).eq.'LMAX') Go To 102
      If (KWord(1:4).eq.'RQUA') Go To 103
      If (KWord(1:4).eq.'NR  ') Go To 104
      If (KWord(1:4).eq.'NGRI') Go To 105
      If (KWord(1:4).eq.'LOBA') Go To 106
      If (KWord(1:4).eq.'GGL ') Go To 107
      If (KWord(1:4).eq.'WHOL') Go To 108
      If (KWord(1:4).eq.'GLOB') Go To 109
      If (KWord(1:4).eq.'DIAT') Go To 110
      If (KWord(1:4).eq.'NOPR') Go To 111
      If (KWord(1:4).eq.'CROW') Go To 112
      If (KWord(1:4).eq.'LEBE') Go To 113
      If (KWord(1:4).eq.'FIXE') Go To 114
      If (KWord(1:4).eq.'MOVI') Go To 115
      If (KWord(1:4).eq.'NORO') Go To 116
      If (KWord(1:4).eq.'RHOT') Go To 117
      If (KWord(1:4).eq.'T_X ') Go To 118
      If (KWord(1:4).eq.'NOSC') Go To 119
      If (KWord(1:4).eq.'T_Y ') Go To 120
      If (KWord(1:4).eq.'NQDI') Go To 121
      If (KWord(1:4).eq.'FADE') Go To 122
      If (KWord(1:4).eq.'MOSS') Go To 123
*
      If (KWord(1:4).eq.'END ') Go To 997
      iChrct=Len(KWord)
      Last=iCLast(KWord,iChrct)
      Write (6,*)
      Call WarningMessage(2,'Error in FUNI_input')
      Write (6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      Write (6,*) ' Error in keyword.'
      Call Quit_OnUserError()
*                                                                      *
****** RTHR ************************************************************
*                                                                      *
*     Read the radial threshold
*
 100  KWord = Get_Ln(LuRd)
      Call Get_F1(1,Threshold)
      Threshold = Abs(Threshold)
      Go To 999
*                                                                      *
****** GRID ************************************************************
*                                                                      *
*     Read quadrature quality
*
 101  KWord = Get_Ln(LuRd)
      Call UpCase(KWord)
      If (Index(KWord,'COARSE').ne.0) Then
*------- a la Gaussian
         nR=35
         L_Quad=17
         Crowding=0.90D0
         Fade=3.0D0
         Quadrature='MHL'
      Else If (Index(KWord,'ULTRAFINE').ne.0) Then
*------- a la Gaussian
         nR=99
         L_Quad=41
         Crowding=1.0D10
         Fade=10.0D0
         Quadrature='MHL'
      Else If (Index(KWord,'FINE').ne.0) Then
*------- a la Gaussian
         nR=75
         L_Quad=29
         Crowding=3.0D0
         Fade=6.0D0
         Quadrature='MHL'
      Else If (Index(KWord,'SG1GRID').ne.0) Then
*------- a la Gaussian
         nR=50
         L_Quad=23
         Crowding=1.0D0
         Fade=5.0D0
         Quadrature='MHL'
      Else
         Call WarningMessage(2,'Funi_Input: Illegal grid')
         Write (6,*) 'Type=',KWord
         Call Abend()
      End If
      Go To 999
*                                                                      *
****** LMAX ************************************************************
*                                                                      *
*     Read angular grid size
*
 102  KWord = Get_Ln(LuRd)
      Call Get_I1(1,L_Quad)
      Go To 999
*                                                                      *
****** RQUA ************************************************************
*                                                                      *
*     Read radial quadrature scheme
*
 103  KWord = Get_Ln(LuRd)
      Quadrature = KWord(1:10)
      Call Upcase(Quadrature)
      Go To 999
*                                                                      *
****** NR   ************************************************************
*                                                                      *
*     Read number of radial grid points
*
 104  KWord = Get_Ln(LuRd)
      Call Get_I1(1,nR)
      Go To 999
*                                                                      *
****** NGRI ************************************************************
*                                                                      *
*     Read max number of grid points to process at one instance
*
 105  KWord = Get_Ln(LuRd)
      Call Get_I1(1,nGridMax)
      Go To 999
*                                                                      *
****** LOBA ************************************************************
*                                                                      *
*     Activate use of Lobatto angular quadrature
*
 106  iOpt_Angular=iOr(iAnd(iOpt_Angular,mask_111010),1)
      Go To 999
*                                                                      *
****** NGRI ************************************************************
*                                                                      *
*     Activate use of Gauss and Gauss-Legendre angular quadrature
*
 107  iOpt_Angular=iAnd(iOpt_Angular,mask_111010)
      Go To 999
*                                                                      *
****** WHOL ************************************************************
*                                                                      *
*     Activate use of routines which scan the whole atomic grid for
*     each sub block.
*
 108  iOpt_Angular=iOr(iAnd(iOpt_Angular,mask_111101),2)
      Go To 999
*                                                                      *
****** GLOB ************************************************************
*                                                                      *
*     Activate use of global partitioning technique.
*
 109  Write (6,*) 'The Global option is redundant!'
      Go To 999
*                                                                      *
****** DIAT ************************************************************
*                                                                      *
*     Activate use of diatomic partitioning technique.
*
 110  Write (6,*) 'The Diatomic option is redundant!'
      Go To 999
*                                                                      *
****** NOPR ************************************************************
*                                                                      *
*     Turn off the the angular prunning
*
 111  Angular_Prunning = Off
      Go To 999
*                                                                      *
****** CROW ************************************************************
*                                                                      *
*     Read the crowding factor
*
 112  KWord = Get_Ln(LuRd)
      Call Get_F1(1,Crowding)
      Go To 999
*                                                                      *
****** LEBE ************************************************************
*                                                                      *
*     Turn off the Lebedev angular grid
*
 113  iOpt_Angular=iOr(iAnd(iOpt_Angular,mask_111011),4)
      Go To 999
*                                                                      *
****** FIXE ************************************************************
*                                                                      *
*     Turn on grid type = fixed
*
 114  Grid_Type=Fixed_Grid
      Go To 999
*                                                                      *
****** MOVE ************************************************************
*                                                                      *
*     Turn on grid type = moving
*
 115  Grid_Type=Moving_Grid
      Go To 999
*                                                                      *
****** NORO ************************************************************
*                                                                      *
*     Turn of rotational invariant energy
*
 116  Rotational_Invariance = Off
      Go To 999
*                                                                      *
****** RHOT ************************************************************
*                                                                      *
*     Threshold for density when grid points are ignored.
*
*     Obsolete command!
*
 117  KWord = Get_Ln(LuRd)
      Call Get_F1(1,Dummy)
      Go To 999
*                                                                      *
****** T_X  ************************************************************
*                                                                      *
*     Screening threshold for density computation.
*
 118  KWord = Get_Ln(LuRd)
      Call Get_F1(1,T_X)
      Go To 999
*                                                                      *
****** NOSC ************************************************************
*                                                                      *
*     Turn of the screening and the prunning.
*
 119  T_X=0.0D0
      T_y=0.0D0
      Crowding=1.0D10
      Angular_Prunning = Off
      Go To 999
*                                                                      *
****** T_Y  ************************************************************
*                                                                      *
*     Screening threshold for integral computation.
*
 120  KWord = Get_Ln(LuRd)
      Call Get_F1(1,T_Y)
      Go To 999
*                                                                      *
****** NQDI ************************************************************
*                                                                      *
*     Recompute the AO values
*
 121  NQ_Direct=On
      Go To 999
*                                                                      *
****** T_Y  ************************************************************
*                                                                      *
*     Fading factor for angular pruning.
*
 122  KWord = Get_Ln(LuRd)
      Call Get_F1(1,Fade)
      Go To 999
*                                                                      *
****** MOSS ************************************************************
*                                                                      *
*     Assign Mossbauer center
*
 123  KWord = Get_Ln(LuRd)
      MBC=KWord(1:8)
      Call UpCase(MBC)
      Go To 999
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
 997  Continue
*
      If (Check(iOpt_Angular,3)) Then
         If (L_Quad.ne. 5 .and.
     &       L_Quad.ne. 7 .and.
     &       L_Quad.ne.11 .and.
     &       L_Quad.ne.17 .and.
     &       L_Quad.ne.23 .and.
     &       L_Quad.ne.29 .and.
     &       L_Quad.ne.35 .and.
     &       L_Quad.ne.41 .and.
     &       L_Quad.ne.47 .and.
     &       L_Quad.ne.53 .and.
     &       L_Quad.ne.59       ) Then
            Write (6,*) 'L_Quad does not comply with Lebedev grid.'
            iOpt_Angular=iAnd(iOpt_Angular,mask_111011)
            Write (6,*) 'Lobatto grid activated!'
            iOpt_Angular=iOr(iAnd(iOpt_Angular,mask_111110),1)
         End If
      End If
*
      Return
      End
