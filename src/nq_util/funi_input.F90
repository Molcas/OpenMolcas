!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Funi_Input(LuRd)
use nq_Grid, only: nGridMax
use nq_Info
implicit real*8(a-h,o-z)
#include "real.fh"
character*180 Get_Ln, Key, KWord
external Get_Ln
logical Check
! Statement function
Check(i,j) = iand(i,2**(j-1)) /= 0

!                                                                      *
!***********************************************************************
!                                                                      *
mask_111110 = 62
mask_111101 = 61
mask_111011 = 59
mask_111010 = 58

! KeyWord directed input

999 continue
Key = Get_Ln(LuRd)
!write(6,*) ' Processing:',Key
KWord = Key
call UpCase(KWord)
if (KWord(1:4) == 'RTHR') Go To 100
if (KWord(1:4) == 'GRID') Go To 101
if (KWord(1:4) == 'LMAX') Go To 102
if (KWord(1:4) == 'RQUA') Go To 103
if (KWord(1:4) == 'NR  ') Go To 104
if (KWord(1:4) == 'NGRI') Go To 105
if (KWord(1:4) == 'LOBA') Go To 106
if (KWord(1:4) == 'GGL ') Go To 107
if (KWord(1:4) == 'WHOL') Go To 108
if (KWord(1:4) == 'GLOB') Go To 109
if (KWord(1:4) == 'DIAT') Go To 110
if (KWord(1:4) == 'NOPR') Go To 111
if (KWord(1:4) == 'CROW') Go To 112
if (KWord(1:4) == 'LEBE') Go To 113
if (KWord(1:4) == 'FIXE') Go To 114
if (KWord(1:4) == 'MOVI') Go To 115
if (KWord(1:4) == 'NORO') Go To 116
if (KWord(1:4) == 'RHOT') Go To 117
if (KWord(1:4) == 'NOSC') Go To 119
if (KWord(1:4) == 'T_Y ') Go To 120
if (KWord(1:4) == 'NQDI') Go To 121
if (KWord(1:4) == 'FADE') Go To 122
if (KWord(1:4) == 'MOSS') Go To 123

if (KWord(1:4) == 'END ') Go To 997
iChrct = len(KWord)
Last = iCLast(KWord,iChrct)
write(6,*)
call WarningMessage(2,'Error in FUNI_input')
write(6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
write(6,*) ' Error in keyword.'
call Quit_OnUserError()
!                                                                      *
!***** RTHR ************************************************************
!                                                                      *
! Read the radial threshold

100 KWord = Get_Ln(LuRd)
call Get_F1(1,Threshold)
Threshold = abs(Threshold)
Go To 999
!                                                                      *
!***** GRID ************************************************************
!                                                                      *
! Read quadrature quality

101 KWord = Get_Ln(LuRd)
call UpCase(KWord)
if (index(KWord,'COARSE') /= 0) then
  ! a la Gaussian
  nR = 35
  L_Quad = 17
  Crowding = 0.90d0
  Fade = 3.0d0
  Quadrature = 'MHL'
else if (index(KWord,'ULTRAFINE') /= 0) then
  ! a la Gaussian
  nR = 99
  L_Quad = 41
  Crowding = 1.0d10
  Fade = 10.0d0
  Quadrature = 'MHL'
else if (index(KWord,'FINE') /= 0) then
  ! a la Gaussian
  nR = 75
  L_Quad = 29
  Crowding = 3.0d0
  Fade = 6.0d0
  Quadrature = 'MHL'
else if (index(KWord,'SG1GRID') /= 0) then
  ! a la Gaussian
  nR = 50
  L_Quad = 23
  Crowding = 1.0d0
  Fade = 5.0d0
  Quadrature = 'MHL'
else
  call WarningMessage(2,'Funi_Input: Illegal grid')
  write(6,*) 'Type=',KWord
  call Abend()
end if
Go To 999
!                                                                      *
!***** LMAX ************************************************************
!                                                                      *
! Read angular grid size

102 KWord = Get_Ln(LuRd)
call Get_I1(1,L_Quad)
Go To 999
!                                                                      *
!***** RQUA ************************************************************
!                                                                      *
! Read radial quadrature scheme

103 KWord = Get_Ln(LuRd)
Quadrature = KWord(1:10)
call Upcase(Quadrature)
Go To 999
!                                                                      *
!***** NR   ************************************************************
!                                                                      *
! Read number of radial grid points

104 KWord = Get_Ln(LuRd)
call Get_I1(1,nR)
Go To 999
!                                                                      *
!***** NGRI ************************************************************
!                                                                      *
! Read max number of grid points to process at one instance

105 KWord = Get_Ln(LuRd)
call Get_I1(1,nGridMax)
Go To 999
!                                                                      *
!***** LOBA ************************************************************
!                                                                      *
! Activate use of Lobatto angular quadrature

106 iOpt_Angular = ior(iand(iOpt_Angular,mask_111010),1)
Go To 999
!                                                                      *
!***** NGRI ************************************************************
!                                                                      *
! Activate use of Gauss and Gauss-Legendre angular quadrature

107 iOpt_Angular = iand(iOpt_Angular,mask_111010)
Go To 999
!                                                                      *
!***** WHOL ************************************************************
!                                                                      *
! Activate use of routines which scan the whole atomic grid for
! each sub block.

108 iOpt_Angular = ior(iand(iOpt_Angular,mask_111101),2)
Go To 999
!                                                                      *
!***** GLOB ************************************************************
!                                                                      *
! Activate use of global partitioning technique.

109 write(6,*) 'The Global option is redundant!'
Go To 999
!                                                                      *
!***** DIAT ************************************************************
!                                                                      *
! Activate use of diatomic partitioning technique.

110 write(6,*) 'The Diatomic option is redundant!'
Go To 999
!                                                                      *
!***** NOPR ************************************************************
!                                                                      *
! Turn off the the angular prunning

111 Angular_Prunning = Off
Go To 999
!                                                                      *
!***** CROW ************************************************************
!                                                                      *
! Read the crowding factor

112 KWord = Get_Ln(LuRd)
call Get_F1(1,Crowding)
Go To 999
!                                                                      *
!***** LEBE ************************************************************
!                                                                      *
! Turn off the Lebedev angular grid

113 iOpt_Angular = ior(iand(iOpt_Angular,mask_111011),4)
Go To 999
!                                                                      *
!***** FIXE ************************************************************
!                                                                      *
! Turn on grid type = fixed

114 Grid_Type = Fixed_Grid
Go To 999
!                                                                      *
!***** MOVE ************************************************************
!                                                                      *
! Turn on grid type = moving

115 Grid_Type = Moving_Grid
Go To 999
!                                                                      *
!***** NORO ************************************************************
!                                                                      *
! Turn of rotational invariant energy

116 Rotational_Invariance = Off
Go To 999
!                                                                      *
!***** RHOT ************************************************************
!                                                                      *
! Threshold for density when grid points are ignored.
!
! Obsolete command!

117 KWord = Get_Ln(LuRd)
call Get_F1(1,Dummy)
Go To 999
!                                                                      *
!***** NOSC ************************************************************
!                                                                      *
! Turn of the screening and the prunning.

119 T_y = 0.0d0
Crowding = 1.0d10
Angular_Prunning = Off
Go To 999
!                                                                      *
!***** T_Y  ************************************************************
!                                                                      *
! Screening threshold for integral computation.

120 KWord = Get_Ln(LuRd)
call Get_F1(1,T_Y)
Go To 999
!                                                                      *
!***** NQDI ************************************************************
!                                                                      *
! Recompute the AO values

121 NQ_Direct = On
Go To 999
!                                                                      *
!***** T_Y  ************************************************************
!                                                                      *
! Fading factor for angular pruning.

122 KWord = Get_Ln(LuRd)
call Get_F1(1,Fade)
Go To 999
!                                                                      *
!***** MOSS ************************************************************
!                                                                      *
! Assign Mossbauer center

123 KWord = Get_Ln(LuRd)
MBC = KWord(1:8)
call UpCase(MBC)
Go To 999
!                                                                      *
!***********************************************************************
!                                                                      *
997 continue

if (Check(iOpt_Angular,3)) then
  if ((L_Quad /= 5) .and. (L_Quad /= 7) .and. (L_Quad /= 11) .and. (L_Quad /= 17) .and. (L_Quad /= 23) .and. (L_Quad /= 29) .and. &
      (L_Quad /= 35) .and. (L_Quad /= 41) .and. (L_Quad /= 47) .and. (L_Quad /= 53) .and. (L_Quad /= 59)) then
    write(6,*) 'L_Quad does not comply with Lebedev grid.'
    iOpt_Angular = iand(iOpt_Angular,mask_111011)
    write(6,*) 'Lobatto grid activated!'
    iOpt_Angular = ior(iand(iOpt_Angular,mask_111110),1)
  end if
end if

return

end subroutine Funi_Input
