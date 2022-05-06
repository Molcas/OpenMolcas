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

do
  Key = Get_Ln(LuRd)
  !write(6,*) ' Processing:',Key
  KWord = Key
  call UpCase(KWord)
  select case (KWord(1:4))

    case ('RTHR')
      !                                                                *
      !***** RTHR ******************************************************
      !                                                                *
      ! Read the radial threshold

      KWord = Get_Ln(LuRd)
      call Get_F1(1,Threshold)
      Threshold = abs(Threshold)

    case ('GRID')
      !                                                                *
      !***** GRID ******************************************************
      !                                                                *
      ! Read quadrature quality

      KWord = Get_Ln(LuRd)
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

    case ('LMAX')
      !                                                                *
      !***** LMAX ******************************************************
      !                                                                *
      ! Read angular grid size

      KWord = Get_Ln(LuRd)
      call Get_I1(1,L_Quad)

    case ('RQUA')
      !                                                                *
      !***** RQUA ******************************************************
      !                                                                *
      ! Read radial quadrature scheme

      KWord = Get_Ln(LuRd)
      Quadrature = KWord(1:10)
      call Upcase(Quadrature)

    case ('NR  ')
      !                                                                *
      !***** NR   ******************************************************
      !                                                                *
      ! Read number of radial grid points

      KWord = Get_Ln(LuRd)
      call Get_I1(1,nR)

    case ('NGRI')
      !                                                                *
      !***** NGRI ******************************************************
      !                                                                *
      ! Read max number of grid points to process at one instance

      KWord = Get_Ln(LuRd)
      call Get_I1(1,nGridMax)

    case ('LOBA')
      !                                                                *
      !***** LOBA ******************************************************
      !                                                                *
      ! Activate use of Lobatto angular quadrature

      iOpt_Angular = ior(iand(iOpt_Angular,mask_111010),1)

    case ('GGL ')
      !                                                                *
      !***** GGL  ******************************************************
      !                                                                *
      ! Activate use of Gauss and Gauss-Legendre angular quadrature

      iOpt_Angular = iand(iOpt_Angular,mask_111010)

    case ('WHOL')
      !                                                                *
      !***** WHOL ******************************************************
      !                                                                *
      ! Activate use of routines which scan the whole atomic grid for
      ! each sub block.

      iOpt_Angular = ior(iand(iOpt_Angular,mask_111101),2)

    case ('GLOB')
      !                                                                *
      !***** GLOB ******************************************************
      !                                                                *
      ! Activate use of global partitioning technique.

      write(6,*) 'The Global option is redundant!'

    case ('DIAT')
      !                                                                *
      !***** DIAT ******************************************************
      !                                                                *
      ! Activate use of diatomic partitioning technique.

      write(6,*) 'The Diatomic option is redundant!'

    case ('NOPR')
      !                                                                *
      !***** NOPR ******************************************************
      !                                                                *
      ! Turn off the the angular prunning

      Angular_Prunning = Off

    case ('CROW')
      !                                                                *
      !***** CROW ******************************************************
      !                                                                *
      ! Read the crowding factor

      KWord = Get_Ln(LuRd)
      call Get_F1(1,Crowding)

    case ('LEBE')
      !                                                                *
      !***** LEBE ******************************************************
      !                                                                *
      ! Turn off the Lebedev angular grid

      iOpt_Angular = ior(iand(iOpt_Angular,mask_111011),4)

    case ('FIXE')
      !                                                                *
      !***** FIXE ******************************************************
      !                                                                *
      ! Turn on grid type = fixed

      Grid_Type = Fixed_Grid

    case ('MOVI')
      !                                                                *
      !***** MOVE ******************************************************
      !                                                                *
      ! Turn on grid type = moving

      Grid_Type = Moving_Grid

    case ('NORO')
      !                                                                *
      !***** NORO ******************************************************
      !                                                                *
      ! Turn of rotational invariant energy

      Rotational_Invariance = Off

    case ('RHOT')
      !                                                                *
      !***** RHOT ******************************************************
      !                                                                *
      ! Threshold for density when grid points are ignored.
      !
      ! Obsolete command!

      KWord = Get_Ln(LuRd)
      call Get_F1(1,Dummy)

    case ('NOSC')
      !                                                                *
      !***** NOSC ******************************************************
      !                                                                *
      ! Turn of the screening and the prunning.

      T_y = 0.0d0
      Crowding = 1.0d10
      Angular_Prunning = Off

    case ('T_Y ')
      !                                                                *
      !***** T_Y  ******************************************************
      !                                                                *
      ! Screening threshold for integral computation.

      KWord = Get_Ln(LuRd)
      call Get_F1(1,T_Y)

    case ('NQDI')
      !                                                                *
      !***** NQDI ******************************************************
      !                                                                *
      ! Recompute the AO values

      NQ_Direct = On

    case ('FADE')
      !                                                                *
      !***** T_Y  ******************************************************
      !                                                                *
      ! Fading factor for angular pruning.

      KWord = Get_Ln(LuRd)
      call Get_F1(1,Fade)

    case ('MOSS')
      !                                                                *
      !***** MOSS ******************************************************
      !                                                                *
      ! Assign Mossbauer center

      KWord = Get_Ln(LuRd)
      MBC = KWord(1:8)
      call UpCase(MBC)

    case ('END ')
      !                                                                *
      !***** END  ******************************************************
      !                                                                *
      exit

    case default
      iChrct = len(KWord)
      Last = iCLast(KWord,iChrct)
      write(6,*)
      call WarningMessage(2,'Error in FUNI_input')
      write(6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      write(6,*) ' Error in keyword.'
      call Quit_OnUserError()

  end select
end do
!                                                                      *
!***********************************************************************
!                                                                      *

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
