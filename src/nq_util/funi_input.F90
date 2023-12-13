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
use nq_Info, only: Angular_Pruning, Crowding, Fade, Fixed_Grid, Grid_Type, iOpt_Angular, L_Quad, MBC, Moving_Grid, NQ_Direct, nR, &
                   Off, On, Quadrature, Rotational_Invariance, T_Y, Threshold
use Constants, only: Zero, One, Three, Five, Six, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuRd
integer(kind=iwp) :: iChrct, Last
real(kind=wp) :: Dummy
character(len=180) :: Key, KWord
integer(kind=iwp), external :: iCLast
character(len=180), external :: Get_Ln

!                                                                      *
!***********************************************************************
!                                                                      *
! KeyWord directed input

do
  Key = Get_Ln(LuRd)
  !write(u6,*) ' Processing:',Key
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
        Crowding = 0.9_wp
        Fade = Three
        Quadrature = 'MHL'
      else if (index(KWord,'ULTRAFINE') /= 0) then
        ! a la Gaussian
        nR = 99
        L_Quad = 41
        Crowding = 1.0e10_wp
        Fade = Ten
        Quadrature = 'MHL'
      else if (index(KWord,'FINE') /= 0) then
        ! a la Gaussian
        nR = 75
        L_Quad = 29
        Crowding = Three
        Fade = Six
        Quadrature = 'MHL'
      else if (index(KWord,'SG1GRID') /= 0) then
        ! a la Gaussian
        nR = 50
        L_Quad = 23
        Crowding = One
        Fade = Five
        Quadrature = 'MHL'
      else
        call WarningMessage(2,'Funi_Input: Illegal grid')
        write(u6,*) 'Type=',KWord
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

      iOpt_Angular = ibset(ibclr(iOpt_Angular,2),0)

    case ('GGL ')
      !                                                                *
      !***** GGL  ******************************************************
      !                                                                *
      ! Activate use of Gauss and Gauss-Legendre angular quadrature

      iOpt_Angular = ibclr(ibclr(iOpt_Angular,2),0)

    case ('WHOL')
      !                                                                *
      !***** WHOL ******************************************************
      !                                                                *
      ! Activate use of routines which scan the whole atomic grid for
      ! each sub block.

      iOpt_Angular = ibset(iOpt_Angular,1)

    case ('GLOB')
      !                                                                *
      !***** GLOB ******************************************************
      !                                                                *
      ! Activate use of global partitioning technique.

      write(u6,*) 'The Global option is redundant!'

    case ('DIAT')
      !                                                                *
      !***** DIAT ******************************************************
      !                                                                *
      ! Activate use of diatomic partitioning technique.

      write(u6,*) 'The Diatomic option is redundant!'

    case ('NOPR')
      !                                                                *
      !***** NOPR ******************************************************
      !                                                                *
      ! Turn off the the angular pruning

      Angular_Pruning = Off

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

      iOpt_Angular = ibset(iOpt_Angular,2)

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
      ! Turn off the screening and the pruning.

      T_y = Zero
      Crowding = 1.0e10_wp
      Angular_Pruning = Off

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
      write(u6,*)
      call WarningMessage(2,'Error in FUNI_input')
      write(u6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      write(u6,*) ' Error in keyword.'
      call Quit_OnUserError()

  end select
end do
!                                                                      *
!***********************************************************************
!                                                                      *

if (btest(iOpt_Angular,2)) then
  if ((L_Quad /= 5) .and. (L_Quad /= 7) .and. (L_Quad /= 11) .and. (L_Quad /= 17) .and. (L_Quad /= 23) .and. (L_Quad /= 29) .and. &
      (L_Quad /= 35) .and. (L_Quad /= 41) .and. (L_Quad /= 47) .and. (L_Quad /= 53) .and. (L_Quad /= 59)) then
    write(u6,*) 'L_Quad does not comply with Lebedev grid.'
    iOpt_Angular = ibclr(iOpt_Angular,2)
    write(u6,*) 'Lobatto grid activated!'
    iOpt_Angular = ibset(iOpt_Angular,0)
  end if
end if

return

end subroutine Funi_Input
