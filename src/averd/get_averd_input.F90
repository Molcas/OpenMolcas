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
!---------------------------------------------------------------------
! Read input for averd.
!
! Title         -Title
! Wset          -Set of weights. Dimension is nSet.
! iPrint        -How much print.(1=minimal,2=print average orbitals,
!                5=print all orbitals,99=wacko!)
! nSet          -Number of orbitals to read
! DensityBased  -Is the procedure density or orbital based
! ThrOcc        -Print which orbitals have occupation number below
!                this threshold.
!-----------------------------------------------------------------------

subroutine Get_Averd_input(Title,iPrint,nSet,DensityBased,ThrOcc)

use Averd_global, only: Wset
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp, u6

implicit none
character(len=72), intent(inout) :: Title
real(kind=wp), intent(inout) :: ThrOcc
integer(kind=iwp), intent(inout) :: iPrint, nSet
logical(kind=iwp), intent(inout) :: DensityBased
integer(kind=iwp) :: iChrct, Last, LuRd
character(len=180) :: Key
character(len=4) :: Kword
character(len=180), external :: Get_Ln
integer(kind=iwp), external :: iCLast
#include "warnings.h"

!-- Call subroutines that handle the input.

LuRd = 21
call SpoolInp(LuRd)
rewind(LuRd)
call RdNLst(LuRd,'AVERD')

do

  !-- Get_Ln read the keyword and skips line starting with *
  !   or is empty.

  Key = Get_Ln(LuRd)
  Kword = trim(Key)
  call UpCase(Kword)

  !-- The keywords...

  select case (Kword(1:4))

    case ('WSET')

      !-- Read weights.

      Key = Get_Ln(LuRd)
      call Get_I1(1,nSet)
      call mma_allocate(Wset,nSet,label='Wset')
      Key = Get_Ln(LuRd)
      call Get_F(1,Wset,nSet)

    case ('PRIN')

      !-- How much print?

      Key = Get_Ln(LuRd)
      call Get_I1(1,iPrint)

    case ('TITL')

      !-- Title

      Key = Get_Ln(LuRd)
      Title = Key(1:len(Title))

    case ('ORBI')

      !-- Should it be density based, or orbital based.

      DensityBased = .false.

    case ('OCCU')

      !-- I want to be told which orbitals have occ.num. below threshold.

      Key = Get_Ln(LuRd)
      call Get_F1(1,ThrOcc)

    case ('END ')

      exit

    case default

      !-- ...and what happens if something else is encountered.

      iChrct = len(KWord)
      Last = iCLast(KWord,iChrct)
      write(u6,*) ' '
      write(u6,'(1x,a,a)') Kword(1:Last),' is not a valid keyword!'
      write(u6,*) ' ERROR!'
      call Quit(_RC_INPUT_ERROR_)

  end select

end do

!-- A most Graceful Exit.

return

end subroutine Get_Averd_input
