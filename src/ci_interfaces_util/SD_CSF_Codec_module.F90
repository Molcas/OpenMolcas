!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

module SD_CSF_Codec_Module

  use constants, only: Zero, One
  use sguga, only: SGS
  use Definitions, only: wp, iwp, u6

  implicit none
  private

  public :: SD_CSF_Codec
  public :: CSF_TO_SD
  public :: SD_TO_CSF

  integer(iwp), parameter :: CSF_TO_SD = 1
  integer(iwp), parameter :: SD_TO_CSF = 2

contains

  subroutine SD_CSF_Codec(SD,CSF,Direction,IRefSm)

  real(wp), intent(inout) :: SD(:)
  real(wp), intent(inout) :: CSF(:)

  integer(iwp), intent(in) :: Direction
  integer(iwp), intent(in) :: IRefSm

  integer(iwp), Parameter :: iState=1
  integer(iwp) :: Mult, iCSF, iSD, nSD, nCSF, MinOpen, MaxOpen, iOpen, nOC
  real(wp), pointer:: T(:,:)

  Mult = SGS(iState)%iSpin
  MinOpen = mult - 1
  MaxOpen = Min(SGS(iState)%nActEl, 2*SGS(iState)%nLev-SGS(iState)%nActel)
  if (mod(MaxOpen-MinOpen,2) /= 0) MaxOpen = MaxOpen - 1

! For each sector of open-shells the data is stored as
!
! CDC (nSD x nCSF)
! SD (nSD x nOC)
! CSF (nCSF x nOC)

  iCSF = 1
  iSD  = 1
  do iOpen = MinOpen, MaxOpen, 2

#ifdef _SANITY_CHECK_
     if (.not. allocated(SGS(iState)%CDCInfo(iOpen,Mult)%CDC)) then
        write(u6,*) 'SD_CSF_Codec: CDC not allocated'
        call Abend()
     endif
     if (.not. allocated(SGS(iState)%CDCInfo(iOpen,Mult)%DTOC) .and. Mult == 1) then
        write(u6,*) 'SD_CSF_Codec: DTOC not allocated'
        call Abend()
     endif
#endif

!    This code is only intended for conversion between the CASVB format and
!    the format of FAROALD and SGUGA.
!    if (Mult == 1) then
!       T => SGS(iState)%CDCInfo(iOpen,Mult)%DTOC
!    else
        T => SGS(iState)%CDCInfo(iOpen,Mult)%CDC
!    endif

     nSD  = size(T,1)
     nCSF = size(T,2)
     nOC = SGS(iState)%CDCInfo(iOpen,Mult)%nOCSym(IRefSm)

     select case(Direction)

     case(CSF_TO_SD)

        call DGEMM_('N','N', &
                    nSD,nOC,nCSF, &
                    One,T,nSD, &
                        CSF(iCSF),nCSF, &
                    Zero,SD(iSD),nSD)

     case(SD_TO_CSF)

        call DGEMM_('T','N',   &
                    nCSF,nOC,nSD,   &
                    One,T,nSD,   &
                        SD(iSD),nSD,   &
                    Zero,CSF(iCSF),nCSF)

     end select

     iCSF = iCSF + nCSF*nOC
     iSD  = iSD  + nSD *nOC

  end do

  end subroutine SD_CSF_Codec

end module SD_CSF_Codec_Module
