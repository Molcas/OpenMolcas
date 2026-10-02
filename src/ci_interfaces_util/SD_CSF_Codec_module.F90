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

#define _SANITY_CHECK_
module SD_CSF_Codec_Module
  use stdalloc, only: mma_allocate, mma_deallocate
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
  real(wp), allocatable :: SD_Work(:)
  real(wp), pointer:: T(:,:)
  integer(iwp) :: i, iSD_REO
  real(wp) :: Phase

  Mult = SGS(iState)%iSpin
  MinOpen = mult - 1
  MaxOpen = Min(SGS(iState)%nActEl, 2*SGS(iState)%nLev-SGS(iState)%nActel)
  if (mod(MaxOpen-MinOpen,2) /= 0) MaxOpen = MaxOpen - 1

#ifdef _SANITY_CHECK_
  if (.not. allocated(SGS(iState)%SDREO)) then
     write(u6,*) 'SD_CSF_Codec: SDREO not allocated'
     call Abend()
  endif
  if (.not. allocated(SGS(iState)%SDREO(IRefSm)%List)) then
     write(u6,*) 'SD_CSF_Codec: SDREO List not allocated for IRefSm=', IRefSm
     call Abend()
  endif
  if (size(SGS(iState)%SDREO(IRefSm)%List) /= size(SD)) then
     write(u6,*) 'SD_CSF_Codec: inconsistent SDREO/SD sizes', &
                 size(SGS(iState)%SDREO(IRefSm)%List), size(SD)
     call Abend()
  endif
  if (any(abs(SGS(iState)%SDREO(IRefSm)%List) < 1_iwp) .or. &
      any(abs(SGS(iState)%SDREO(IRefSm)%List) > size(SD))) then
     write(u6,*) 'SD_CSF_Codec: SDREO List contains an invalid address'
     call Abend()
  endif
#endif

  call mma_allocate(SD_Work,size(SD),Label='SD_Work')

! SDREO(iRefSm)%List maps a determinant in the CDC/LUCIA ordering used by
! the codec to its position in the externally stored, resorted SD vector:
!
!     SDREO(iRefSm)%List(iCDC) = iResorted
!
! DGEMM therefore always operates on SD_Work in CDC/LUCIA order.  The
! permutation is applied after CSF_TO_SD and undone before SD_TO_CSF.

  select case(Direction)
  case(CSF_TO_SD)
     SD_Work(:) = Zero
  case(SD_TO_CSF)
     do i = 1, size(SD)
        iSD_REO = abs(SGS(iState)%SDREO(IRefSm)%List(i))
        Phase = merge(One,-One,SGS(iState)%SDREO(IRefSm)%List(i) > 0_iwp)
        SD_Work(i) = Phase * SD(iSD_REO)
     end do
  case default
     write(u6,*) 'SD_CSF_Codec: invalid Direction=', Direction
     call Abend()
  end select

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
                        CSF(iCSF:iCSF+nCSF*nOC-1),nCSF, &
                    Zero,SD_Work(iSD:iSD+nSD*nOC-1),nSD)
     case(SD_TO_CSF)
        call DGEMM_('T','N', &
                    nCSF,nOC,nSD, &
                    One,T,nSD, &
                        SD_Work(iSD:iSD+nSD*nOC-1),nSD, &
                    Zero,CSF(iCSF:iCSF+nCSF*nOC-1),nCSF)
     end select

     iCSF = iCSF + nCSF*nOC
     iSD  = iSD  + nSD*nOC
  end do

#ifdef _SANITY_CHECK_
  if (iSD-1 /= size(SD_Work)) then
     write(u6,*) 'SD_CSF_Codec: codec SD size mismatch', iSD-1, size(SD_Work)
     call Abend()
  endif
#endif

  if (Direction == CSF_TO_SD) then
     do i = 1, size(SD)
        iSD_REO = abs(SGS(iState)%SDREO(IRefSm)%List(i))
        Phase = merge(One,-One, SGS(iState)%SDREO(IRefSm)%List(i) > 0_iwp)

        SD(iSD_REO) = Phase * SD_Work(i)
     end do
  endif


  Call mma_deallocate(SD_Work)
  end subroutine SD_CSF_Codec
end module SD_CSF_Codec_Module
