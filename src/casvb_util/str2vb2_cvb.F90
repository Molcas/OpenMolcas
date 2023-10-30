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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine str2vb2_cvb(bikcof,cvb,cvbdet,iway,idetvb,i2s,nS,nalf1,nMs,absym,ndetvb,nvb,kbasis,nel,nconfion)

use casvb_global, only: ifnss1, ifnss2, ikcoff, ndetvbs
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: bikcof(*)
integer(kind=iwp), intent(in) :: iway, ndetvb, idetvb(ndetvb), nS, i2S(nS), nMs, nalf1(nMs), nvb, kbasis, nel, nconfion(0:*)
real(kind=wp), intent(inout) :: cvb(nvb), cvbdet(ndetvb)
logical(kind=iwp), intent(in) :: absym
integer(kind=iwp) :: i2s_keep, i_det, i_spin, iconfadd, idadd, idet, iMs, ioff, ioff_bikcof, ion, iS, isadd, j_spin, n_det, &
                     n_det_values, n_spin, n_spin_values, nalfsing, nalfsing_det, nalfsing_keep, nelsing
integer(kind=iwp), pointer :: ifnss(:,:)
real(kind=wp), allocatable :: tmp(:,:), w(:)
real(kind=wp), parameter :: sq2 = sqrt(Two), sqp5 = sqrt(Half)

if (kbasis == 6) then
  ifnss => ifnss2
else
  ifnss => ifnss1
end if

i2s_keep = 0 ! dummy initialize
nalfsing_keep = 0 ! dummy initialize
! Determinant to structure transformation
call mma_allocate(w,ndetvb,label='w')
if (iway == 1) then
  cvb(:) = Zero
  do idet=1,ndetvb
    w(idet) = cvbdet(idetvb(idet))
  end do
else if (iway == 2) then
  w(:) = Zero
end if
idadd = 0
isadd = 0
iconfadd = 0
do ion=0,nel/2
  if (nconfion(ion) == 0) cycle
  nelsing = nel-2*ion
  ! Investigate different S and Ms possibilities and
  ! prepare to collect different BIKCOF matrices if necessary:
  n_spin = 0
  n_spin_values = 0
  do iS=1,nS
    if (i2s(iS) <= nelsing) then
      n_spin = n_spin+ifnss(nelsing,i2s(iS))
      n_spin_values = n_spin_values+1
      i2s_keep = i2s(iS)
    end if
  end do
  n_det = 0
  n_det_values = 0
  do iMs=1,nMs
    nalfsing = nalf1(iMs)-ion
    if (nalfsing >= 0) then
      n_det = n_det+ndetvbs(nelsing,nalfsing)
      n_det_values = n_det_values+1
      nalfsing_keep = nalfsing
    end if
  end do
  if (kbasis == 6) then
    do iS=1,nS
      nalfsing_det = (nelsing+i2s(iS))/2
      if (i2s(iS) <= nelsing) then
        do iMs=1,nMs
          nalfsing = nalf1(iMs)-ion
          if (nalfsing >= 0) then
            if (nalfsing /= nalfsing_det) cycle
            if (iway == 1) then
              do idet=1,ifnss(nelsing,i2s(iS))
                if ((i2s(iS) == 0) .and. absym .and. (ndetvbs(nelsing,nalfsing) /= 1)) then
                  call daxpy_(nconfion(ion),sq2,w(idet+idadd),n_det,cvb(idet+isadd),n_spin)
                else
                  call daxpy_(nconfion(ion),One,w(idet+idadd),n_det,cvb(idet+isadd),n_spin)
                end if
              end do
            else if (iway == 2) then
              do idet=1,ifnss(nelsing,i2s(iS))
                if ((i2s(iS) == 0) .and. absym .and. (ndetvbs(nelsing,nalfsing) /= 1)) then
                  call daxpy_(nconfion(ion),sqp5,cvb(idet+isadd),n_spin,w(idet+idadd),n_det)
                  call daxpy_(nconfion(ion),sqp5,cvb(idet+isadd),n_spin,w(ndetvbs(nelsing,nalfsing)-idet+1+idadd),n_det)
                else
                  call daxpy_(nconfion(ion),One,cvb(idet+isadd),n_spin,w(idet+idadd),n_det)
                end if
              end do
            end if
          end if
        end do
      end if
    end do
    ! Skip collection if not necessary ...
  else if ((n_spin_values == 1) .and. (n_det_values == 1)) then
    if (iway == 1) then
      call mxattbp_cvb(bikcof(1+ikcoff(nelsing,nalfsing_keep,i2s_keep)),w(1+idadd),n_spin,n_det,nconfion(ion),cvb(1+isadd))
    else if (iway == 2) then
      call mxatbp_cvb(bikcof(1+ikcoff(nelsing,nalfsing_keep,i2s_keep)),cvb(1+isadd),n_det,n_spin,nconfion(ion),w(1+idadd))
    end if
  else
    call mma_allocate(tmp,n_det,n_spin,label='tmp')
    tmp(:,:) = Zero
    i_spin = 1
    do iS=1,nS
      if (i2s(iS) <= nelsing) then
        i_det = 1
        do iMs=1,nMs
          nalfsing = nalf1(iMs)-ion
          if (nalfsing >= 0) then
            if (ikcoff(nelsing,nalfsing,i2s(iS)) /= -1) then
              ioff_bikcof = 1+ikcoff(nelsing,nalfsing,i2s(iS))
              ioff = i_spin
              do j_spin=1,ifnss(nelsing,i2s(iS))
                tmp(i_det:idet+ndetvbs(nelsing,nalfsing)-1,ioff) = bikcof(ioff_bikcof:ioff_bikcof+ndetvbs(nelsing,nalfsing)-1)
                ioff_bikcof = ioff_bikcof+ndetvbs(nelsing,nalfsing)
                ioff = ioff+1
              end do
            end if
            i_det = i_det+ndetvbs(nelsing,nalfsing)
          end if
        end do
        i_spin = i_spin+ifnss(nelsing,i2s(iS))
      end if
    end do

    if (iway == 1) then
      call mxattbp_cvb(tmp,w(1+idadd),n_spin,n_det,nconfion(ion),cvb(1+isadd))
    else if (iway == 2) then
      call mxatbp_cvb(tmp,cvb(1+isadd),n_det,n_spin,nconfion(ion),w(1+idadd))
    end if
    call mma_deallocate(tmp)
  end if
  isadd = isadd+nconfion(ion)*n_spin
  idadd = idadd+nconfion(ion)*n_det
  iconfadd = iconfadd+nconfion(ion)
end do
if (iway == 2) then
  do idet=1,ndetvb
    cvbdet(idetvb(idet)) = w(idet)
  end do
end if
call mma_deallocate(w)
nullify(ifnss)

return

end subroutine str2vb2_cvb
