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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                  Fail)
!***********************************************************************
!                                                                      *
! Object: to change the length of the primitive and basis functions    *
!         vectors of center D and B, and center  D, B, C, and A,       *
!         respectively.                                                *
!                                                                      *
! Called from: PSOAOx                                                  *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!***********************************************************************

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iBas, kBas, jBas, lBas, jPrim, lPrim
integer(kind=iwp), intent(inout) :: iBsInc, kBsInc, jBsInc, lBsInc, jPrInc, lPrInc
logical(kind=iwp), intent(inout) :: QiBas, QkBas, QjBas, QlBas, QjPrim, QlPrim
logical(kind=iwp), intent(out) :: Fail
integer(kind=iwp) :: i

Fail = .false.
if (QlPrim) then
  if (lPrInc == 1) then
    QlPrim = .false.
    QjPrim = .true.
  else
    do i=2,lPrim
      if ((lPrim+1)/i < lPrInc) then
        lPrInc = max(1,(lPrim+1)/i)
        return
      end if
    end do
  end if
end if
if (QjPrim) then
  lPrInc = lPrim
  if (jPrInc == 1) then
    QjPrim = .false.
    QlBas = .true.
  else
    do i=2,jPrim
      if ((jPrim+1)/i < jPrInc) then
        jPrInc = max(1,(jPrim+1)/i)
        return
      end if
    end do
  end if
end if
lPrInc = lPrim
jPrInc = jPrim
if (QlBas) then
  if (lBsInc == 1) then
    QlBas = .false.
    QjBas = .true.
  else
    do i=2,lBas
      if ((lBas+1)/i < lBsInc) then
        lBsInc = max(1,(lBas+1)/i)
        QlPrim = .true.
        return
      end if
    end do
  end if
end if
if (QjBas) then
  lBsInc = lBas
  if (jBsInc == 1) then
    QjBas = .false.
    QkBas = .true.
  else
    do i=2,jBas
      if ((jBas+1)/i < jBsInc) then
        jBsInc = max((jBas+1)/i,1)
        QlPrim = .true.
        return
      end if
    end do
  end if
end if
if (QkBas) then
  lBsInc = lBas
  jBsInc = jBas
  if (kBsInc == 1) then
    QkBas = .false.
    QiBas = .true.
  else
    do i=2,kBas
      if ((kBas+1)/i < kBsInc) then
        kBsInc = max((kBas+1)/i,1)
        QlPrim = .true.
        return
      end if
    end do
  end if
end if
if (QiBas) then
  lBsInc = lBas
  jBsInc = jBas
  kBsInc = kBas
  if (iBsInc == 1) then
    Fail = .true.
    return
  end if
  do i=2,iBas
    if ((iBas+1)/i < iBsInc) then
      iBsInc = max(1,(iBas+1)/i)
      QlPrim = .true.
      return
    end if
  end do
end if

return

end subroutine Change
