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

module Sizes_of_Seward

implicit none
private

public :: S, Size_Dmp, Size_Get

#include  "itmax.fh"
integer i
type Sizes_of_Stuff
  sequence
  integer :: Low_Anchor
  integer :: m2Max = 0
  integer :: mCentr = 0
  integer :: mCentr_Aux = 0
  integer :: mCentr_Frag = 0
  integer :: Mx_mdc = 0
  integer :: Mx_Shll = 0
  integer :: n2Tot = 0
  integer :: jMax = 5
  integer :: MaxPrm(0:iTabMx) = [(0,i=0,iTabMx)]
  integer :: MaxBas(0:iTabMx) = [(0,i=0,iTabMx)]
  integer :: nDim = 0
  integer :: nShlls = 0
  integer :: Max_Center = 15
  integer :: kCentr = 0
  integer :: nMltpl = 2
  integer :: iAngMx = -1
  integer :: High_Anchor
end type Sizes_of_Stuff

type(Sizes_of_Stuff), target :: S
integer, pointer :: p_ix(:)
integer Len, Len2
logical Found

interface
  subroutine Abend()
  end subroutine Abend
  subroutine Put_iArray(Label,data,nData)
    character*(*) Label
    integer nData
    integer data(nData)
  end subroutine Put_iArray
  subroutine Get_iArray(Label,data,nData)
    character*(*) Label
    integer nData
    integer data(nData)
  end subroutine Get_iArray
  subroutine Qpg_iArray(Label,Found,nData)
    character*(*) Label
    logical Found
    integer nData
  end subroutine Qpg_iArray
end interface

contains

subroutine Size_Init()

  use iso_c_binding

  integer, external :: ip_of_iWork

  Len = ip_of_iWork(S%High_Anchor)-ip_of_iWork(S%Low_Anchor)+1
  call c_f_pointer(c_loc(S%Low_Anchor),p_ix,[Len])

end subroutine Size_Init

subroutine Size_Dmp()

  call Size_Init()
  call Put_iArray('Sizes',p_ix,Len)
  nullify(p_ix)

end subroutine Size_Dmp

subroutine Size_Get()

  call Qpg_iArray('Sizes',Found,Len2)
  if (.not. Found) then
    write(6,*) 'Size_Get: Sizes not found.'
    call Abend()
  end if
  call Size_Init()
  if (Len /= Len2) then
    write(6,*) 'Size_Get: Len/=Len2.'
    call Abend()
  end if
  call Get_iArray('Sizes',p_ix,Len)
  nullify(p_ix)

end subroutine Size_Get

end module Sizes_of_Seward
