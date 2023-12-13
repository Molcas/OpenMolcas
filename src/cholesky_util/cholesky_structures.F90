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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

module Cholesky_Structures

! An extension to Data_Structures
! This is in a separate module because it needs the Cholesky module,
! which needs Data_Structures, so there would be a circular dependency.

use Data_Structures, only: Allocate_DT, Deallocate_DT, V1
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
private

public :: Allocate_DT, Deallocate_DT, L_Full_Type, Lab_Type

type L_Full_Pointers
  real(kind=wp), contiguous, pointer :: A3(:,:,:) => null()
  real(kind=wp), contiguous, pointer :: A21(:,:) => null()
  real(kind=wp), contiguous, pointer :: A12(:,:) => null()
end type L_Full_Pointers

type L_Full_Type
  integer(kind=iwp) :: iCase = 0
  integer(kind=iwp) :: iSym = 0
  integer(kind=iwp) :: nSym = 0
  integer(kind=iwp) :: nShell = 0
  real(kind=wp), allocatable :: A0(:)
  type(L_Full_Pointers), allocatable :: SPB(:,:,:)
end type L_Full_Type

type Lab_Type
  integer(kind=iwp) :: nSym = 0
  integer(kind=iwp) :: nDen = 0
  integer(kind=iwp) :: nShell = 0
  real(kind=wp), allocatable :: A0(:)
  logical(kind=iwp), allocatable :: Keep(:,:)
  type(V1), allocatable :: SB(:,:,:)
end type Lab_Type

! Extend allocate/deallocate data types
interface Allocate_DT
  module procedure :: Allocate_L_Full, Allocate_Lab
end interface Allocate_DT
interface Deallocate_DT
  module procedure :: Deallocate_L_Full, Deallocate_Lab
end interface Deallocate_DT

! Private extensions to mma interfaces
interface cptr2loff
  module procedure :: lfp_cptr2loff, v1_cptr2loff
end interface
interface mma_allocate
  module procedure :: lfp_mma_allo_3D, lfp_mma_allo_3D_lim, v1_mma_allo_3D, v1_mma_allo_3D_lim
end interface
interface mma_deallocate
  module procedure :: lfp_mma_free_3D, v1_mma_free_3D
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                      L F u l l - T Y P E   S E C T I O N             !
!                                                                      !
!                      L  full storage shell-pair blocked              !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Allocate_L_Full(Adam,nShell,iShp_rs,JNUM,JSYM,nSym,Memory)

  use Index_Functions, only: iTri, nTri_Elem
  use Cholesky, only: nBasSh, nnBstRSh

  type(L_Full_Type), target, intent(out) :: Adam
  integer(kind=iwp), intent(in) :: nShell, iShp_rs(nTri_Elem(nShell)), JNUM, JSYM, nSym
  integer(kind=iwp), optional, intent(out) :: Memory(2)
  integer(kind=iwp) :: iaSh, ibSh, iShp, iSyma, iSymb, LFULL, iS, iE, MemSPB, n1, n2

  LFULL = 0
  do iaSh=1,nShell
    do ibSh=1,iaSh
      iShp = iTri(iaSh,ibSh)

      if (iShp_rs(iShp) <= 0) cycle

      if (nnBstRSh(Jsym,iShp_rs(iShp),1) <= 0) cycle

      do iSymb=1,nSym
        iSyma = Mul(iSymb,Jsym)
        if (iSyma < iSymb) cycle

        LFULL = LFULL+nBasSh(iSyma,iaSh)*nBasSh(iSymb,ibSh)
        if (iaSh == ibSh) cycle

        LFULL = LFULL+nBasSh(iSyma,ibSh)*nBasSh(iSymb,iaSh)

      end do

    end do
  end do
  LFULL = LFULL*JNUM

  if (present(Memory)) then
    MemSPB = nSym*nShell*(nShell+1)
    MemSPB = (MemSPB*storage_size(Adam%SPB)-1)/storage_size(Adam%A0)+1
    Memory = [LFULL,MemSPB]
    return
  end if

  Adam%iCase = 1
  Adam%nSym = nSym
  Adam%iSym = JSYM
  Adam%nShell = nShell

  call mma_allocate(Adam%A0,LFULL,Label='Adam%A0')

  call mma_allocate(Adam%SPB,nSym,nTri_Elem(nShell),2,label='Adam%SPB')
# include "macros.fh"
  unused_proc(mma_allocate(Adam%SPB,[0,0],[0,0],[0,0]))

  iE = 0
  do iaSh=1,nShell
    do ibSh=1,iaSh
      iShp = iTri(iaSh,ibSh)

      if (iShp_rs(iShp) <= 0) cycle

      if (nnBstRSh(Jsym,iShp_rs(iShp),1) <= 0) cycle

      do iSymb=1,nSym
        iSyma = Mul(iSymb,Jsym)
        if (iSyma < iSymb) cycle

        iS = iE+1

        n1 = nBasSh(iSyma,iaSh)
        n2 = nBasSh(iSymb,ibSh)

        iE = iE+n1*JNUM*n2

        Adam%SPB(iSyma,iShp_rs(iShp),1)%A3(1:n1,1:JNUM,1:n2) => Adam%A0(iS:iE)
        Adam%SPB(iSyma,iShp_rs(iShp),1)%A21(1:n1*JNUM,1:n2) => Adam%A0(iS:iE)
        Adam%SPB(iSyma,iShp_rs(iShp),1)%A12(1:n1,1:JNUM*n2) => Adam%A0(iS:iE)

        if (iaSh == ibSh) cycle

        iS = iE+1

        n1 = nBasSh(iSyma,ibSh)
        n2 = nBasSh(iSymb,iaSh)

        iE = iE+n1*JNUM*n2

        Adam%SPB(iSyma,iShp_rs(iShp),2)%A3(1:n1,1:JNUM,1:n2) => Adam%A0(iS:iE)
        Adam%SPB(iSyma,iShp_rs(iShp),2)%A21(1:n1*JNUM,1:n2) => Adam%A0(iS:iE)
        Adam%SPB(iSyma,iShp_rs(iShp),2)%A12(1:n1,1:JNUM*n2) => Adam%A0(iS:iE)

      end do

    end do
  end do

end subroutine Allocate_L_Full

subroutine Deallocate_L_Full(Adam)

  use Index_Functions, only: iTri

  type(L_Full_Type), intent(inout) :: Adam
  integer(kind=iwp) :: iaSh, ibSh, iShp, iSyma

  do iaSh=1,Adam%nShell
    do ibSh=1,iaSh
      iShp = iTri(iaSh,ibSh)

      do iSyma=1,Adam%nSym

        Adam%SPB(iSyma,iShp,1)%A3 => null()
        Adam%SPB(iSyma,iShp,1)%A21 => null()
        Adam%SPB(iSyma,iShp,1)%A12 => null()
        Adam%SPB(iSyma,iShp,2)%A3 => null()
        Adam%SPB(iSyma,iShp,2)%A21 => null()
        Adam%SPB(iSyma,iShp,2)%A12 => null()

      end do

    end do
  end do

  call mma_deallocate(Adam%SPB)
  call mma_deallocate(Adam%A0)
  Adam%iCase = 0
  Adam%nSym = 0
  Adam%iSym = 0
  Adam%nShell = 0

end subroutine Deallocate_L_Full

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                      L a b - T Y P E   S E C T I O N                 !
!                                                                      !
!                      Lab storaged                                    !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Allocate_Lab(Lab,JNUM,nBasSh,nBas,nShell,nSym,nDen,Memory)

  type(Lab_Type), target, intent(out) :: Lab
  integer(kind=iwp), intent(in) :: JNUM, nShell, nSym, nBasSh(nSym,nShell), nBas(nSym), nDen
  integer(kind=iwp), optional, intent(out) :: Memory(2)
  integer(kind=iwp) :: iSym, iDen, Lab_Memory, iE, iS, iSh, MemKeep, MemSB

  Lab_Memory = 0
  do iSym=1,nSym
    Lab_Memory = max(nBas(iSym),Lab_Memory)
  end do
  Lab_Memory = Lab_Memory*JNUM*nDen

  if (present(Memory)) then
    MemKeep = nShell*nDen
    MemKeep = (MemKeep*storage_size(Lab%Keep)-1)/storage_size(Lab%A0)+1
    MemSB = nShell*nSym*nDen
    MemSB = (MemSB*storage_size(Lab%SB)-1)/storage_size(Lab%A0)+1
    Memory = [Lab_Memory,MemKeep+MemSB]
    return
  end if

  Lab%nSym = nSym
  Lab%nDen = nDen
  Lab%nShell = nShell
  call mma_allocate(Lab%A0,Lab_Memory,Label='Lab%A0')
  call mma_allocate(Lab%Keep,nShell,nDen,Label='Lab%Keep')
  call mma_allocate(Lab%SB,nShell,nSym,nDen,Label='Lab%SB')
# include "macros.fh"
  unused_proc(mma_allocate(Lab%SB,[0,0],[0,0],[0,0]))

  do iSym=1,nSym
    iE = 0
    do iDen=1,nDen
      do iSh=1,nShell

        iS = iE+1
        iE = iE+nBasSh(iSym,iSh)*JNUM

        Lab%SB(iSh,iSym,iDen)%A(1:nBasSh(iSym,iSh)*JNUM) => Lab%A0(iS:iE)

      end do
    end do
  end do

end subroutine Allocate_Lab

subroutine Deallocate_Lab(Lab)

  type(Lab_Type), intent(inout) :: Lab
  integer(kind=iwp) :: iSym, iDen, iSh

  do iSym=1,Lab%nSym
    do iDen=1,Lab%nDen
      do iSh=1,Lab%nShell

        Lab%SB(iSh,iSym,iDen)%A => null()

      end do
    end do
  end do

  Lab%nSym = 0
  Lab%nDen = 0
  Lab%nShell = 0
  call mma_deallocate(Lab%A0)
  call mma_deallocate(Lab%Keep)
  call mma_deallocate(Lab%SB)

end subroutine Deallocate_Lab

! Define lfp_cptr2loff, lfp_mma_allo_3D, lfp_mma_allo_3D_lim, lfp_mma_free_3D
!        v1_cptr2loff, v1_mma_allo_3D, v1_mma_allo_3D_lim, v1_mma_free_3D
#define _TYPE_ type(L_Full_Pointers)
#  define _FUNC_NAME_ lfp_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ lfp_mma
#  define _DIMENSIONS_ 3
#  define _DEF_LABEL_ 'lfp_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

#define _TYPE_ type(V1)
#  define _FUNC_NAME_ v1_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ v1_mma
#  define _DIMENSIONS_ 3
#  define _DEF_LABEL_ 'v1_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module Cholesky_Structures
