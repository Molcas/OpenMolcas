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

subroutine TMOMInt(wavevector,iOpt)
!***********************************************************************
!                                                                      *
! Object: driver for computation of TMOM integrals                     *
!                                                                      *
!***********************************************************************

!#define _DEBUGPRINT_
use Integral_interfaces, only: int_kernel, int_mem

#ifdef _DEBUGPRINT_
use Index_Functions, only: nTri_Elem, nTri_Elem1
use OneDat, only: sOpSiz
use Sizes_of_Seward, only: S
use Basis_Info, only: nBas
use Symmetry_Info, only: Mul, nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: u6
#endif
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: wavevector(3)
integer(kind=iwp), intent(in) :: iOpt
integer(kind=iwp) :: nComp, nOrdOp
real(kind=wp) :: dum(1), rHrmt
character(len=8) :: Label
integer(kind=iwp), allocatable :: ipList(:), OperC(:), OperI(:)
real(kind=wp), allocatable :: CoorO(:), Nuc(:)
procedure(int_kernel) :: EMFInt
procedure(int_mem) :: EMFMem
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, iCase, iComp, idum(1), ij, ilen, iMltpl, IOFF(8,8), iOpt0, iOpt1, iRc, iSyLbl, iSyLbl_TMOM, ix, iy, iz, j, &
                     jOff, Len_, Length, nInts, nInts_TMOM, Phase
real(kind=wp), allocatable :: Fact, Int_I(:), Int_I_O(:), Int_R(:), Int_R_O(:), Temp, Temp_Int(:), x, xy, xyz
#endif

#include "warnings.h"
! ipList: list of pointers to the integrals of each component
!         of the operator
! OperI: list which irreps a particular component of the operator
!        belongs to
! OperC: list the character of each component of the operator
! CoorO: list of origins of the operator, one for each component

call Set_Basis_Mode('Valence')
call Setup_iSD()

!***********************************************************************
!***********************************************************************
!                                                                      *
!     Electromagnetic field radiation integrals.                       *
!                                                                      *
!     Note that the integral is neither symmetric nor antisymmetric!   *
!                                                                      *
!***********************************************************************
!***********************************************************************
rHrmt = -One ! Not used

! B*s Magnetic * Spin: not that this boils down to just integrals
! over A.

if (iOpt == 2) then
  rHrmt = One
  nOrdOp = 0
  Label = 'TMOM0'
  nComp = 2
  call Allocate_Aux()
  ! Here we put in the k-vector
  CoorO(1:3) = wavevector(:)
  CoorO(4:) = Zero

  ! The electromagnetic field operator contributes to all
  ! irreducible irreps, hence OperI=255. Since the operator
  ! itself is not symmetry-adapted OperC is set to a dummy value.

  OperI(:) = 255
  OperC(:) = 0 ! Dummy

  Nuc(:) = Zero
  call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,0)

  call Deallocate_Aux()
# ifdef _DEBUGPRINT_

  call mma_allocate(CoorO,6,Label='CoorO')
  CoorO(:) = Zero
  CoorO(1:3) = wavevector
  write(u6,*) 'Wavevector=',Wavevector

  ! This section of the code is for pure debugging and will replace
  ! exact operator with truncated expansions of the operator in
  ! terms of multipole integrals

  iOpt0 = 0 ! Write
  iOpt1 = ibset(0,sOpSiz) ! Read just data size and symmetry
  iRc = -1
  Label = 'TMOM0  R'
  iComp = 1
  ! Pick up the size and the symmetry label.
  call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl_TMOM)
  nInts_TMOM = idum(1)
  call mma_allocate(Int_R,nInts_TMOM+4,Label='Int_R')
  call mma_allocate(Int_I,nInts_TMOM+4,Label='Int_I')
  call mma_allocate(Int_R_O,nInts_TMOM+4,Label='Int_R_O')
  call mma_allocate(Int_I_O,nInts_TMOM+4,Label='Int_I_O')

  call RdOne(iRc,iOpt0,Label,iComp,Int_R_O,iSyLbl_TMOM)
  Label = 'TMOM0  I'
  iComp = 1
  call RdOne(iRc,iOpt0,Label,iComp,Int_I_O,iSyLbl_TMOM)
  Length = 0
  IOFF(:,:) = -1
  do i=1,nIrrep
    do j=1,i
      ij = Mul(i,j)-1
      if (btest(iSyLbl_TMOM,ij)) then
        IOFF(i,j) = Length+1
        if (i == j) then
          Len_ = nTri_Elem(nBas(i-1))
        else
          Len_ = nBas(i-1)*nBas(j-1)
        end if
        Length = Length+Len_
      end if
    end do
  end do

  Int_R(:) = Zero
  Int_R(nInts_TMOM+1:nInts_TMOM+3) = CoorO
  Int_I(:) = Zero
  Int_I(nInts_TMOM+1:nInts_TMOM+3) = CoorO

  S%nMltpl = 9
  iCase = 1
  Phase = One
  do iMltpl=0,S%nMltpl
    write(Label,'(A,I2)') 'Mltpl ',iMltpl
    nComp = nTri_Elem1(iMltpl)
    iComp = 0
    do ix=iMltpl,0,-1
      x = CoorO(1)**ix
      do iy=iMltpl-ix,0,-1
        xy = x*CoorO(2)**iy
        iz = iMltpl-ix-iy
        xyz = xy*CoorO(3)**iz

        Fact = Phase*xyz/(gamma(real(ix+1,kind=wp))*gamma(real(iy+1,kind=wp))*gamma(real(iz+1)))

        iComp = iComp+1
        if (Fact == Zero) cycle
        call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
        if (iRC /= 0) then
          write(u6,*) 'TMOMINT: Error reading ',Label
          call Abend()
        end if
        nInts = idum(1)
        call mma_allocate(Temp_Int,nInts+4,Label='Temp_Int')
        call RdOne(iRc,iOpt0,Label,iComp,Temp_Int,iSyLbl)

        Length = 0
        do i=1,nIrrep
          do j=1,i
            ij = Mul(i,j)-1
            if (btest(iSyLbl,ij)) then
              jOff = Length+1
              if (i == j) then
                Len_ = nTri_Elem(nBas(i-1))
              else
                Len_ = nBas(i-1)*nBas(j-1)
              end if
              if (iCase == 1) then
                ! Contribution to the real part
                Int_R(IOFF(i,j):IOFF(i,j)+Len_-1) = Int_R(IOFF(i,j):IOFF(i,j)+Len_-1)+Fact*Temp_Int(jOff:jOff+Len_-1)
              else
                ! Contribution to the imaginary part
                Int_I(IOFF(i,j):IOFF(i,j)+Len_-1) = Int_I(IOFF(i,j):IOFF(i,j)+Len_-1)+Fact*Temp_Int(jOff:jOff+Len_-1)
              end if
              Length = Length+Len_
            end if
          end do
        end do
        call mma_deallocate(Temp_Int)
      end do
    end do

    if (iCase == 1) then
      iCase = 2
    else
      iCase = 1
      Phase = -Phase
    end if
  end do

  ! Compare exact integrals with approximated.

# define _COMPARE_
# ifdef _COMPARE_
  Length = 0
  do i=1,nIrrep
    do j=1,i
      ij = Mul(i,j)-1
      if (btest(iSyLbl_TMOM,ij)) then
        if (i == j) then
          Len_ = nTri_Elem(nBas(i-1))
        else
          Len_ = nBas(i-1)*nBas(j-1)
        end if
        do iLen=1,Len_
          temp = abs(Int_R_O(Length+iLen)-Int_R(Length+iLen))/max(abs(Int_R_O(Length+iLen)),abs(Int_R(Length+iLen)),1.0e-8_wp)
          if (temp > 1.0e-2_wp) then
            write(u6,*) 'isym,jsym,iLen=',i,j,iLen
            write(u6,*) 'Int_R,Int_Q=',Int_R_O(Length+iLen),Int_R(Length+iLen)
          end if
          temp = abs(Int_I_O(Length+iLen)-Int_I(Length+iLen))/max(abs(Int_I_O(Length+iLen)),abs(Int_I(Length+iLen)),1.0e-8_wp)
          if (temp > 1.0e-2_wp) then
            write(u6,*) 'isym,jsym,iLen=',i,j,iLen
            write(u6,*) 'Int_I,Int_J=',Int_I_O(Length+iLen),Int_I(Length+iLen)
          end if
        end do
        Length = Length+Len_
      end if
    end do
  end do
# endif

  ! Overwrite the integrals with a truncated expansion.

  Label = 'TMOM0  R'
  iComp = 1
  call WrOne(iRc,iOpt0,Label,iComp,Int_R,iSyLbl_TMOM)
  Label = 'TMOM0  I'
  iComp = 1
  call WrOne(iRc,iOpt0,Label,iComp,Int_I,iSyLbl_TMOM)

  call mma_deallocate(Int_R_O)
  call mma_deallocate(Int_I_O)
  call mma_deallocate(Int_R)
  call mma_deallocate(Int_I)
  call mma_deallocate(CoorO)

#endif
end if

! A*nabla. Note that when used the numbers are multiplied with -i to
! generate A*p.

if (iOpt <= 2) then
  nOrdOp = 1
  Label = 'TMOM'
  nComp = 12
  call Allocate_Aux()
  ! Here we put in the k-vector
  CoorO(1:3) = wavevector(:)
  CoorO(4:) = Zero

  ! The electromagnetic field operator contributes to all
  ! irreducible irreps, hence OperI=255. Since the operator
  ! itself is not symmetry-adapted OperC is set to a dummy value.

  OperI(:) = 255
  OperC(:) = 0 ! Dummy

  Nuc(:) = Zero
  call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,0)

  call Deallocate_Aux()
end if

! The A^2 term

if (iOpt > 2) then
  nOrdOp = 0
  Label = 'TMOM2'
  nComp = 2
  call Allocate_Aux()
  ! Here we put in the k-vector
  ! Change the argument to 2xA
  CoorO(1:3) = Two*wavevector(:)
  CoorO(4:) = Zero

  ! The electromagnetic field operator contributes to all
  ! irreducible irreps, hence OperI=255. Since the operator
  ! itself is not symmetry-adapted OperC is set to a dummy value.

  OperI(:) = 255
  OperC(:) = 0 ! Dummy

  Nuc(:) = Zero
  call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,CoorO,nOrdOp,Nuc,rHrmt,OperC,dum,1,0)

  call Deallocate_Aux()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *

contains

subroutine Allocate_Aux()

  use stdalloc, only: mma_allocate

  call mma_Allocate(ipList,nComp,Label='ipList')
  call mma_Allocate(OperI,nComp,Label='OperI')
  call mma_Allocate(OperC,nComp,Label='OperC')
  call mma_Allocate(CoorO,3*nComp,Label='CoorO')
  call mma_Allocate(Nuc,nComp,Label='Nuc')

end subroutine Allocate_Aux

subroutine Deallocate_Aux()

  use stdalloc, only: mma_deallocate

  call mma_Deallocate(OperC)
  call mma_Deallocate(OperI)
  call mma_Deallocate(ipList)
  call mma_Deallocate(CoorO)
  call mma_Deallocate(Nuc)

end subroutine Deallocate_Aux

end subroutine TMOMInt
