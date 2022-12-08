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
! Copyright (C) 1991, Roland Lindh                                     *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Drvh1_mck(Nona)
!***********************************************************************
!                                                                      *
! Object: driver for computation of gradients of one-electron matrices.*
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '91                                              *
!             Modified by Anders Bernhardsson for Gradients            *
!             May 95                                                   *
!***********************************************************************

use mck_interface, only: grd_mck_kernel, mck_mem
use Index_Functions, only: nTri_Elem
use Basis_Info, only: dbsc, nBas, nCnttp
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: Nona
#include "print.fh"
integer(kind=iwp) :: i, iCnt, iCnttp, idCar, idcnt, iIrrep, loper, nDens, nFock
character(len=8) :: Label
logical(kind=iwp) :: lECP
real(kind=wp), allocatable :: D0(:), Fock(:)
procedure(grd_mck_kernel) :: KneGrd_mck, m1grd_mck, nagrd_mck, nonatwo, OvrGrd_mck, prjgrd_mck, srogrd_mck
procedure(mck_mem) :: KneMem_mck, m1mm1, na2mem, namem_mck, OvrMem_mck, prjmm1, sromm1

if (show) then
  nFock = 0
  nDens = 0
  do iIrrep=0,nIrrep-1
    nFock = nFock+nTri_Elem(nBas(iIrrep))
    nDens = nDens+nTri_Elem(nBas(iIrrep))
  end do

  ! Read the variational 1st order density matrix
  ! density matrix in AO/SO basis
  call mma_allocate(D0,nDens,Label='D0')
  call Get_D1ao_Var(D0,nDens)
  ! Read the generalized Fock matrix
  ! Fock matrix in AO/SO basis
  call mma_allocate(Fock,nFock,Label='Fock')
  call Get_dArray_chk('FockOcc',Fock,nFock)
else
  nFock = 1
  nDens = 1
  call mma_allocate(Fock,nFock,Label='Fock')
  call mma_allocate(D0,nDens,Label='D0')
  Fock(1) = Zero
  D0(1) = Zero
end if
if (Nona) then
  !*********************************************************************
  !0a)                                                                 *
  !     Antisymmetric gradient of Overlap matrix                       *
  !                                                                    *
  !*********************************************************************
  Label = 'OVRGRDA'
  idcnt = 0
  do iCnttp=1,nCnttp
    do iCnt=1,dbsc(iCnttp)%nCntr
      idcnt = idcnt+1
      do idCar=1,3
        call Cnt1El(OvrGrd_mck,OvrMem_mck,Label,idcnt,idcar,loper,-One,.false.,Fock,'OVRGRDA ',0)
      end do
    end do
  end do
  !*********************************************************************
  !0b)                                                                 *
  !     Non-adiabatic second derivative integrals                      *
  !                                                                    *
  !*********************************************************************
  Label = 'NONA2'
  idcnt = 0
  do iCnttp=1,nCnttp
    do iCnt=1,dbsc(iCnttp)%nCntr
      idcnt = idcnt+1
      do idCar=1,3
        call Cnt1El(NONATWO,NA2Mem,Label,idcnt,idcar,loper,One,.false.,Fock,'NONA2   ',0)
      end do
    end do
  end do

end if

!***********************************************************************
!1)                                                                    *
!     Gradient of Overlap matrix                                       *
!                                                                      *
!***********************************************************************
Label = 'OVRGRD'
idcnt = 0
do iCnttp=1,nCnttp
  do iCnt=1,dbsc(iCnttp)%nCntr
    idcnt = idcnt+1
    do idCar=1,3
      call Cnt1El(OvrGrd_mck,OvrMem_mck,Label,idcnt,idcar,loper,One,.false.,Fock,'OVRGRD  ',0)
    end do
  end do
end do
!
!***********************************************************************
!2)                                                                    *
!     Gradient of Kinetic operator                                     *
!                                                                      *
!***********************************************************************
Label = 'KNEGRD'
idcnt = 0
do iCnttp=1,nCnttp
  do iCnt=1,dbsc(iCnttp)%nCntr
    idcnt = idcnt+1
    do idCar=1,3
      call Cnt1El(KneGrd_mck,KneMem_mck,Label,idcnt,idcar,loper,One,.false.,D0,'ONEGRD  ',0)
    end do
  end do
end do
!
!***********************************************************************
!3)                                                                    *
!     Gradient of Nuclear attraction Operator                          *
!                                                                      *
!***********************************************************************
Label = 'NAGRD'
idcnt = 0
do iCnttp=1,nCnttp
  do iCnt=1,dbsc(iCnttp)%nCntr
    idcnt = idcnt+1
    do idCar=1,3
      call Cnt1El(NaGrd_mck,NaMem_mck,Label,idcnt,idcar,loper,One,.true.,D0,'ONEGRD  ',1)
    end do
  end do
end do
!
!
!***********************************************************************
!3)                                                                    *
!     Gradient of Nuclear attraction Operator  ECP-part                *
!                                                                      *
!***********************************************************************
lECP = .false.
do i=1,nCnttp
  lECP = lECP .or. dbsc(i)%ECP
end do
if (lecp) then
  idcnt = 0
  do iCnttp=1,nCnttp
    do iCnt=1,dbsc(iCnttp)%nCntr
      idcnt = idcnt+1
      do idCar=1,3
        Label = 'PRJGRD'
        call Cnt1El(Prjgrd_mck,PrjMm1,Label,idcnt,idcar,loper,One,.true.,D0,'ONEGRD  ',1)
        Label = 'M1GRD'
        call Cnt1El(m1grd_mck,m1Mm1,Label,idcnt,idcar,loper,One,.true.,D0,'ONEGRD  ',1)
        Label = 'SROGRD'
        call Cnt1El(Srogrd_mck,sroMm1,Label,idcnt,idcar,loper,One,.true.,D0,'ONEGRD  ',1)
      end do
    end do
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(D0)
call mma_deallocate(Fock)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Drvh1_mck
