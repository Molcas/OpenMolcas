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
! Copyright (C) 1998, Anders Bernhardsson                              *
!***********************************************************************

subroutine RHS_PT2(rkappa,CLag,SLag)
! The RHS array for CASPT2 has been already calculated in the
! CASPT2 module, so here we only need to read it from file

use MCLR_Data, only: LuPT2, nDens
use input_mclr, only: nSym, nRoots, nCSF, nOrb

implicit none
real*8 rKappa(nDens), CLag(*), SLag(nRoots,nRoots)
integer nOLag, nCLag, i, j
real*8 Tmp
integer istatus

! Read in a and b part of effective gradient from CASPT2

nOLag = 0
nCLag = 0
do i=1,nSym
  nOLag = nOLag+nOrb(i)*nOrb(i)
  nCLag = nCLag+nRoots*nCSF(i)
end do

do i=1,nCLag
  read(LUPT2,*,iostat=istatus) CLag(i)
  if (istatus < 0) call Error()
end do
do i=1,nOLag
  read(LUPT2,*,iostat=istatus) tmp ! rKappa(i)
  if (istatus < 0) call Error()
  rKappa(i) = rKappa(i)+tmp
end do
do j=1,nRoots
  do i=1,nRoots
    read(LUPT2,*,iostat=istatus) SLag(i,j)
    if (istatus < 0) call Error()
  end do
end do

return

contains

subroutine Error()

  use Definitions, only: u6

  write(u6,*)
  write(u6,'(1x,A)') 'The file which has to be written in CASPT2 module does not exist in RHS_PT2.'
  write(u6,'(1x,A)') 'For single-point gradient calculation, you need GRAD or GRDT keyword in &CASPT2.'
  write(u6,'(1x,A)') 'For geometry optimization, you do not need anything, so something is wrong with the code.'
  write(u6,*)
  call abend()

end subroutine Error

end subroutine RHS_PT2
