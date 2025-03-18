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

use MCLR_Data, only: LuPT2
use input_mclr, only: nSym, nRoots, nCSF, nOrb

implicit none
real*8 rKappa(*), CLag(*), SLag(*)
integer nOLag, nCLag, i, nSLag
real*8 Tmp

! Read in a and b part of effective gradient from CASPT2

nOLag = 0
nCLag = 0
do i=1,nSym
  nOLag = nOLag+nOrb(i)*nOrb(i)
  nCLag = nCLag+nRoots*nCSF(i)
end do
nSLag = nRoots*nRoots

do i=1,nCLag
  read(LUPT2,*,end=200) CLag(i)
end do
do i=1,nOLag
  read(LUPT2,*,end=200) tmp ! rKappa(i)
  rKappa(i) = rKappa(i)+tmp
end do
do i=1,nSLag
  read(LUPT2,*,end=200) SLag(i)
end do

return

200 continue
write(6,*)
write(6,'(1x,A)') 'The file which has to be written in CASPT2 module does not exist in RHS_PT2.'
write(6,'(1x,A)') 'For single-point gradient calculation, you need GRAD or GRDT keyword in &CASPT2.'
write(6,'(1x,A)') 'For geometry optimization, you do not need anything, so something is wrong with the code.'
write(6,*)
call abend()

end subroutine RHS_PT2
