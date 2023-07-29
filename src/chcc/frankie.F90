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

subroutine frankie(nfro,no,nv,printkey)

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use Cholesky, only: timings
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nfro, no, nv, printkey
integer(kind=iwp) :: idum(1), nbas, ndel, nocc, norb, rc
real(kind=wp) :: FracMem
type(DSBA_Type) :: CMO

!.1 - get the info on  nBas, nOrb, nOcc. Use nFro from input

!# nbas = nfro + nocc + nvirt + ndel

call Get_iArray('nBas',idum,1)
nBas = idum(1)
call Get_iArray('nOrb',idum,1)
nOrb = idum(1)
call Get_iArray('nIsh',idum,1) ! in general > no
nOcc = idum(1)

ndel = nbas-no-nv-nfro

if (printkey >= 10) then
  write(u6,*) 'nbas = ',nbas
  write(u6,*) 'norb = ',norb
  write(u6,*) 'nocc = ',nocc
  write(u6,*) 'nfro = ',nfro
  write(u6,*) 'no   = ',no,' (nocc-nfro)'
  write(u6,*)
  write(u6,*) 'ndel = ',ndel
end if

if ((no+nfro+nv+ndel) /= nbas) then
  write(u6,*) 'Problem '
  write(u6,*) 'nbas from Runfile : ',nbas
  write(u6,*) 'nbas control      : ',nfro+no+nv+ndel
  call abend()
end if

timings = .false.
if (printkey > 1) timings = .true.

!.2 - allocate space for CMO with removed SCF deleted and frozen orbitals
!     final ordering of indices : (o+v,nbas)

call Allocate_DT(CMO,[no+nv],[nbas],1)
if (printkey >= 10) write(u6,*) 'Dopice 1 - Allo'

!.3 - read CMO
call read_mo(Cmo,nfro,no,nv,ndel,nbas,nOrb)
!.3 - invert the CMO matrix
FracMem = Zero ! in a parallel run set it to a sensible value
rc = 0
call Cho_X_init(rc,FracMem) ! initialize cholesky info
if (printkey >= 10) write(u6,*) 'Dopice 2 ',rc

call CHO_CC_drv(rc,CMO)
if (printkey >= 10) write(u6,*) 'Dopice 3 '

call Cho_X_final(rc)
if (printkey >= 10) write(u6,*) 'Dopice 4 '

if (rc /= 0) then
  write(u6,*) 'cho_cc_drv failed'
  call abend()
end if

!. - deallocate CMO
call Deallocate_DT(CMO)

return

end subroutine frankie
