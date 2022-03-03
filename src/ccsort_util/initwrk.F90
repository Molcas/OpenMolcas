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

subroutine initwrk(length)
! this routine calculates required size of work space and
! defines initial positions of work vectors

use ccsort_global, only: fullprint, mapdri, mapiri, noa, NORB, NSYM, pos10, pos20, pos30, posri0, t3key
use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: length
integer(kind=iwp) :: n, norbmax, sizempq, sizeri, sizev1, sizev2, sizevint, syma, symi, symj, symm, symmp, symp, sympq, sympqi, symq

!1 def maxsize of vint

norbmax = norb(1)
do n=1,nsym
  if (norb(n) > norbmax) norbmax = norb(n)
end do

sizevint = norbmax*norbmax*norbmax

!2 def size of <pq|i>=j>, <pq|i,j>

sizev1 = 0
sizev2 = 0
do symp=1,nsym
  do symq=1,nsym
    sympq = mul(symp,symq)
    do symi=1,nsym
      sympqi = mul(sympq,symi)
      symj = sympqi
      ! calc. length
      if (symj > symi) then
        sizev2 = sizev2+noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
      else
        sizev1 = sizev1+noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
        sizev2 = sizev2+noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
      end if
    end do
  end do
end do

!3 def maxsize of <_am|pq>

sizempq = 0
do syma=1,nsym

  length = 0
  do symm=1,nsym
    do symp=1,nsym
      symmp = mul(symm,symp)
      symq = mul(syma,symmp)
      ! calc. length
      length = length+noa(symm)*NORB(symp)*NORB(symq)
    end do
  end do

  if (sizempq < length) sizempq = length

end do

!4 def maxsize of R_i if needed

sizeri = 0

if (t3key == 1) then
  do symi=1,nsym
    call ccsort_t3grc0(3,8,4,4,4,0,symi,1,length,mapdri,mapiri)
    length = length-1
    if (length > sizeri) sizeri = length
  end do
end if

! ******* distribution of memory ******

pos10 = 1+sizevint
pos20 = pos10+sizev1
pos30 = pos20+sizev2
posri0 = pos30+sizempq
length = posri0+sizeri-1

if (fullprint > 1) then
  write(u6,*)
  write(u6,'(6X,A)') 'size of help (work) vectors:'
  write(u6,'(6X,A)') '----------------------------'
  write(u6,*)
  write(u6,'(6X,A,I8)') 'Vints     V0 required : ',sizevint
  write(u6,'(6X,A,I8)') 'PQIJ ints V1 required : ',sizev1
  write(u6,'(6X,A,I8)') '          V2 required : ',sizev2
  write(u6,'(6X,A,I8)') 'AMIJ ints V3 required : ',sizempq
  write(u6,'(6X,A,I8)') 'R_i mtx   Ri required : ',sizeri
end if

if (fullprint >= 0) write(u6,'(6X,A,I20)') 'Required WRK size-sum : ',length

return

end subroutine initwrk
