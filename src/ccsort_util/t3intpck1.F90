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

subroutine t3intpck1(vint,r,dimv1,dimv2,dimv3,dima,dimbc,symq,symr,syms,nob,nvb)
! this routines pack integral block symi,symq,symr,syms
! R_i(a,bc) = V_i(b,a,c)
! for symq(b)=syms(c)
! and writes R block onto proper place of open t3nam file - lunt3
!
! vint  - integrals for given symmetries for given i (I)
! r     - final R_i matrix (O)
! dimv1 - 1st dimension of V (I)
! dimv2 - 2nd dimension of V (I)
! dimv3 - 3rd dimension of V (I)
! dima  - dimension of a in R (I)
! dimbc - dimension of bc in R (I)
! symq  - symmetry of q (b) (I)
! symr  - symmetry of r (a) (I)
! syms  - symmetry of s (c) (I)
! nob   - number of beta occupied in each irrep (I)
! nvb   - number of beta virtuals in each irrep (I)

implicit none
#include "reorg.fh"
#include "files_ccsd.fh"
integer symq, symr, syms
integer dimv1, dimv2, dimv3, dima, dimbc
integer nob(1:8)
integer nvb(1:8)
real*8 vint(1:dimv1,1:dimv2,1:dimv3)
real*8 r(1:dima,1:dimbc)
! help variables
integer a, b, c, bc, adda, length, iaddr

! if there are no beta virtuals - skip to write section
if (nvb(symq)*nvb(symr)*nvb(syms) /= 0) then

  ! calc additional constant for a
  adda = nob(symr)

  ! do packing

  bc = 0
  do b=nob(symq)+1,nob(symq)+nvb(symq)
    do c=nob(syms)+1,b
      bc = bc+1
      do a=1,nvb(symr)
        r(a,bc) = vint(b,a+adda,c)
      end do
    end do
  end do

end if

! write section

length = dima*dimbc
if (length > 0) then
  iaddr = daddr(lunt3)
  call ddafile(lunt3,1,r(1,1),length,iaddr)
end if

return

end subroutine t3intpck1
