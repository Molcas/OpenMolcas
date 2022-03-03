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

subroutine t3intpck2(vint,r,dimv1,dimv2,dimv3,dima,dimb,dimc,symq,symr,syms,nob,nvb)
! this routine packs integral block symi,symq,symr,syms
! R_i(a,b,c) = V_i(b,a,c)
! for symq(b)>syms(c)
! and writes R block onto proper place of open t3nam file - lunt3
!
! vint  - integrals for given symmetries for given i (I)
! r     - final R_i matrix (O)
! dimv1 - 1st dimension of V (I)
! dimv2 - 2nd dimension of V (I)
! dimv3 - 3rd dimension of V (I)
! dima  - dimension of a in R (I)
! dimb  - dimension of b in R (I)
! dimc  - dimension of c in R (I)
! symq  - symmetry of q (b) (I)
! symr  - symmetry of r (a) (I)
! syms  - symmetry of s (c) (I)
! nob   - number of beta occupied in each irrep (I)
! nvb   - number of beta virtuals in each irrep (I)

use ccsort_global, only: daddr, lunt3
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimv1, dimv2, dimv3, dima, dimb, dimc, symq, symr, syms, nob(8), nvb(8)
real(kind=wp), intent(in) :: vint(dimv1,dimv2,dimv3)
real(kind=wp), intent(out) :: r(dima,dimb,dimc)
integer(kind=iwp) :: a, adda, addb, addc, b, c, iaddr, length

! if there are no beta virtuals - skip
if (nvb(symq)*nvb(symr)*nvb(syms) == 0) return

! calc additional constants for a,b,c
adda = nob(symr)
addb = nob(symq)
addc = nob(syms)

! do packing

do c=1,nvb(syms)
  do b=1,nvb(symq)
    do a=1,nvb(symr)
      r(a,b,c) = vint(b+addb,a+adda,c+addc)
    end do
  end do
end do

! write section

length = dima*dimb*dimc
if (length > 0) then
  iaddr = daddr(lunt3)
  call ddafile(lunt3,1,r(1,1,1),length,iaddr)
end if

return

end subroutine t3intpck2
