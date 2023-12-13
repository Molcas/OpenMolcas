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

subroutine gen_vvoo(w,l1,tmp,l2)
! this routine does:
!
! regenerate (ab,ij) integrals from blocked
! MO cholesky vectors
!
! <vv|oo> = (vo|vo)
!
! --------
!
!       L1(m,I,A')

use ChT3_global, only: DimGrpaR, L1Name, LunAux, nc, no, nv, NvGrp
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: w(*)
real(kind=wp), intent(_OUT_) :: l1(*), tmp(*), l2(*)
integer(kind=iwp) :: a, a_tmp, b, b_tmp, dima, dimb, lasta, lastb, length

do a=1,NvGrp

  !1 read tmp(m,I,A')

  dima = DimGrpaR(a)
  length = nc*no*dima
  call GetX_t3(tmp,length,LunAux,L1Name(a),1,1)

  !5 map l1 (A',I,m) <- tmp (m,I,A')

  call Map3_321_t3(tmp,l1,nc,no,dima)

  ! ----- read tmp(m,I,B')

  do b=1,a

    dimb = DimGrpaR(b)
    length = nc*no*dimb
    call GetX_t3(tmp,length,LunAux,L1Name(b),1,1)

    !4 map l2 (m,B',I) <- tmp (m,I,B')

    call Map3_132_t3(tmp,l2,nc,no,dimb)

    ! zero tmp

    tmp(1:dima*no*dimb*no) = Zero

    !7 mult tmp(A',I,B',J) <- l1 (A',I,m) . l2(m,B',J)

    call mc0c1a3b(dima*no,nc,nc,dimb*no,dima*no,dimb*no,dima*no,nc,dimb*no,l1,l2,tmp)

    lasta = 0
    if (a > 1) then
      do a_tmp=1,a-1
        lasta = lasta+DimGrpaR(a_tmp)
      end do
    end if

    lastb = 0
    if (b > 1) then
      do b_tmp=1,b-1
        lastb = lastb+DimGrpaR(b_tmp)
      end do
    end if

    !write(u6,'(A,4(i4,x))') 'BB1 dima, dimb, lasta, lastb ',dima,dimb,lasta,lastb
    call grow_vvoo(w,tmp,no,nv,dima,dimb,lasta,lastb)

    if (a /= b) then

      call Map4_3412_t3(tmp,l2,dima,no,dimb,no)

      !write(u6,'(A,4(i4,x))') 'BB2 dima, dimb, lasta, lastb ',dimb,dima,lastb,lasta
      call grow_vvoo(w,l2,no,nv,dimb,dima,lastb,lasta)

    end if

    !3 end loop over B'

  end do

  !3 end loop over A'

end do

return

end subroutine gen_vvoo
