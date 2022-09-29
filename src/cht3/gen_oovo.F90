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

subroutine gen_oovo(w,l0,l1,tmp)
! this routine genetates (ij,a,k) integrals from
! blocked MO cholesky vectors
!
! --------
!
!       L0(m,IJ)    L0vctr  I>=J
!       L1(m,I ,A') L1vcxx xx - Group of A'

use ChT3_global, only: DimGrpaR, L1Name, LunAux, nc, no, nv, NvGrp
use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: w(*), l0(*), tmp(*)
real(kind=wp), intent(inout) :: l1(*)
integer(kind=iwp) :: a, a_tmp, dima, last, length

!1 read tmp(m,IJ)

length = nc*nTri_Elem(no)

!mp !write(u6,'(A,A6)') 'L0vcrt ','L0vcrt'
!mp !write(u6,*) 'length = ',length
!mp !write(u6,*) 'file size (ifort) = ',8+8*length

call GetX_t3(tmp,length,LunAux,'L0vctr',1,1)

!2 map l0(IJ,m)   <- tmp(m,IJ)

call Map2_21_t3(tmp,l0,nc,nTri_Elem(no))

!3 loop over A'

do a=1,NvGrp
  !mp@@ dima = nv/NvGrp
  dima = DimGrpaR(a)

  !4 read tmp(m,I,A')

  !mp !write(u6,'(A,i3,2x,A6)') 'a,L1Name(a) ',a,L1Name(a)

  !mp@@ if (a == NvGrp) dima = nv-((NvGrp-1)*dima)

  !mp !write(u6,*) 'dima = ',dima
  length = nc*no*dima
  !mp !write(u6,*) 'length = ',length
  !mp !write(u6,*) 'file size (ifort) = ',8+8*length

  call GetX_t3(tmp,length,LunAux,L1Name(a),1,1)

  !5 grow l1(m,I,A)

  !mp last = (a-1)*(nv/NvGrp)

  last = 0
  if (a > 1) then
    do a_tmp=1,a-1
      last = last+DimGrpaR(a_tmp)
    end do
  end if

  call grow_l1(l1,tmp,dima,nc,no,nv,last)

  !6 end loop over A'

end do

!7 map tmp(m,A,I) <- l1(m,I,A)

call Map3_132_t3(l1,tmp,nc,no,nv)

!7.1 zero w

w(1:nTri_Elem(no)*nv*no) = Zero

!8  mult w(IJ,A,I)  <- l0(IJ,m) . tmp(m,A,I)

call mc0c1a3b(nTri_Elem(no),nc,nc,nv*no,nTri_Elem(no),nv*no,nTri_Elem(no),nc,nv*no,l0,tmp,w)

return

end subroutine gen_oovo
