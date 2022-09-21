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

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: w(*), l0(*), l1(*), tmp(*)
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
integer(kind=iwp) :: a, a_tmp, dima, last, length

!1 read tmp(m,IJ)

length = nc*(no*(no+1))/2

!mp !write(u6,'(A,A6)') 'L0vcrt ','L0vcrt'
!mp !write(u6,*) 'length = ',length
!mp !write(u6,*) 'file size (ifort) = ',8+8*length

call GetX_t3(tmp,length,LunAux,'L0vctr',1,1)

!2 map l0(IJ,m)   <- tmp(m,IJ)

call Map2_21_t3(tmp,l0,nc,(no*(no+1)/2))

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

call zeroma(w,1,((no*(no+1))/2)*nv*no)

!8  mult w(IJ,A,I)  <- l0(IJ,m) . tmp(m,A,I)

call mc0c1a3b((no*(no+1))/2,nc,nc,nv*no,(no*(no+1))/2,nv*no,(no*(no+1))/2,nc,nv*no,l0,tmp,w)

return

end subroutine gen_oovo
