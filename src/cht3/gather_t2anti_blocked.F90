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

subroutine gather_t2anti_blocked(length1,length2,ngbf,ngbl,t2,t2_tmp,tmp)
! length1 = length of the 1st VO index (nv)
! length2 = length of the 2nd VO index (=< vblock)
!
! This routine generates T2 amplitudes in this form :
!
! T2 = t2(a,b,j<i) - t2(b,a,j<i) finaly : <a,b,(i<j) >  ;  a in nv, b in vblock

use ChT3_global, only: DimGrpaR, LunAux, no, NvGrp, T2Name
use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: length1, length2, ngbf, ngbl
real(kind=wp), intent(inout) :: t2(*)
real(kind=wp), intent(_OUT_) :: t2_tmp(*), tmp(*)
integer(kind=iwp) :: a, a_tmp, aa, b, b_tmp, bb, dima, dimb, lasta, lastb, length
logical(kind=iwp) :: switch

!mp write(u6,*)
!mp write(u6,*) '------ DimGrpaR ------'
!mp write(u6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
!mp write(u6,*)

do a=1,NvGrp
  do b=ngbf,ngbl

    switch = .false.
    if (a >= b) then
      aa = a
      bb = b
    else
      aa = b
      bb = a
      switch = .true.
    end if

    dima = DimGrpaR(aa)
    dimb = DimGrpaR(bb)

    if (aa == bb) then  ! aa=bb
      length = nTri_Elem(dima)*no*no
      call GetX_t3(tmp,length,LunAux,T2Name(aa,bb),1,1)
    else ! aa>bb
      length = dima*dimb*no*no
      call GetX_t3(tmp,length,LunAux,T2Name(aa,bb),1,1)
    end if

    lasta = 0
    if (a > 1) then
      do a_tmp=1,a-1
        lasta = lasta+DimGrpaR(a_tmp)
      end do
    end if

    lastb = 0
    if (b > ngbf) then
      do b_tmp=ngbf,b-1
        lastb = lastb+DimGrpaR(b_tmp)
      end do
    end if

    !mp write(u6,'(A,2(i5,2x),L,2(i3,x))') 'lasta, lastb, switch, a, b = ',lasta,lastb,switch,a,b

    if (aa == bb) then ! expand tmp
      call expand4_12(tmp,t2_tmp,dima,no,no)
      call grow_t2anti_blocked1(t2,t2_tmp,dima,dimb,no,lasta,lastb,length1,length2)
    else
      if (.not. switch) then
        call grow_t2anti_blocked1(t2,tmp,dima,dimb,no,lasta,lastb,length1,length2)
      else
        call grow_t2anti_blocked2(t2,tmp,dima,dimb,no,lasta,lastb,length1,length2)
      end if
    end if

  end do
end do

return

end subroutine gather_t2anti_blocked
