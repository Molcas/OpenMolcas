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

subroutine gather_t2_blocked(length1,length2,ngaf,ngal,ngbf,ngbl,t2,t2_tmp,tmp)
! length1 = length of the 1st VO index (=< vblock)
! length2 = length of the 2nd VO index (=< vblock)

use ChT3_global, only: DimGrpaR, LunAux, maxdim, no, T2Name
use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: length1, length2, ngaf, ngal, ngbf, ngbl
real(kind=wp), intent(inout) :: t2(length1*length2*no*no)
real(kind=wp), intent(out) :: t2_tmp(maxdim*maxdim*no*no), tmp(maxdim*maxdim*no*no)
integer(kind=iwp) :: a, a_tmp, b, b_tmp, dima, dimb, lasta, lastb, length
logical(kind=iwp) :: sym

sym = .false.

!mp if (ngaf == ngbf) sym = .true.
if ((ngaf == ngbf) .and. (ngal == ngbl)) sym = .true.

!mp write(u6,*)
!mp write(u6,*) '------ DimGrpaR ------'
!mp write(u6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
!mp write(u6,*)

do a=ngaf,ngal
  do b=ngbf,min(a,ngbl)

    dima = DimGrpaR(a)
    dimb = DimGrpaR(b)

    if (a == b) then  ! a=b
      length = nTri_Elem(dima)*no*no
      call GetX_t3(tmp,length,LunAux,T2Name(a,b),1,1)
    else ! a>b
      length = dima*dimb*no*no
      call GetX_t3(tmp,length,LunAux,T2Name(a,b),1,1)
    end if

    lasta = 0
    if (a > ngaf) then
      do a_tmp=ngaf,a-1
        lasta = lasta+DimGrpaR(a_tmp)
      end do
    end if

    lastb = 0
    if (b > ngbf) then
      do b_tmp=ngbf,b-1
        lastb = lastb+DimGrpaR(b_tmp)
      end do
    end if

    !mp write(u6,'(A,2(i5,2x),2(i3,x))') 'lasta, lastb, a, b = ',lasta,lastb,a,b

    if (a == b) then ! expand and map
      call expand4_12(tmp,t2_tmp,dima,no,no)
      call grow_t2_blocked(t2,t2_tmp,dima,dimb,no,lasta,lastb,length1,length2,sym)
    else
      call grow_t2_blocked(t2,tmp,dima,dimb,no,lasta,lastb,length1,length2,sym)
    end if

  end do
end do

return

end subroutine gather_t2_blocked
