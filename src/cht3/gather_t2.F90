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

subroutine gather_t2(t2,t2_tmp,tmp)
! temporary routine. In future T2 block structure will be merged
!                    with the block structure of the (T) code
!
! this routine does:
!
! cycle through T2XY files and gather them into on T2 array
! T2(nv_beta,nv_alpha,no_beta,no_alpha)

use ChT3_global, only: DimGrpaR, LunAux, no, nv, NvGrp, T2Name
use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: t2(*)
real(kind=wp), intent(_OUT_) :: t2_tmp(*), tmp(*)
integer(kind=iwp) :: a, a_tmp, b, b_tmp, dima, dimb, lasta, lastb, length

!mp write(u6,*)
!mp write(u6,*) '------ DimGrpaR ------'
!mp write(u6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
!mp write(u6,*)

do a=1,NvGrp
  do b=1,a

    !mp@@ dima = nv/NvGrp
    dima = DimGrpaR(a)
    dimb = DimGrpaR(b)

    !mp write(u6,'(A,i3,i3,2x,A6)') 'a,b,T2Name(a,b) ',a,b,T2Name(a,b)

    if (a == b) then  ! a=b
      ! open the pertinent file

      !mp@@ if (a == NvGrp) dima = nv-((NvGrp-1)*dima)
      dimb = dima

      !mp write(u6,*) 'dima = ',dima
      length = nTri_Elem(dima)*no*no
      !mp write(u6,*) 'length = ',length
      !mp !write(u6,*) 'file size (g77) = ',16+length*8
      !mp write(u6,*) 'file size (ifort) = ',8+length*8

      call GetX_t3(tmp,length,LunAux,T2Name(a,b),1,1)

    else ! a>b
      ! open the pertinent file

      !mp@@ if (a == NvGrp) dima = nv-((NvGrp-1)*dima)
      !mp@@ if (b == NvGrp) dimb = nv-((NvGrp-1)*dimb)

      !mp write(u6,*) 'dima = ',dima
      !mp write(u6,*) 'dimb = ',dimb
      length = dima*dimb*no*no
      !mp write(u6,*) 'length = ',length
      !mp !write(u6,*) 'file size (g77) = ',16+length*8
      !mp write(u6,*) 'file size (ifort) = ',8+length*8

      call GetX_t3(tmp,length,LunAux,T2Name(a,b),1,1)

    end if

    ! add its contents to the t2 array

    !mp@@ lasta = (a-1)*(nv/NvGrp)
    !mp@@ lastb = (b-1)*(nv/NvGrp)
    !mp@@
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
    !mp@@

    !mp write(u6,'(A,2(i5,2x),2(i3,x))') 'lasta, lastb, a, b = ',lasta,lastb,a,b

    if (a == b) then ! expand and map
      ! expand and map l2_1 (a',b',m) <- tmp (m,ab')
      call expand4_12(tmp,t2_tmp,dima,no,no)
      call grow_t2neq(t2,t2_tmp,dima,dimb,nv,no,lasta,lastb)
    else
      call grow_t2neq(t2,tmp,dima,dimb,nv,no,lasta,lastb)
    end if

  end do
end do

return

end subroutine gather_t2
