************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine small_svd(m,n,amat,umat,vmat,svals)
      implicit real*8 (a-h,o-z)
      dimension amat(m,*)
      dimension umat(m,*)
      dimension vmat(n,*)
      dimension svals(*)

      if (m.ge.n) then
       call svd_mgen(m,n,amat,umat,vmat,svals)
      else
       call svd_mlen(m,n,amat,umat,vmat,svals)
      end if
      return
      end

      subroutine svd_mgen(m,n,amat,umat,vmat,svals)
      implicit real*8 (a-h,o-z)
      dimension amat(m,n)
      dimension umat(m,n)
      dimension vmat(n,n)
      dimension svals(n)
#include "WrkSpc.fh"

      sq2=sqrt(2.0D0)
      mn=m+n
      mn3=(mn*(mn+1))/2
      mnsq=mn**2
      call getmem('bigA','allo','real',labig,mn3)
      call getmem('bigU','allo','real',lubig,mnsq)
      call dcopy_(mn3,0.0D0,0,work(labig),1)
      call dcopy_(mnsq,0.0D0,0,work(lubig),1)
      call dcopy_(mn,sq2,0,work(lubig),mn+1)
      do i=1,m
       do j=1,n
        work(labig-1+((n+i)*(n+i-1))/2+j)=amat(i,j)
       end do
      end do
      call jacob(work(labig),work(lubig),mn,mn)
      call jacord(work(labig),work(lubig),mn,mn)
      do i=1,n
       irev=mn+1-i
       svals(i)=work(labig-1+(irev*(irev+1))/2)
       do j=1,n
        vmat(j,i)=work(lubig-1+j+mn*(irev-1))
       end do
       do j=1,m
        umat(j,i)=work(lubig-1+n+j+mn*(irev-1))
       end do
      end do
      call getmem('bigA','free','real',labig,mn3)
      call getmem('bigU','free','real',lubig,mnsq)
      return
      end

      subroutine svd_mlen(m,n,amat,umat,vmat,svals)
      implicit real*8 (a-h,o-z)
      dimension amat(m,n)
      dimension umat(m,m)
      dimension vmat(n,m)
      dimension svals(m)
#include "WrkSpc.fh"

      mn=m+n
      mn3=(mn*(mn+1))/2
      mnsq=mn**2
      call getmem('bigA','allo','real',labig,mn3)
      call getmem('bigU','allo','real',lubig,mnsq)
      call dcopy_(mn3,0.0D0,0,work(labig),1)
      call dcopy_(mnsq,0.0D0,0,work(lubig),1)
      call dcopy_(mn,sq2,0,work(lubig),mn+1)
      do i=1,m
       do j=1,n
        work(labig-1+((n+i)*(n+i-1))/2+j)=amat(i,j)
       end do
      end do
      call jacob(work(labig),work(lubig),mn,mn)
      call jacord(work(labig),work(lubig),mn,mn)
      do i=1,m
       irev=mn+1-i
       svals(i)=work(labig-1+(irev*(irev+1))/2)
       do j=1,n
        vmat(j,i)=work(lubig-1+j+mn*(irev-1))
       end do
       do j=1,m
        umat(j,i)=work(lubig-1+n+j+mn*(irev-1))
       end do
      end do
      call getmem('bigA','free','real',labig,mn3)
      call getmem('bigU','free','real',lubig,mnsq)
      return
      end
