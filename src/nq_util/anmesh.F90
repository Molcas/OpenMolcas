!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2001, Roland Lindh                                     *
!               2001, Laura Gagliardi                                  *
!***********************************************************************
      Subroutine AnMesh(  nscheme,    pa,       rPt,   wPt)
      use nq_Info
      Implicit Real*8 (a-h,o-z)
      Implicit Integer(i-n)
#include "real.fh"
#include "debug.fh"
#include "WrkSpc.fh"
      Integer nscheme(8)
      Real*8  pa(*),    rPt(3,*),    wPt(*)
!***********************************************************************
!                                                                      *
!
         If (Debug) Then
        Write (6,*)
        Write (6,*) ' ******** The Angular Lebedev Grid ********'
        Write (6,*)
        Write (6,*)
         End If
!
      i = 0
      ip = 0
!
!     nscheme(2) -> 6 points
!
!lg      write (*,*) 'nscheme',  (nscheme(i), i=1,8)
      If(nscheme(2) .gt. 0) Then
!lg          write (*,*) 'nscheme(2)', nscheme(2)
         ip = ip + 1
         Do ix=1,3
            Do iy=1,-1,-2
               i = i + 1
               wPt(i) = pa(ip)
               Do j=1,3
                  rPt(j,i) = Zero
               Enddo
               rPt(ix,i) = DBLE(iy)
!lg              write (*,*) rPt(ix,i), wPt(i)
            Enddo
         Enddo
      Endif
!
!     nscheme(3) -> 8 points
!
      If(nscheme(3) .gt. 0) Then   !
         c = One/sqrt(Three)
         ip = ip + 1
         Do ix=1,-1,-2
            Do iy=1,-1,-2
               Do iz=1,-1,-2
                  i = i + 1
                  wPt(i) = pa(ip)
                  rPt(1,i) = DBLE(ix)*c
                  rPt(2,i) = DBLE(iy)*c
                  rPt(3,i) = DBLE(iz)*c
               Enddo
            Enddo
         Enddo
      Endif
!
!      nscheme(4) -> 12 points
!
      If(nscheme(4) .gt. 0) Then
         c = One/sqrt(Two)
         ip = ip + 1
         Do ix=1,-1,-2
            Do iy=1,-1,-2
               Do iz=1,3
                  i = i + 1
                  wPt(i) = pa(ip)
                  rPt(iz,i) = DBLE(ix)*c
                  j = mod(iz,3) + 1
                  rPt(j,i) = DBLE(iy)*c
                  j = 6 - iz - j
                  rPt(j,i) = Zero
               Enddo
            Enddo
         Enddo
      endif
!
!     24a points
!
      n1 = nscheme(5)
      Do jj=1,n1
         ip = ip + 1
         uu = pa(ip)
         vv = sqrt(One - Two*uu*uu)
         ip = ip + 1
         Do ix=1,-1,-2
            Do iy=1,-1,-2
               Do iz=1,-1,-2
                  Do j=1,3
                     i = i + 1
                     wPt(i) = pa(ip)
                     Do j1=1,3
                        rPt(j1,i) = uu
                     Enddo
                     rPt(j,i) = vv
                     rPt(1,i) = rPt(1,i)*DBLE(ix)
                     rPt(2,i) = rPt(2,i)*DBLE(iy)
                     rPt(3,i) = rPt(3,i)*DBLE(iz)
                  Enddo
               Enddo
            Enddo
         Enddo
      Enddo
!
!     24b points
!
      n1 = nscheme(6)
      Do jj=1,n1
         ip = ip + 1
         pp = pa(ip)
         qq = sqrt(One - pp*pp)
         ip = ip + 1
         Do ix=1,-1,-2
            Do iy=1,-1,-2
               Do ii=0,1
                  Do j=1,3
                     i = i + 1
                     wPt(i) = pa(ip)
                     j1 = mod(j+ii,3)+1
                     rPt(j1,i) = pp*DBLE(ix)
                     j1 = mod(j+1-ii,3) + 1
                     rPt(j1,i) = qq*DBLE(iy)
                     rPt(j,i) = zero
                  Enddo
               Enddo
            Enddo
         Enddo
      Enddo
!
!     48 points
!
      n1 = nscheme(7)
!lg      write (*,*) 'i, n1 =', i, n1
      Do jj=1,n1
         ip = ip + 1
         rr = pa(ip)
         ip = ip + 1
         ss = pa(ip)
         tt = sqrt(One - rr*rr - ss*ss)
         ip = ip + 1
         Do ix=1,-1,-2
            Do iy=1,-1,-2
               Do iz=1,-1,-2
                  Do j=1,3
                     Do ii=0,1
                        i = i + 1
                        wPt(i) = pa(ip)
                        rPt(j,i) = rr*DBLE(ix)
                        j1 = mod(j+ii,3) + 1
                        rPt(j1,i) = ss*DBLE(iy)
                        j1 = mod(j+1-ii,3) + 1
                        rPt(j1,i) = tt*DBLE(iz)
!lg                        write (*,*) rPt(j1,i), wPt(i),j1,i
                     Enddo
!lg      write (*,*) 'Enddo1', i
                  Enddo
!lg      write (*,*) 'Enddo2', i
               Enddo
!lg      write (*,*) 'Enddo3', i
            Enddo
!lg      write (*,*) 'Enddo4', i
         Enddo
!lg      write (*,*) 'Enddo5', i
      Enddo
!lg                        write (*,*) 'enddo', n1
!                                                                      *
!***********************************************************************
!                                                                      *
!
!lg       write (*,*) 'End of AnMesh'
      Return
      End
