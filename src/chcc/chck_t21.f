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
        subroutine Chck_T21 (T21,beSGrp,gaSGrp)
c
c        test T2n+
c
        implicit none
#include "chcc1.fh"
        integer beSGrp,gaSGrp
        real*8 T21(1:16*31,1:no*(no-1)/2)
        integer a,b,u,v,be,ga,bega,uv,bad,gap,bep
        real*8 s
c
        if (beSGrp.eq.2) then
          bep=nv/2
        else
          bep=0
        end if
c
        if (gaSGrp.eq.2) then
          gap=nv/2
        else
          gap=0
        end if

        bad=0
c
        uv=0
        do u=2,no
        do v=1,u-1
        uv=uv+1
c
          bega=0
          do be=2,nv/2
          do ga=1,be-1
          bega=bega+1
c
            s=0.0d0
            do a=1,nv
            b=a
            s=s+(Q4(b,gap+ga,a,bep+be)+Q4(b,bep+be,a,gap+ga))*
     c          (T2c(b,a,v,u)+T2c(b,a,u,v))/4
            end do
c
            s=0.0d0
            do a=2,nv
            do b=1,a-1
            s=s+(Q4(b,gap+ga,a,bep+be)-Q4(b,bep+be,a,gap+ga))*
     c          (T2c(b,a,v,u)-T2c(b,a,u,v))/2
            end do
            end do
c
          if (abs(T21(bega,uv)-s).gt.1.0d-10) then
            bad=bad+1
c        write (6,99) be,ga,u,v
c99        format (4(i3,1x))
          end if
          T21(bega,uv)=s
c
          end do
          end do
c
        end do
        end do
c
        if (bad.eq.0) then
        write (6,*) ' Chck T2 OK ', bad
        else
        write (6,*) ' Chck T2 Bug !!!!!!! ', bad
        end if
c
        return
        end
