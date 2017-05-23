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
c
c        this file contains following routines:
c
c        T1_div
c        T2d_div
c        T2od_div
c        MkTau_chcc
c        ExpT2
c        GetRest
c        SaveRest
c        VanishT1
c
c        ----------------------------
c
        subroutine T1_div (T1n,OE,no,nv)
c
c        this routine do:
c        T1n(a,i) = T1n(a,i)/(e(i)-e(a))
c
c        divison of T1n amplitides by denominator
c
        implicit none
        integer no,nv
        real*8 T1n(1:nv,1:no)
        real*8 OE(1:no+nv)
c
c        help variables
        integer i,a
        real*8 ei
c
        do i=1,no
        ei=OE(i)
          do a=1,nv
            t1n(a,i)=t1n(a,i)/(ei-OE(no+a))
          end do
        end do
c
c
        return
        end
c
c        ----------------------------
c
        subroutine T2d_div (T2,OE,dima,dimb,adda,addb,no,nv)
c
c        this routine do:
c        T2(a',b',i,j) = T2(a',b',i,j/(e(i)+e(j)-e(a)-e(b))
c        for aGrp=beGrp, where only a'>=b',i,j are valid,
c        and completed also cases a'<b',i,j
c
c        divison of T2n amplitides by denominator
c
        implicit none
        integer dima,dimb,adda,addb,no,nv
        real*8 T2(1:dima,1:dimb,1:no,1:no)
        real*8 OE(1:no+nv)
c
c        help variables
        integer i,j,a,av,b,bv
        real*8 eija
c
        av=no+adda
        bv=no+addb
c
c1        divison by denominators
c
        do j=1,no
        do i=1,no
          do a=1,dima
          eija=OE(i)+OE(j)-OE(av+a)
          do b=1,a
            T2(a,b,i,j)=T2(a,b,i,j)/(eija-OE(bv+b))
          end do
        end do
        end do
        end do
c
c2        completing upper triangle
c
        do j=1,no
        do i=1,no
          do b=2,dima
          do a=1,b-1
             T2(a,b,i,j)=T2(b,a,j,i)
          end do
        end do
        end do
        end do
c
c
        return
        end
c
c        ----------------------------
c
        subroutine T2od_div (T2,OE,dima,dimb,adda,addb,no,nv)
c
c        this routine do:
c        T2(a',b',i,j) = T2(a',b',i,j/(e(i)+e(j)-e(a)-e(b))
c
c        divison of T1n amplitides by denominator
c
        implicit none
        integer dima,dimb,adda,addb,no,nv
        real*8 T2(1:dima,1:dimb,1:no,1:no)
        real*8 OE(1:no+nv)
c
c        help variables
        integer i,j,a,av,b,bv
        real*8 eijb
c
        av=no+adda
        bv=no+addb
c
        do j=1,no
        do i=1,no
        do b=1,dimb
          eijb=OE(i)+OE(j)-OE(bv+b)
          do a=1,dima
            T2(a,b,i,j)=T2(a,b,i,j)/(eijb-OE(av+a))
          end do
        end do
        end do
        end do
c
        return
        end
c
c        ----------------------------
c
        subroutine MkTau_chcc (T2,T11,T12,dima,dimb,no,f1,f2)
c
c        this routine do:
c        T2(a',b',i,j) = f1 . T2(a',b',i,j) + f2 . T11(a',i) . T12(b',j)
c
c        N.B. Kvajt odflaknute
c
        implicit none
        integer dima,dimb,no
        real*8 f1,f2
        real*8 T2(1:dima,1:dimb,1:no,1:no)
        real*8 T11(1:dima,1:no)
        real*8 T12(1:dimb,1:no)
c
c        help variables
        integer i,j,a,b
        real*8 c
c
c
        do j=1,no
          do i=1,no
            do b=1,dimb
              c=t12(b,j)*f2
              do a=1,dima
                t2(a,b,i,j)=f1*t2(a,b,i,j)+c*t11(a,i)
              end do
            end do
          end do
        end do
c
c
        return
        end
c
c        ----------------------------
c
        subroutine ExpT2 (T2p,T2u,dima,dimab,no)
c
c        this routine do:
c        Make T2u(a',b',i,j) from T2p(a'b',i,j)
c
c        N.B. Kvajt odflaknute, vypocet ab sa da dat zefektivnit
c
        implicit none
        integer dima,dimab,no
        real*8 T2p(1:dimab,1:no,1:no)
        real*8 T2u(1:dima,1:dima,1:no,1:no)
c
c        help variables
        integer i,j,a,b,ab,ba0
c
c
        do j=1,no
        do i=1,no
          do b=1,dima
            ba0=b*(b-1)/2
            do a=1,b
              T2u(a,b,i,j)=T2p(ba0+a,j,i)
            end do
            do a=1+b,dima
              ab=a*(a-1)/2+b
              T2u(a,b,i,j)=T2p(ab,i,j)
            end do
          end do
        end do
        end do
c
c
        return
        end
c
c        ----------------------------
c
        subroutine GetRest (wrk,wrksize,LunAux,niter,E1old,E2old)
c
c        this file read 1) T1o
c                      2) E1old,E2old,niter
c        from RstFil file
c
        implicit none
#include "chcc1.fh"
#include "wrk.fh"
        integer LunAux,niter
        real*8 E1old,E2old
c
c        help variables
        integer len
c
*       open (unit=LunAux,File='RstFil',form='unformatted')
        Call MOLCAS_BinaryOpen_Vanilla(LunAux,'RstFil')
        len=no*nv
        call rea_chcc (LunAux,len,wrk(PossT1o))
        read (LunAux) E1old,E2old,niter
        close (LunAux)
c
c
        return
        end
c
c        ----------------------------
c
        subroutine SaveRest (wrk,wrksize,LunAux,niter,E1old,E2old)
c
c        this file save 1) T1o,OE
c                      2) E1old,E2old,niter
c        into RstFil file
c
        implicit none
#include "chcc1.fh"
#include "wrk.fh"
        integer LunAux,niter
        real*8 E1old,E2old
c
c        help variables
        integer len
c
*       open (unit=LunAux,File='RstFil',form='unformatted')
        Call MOLCAS_BinaryOpen_Vanilla(LunAux,'RstFil')
        len=no*nv
        call wri_chcc (LunAux,len,wrk(PossT1o))
        write (LunAux) E1old,E2old,niter
        close (LunAux)
c
c
        return
        end
c
c        ----------------------------
c
        subroutine VanishT1 (wrk,wrksize)
c
c        this routine do:
c        Vanish T1o
c
        implicit none
#include "chcc1.fh"
#include "wrk.fh"
c
c        help variables
        integer len
c
        len=no*nv
        call mv0zero (len,len,wrk(PossT1o))
c
        return
        end
