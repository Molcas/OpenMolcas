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
c        MkT_CAlld
c        MkT_C136od
c        MkT_C245od
c        MkT_C78d
c        MkT_C78od
c        Energy_E2od
c        MkT_red
c        MkT_exp
c        Energy_E1
c        DistMemSum
c
c        -----------------------
c
        subroutine MkT_CAlld (T2,X,Y,dimbe,no)
c
c        this routine do:
c       T2n(be',ga',u,v) <-
c           C1                + 1/2 X(be',u,ga',v)
c           C2                + 1/2 X(ga',v,be',u)
c           C3                - 1/2 Y(be',u,ga',v)
c           C4                - 1/2 Y(ga',v,be',u)
c           C5                - 1   Y(ga',u,be',v)
c           C6                - 1   Y(be',v,ga',u)
c        for beGrp=gaGrp
c
        implicit none
        integer dimbe,no
        real*8 T2(1:dimbe,1:dimbe,1:no,1:no)
        real*8 X(1:dimbe,1:no,1:dimbe,1:no)
        real*8 Y(1:dimbe,1:no,1:dimbe,1:no)
c
c        help variables
        integer u,v,be,ga
c
        do v=1,no
          do u=1,no
            do ga=1,dimbe
c              do be=1,dimbe  - povodne, stacia iba cleny be>=ga
              do be=ga,dimbe
        T2(be,ga,u,v)=(X(be,u,ga,v)-Y(be,u,ga,v))/2-Y(be,v,ga,u)
              end do
            end do
          end do
        end do
c
        do v=1,no
          do u=1,no
            do be=1,dimbe
c              do ga=1,dimbe  - povodne, stacia iba cleny be>=ga
              do ga=1,be
                 T2(be,ga,u,v)=T2(be,ga,u,v)
     c                       +(X(ga,v,be,u)-Y(ga,v,be,u))/2
     c                       -Y(ga,u,be,v)
              end do
            end do
          end do
        end do
c
        return
        end
c
c        -----------------------
c
        subroutine MkT_C136od (T2,X,Y,dimbe,dimga,no)
c
c        this routine do:
c       T2n(be',ga',u,v) <-
c        C1                + 1/2 X(be',u,ga',v)
c        C3                - 1/2 Y(be',u,ga',v)
c        C6                - 1   Y(be',v,ga',u)
c        for beGrp>gaGrp
c
        implicit none
        integer dimbe,dimga,no
        real*8 T2(1:dimbe,1:dimga,1:no,1:no)
        real*8 X(1:dimbe,1:no,1:dimga,1:no)
        real*8 Y(1:dimbe,1:no,1:dimga,1:no)
c
c        help variables
        integer u,v,be,ga
c
        do v=1,no
          do u=1,no
            do ga=1,dimga
              do be=1,dimbe
        T2(be,ga,u,v)=(X(be,u,ga,v)-Y(be,u,ga,v))/2-Y(be,v,ga,u)
              end do
            end do
          end do
        end do
c
        return
        end
c
c        -----------------------
c
        subroutine MkT_C245od (T2,X,Y,dimbe,dimga,no)
c
c        this routine do:
c       T2n(be',ga',u,v) <<-
c        C2                + 1/2 X(ga',v,be',u)
c        C4                - 1/2 Y(ga',v,be',u)
c        C5                - 1   Y(ga',u,be',v)
c        for beGrp>gaGrp
c
        implicit none
        integer dimbe,dimga,no
        real*8 T2(1:dimbe,1:dimga,1:no,1:no)
        real*8 X(1:dimga,1:no,1:dimbe,1:no)
        real*8 Y(1:dimga,1:no,1:dimbe,1:no)
c
c        help variables
        integer u,v,be,ga
c
        do v=1,no
          do u=1,no
            do be=1,dimbe
              do ga=1,dimga
                 T2(be,ga,u,v)=T2(be,ga,u,v)
     c                       +(X(ga,v,be,u)-Y(ga,v,be,u))/2
     c                       -Y(ga,u,be,v)
              end do
            end do
          end do
        end do
c
        return
        end
c
c        -----------------------
c
        subroutine MkT_C78d (T2,Tp,Tm,dimbe,dimbepp,addbepp,no)
c
c        this routine do:
c       T2n(be',ga',u,v) <<-
c        C7                    T2+(bega",uv)
c        C8                    T2-(bega",uv)
c        for beGrp=gaGrp, beSGrp=gaSGrp
c        N.B. calc only contributions to be',ga' (not ga',be')
c
        implicit none
        integer dimbe,dimbepp,addbepp,no
        real*8 T2(1:dimbe,1:dimbe,1:no,1:no)
        real*8 Tp(1:dimbepp*(dimbepp+1)/2,1:no*(no+1)/2)
        real*8 Tm(1:dimbepp*(dimbepp-1)/2,1:no*(no-1)/2)
c
c        help variables
        integer u,v,be,ga,uv,bega,bep,gap
        real*8 fact
c
c
c1        Distribute symmetric T2+ on proper possitions
c
        uv=0
        do u=1,no
        do v=1,u
        uv=uv+1
        if (u.eq.v) then
          fact=0.5d0
        else
          fact=1.0d0
        end if
c
c          case be".ne.ga"
          bep=addbepp+1
          do be=2,dimbepp
          bega=be*(be-1)/2
          bep=bep+1
          gap=addbepp
          do ga=1,be-1
          bega=bega+1
          gap=gap+1
c
            T2(bep,gap,u,v)=T2(bep,gap,u,v)+Tp(bega,uv)*fact
            T2(bep,gap,v,u)=T2(bep,gap,v,u)+Tp(bega,uv)*fact
c             T2(gap,bep,u,v)=T2(gap,bep,u,v)+Tp(bega,uv)*fact
c            T2(gap,bep,v,u)=T2(gap,bep,v,u)+Tp(bega,uv)*fact
c
          end do
          end do
c
c          case be=ga
          bep=addbepp
          do be=1,dimbepp
          bega=be*(be+1)/2
          bep=bep+1
c
            T2(bep,bep,u,v)=T2(bep,bep,u,v)+Tp(bega,uv)*fact
            T2(bep,bep,v,u)=T2(bep,bep,v,u)+Tp(bega,uv)*fact
c
          end do
c
        end do
        end do
c
c
c2        Distribute anti-symmetric T2- on proper possitions
c
        uv=0
        do u=2,no
        do v=1,u-1
        uv=uv+1
c
          bep=addbepp+1
          bega=0
          do be=2,dimbepp
          bep=bep+1
          gap=addbepp
          do ga=1,be-1
          bega=bega+1
          gap=gap+1
c
            T2(bep,gap,u,v)=T2(bep,gap,u,v)+Tm(bega,uv)
            T2(bep,gap,v,u)=T2(bep,gap,v,u)-Tm(bega,uv)
c            T2(gap,bep,u,v)=T2(gap,bep,u,v)-Tm(bega,uv)
c            T2(gap,bep,v,u)=T2(gap,bep,v,u)+Tm(bega,uv)
c
          end do
          end do
c
        end do
        end do
c
c
        return
        end
c
c        -----------------------
c
        subroutine MkT_C78od (T2,Tp,Tm,
     c             dimbe,dimga,dimbepp,dimgapp,addbepp,addgapp,no)
c
c        this routine do:
c       T2n(be',ga',u,v) <<-
c        C7                    T2+(be',ga',uv)
c        C8                    T2-(be',ga',uv)
c        for beSGrp.ne.gaSGrp
c
        implicit none
        integer dimbe,dimga,dimbepp,dimgapp,addbepp,addgapp,no
        real*8 T2(1:dimbe,1:dimga,1:no,1:no)
        real*8 Tp(1:dimbepp,1:dimgapp,1:no*(no+1)/2)
        real*8 Tm(1:dimbepp,1:dimgapp,1:no*(no-1)/2)
c
c        help variables
        integer u,v,be,ga,uvp,uvm,uup,bep,gap
c
c        diagonal part - contributoion from T+ only
        do u=1,no
        uup=u*(u+1)/2
c
c          cycle over be",ga"
          gap=addgapp
          do ga=1,dimgapp
            gap=gap+1
            bep=addbepp
            do be=1,dimbepp
              bep=bep+1
              T2(bep,gap,u,u)=T2(bep,gap,u,u)+Tp(be,ga,uup)
            end do
          end do
c
        end do
c
c
c        off diagonal - both T+,T-
        uvm=0
        do u=2,no
        uvp=u*(u-1)/2
        do v=1,u-1
          uvm=uvm+1
          uvp=uvp+1
c
c          cycle over be",ga"
          gap=addgapp
          do ga=1,dimgapp
            gap=gap+1
            bep=addbepp
            do be=1,dimbepp
              bep=bep+1
              T2(bep,gap,u,v)=T2(bep,gap,u,v)
     c                       +Tp(be,ga,uvp)+Tm(be,ga,uvm)
              T2(bep,gap,v,u)=T2(bep,gap,v,u)
     c                       +Tp(be,ga,uvp)-Tm(be,ga,uvm)
            end do
          end do
c
c          cycle over be",ga"
          gap=addgapp
          do ga=1,dimgapp
            gap=gap+1
            bep=addbepp
            do be=1,dimbepp
              bep=bep+1
c              T2(bep,gap,v,u)=T2(bep,gap,v,u)
c    c                       +Tp(be,ga,uvp)-Tm(be,ga,uvm)
            end do
          end do
c
        end do
        end do
c
c
        return
        end
c
c        -----------------------
c
        subroutine Energy_E2d (V,Tau,e,eos,dima,no)
c
c        this routine calc:
c        E2 = sum(a,b,i,j) (2V(ai|bj)-V(aj|bi)) . Tau(a,b,i,j
c        E2os=sum(a,b,i,j) (  (ai|bj)         ) . Tau(a,b,i,j
c        where
c           E2    - complete E2 component of energy
c           E2os  - other spin E2 component of energy
c        for aGrp=bGrp (Tau array is full, but only a>=b values are completed)
c
        implicit none
        integer dima,no
        real*8 V(1:dima,1:no,1:dima,1:no)
        real*8 Tau(1:dima,1:dima,1:no,1:no)
        real*8 e,eos
c
c        help variables
        integer a,b,i,j
        real*8 ehlp
c
        e=0.0d0
        eos=0.0d0
        ehlp=0.0d0
c
c        off diagonal
c
        do j=1,no
        do i=1,no
        do b=1,dima-1
           ehlp=ehlp+V(b,i,b,j)*Tau(b,b,i,j)
        do a=b+1,dima
           e=e+(2.0d0*V(a,i,b,j)-V(a,j,b,i))*Tau(a,b,i,j)
           eos=eos+V(a,i,b,j)*Tau(a,b,i,j)
        end do
        end do
           ehlp=ehlp+V(dima,i,dima,j)*Tau(dima,dima,i,j)
        end do
        end do
c
        e=2.0d0*e
        eos=2.0d0*eos
c
        e=e+ehlp
        eos=eos+ehlp
c
        return
        end
c
c        -----------------------
c
        subroutine Energy_E2od (V,Tau,e,eos,dima,dimb,no)
c
c        this routine calc:
c        E2 = sum(a,b,i,j) (2V(ai|bj)-V(aj|bi)) . Tau(a,b,i,j)
c        E2os=sum(a,b,i,j)  V(ai|bj)            . Tau(a,b,i,j)
c        where
c           E2    - complete E2 component of energy
c           E2os  - other spin E2 component of energy
c
        implicit none
        integer dima,dimb,no
        real*8 V(1:dima,1:no,1:dimb,1:no)
        real*8 Tau(1:dima,1:dimb,1:no,1:no)
        real*8 e,eos
c
c        help variables
        integer a,b,i,j
c
        e=0.0d0
        eos=0.0d0
c
        do j=1,no
        do i=1,no
        do b=1,dimb
        do a=1,dima
           e=e+(2.0d0*V(a,i,b,j)-V(a,j,b,i))*Tau(a,b,i,j)
           eos=eos+V(a,i,b,j)*Tau(a,b,i,j)
        end do
        end do
        end do
        end do
c
        return
        end
c
c        -----------------------
c
        subroutine MkT_red (T2red,T2full,dimbe,no)
c
c        this routine produce reduced set of ampitudes for storing:
c       T2red(be'ga',u,v) <- T2full(be',ga',u,v)
c        for beGrp=gaGrp
c
        implicit none
        integer dimbe,no
        real*8 T2full(1:dimbe,1:dimbe,1:no,1:no)
        real*8 T2red(1:dimbe*(dimbe+1)/2,1:no,1:no)
c
c        help variables
        integer i,j,be,ga,bega
c
        do j=1,no
        do i=1,no
c
        bega=0
        do be=1,dimbe
        do ga=1,be
        bega=bega+1
c
          T2red(bega,i,j)=T2Full(be,ga,i,j)
c
        end do
        end do
c
        end do
        end do
c
        return
        end
c
c        -----------------------
c
        subroutine MkT_exp (T2red,T2full,dimbe,no)
c
c        this routine produce expanded set of amplitudes from reduced set
c        (used in parallel case)
c       T2red(be'ga',u,v) -> T2full(be',ga',u,v)
c        for beGrp=gaGrp
c
        implicit none
        integer dimbe,no
        real*8 T2full(1:dimbe,1:dimbe,1:no,1:no)
        real*8 T2red(1:dimbe*(dimbe+1)/2,1:no,1:no)
c
c        help variables
        integer i,j,be,ga,bega
c
        do j=1,no
        do i=1,no
c
        bega=0
        do be=1,dimbe
        do ga=1,be
        bega=bega+1
c
cdir          T2red(bega,i,j)=T2Full(be,ga,i,j)
cinv
          T2Full(be,ga,i,j)=T2red(bega,i,j)
          T2Full(ga,be,j,i)=T2red(bega,i,j)
c
        end do
        end do
c
        end do
        end do
c
        return
        end
c
c        ----------------------------
c
        subroutine Energy_E1 (T1n,Fvo,no,nv,E1)
c
c        this routine do:
c        E1 = sum(a,i) T1n(a,i) . Fvo(a,i)
c
c        calculate E1 component of energy
c
        implicit none
        integer no,nv
        real*8 T1n(1)
        real*8 Fvo(1)
        real*8 e1
c
c        help variables
        integer dim1
c
        e1=0.0d0
        dim1=nv*no
        call mr0u3wt (dim1,dim1,dim1,1,1,T1n(1),Fvo(1),e1)
c
        return
        end
c
c        ------------------------
c
        subroutine DistMemSum (NvGrp,maxdim,
     c        PossV1,PossV2,PossV3,
     c        PossH1,PossH2,
     c        PossT)

c
c       This routine do:
c       define initial possitions of V and H
c       used in summary routine
c

        implicit none
#include "chcc1.fh"
c
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3
        integer PossH1,PossH2
        integer PossT
c
c
c       help variables
        integer length
c
c
c1        V files
c
        length=no*no*maxdim*maxdim
        PossV1=possT
        PossT=PossT+length
        PossV2=possT
        PossT=PossT+length
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,99) 'DM V ',PossV1,PossV2,PossV3,length
        end if
c
c
c2        H files
c
        length=no*maxdim
        PossH1=possT
        PossT=PossT+length
        PossH2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,99) 'DM H ',PossH1,PossH2,PossV3,length
        end if
c
c
99      format (a7,10(i10,1x))
c
        if (printkey.ge.10) then
        write (6,99) 'PossT ',PossT
        end if
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(NvGrp)
        end
c
c        ------------------------
c
        subroutine UpG_T1 (T1)
c
c        upgrade T1
c
        implicit none
#include "chcc1.fh"
        real*8 T1(1:nv,1:no)
c
c        help var
        integer i,a
c
        do i=1,no
        do a=1,nv
          T1c(a,i)=T1(a,i)
        end do
        end do
c
        return
        end
c
c        ------------------------
c
        subroutine UpG_T2d (T2,dima,adda)
c
c        upgrade T2 (diagonal)
c
        implicit none
#include "chcc1.fh"
        integer dima,adda
        real*8 T2(1:dima*(dima+1)/2,1:no,1:no)
c
c        help var
        integer i,j,a,b,ab
c
        do j=1,no
        do i=1,no
          ab=0
          do a=1,dima
          do b=1,a
          ab=ab+1
            T2c(a+adda,b+adda,i,j)=T2(ab,i,j)
            T2c(b+adda,a+adda,j,i)=T2(ab,i,j)
          end do
          end do
        end do
        end do
c
        return
        end
c
c        ------------------------
c
        subroutine UpG_T2od (T2,dima,adda,dimb,addb)
c
c        upgrade T2 (off-diagonal)
c
        implicit none
#include "chcc1.fh"
        integer dima,adda,dimb,addb
        real*8 T2(1:dima,1:dimb,1:no,1:no)
c
c        help var
        integer i,j,a,b
c
        do j=1,no
        do i=1,no
          do b=1,dimb
          do a=1,dima
            T2c(a+adda,b+addb,i,j)=T2(a,b,i,j)
            T2c(b+addb,a+adda,j,i)=T2(a,b,i,j)
          end do
          end do
        end do
        end do
c
        return
        end
c
c        ------------------------
c
        subroutine Chck_t2sym
c
c        chek T2c symmetry abij = baji
c
        implicit none
#include "chcc1.fh"
c
        integer i,j,a,b,bad
c
        bad=0
        do j=1,no
        do i=1,no
        do b=1,nv
        do a=1,nv
          if (abs(T2c(a,b,i,j)-T2c(b,a,j,i)).gt.1.0d-10) then
          bad=bad+1
          end if
        end do
        end do
        end do
        end do
c
        write (6,*) ' T2 Symm Check: ',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_Energ
c
c        calc energy from T1c and t2c
c
        implicit none
#include "chcc1.fh"
c
c        help var
        integer i,j,a,b
        real*8 e
c
        e=0.0d0
c
        do j=1,no
        do i=1,no
        do b=1,nv
        do a=1,nv
          e=e+(2.0d0*Q21(a,i,b,j)-Q21(a,j,b,i))*
     c        (T2c(a,b,i,j)+T1c(a,i)*T1c(b,j))
        end do
        end do
        end do
        end do
c
        write (6,*) ' Energia Checkeroo',e
c
        return
        end
