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
c        files checkeroo
c        MakeChckData
c          MkL0
c          MkL2_chcc
c          MkQ0
c          MkQ1
c          MkQ22
c          MkQ4
c          MkQ3
c          MkOE
c          MkT1T2
c        SaveChckData
c        GetChckData
c        DistMemChck
c
c        ------------------------
c
        subroutine MakeChckData (wrk,wrksize,LunAux)
c
c        this routine generate checkeroo data
c        Use id possible only when NaGrp=NaSGRp=1
c
c        assumption nc>=nv>=no (inac preverit dimenzovanie)
c        dimensions: V1    - no*nv*nv2, nv2*nv2
c                    V2    - nc*nv2
c                    V3    - nc*nv*no
c
        implicit none
#include "wrk.fh"
#include "chcc1.fh"
#include "chcc_files.fh"
c
        integer LunAux
c
c        help variables
        integer dim1
        character*6 LunName
        integer PossV1,PossV2,PossV3,PossT
c
        call DistMemChck (PossV1,PossV2,PossV3,PossT)
c
c1        make Q0
        LunName=I0name
        dim1=no*(no+1)*no*(no+1)/4
        call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
        call MkQ0 (wrk(PossV1))
c
c2        make Q1
        LunName=I1name(1)
        dim1=no*(no+1)*no*nv/2
        call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
        call MkQ1 (wrk(PossV1))
c
c3        Read Q21
        LunName=I2name(1,1)
        dim1=no*no*nv*nv
        call GetX (Q21(1,1,1,1),dim1,LunAux,LunName,1,1)
c
c
c4        make Q22
        LunName=I3name(1,1)
        dim1=no*(no+1)*nv*(nv+1)/4
        call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
        call MkQ22 (wrk(PossV1))
c
c5        make Q4, L2k
        LunName=L2name(1,1)
        dim1=nc*nv*(nv+1)/2
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
         call MkL2_chcc (wrk(PossV2))
        dim1=nv*(nv+1)*nv*(nv+1)/4
        call mv0zero (dim1,dim1,wrk(PossV1))
        dim1=nv*(nv+1)/2
        call mc0c1at3b (nc,dim1,nc,dim1,dim1,dim1,
     c                 dim1,nc,dim1,
     c                 wrk(PossV2),wrk(PossV2),wrk(PossV1))
        call MkQ4 (wrk(PossV1))

c
c        make Q3, L1k
        LunName=L1name(1)
        dim1=nc*no*nv
        call GetX (wrk(PossV3),dim1,LunAux,LunName,1,1)
         call mv0u (dim1,wrk(PossV3),1,L1k(1,1,1),1)
        dim1=no*nv*nv*(nv+1)/2
        call mv0zero (dim1,dim1,wrk(PossV1))
        dim1=nv*(nv+1)/2
        call mc0c1at3b (nc,dim1,nc,no*nv,dim1,no*nv,dim1,nc,no*nv,
     c                 wrk(PossV2),wrk(PossV3),wrk(PossV1))
        call MkQ3 (wrk(PossV1))
c
c
c        make L0
         LunName=L0name
         dim1=nc*no*(no+1)/2
         call GetX (wrk(PossV3),dim1,LunAux,LunName,1,1)
         call MkL0 (wrk(PossV3))
c
c
c        make OEo,OEv
         call MkOE (wrk(PossOE))
c
        call MkT1T2
c
        return
        end
c
c        --------------------
c
        subroutine MkL0 (V)
c
c        L0(m,i,j) <- V(m,ij)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nc,1:no*(no+1)/2)
c
c        help var
        integer i,j,ij,m
c
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
          do m=1,nc
             L0k(m,i,j)=V(m,ij)
             L0k(m,j,i)=V(m,ij)
          end do
        end do
        end do
c
        return
        end
c
c        --------------------
c
        subroutine MkL2_chcc (V)
c
c        L2(m,a,b) <- V(m,ab)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nc,1:nv*(nv+1)/2)
c
c        help var
        integer a,b,ab,m
c
        ab=0
        do a=1,nv
        do b=1,a
        ab=ab+1
          do m=1,nc
             L2k(m,a,b)=V(m,ab)
             L2k(m,b,a)=V(m,ab)
          end do
        end do
        end do
c
        return
        end
c
c        --------------------
c
        subroutine MkQ0 (V)
c
c        Q0(i,j,k,l) <- V(ij,kl)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:no*(no+1)/2,1:no*(no+1)/2)
c
c        help variables
        integer i,j,k,l,ij,kl
c
        kl=0
        do k=1,no
        do l=1,k
        kl=kl+1
          ij=0
          do i=1,no
          do j=1,i
          ij=ij+1
            Q0(i,j,k,l)=V(ij,kl)
            Q0(i,j,l,k)=V(ij,kl)
            Q0(j,i,k,l)=V(ij,kl)
            Q0(j,i,l,k)=V(ij,kl)
          end do
          end do
        end do
        end do
c
        return
        end
c
c        --------------------
c
        subroutine MkQ1 (V)
c
c        Q1(a,j,k,l) <- V(aj,kl)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nv,1:no,1:no*(no+1)/2)
c
c        help variables
        integer a,j,k,l,kl
c
        kl=0
        do k=1,no
        do l=1,k
        kl=kl+1
          do j=1,no
          do a=1,nv
            Q1(a,j,k,l)=V(a,j,kl)
            Q1(a,j,l,k)=V(a,j,kl)
          end do
          end do
        end do
        end do
c
        return
        end
c
c        --------------------
c
        subroutine MkQ22 (V)
c
c        Q22(a,b,k,l) <- V(ij,kl)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nv*(nv+1)/2,1:no*(no+1)/2)
c
c        help variables
        integer a,b,k,l,ab,kl
c
        kl=0
        do k=1,no
        do l=1,k
        kl=kl+1
          ab=0
          do a=1,nv
          do b=1,a
          ab=ab+1
            Q22(a,b,k,l)=V(ab,kl)
            Q22(a,b,l,k)=V(ab,kl)
            Q22(b,a,k,l)=V(ab,kl)
            Q22(b,a,l,k)=V(ab,kl)
          end do
          end do
        end do
        end do
c
        return
        end
c
c        --------------------
c
        subroutine MkQ4 (V)
c
c        Q4(a,b,c,d) <- V(ab,cd)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nv*(nv+1)/2,1:nv*(nv+1)/2)
c
c        help variables
        integer a,b,c,d,ab,cd
c
        cd=0
        do c=1,nv
        do d=1,c
        cd=cd+1
          ab=0
          do a=1,nv
          do b=1,a
          ab=ab+1
            Q4(a,b,c,d)=V(ab,cd)
            Q4(a,b,d,c)=V(ab,cd)
            Q4(b,a,c,d)=V(ab,cd)
            Q4(b,a,d,c)=V(ab,cd)
          end do
          end do
        end do
        end do
c
        return
        end
c
c        --------------------
c
        subroutine MkQ3 (V)
c
c        Q1(a,b,c,l) <- V(ab,cl)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nv*(nv+1)/2,1:no,1:nv)
c
c        help variables
        integer a,b,ab,c,l
c
        do l=1,no
        do c=1,nv
          ab=0
          do a=1,nv
          do b=1,a
          ab=ab+1
            Q3(a,b,c,l)=V(ab,l,c)
            Q3(b,a,c,l)=V(ab,l,c)
          end do
          end do
        end do
        end do
c
        return
        end
c
c        --------------------
c
        subroutine MkOE (OE)
c
c        OEo(i) <- OE(i)
c        OEv(a) <- OE(a)
c
        implicit none
#include "chcc1.fh"
        real*8 OE(1:(nv+no))
c
c        help variables
        integer a,i
c
        do i=1,no
          OEo(i)=OE(i)
        end do
c
        do a=1,nv
          OEv(a)=OE(no+a)
        end do
c
        return
        end
c
c        --------------------
c
        subroutine MkT1T2
c
c        T1(a,i) = 0 (mozno neskor ine)
c        T2(a,b,i,j) = (ai|bj)/Dabij
c
        implicit none
#include "chcc1.fh"
c
c        help variables
        integer a,b,i,j
c
        do i=1,no
        do a=1,nv
          T1c(a,i)=0.0d0
        end do
        end do
c
        do j=1,no
        do i=1,no
          do b=1,nv
          do a=1,nv
            T2c(a,b,i,j)=Q21(a,i,b,j)/(OEo(i)+OEo(j)-OEv(a)-OEv(b))
          end do
          end do
        end do
        end do
c
        return
        end
c
c        --------------------
c
        subroutine SaveChckData (LunAux)
c
c        this routine do:
c        Save Chck Data to ChkDat file
c
        implicit none
#include "chcc1.fh"
        integer LunAux
c
C       open (unit=LunAux,file='ChKDat',form='unformatted')
        Call Molcas_BinaryOpen_Vanilla(LunAux,'ChKDat')
          write(LunAux) T1c,T2c,OEo,OEv,Q0,Q1,Q21,Q22,Q3,Q4
     c                 ,L0k,L1k,L2k
        close (LunAux)
c
        return
        end
c
c        --------------------
c
        subroutine GetChckData (LunAux)
c
c        this routine do:
c        Get Chck Data from ChkDat file
c
        implicit none
#include "chcc1.fh"
        integer LunAux
c
*       open (unit=LunAux,file='ChKDat',form='unformatted')
        Call Molcas_BinaryOpen_Vanilla(LunAux,'ChKDat')
          read(LunAux) T1c,T2c,OEo,OEv,Q0,Q1,Q21,Q22,Q3,Q4
     c                 ,L0k,L1k,L2k
        close (LunAux)
c
        return
        end
c
c        --------------------
c
        subroutine  DistMemChck (PossV1,PossV2,PossV3,PossT)
c
c        Dimensions for DistMemCheck
c        dimensions: V1    - no*nv*nv2, nv2*nv2
c                    V2    - nc*nv2
c                    V3    - nc*nv*no

        implicit none
#include "chcc1.fh"
        integer PossV1,PossV2,PossV3,PossT
c
c        help variables
        integer len
c
        PossT=PossFree
c
        PossV1=PossT
        len=no*nv*nv*(nv+1)/2
        if ((nv*(nv+1)*nv*(nv+1)/4).gt.len) then
        len=nv*(nv+1)*nv*(nv+1)/4
        end if
        PossT=PossT+len
c
        PossV2=PossT
        len=nc*nv*(nv+1)/2
        PossT=PossT+len
c
        PossV3=PossT
        len=nc*no*nv
c
        PossT=PossT+len
c
        write (6,*) ' Poss ChCk',PossV1,PossV2,PossV3,PossT
c
        return
        end
