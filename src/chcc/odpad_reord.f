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
c       this file contains:
c
c       DefParReord
c       DistMemReord
c       DistMemPerm
c       Ext_L2s
c       Ext_L2u
c       Ext_L1
c       Ext_L0
c       CVE2
c        MkT20p
c        MkT20u
c        InsL
c
c        integral based approach routines
c        Ext_W4
c          Ext_W4hlp1
c          Ext_W4hlp2
c          Ext_W4hlp3
c        Ext_W3
c
c        JoinLvec
c
c        DefW34y
c        DefW3y
c        DefW4y
c        DefW4y2
c
c       --------------------------------
c
        subroutine DefParReord (NaGrpR,maxdim)

c
c       This routine do:
c       define parameters in chcc_reord.fh using NaGrpR,maxdim
c
c       I/O parameter description:
c       NxGrpR   - # of groups in a (=b) set (I)
c       maxdim   - # maximal dimension of a (=b) Groups(O)
c
        implicit none
#include "chcc1.fh"
#include "chcc_reord.fh"
#include "chcc_files.fh"
c
        integer NaGrpR,maxdim
c
c       help variables
c
        real*8 rdim
        integer i,j
        integer Up(1:MaxGrp),Low(1:MaxGrp)
c
c
c1      define parameters of Groups of a set
c
        rdim=1.0d0*nv/(1.0d0*NaGrpR)
c
        do i=1,NaGrpR
c
           if (i.eq.1) then
             Up(i)=int(rdim*i)
             Low(i)=1
           else if (i.eq.NaGrpR) then
             Up(i)=nv
             Low(i)=Up(i-1)+1
           else
             Up(i)=int(rdim*i)
             Low(i)=Up(i-1)+1
           end if
c
           DimGrpaR(i)=(Up(i)-Low(i))+1
c
        end do
c
c
c2      find maximal dimensions of a'
c
        maxdim=DimGrpaR(1)
        do i=1,NaGrpR
          if (DimGrpaR(i).gt.maxdim) then
          maxdim=DimGrpaR(i)
          end if
        end do
c
c
c3.1    def L2Name, T2Name, I2Name,I3Name
c
        do i=1,MaxGrp
        do j=1,MaxGrp
          call DefParo3v3Hlp1(i,j,'L2',L2Name(i,j))
          call DefParo3v3Hlp1(i,j,'T2',T2Name(i,j))
          call DefParo3v3Hlp1(i,j,'I2',I2Name(i,j))
          call DefParo3v3Hlp1(i,j,'I3',I3Name(i,j))
        end do
        end do
c
c3.2    def L1Name,I1Name
c
        do i=1,MaxGrp
          call DefParo3v3Hlp2 (i,'L1vc',L1Name(i))
          call DefParo3v3Hlp2 (i,'I1in',I1Name(i))
        end do
c
c3.3    def L0Name,I0Name
c
        L0Name='L0vctr'
        I0Name='I0intg'
c
        return
        end
c
c       ------------------------------------
c
        subroutine DistMemReord (NaGrpR,maxdim,maxdimSG,NchBlk,
     c        PossV1,PossV2,PossV3,PossV4,PossM1,PossM2,
     c        PossT)

c
c       This routine do:
c       define initial possitions of OE,V1-V4,M1,2 arrays
c       described in Reord routine
c
c       Memory requirements:
c        intkey=0
c       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom}
c       V2   - max {ov'ov'; v'v'm, ov'm; oom}
c       V3   - max {ov'm; oom}
c       V4   - oom
c        intkey=0
c       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm; oom; V"V"V"V"}
c       V2   - max {ov'ov'; v'v'm, ov'm; oom}
c       V3   - max {ov'm; oom; V'V'M}
c       V4   - oom
c        M1   - V"V"m
c        M2   - max {V"V"M; OV"M)
c
c       I/O parameter description:
c       NxGrp    - # of groups in a,b,be,ga set (I)
c       maxdim   - maximal dimension of V' (I)
c       NChBlk   - # of Choleski vectors in one Block - m' (I)
c       Possx    - initial possitinos of arrays (O-all)
c       PossT    - initial and last possition (I/O)
c
        implicit none
#include "chcc1.fh"
c
        integer NaGrpR,maxdim,maxdimSG,NchBlk
        integer PossV1,PossV2,PossV3,PossV4,PossM1,PossM2
        integer PossT
c
c       help variables
        integer length,nbas
c
c
c2      V1 file
c       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom, V"V"V"V"}
c
        nbas=no+nv
c
        length=maxdim*maxdim*no*no
        if ((nbas*nbas*NChBlk).gt.length) then
          length=nbas*nbas*NChBlk
        end if
        if ((no*maxdim*nc).gt.length) then
          length=no*maxdim*nc
        end if
        if ((maxdim*maxdim*nc).gt.length) then
          length=maxdim*maxdim*nc
        end if
        if ((no*no*nc).gt.length) then
          length=no*no*nc
        end if
        if ((intkey.eq.1).and.(length.le.maxdimSG**4)) then
          length=maxdimSG**4
        end if
c
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
c
c3      V2 file
c       V2   - max {ov'ov'; v'v'm ; ov'm; oom}
c
        length=maxdim*maxdim*no*no
        if ((maxdim*maxdim*nc).gt.length) then
          length=maxdim*maxdim*nc
        end if
        if ((no*maxdim*nc).gt.length) then
          length=no*maxdim*nc
        end if
        if ((no*no*nc).gt.length) then
          length=no*no*nc
        end if
c
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
c
c
c4      V3 file
c       V3   - max {ov'm; oom, V'V'M}
c
        length=no*maxdim*nc
        if ((no*no*nc).gt.length) then
          length=no*no*nc
        end if
        if ((intkey.eq.1).and.(length.le.maxdim*maxdim*nc)) then
          length=maxdim*maxdim*nc
        end if
c
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
c
c5      V4 file
c       V4   - oom
c
        length=no*no*nc
c
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
c
c6        M1   - V"V"m
c
        length=maxdimSG*maxdimSG*nc
        if (intkey.eq.0) then
          length=0
        end if
        PossM1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M1 ',PossM1,length
        end if
c
c7        M2   - max {V"V"M; OV"M)
c
        length=maxdimSG*maxdimSG*nc
        if (length.lt.no*nc*maxdimSG) then
          length=no*nc*maxdimSG
        end if
        if (intkey.eq.0) then
          length=0
        end if
        PossM2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M2 ',PossM2,length
        end if
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(NaGrpR)
        end
c
c
c       ------------------------------------
c
        subroutine DistMemPerm (PossT)
c
c       This routine do:
c       define initial possitions of Permanent
c
c       PossT    - initial and last possition (I/O)
c
        implicit none
#include "chcc1.fh"
c
        integer PossT
c
c       help variables
        integer length,nbas(1)
c
c
c1.1    Foo file
        length=no*no
        PossFoo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Foo ',PossFoo,length
        end if
c
c1.2    Fvo file
        length=no*nv
        PossFvo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Fvo ',PossFvo,length
        end if
c
c1.3    Fvv file
        length=nv*nv
        PossFvv=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Fvv ',PossFvv,length
        end if
c
c
c2      OE file
c
        call Get_iArray('nBas',nBas,1)
        length=nbas(1)
c
        PossOE=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM OE ',PossOE,length
        end if
c
c
c3.1    T1o file
        length=no*nv
        PossT1o=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM T1o ',PossT1o,length
        end if
c
c3.2    T1n file
        length=no*nv
        PossT1n=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM T1n ',PossT1n,length
        end if
c
c
c4.1    Hoo file
        length=no*no
        PossHoo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Hoo ',PossHoo,length
        end if
c
c4.2    Hvo file
        length=no*nv
        PossHvo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Hvo ',PossHvo,length
        end if
c
c4.3    Hvv file
        length=nv*nv
        PossHvv=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Hvv ',PossHvv,length
        end if
c
c
c5.1    Goo file
        length=no*no
        PossGoo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Goo ',PossGoo,length
        end if
c
c5.2    Hvv file
        length=nv*nv
        PossGvv=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Gvv ',PossGvv,length
        end if
c
c6.1    A files @@ A file medzi fixnymi je na zamyslenie (lebo je o4)
        length=no*no*no*(no+1)/2
        PossA=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM A   ',PossA,length
        end if
        if (intkey.eq.0) then
          PossAex=PossT
          PossT=PossT+length
          if (printkey.ge.10) then
          write (6,*) 'DM Aex ',PossAex,length
          end if
        else
          PossAex=PossT
        end if
c
c
comega  PossFree - Possition of the space, where work arrays started
        PossFree=PossT
c
        return
        end
c
c
c       ---------------------
c
        subroutine Ext_L2s (V1,V2,
     c                      dima,dimab,dimc,adda,addb,nbs)
c
c       this routine do:
c       V2(a'b',m') <- V1(p,q,m') for aGrp=bGrp
c
        implicit none
        integer dima,dimab,dimc,adda,addb,nbs
        real*8 V1(1:nbs,1:nbs,1:dimc)
        real*8 V2(1:dimab,1:dimc)
c
c       help variables
        integer a,b,ab,m
c
        do m=1,dimc
        ab=0
        do a=1,dima
        do b=1,a
        ab=ab+1
        V2(ab,m)=V1(adda+a,addb+b,m)
        end do
        end do
        end do
c
        return
        end
c
c       ---------------------
c
        subroutine Ext_L2u (V1,V2,
     c                      dima,dimb,dimc,adda,addb,nbs)
c
c       this routine do:
c       V2(a',b',m') <- V1(p,q,m') for aGrp>bGrp
c
        implicit none
        integer dima,dimb,dimc,adda,addb,nbs
        real*8 V1(1:nbs,1:nbs,1:dimc)
        real*8 V2(1:dima,1:dimb,1:dimc)
c
c       help variables
        integer a,b,m
c
        do m=1,dimc
        do b=1,dimb
        do a=1,dima
        V2(a,b,m)=V1(adda+a,addb+b,m)
        end do
        end do
        end do
c
        return
        end
c
c       ---------------------
c
        subroutine  Ext_L1 (V1,V2,
     c                      no,dima,dimc,adda,nbs)
c
c       this routine do:
c       V2(i,a',m') <- V1(p,q,m')
c
        implicit none
        integer no,dima,dimc,adda,nbs
        real*8 V1(1:nbs,1:nbs,1:dimc)
        real*8 V2(1:no,1:dima,1:dimc)
c
c       help variables
        integer i,a,m
c
        do m=1,dimc
        do a=1,dima
        do i=1,no
        V2(i,a,m)=V1(i,adda+a,m)
        end do
        end do
        end do
c
        return
        end
c
c       ---------------------
c
        subroutine  Ext_L0 (V1,V2,
     c                      no,dimij,dimc,nbs)
c
c       this routine do:
c       V2(i,j,m') <- V1(p,q,m')
c
        implicit none
        integer no,dimij,dimc,nbs
        real*8 V1(1:nbs,1:nbs,1:dimc)
        real*8 V2(1:dimij,1:dimc)
c
c       help variables
        integer i,j,ij,m
c
        do m=1,dimc
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
        V2(ij,m)=V1(i,j,m)
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine CVE2 (V,oe,dima,dimb,adda,addb,no,e2,e2os)
c
c       this routine do:
c       e2  =e2   + sum(a',b',i,j)[(a'i|b'j).{2(a'i|b'j)-(a'j|b'i)}]/Dija'b'
c       e2os=e2os + sum(a',b',i,j)[(a'i|b'j).{ (a'i|b'j)          }]/Dija'b'
c       N.B. Cierny Vypocet E2
c       N.B.2 qvajt odflaknute
c
        implicit none
        integer dima,dimb,adda,addb,no
        real*8 V(1:dima,1:no,1:dimb,1:no)
        real*8 oe(*)
        real*8 e2,e2os
c
c       help variables
        integer i,j,a,b
        real*8 dijab
c
c
        do j=1,no
          do b=1,dimb
            do i=1,no
              do a=1,dima
                dijab=oe(i)+oe(j)-oe(adda+a)-oe(addb+b)
                e2  =e2   +
     &     V(a,i,b,j)*(2.0d0*V(a,i,b,j)-V(a,j,b,i))/dijab
                e2os=e2os +
     &     V(a,i,b,j)*(      V(a,i,b,j)           )/dijab
              end do
            end do
          end do
        end do
c
c
        return
        end
c
c       --------------------------------
c
        subroutine MkT20p (T2,V,oe,dima,adda,no)
c
c       this routine do:
c       T2(a'b',i,j)=(a'i|b'j)/Dija'b'
c        for aGrp.eq.bGrp
c       N.B.2 qvajt odflaknute, neurobene
c
        implicit none
        integer dima,adda,no
        real*8 V(1:dima,1:no,1:dima,1:no)
        real*8 T2(1:dima*(dima+1)/2,1:no,1:no)
        real*8 oe(*)
c
c       help variables
        integer i,j,a,b,ab
        real*8 dijab
c
c
        do j=1,no
          do i=1,no
            ab=0
            do a=1,dima
            do b=1,a
              ab=ab+1
              dijab=oe(i)+oe(j)-oe(adda+a)-oe(adda+b)
              T2(ab,i,j)=V(a,i,b,j)/dijab
            end do
            end do
          end do
        end do
c
c
        return
        end
c
c       --------------------------------
c
        subroutine MkT20u (T2,V,oe,dima,dimb,adda,addb,no)
c
c       this routine do:
c       T2(a',b',i,j)=(a'i|b'j)/Dija'b'
c        for aGrp.ne.bGrp
c       N.B.2 qvajt odflaknute
c
        implicit none
        integer dima,dimb,adda,addb,no
        real*8 T2(1:dima,1:dimb,1:no,1:no)
        real*8 V(1:dima,1:no,1:dimb,1:no)
        real*8 oe(*)
c
c       help variables
        integer i,j,a,b
        real*8 dijab
c
c
        do j=1,no
          do i=1,no
            do b=1,dimb
              do a=1,dima
                dijab=oe(i)+oe(j)-oe(adda+a)-oe(addb+b)
                T2(a,b,i,j)=V(a,i,b,j)/dijab
              end do
            end do
          end do
        end do
c
c
        return
        end
c
c       --------------------------------
c
        subroutine InsL (Ll,Lg,ncLoc,nc,ncOff,dim1)
c
c        this routine do:
c        insert Llocal(ml,dim1) into Lglobal(m,dim1)
c        on a corresponding place
c
        implicit none
        integer ncLoc,nc,ncOff,dim1
        real*8 Ll(1:ncLoc,1:dim1)
        real*8 Lg(1:nc,1:dim1)
c
c        help variables
        integer m,p
c
        do p=1,dim1
          do m=1,ncLoc
            Lg(ncOff+m,p)=Ll(m,p)
          end do
        end do
c
        return
        end
c
c       --------------------------------
c       --------------------------------
c
        subroutine Ext_W4 (V2,M1,
     c                     nc,dima,dimb,dimab,dimapp,dimbpp,dimabpp,
     c                     addapp,addbpp,aGrp,bGrp,aSGrp,bSGrp)
c
c        this routine is a control routine to:
c        Extract M1(m,a"b") <- V2(m,a'b')
c        for all combinations of aGrp,bGrp,aSGrp,bSGrp
c
        implicit none
        integer nc,dima,dimb,dimab,dimapp,dimbpp,dimabpp
        integer addapp,addbpp,aGrp,bGrp,aSGrp,bSGrp
        real*8 M1(1)
        real*8 V2(1)
c
        if (aGrp.eq.bGrp) then
c          case V2(m,a'b')
          if (aSGrp.eq.bSGrp) then
c            subcase M1(m,a"b") <- V2(m,a'b')
            call Ext_W4hlp1 (V2,M1,nc,dima,dimab,dimapp,dimabpp,addapp)
          else
c            subcase M1(m,a",b") <- V2(m,a'b')
            call Ext_W4hlp2 (V2,M1,
     c                       nc,dima,dimab,dimapp,dimbpp,addapp,addbpp)
          end if
        else
c          case M1(m,a",b") <- V2(m,a',b')
          call Ext_W4hlp3 (V2,M1,
     c                     nc,dima,dimb,dimapp,dimbpp,addapp,addbpp)
        end if
c
        return
        end
c
c       --------------------------------
c
        subroutine Ext_W4hlp1 (V2,M1,
     c                         nc,dima,dimab,dimapp,dimabpp,addapp)
c
c        this routine do:
c        Extract M1(m,a"b") <- V2(m,a'b')
c
        implicit none
        integer nc,dima,dimab,dimapp,dimabpp,addapp
        real*8 V2(1:nc,1:dimab)
        real*8 M1(1:nc,1:dimabpp)
c
c        help variables
c
        integer a,ab,app,bpp,abpp,m
c
        abpp=0
        do app=1,dimapp
        a=addapp+app
        ab=a*(a-1)/2+addapp
        do bpp=1,app
        ab=ab+1
        abpp=abpp+1
c
          do m=1,nc
          M1(m,abpp)=V2(m,ab)
          end do
c
        end do
        end do
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(dima)
        end
c
c       --------------------------------
c
        subroutine Ext_W4hlp2 (V2,M1,
     c                       nc,dima,dimab,dimapp,dimbpp,addapp,addbpp)
c
c        this routine do:
c        Extract M1(m,a",b") <- V2(m,a'b')
c
        implicit none
        integer nc,dima,dimab,dimapp,dimbpp,addapp,addbpp
        real*8 V2(1:nc,1:dimab)
        real*8 M1(1:nc,1:dimapp,1:dimbpp)
c
c        help variables
c
        integer a,ab,app,bpp,m
c
        do app=1,dimapp
        a=addapp+app
        ab=a*(a-1)/2+addbpp
        do bpp=1,dimbpp
        ab=ab+1
c
          do m=1,nc
          M1(m,app,bpp)=V2(m,ab)
          end do
c
        end do
        end do
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(dima)
        end
c
c       --------------------------------
c
        subroutine Ext_W4hlp3 (V2,M1,
     c                     nc,dima,dimb,dimapp,dimbpp,addapp,addbpp)
c
c        this routine do:
c        Extract M1(m,a",b") <- V2(m,a',b')
c
        implicit none
        integer nc,dima,dimb,dimapp,dimbpp,addapp,addbpp
        real*8 V2(1:nc,1:dima,1:dimb)
        real*8 M1(1:nc,1:dimapp,1:dimbpp)
c
c        help variables
c
        integer a,b,app,bpp,m
c
        do app=1,dimapp
        a=addapp+app
        do bpp=1,dimbpp
        b=addbpp+bpp
c
          do m=1,nc
          M1(m,app,bpp)=V2(m,a,b)
          end do
c
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine Ext_W3 (V3,M2,nc,no,dimc,dimcpp,addcpp)
c
c        this routine do:
c         Extract M2(m,c",i) <- V3(m,c',i)
c
        implicit none
        integer nc,no,dimc,dimcpp,addcpp
        real*8 V3(1:nc,1:dimc,1:no)
        real*8 M2(1:nc,1:dimcpp,1:no)
c
c        help varibles
        integer m,i,c,cpp
c
        do i=1,no
          c=addcpp
          do cpp=1,dimcpp
          c=c+1
            do m=1,nc
              M2(m,cpp,i)=V3(m,c,i)
            end do
          end do
        end do
c
        return
        end
c
c        ---------------------------
c
        subroutine JoinLvec (wrk,wrksize,
     c             PossV1,PossV2,NaGrpR,LunAux)
c
c       This routine do:
c        Join L0-L2 files (Global, dimensioned as nc)
c        from local L0-L2 files (dimensioned as ncLoc)
c        N.B. This file have sense only for paralell run
c
c
c       Structure of Choleski vector files
c
c       L0(m,IJ)    L0vctr  I>=J
c
c       L1(m,I ,A') L1vcxx xx - Group of A'
c@@        kokot som, ze som to takto urobil, prerobit na L(m,a,i) to treba
c
c       L2(m,A'B')  L2xxyy xx - Group of A', A'>=B'
c                          yy - Group of B'
c
c       Memory requirements:
c        real:
c       V1   - max {ov'm; v'v'm;  oom}
c       V2   - max {v'v'm, ov'm; oom}
c        actual:
c       reord routine requirements are used (DistMemReord)
c
        implicit none
        integer PossV1,PossV2,NaGrpR,LunAux
#include "chcc1.fh"
#include "chcc_reord.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
#include "wrk.fh"
#ifdef _MOLCAS_MPP_
#include "para_info.fh"
#include "chcc_parcc.fh"
c
c       help variables
        character*6 LunName
        integer aGrp,bGrp
        integer dim1,dima,dimb,dimab
        integer i,ncLoc,ncOff
c
c
c7.0        calculate ncOffset for given node
        ncLoc=NChLoc(myRank)
        ncOff=0
        if (myRank.gt.0) then
          do i=0,myRank-1
          ncOff=ncOff+NChLoc(i)
          end do
        end if
c
c
c
c7.1        Make global L0
        LunName=L0Name
c
c7.1.1  Read Local L0 = V2(ml,ij) on proper place in V2
        dim1=ncLoc*no*(no+1)/2
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
c7.1.2        Vanish V1(m,ij)
        dim1=nc*no*(no+1)/2
        call mv0zero (dim1,dim1,wrk(PossV1))
c
c7.1.3        Insert V2(ml,ij) -> V1(m,ij)
        dim1=no*(no+1)/2
        call InsL (wrk(PossV2),wrk(PossV1),ncLoc,nc,ncOff,dim1)
c
c##        Synchronizacny bod
c7.1.4        Allreduce V1
        dim1=nc*no*(no+1)/2
        call gadgop (wrk(PossV1),dim1,'+')
c
c7.1.5        Save L0 (Global)
        dim1=nc*no*(no+1)/2
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c
c
c7.2    make global L1 files
c
        do aGrp=1,NaGrpR
        dima=DimGrpaR(aGrp)
        LunName=L1Name(aGrp)
c
c7.2.1  read L1 = V2(ml,i,a') back into file
        dim1=no*dima*ncLoc
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
c7.2.2        Vanish V1(m,i,a')
        dim1=nc*no*dima
        call mv0zero (dim1,dim1,wrk(PossV1))
c
c7.2.3        Insert V2(ml,i,a') -> V1(m,i,a')
        dim1=no*dima
        call InsL (wrk(PossV2),wrk(PossV1),ncLoc,nc,ncOff,dim1)
c
c##        Synchronizacny bod
c7.2.4        Allreduce V1
        dim1=nc*no*dima
        call gadgop (wrk(PossV1),dim1,'+')
c
c7.2.5        Save V1(m,i,a') (Global)
        dim1=nc*no*dima
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
        end do
c
c
c
c2.3    Make glogal L2(m,a'b') files
c
        do aGrp=1,NaGrpR
        dima=DimGrpaR(aGrp)
        do bGrp=1,aGrp
        dimb=DimGrpaR(bGrp)
        if(aGrp.eq.bGrp) then
          dimab=dima*(dima+1)/2
        else
          dimab=dima*dimb
        end if
        LunName=L2Name(aGrp,bGrp)
c
c
c7.3.1  read L2 = V2(ml,a'b') back into file
        dim1=dimab*ncLoc
        call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
c7.3.2        Vanish V1(m,a'b')
        dim1=nc*dimab
        call mv0zero (dim1,dim1,wrk(PossV1))
c
c7.3.3        Insert V2(ml,a',b') -> V1(m,a',b')
        dim1=dimab
        call InsL (wrk(PossV2),wrk(PossV1),ncLoc,nc,ncOff,dim1)
c
c##        Synchronizacny bod
c7.3.4        Allreduce V1
        dim1=nc*dimab
        call gadgop (wrk(PossV1),dim1,'+')
c
c7.3.5        Save V1(ml,a',b') (Global)
        dim1=nc*dimab
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
        end do
        end do
c
#else
c Avoid unused argument warnings
        if (.false.) then
          call Unused_real_array(wrk)
          call Unused_integer(PossV1)
          call Unused_integer(PossV2)
          call Unused_integer(NaGrpR)
          call Unused_integer(LunAux)
        end if
#endif
c
        return
        end
c
c        ---------------------------
c
        subroutine  DefW34y (aGrp,bGrp,w3y,w4y,NaGrp)
c
c        this routine do:
c        define w3y and w4y keys, which indicates, if atleast one
c        W3/W4 file need to be calculated on this node for given a',b'
c
        implicit none
        integer aGrp,bGrp,w3y,w4y,NaGrp
#include "chcc1.fh"
#include "o2v4.fh"
#ifdef _MOLCAS_MPP_
#include "chcc_parcc.fh"
c
c        help variables
        integer aSGrp,bSGrp,abSGrp,bSGrpUp,cSGrp,cdSGrp
        integer NSGrp
c
        w3y=0
        w4y=0
c
c       cycle over a">=b" subgroups for a',b'
        do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
        if (aGrp.eq.bGrp) then
          bSGrpUp=aSGrp
        else
          bSGrpUp=GrpaUp(bGrp)
        end if
        do bSGrp=GrpaLow(bGrp),bSGrpUp
        abSGrp=aSGrp*(aSGrp-1)/2+bSGrp
c
c          cycle over all c"
          NSGrp=GrpaUp(NaGrp)
          do cSGrp=1,NSGrp
          if (InqW3(abSGrp,cSGrp)) then
            w3y=w3y+1
          end if
          end do
c
c          cycle over all c">=d"
          do cdSGrp=1,NSGrp*(NSGrp+1)/2
          if (InqW4(abSGrp,cdSGrp)) then
            w4y=w4y+1
          end if
          end do
c
c       end cycle over a">=b" subgroups
        end do
        end do
c
#else
        w3y=1
        w4y=1
c Avoid unused argument warnings
        if (.false.) then
          call Unused_integer(aGrp)
          call Unused_integer(bGrp)
          call Unused_integer(NaGrp)
        end if
#endif
c
        return
        end
c
c        ---------------------------
c
        subroutine  DefW3y (aGrp,bGrp,cGrp,w3y)
c
c        this routine do:
c        define w3y key, which indicates, if atleast one
c        W3 file need to be calculated on this node for given
c        a',b',c'
c
        implicit none
        integer aGrp,bGrp,cGrp,w3y
#include "chcc1.fh"
#include "o2v4.fh"
#ifdef _MOLCAS_MPP_
#include "chcc_parcc.fh"
c
c        help variables
        integer aSGrp,bSGrp,abSGrp,bSGrpUp,cSGrp
c
        w3y=0
c
c       cycle over a">=b" subgroups for a',b'
        do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
        if (aGrp.eq.bGrp) then
          bSGrpUp=aSGrp
        else
          bSGrpUp=GrpaUp(bGrp)
        end if
        do bSGrp=GrpaLow(bGrp),bSGrpUp
        abSGrp=aSGrp*(aSGrp-1)/2+bSGrp
c
c          cycle over c" for c'
          do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
          if (InqW3(abSGrp,cSGrp)) then
            w3y=w3y+1
          end if
          end do
c
c       end cycle over a">=b" subgroups
        end do
        end do
c
#else
        w3y=1
c Avoid unused argument warnings
        if (.false.) then
          call Unused_integer(aGrp)
          call Unused_integer(bGrp)
          call Unused_integer(cGrp)
        end if
#endif
c
        return
        end
c
c        ---------------------------
c
        subroutine  DefW4y (aGrp,bGrp,cGrp,dGrp,w4y)
c
c        this routine do:
c        define w4y key, which indicates, if atleast one
c        W4 file need to be calculated on this node for
c        given a',b',c',d'
c
        implicit none
        integer aGrp,bGrp,cGrp,dGrp,w4y
#include "chcc1.fh"
#include "o2v4.fh"
#ifdef _MOLCAS_MPP_
#include "chcc_parcc.fh"
c
c        help variables
        integer aSGrp,bSGrp,abSGrp,bSGrpUp,cSGrp,cdSGrp,dSGrp,dSGrpUp
c
        w4y=0
c
c       cycle over a">=b" subgroups for a',b'
        do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
        if (aGrp.eq.bGrp) then
          bSGrpUp=aSGrp
        else
          bSGrpUp=GrpaUp(bGrp)
        end if
        do bSGrp=GrpaLow(bGrp),bSGrpUp
        abSGrp=aSGrp*(aSGrp-1)/2+bSGrp
c
c          cycle over  c">=d" for c',d'
          do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
          if (cGrp.eq.dGrp) then
            dSGrpUp=cSGrp
          else
            dSGrpUp=GrpaUp(dGrp)
          end if
          do dSGrp=GrpaLow(dGrp),dSGrpUp
          cdSGrp=cSGrp*(cSGrp-1)/2+dSGrp
            if (InqW4(abSGrp,cdSGrp)) then
              w4y=w4y+1
            end if
          end do
          end do
c
c       end cycle over a">=b" subgroups
        end do
        end do
c
#else
        w4y=1
c Avoid unused argument warnings
        if (.false.) then
          call Unused_integer(aGrp)
          call Unused_integer(bGrp)
          call Unused_integer(cGrp)
          call Unused_integer(dGrp)
        end if
#endif
c
        return
        end
c
c        ---------------------------
c
        subroutine  DefW4y2 (aSGrp,bSGrp,cGrp,dGrp,w4y)
c
c        this routine do:
c        define w4y key, which indicates, if atleast one
c        W4 file need to be calculated on this node for
c        given a",b",c',d'
c
        implicit none
        integer aSGrp,bSGrp,cGrp,dGrp,w4y
#include "chcc1.fh"
#include "o2v4.fh"
#ifdef _MOLCAS_MPP_
#include "chcc_parcc.fh"
c
c        help variables
        integer abSGrp,cSGrp,cdSGrp,dSGrp,dSGrpUp
c
        w4y=0
c
c       def abSGrp
        abSGrp=aSGrp*(aSGrp-1)/2+bSGrp
c
c          cycle over  c">=d" for c',d'
          do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
          if (cGrp.eq.dGrp) then
            dSGrpUp=cSGrp
          else
            dSGrpUp=GrpaUp(dGrp)
          end if
          do dSGrp=GrpaLow(dGrp),dSGrpUp
          cdSGrp=cSGrp*(cSGrp-1)/2+dSGrp
            if (InqW4(abSGrp,cdSGrp)) then
              w4y=w4y+1
            end if
          end do
          end do
c
c
#else
        w4y=1
c Avoid unused argument warnings
        if (.false.) then
          call Unused_integer(aSGrp)
          call Unused_integer(bSGrp)
          call Unused_integer(cGrp)
          call Unused_integer(dGrp)
        end if
#endif
c
        return
        end
