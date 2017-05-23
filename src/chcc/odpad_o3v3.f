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
c       MkV_K22
c       MkV_Q22
c       ExtT1
c       MkT_QK42
c       MkD_Q46
c       MkI_Q47
c       MkE_Y3
c       AdV_G2
c       DefParo3v3
c         DefParo3v3Hlp1
c         DefParo3v3Hlp2
c       DistMemo3v3jk
c       DistMemo3v3t2
c       DistMemo3v3chol
c        MkV_Hvv2
c        AdH_Hvv2
c        ExH_X2
c        MkT_T17
c        MkV_Hvo2
c        MkV_Hoo2
c        MkV_Goo3
c        AdV_A23
c        MkV_A1
c        MkV_A4
c        ExV_X41
c        ExV_X43
c        ExV_X42
c        ExA_X4
c        MkV_T18
c        AdT_T17
c        MkT_T15
c        DfH_Hvv1 (suspended)
c        ExH_T13
c        Ext_Aex
c
c
c       --------------------------------
c
        subroutine MkV_K22 (W1,W2,dim)
c
c       this routine do
c       W1(p) = -W2(p)
c
c       N.B. toto je iba plytke copy s opacnym znamienkom,
c            da sa aj s blasmi spravit
c       N.B. Kvajt odflaknute
c
c
        implicit none
        integer dim
        real*8 W1(1:dim)
        real*8 W2(1:dim)
c
c       help variables
        integer p
c
c
        do p=1,dim
          W1(p)=-W2(p)
        end do
c
c
        return
        end
c
c       --------------------------------
c
        subroutine MkV_Q22 (W2,W1,dima)
c
c       this routine do
c       W1(j,u,i,a') = 2W2(i,u,j,a')-W2(j,u,i,a')
c
c       N.B. Kvajt odflaknute
c
        implicit none
#include "chcc1.fh"
        integer dima
        real*8 W1(1:no,1:no,1:no,1:dima)
        real*8 W2(1:no,1:no,1:no,1:dima)
c
c       help variables
        integer i,j,u,a
c
c
        do a=1,dima
        do i=1,no
        do u=1,no
        do j=1,no
          W1(j,u,i,a)=W2(j,u,i,a)-2.0d0*W2(i,u,j,a)
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
        subroutine ExtT1 (H,T1,dima,adda)
c
c       this routine do:
c       Extract H(a',i) <- T1(a,i) for given aGrp
c
c       parameter description:
c       H       - Output file (O)
c       T1      - T1 amplitudes (I)
c       dima    - dimension of given Group (I)
c       adda    - shift of a' in full a set (I)
c
c       N.B. Kvajt odflaknute
c
        implicit none
#include "chcc1.fh"
        integer dima,adda
        real*8 T1(1:nv,1:no)
        real*8 H(1:dima,1:no)
c
c       help variables
        integer a,i
c
        do i=1,no
          do a=1,dima
          H(a,i)=T1(adda+a,i)
          end do
        end do
c
c
        return
        end
c
c       --------------------------------
c
        subroutine MkT_QK42 (T2,T1a,T1b,dima,dimb,no,f1,f2)
c
c       this routine do:
c       T2(a',b',i,j) <- f1 . T2(a',b',i,j) + f2 . T1a(a,i) . T1b(b,j)
c
c       N.B. Kvajt odflaknute
c
        implicit none
        integer dima,dimb,no
        real*8 f1,f2
        real*8 T2(1:dima,1:dimb,1:no,1:no)
        real*8 T1a(1:dima,1:no)
        real*8 T1b(1:dimb,1:no)
c
c       help variables
        integer a,b,i,j
        real*8 c
c

        do j=1,no
          do b=1,dimb
          c=f2*T1b(j,b)
            do i=1,no
              do a=1,dima
                T2(a,b,i,j)=f1*T2(a,b,i,j)+T1a(a,i)*c
              end do
            end do
          end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine MkD_Q46 (D,T2,T1a,T1b,dima,dimb,no)
c
c       this routine do:
c       Create D(a',j,b',i) [NB. D(a,b,i,j) =  t2(a,b,j,i)-T2(a,b,i,j)]
c       from following available arrays (permuted as given):
c       T2(a',j,b',i) [NB. T2abij = 1/2 t2abij + tai . tjb ]
c       T1a(a',p)
c       T1b(b',p)
c
c       N.B. Kvajt odflaknute

c
        implicit none
        integer dima,dimb,no
        real*8 D(1:dima,1:no,1:dimb,1:no)
        real*8 T2(1:dima,1:no,1:dimb,1:no)
        real*8 T1a(1:dima,1:no)
        real*8 T1b(1:dimb,1:no)
c
c       help variables
        integer a,b,i,j
c
        do i=1,no
        do b=1,dimb
        do j=1,no
        do a=1,dima
        D(a,j,b,i)=2.0d0*(T2(a,i,b,j)-T1a(a,j)*T1b(b,i))-T2(a,j,b,i)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine MkI_Q47 (Va,V,dimb,dima,no)
c
c       this routine do:
c       Create Va(B',o_A,o_B,A')  needed in step Q48
c       from following available array (permuted as given):
c       V1(B',o_A,o_B,A') = (A',o_A|B',o_B)
c       N.B.
c
c       N.B. Kvajt odflaknute, aj koment k rutine odflaknuty

c
        implicit none
        integer dima,dimb,no
        real*8 Va(1:dimb,1:no,1:no,1:dima)
        real*8 V(1:dimb,1:no,1:no,1:dima)
c
c       help variables
        integer a,b,i,j
c
        do a=1,dima
        do j=1,no
        do i=1,no
        do b=1,dimb
          Va(b,i,j,a)=2.0d0*V(b,j,i,a)-V(b,i,j,a)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine MkE_Y3 (Va,V,dima,dimb,no)
c
c       this routine do:
c       Va(a',i,b',j) = 2 V(a',j,b',i) - V(a',i,b',j)
c
c       N.B. Kvajt odflaknute

c
        implicit none
        integer dima,dimb,no
        real*8 Va(1:dima,1:no,1:dimb,1:no)
        real*8 V(1:dima,1:no,1:dimb,1:no)
c
c       help variables
        integer a,b,i,j
c
        do j=1,no
        do b=1,dimb
        do i=1,no
        do a=1,dima
          Va(a,i,b,j)=2.0d0*V(a,j,b,i)-V(a,i,b,j)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine AdV_G2 (G2,V,nv,dimbe,dima,no,addbe,adda,fact)
c
c       this routine do:
c       G2(be',a') <- fact sum(i) V(be',i,i,a')
c
c       N.B. Kvajt odflaknute

c
        implicit none
        integer nv,dimbe,dima,no,addbe,adda
        real*8 fact
        real*8 G2(1:nv,1:nv)
        real*8 V(1:dimbe,1:no,1:no,1:dima)
c
c       help variables
c
        integer be,a,i,afull
c
        do a=1,dima
        afull=adda+a
        do i=1,no
        do be=1,dimbe
          G2(addbe+be,afull)=G2(addbe+be,afull)+fact*V(be,i,i,a)
        end do
        end do
        end do
c
        return
        end
c
c       ------------------------------------
c
        subroutine DefParo3v3 (NvGrp,maxdim)
c
c       This routine do:
c       define parameters in o3v3.fh using NvGrp,maxdim
c
c       I/O parameter description:
c       NvGrp    - # of groups in a,b,be,ga set (I)
c       maxdim   - # maximal dimension of (a,b,be,ga)" Groups(O)
c
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NvGrp,maxdim
c
c       help variables
c
        real*8 rdim
        integer i,j
        integer Up(1:MaxGrp),Low(1:MaxGrp)
c
c
c1      define parameters of Groups of v set
c
        rdim=1.0d0*nv/(1.0d0*NvGrp)
c
        do i=1,NvGrp
c
           if (i.eq.1) then
             Up(i)=int(rdim*i)
             Low(i)=1
           else if (i.eq.NvGrp) then
             Up(i)=nv
             Low(i)=Up(i-1)+1
           else
             Up(i)=int(rdim*i)
             Low(i)=Up(i-1)+1
           end if
c
          DimGrpv(i)=(Up(i)-Low(i))+1
c
        end do
c
c
c5      find maximal dimensions of v'
c
        maxdim=DimGrpv(1)
        do i=1,NvGrp
          if (DimGrpv(i).gt.maxdim) then
          maxdim=DimGrpv(i)
          end if
        end do
c
c
c6.1    def L2Name, T2Name, I2Name,I3Name,Tmp1Name,Tmp2Name
c
        do i=1,MaxGrp
        do j=1,MaxGrp
          call DefParo3v3Hlp1(i,j,'L2',L2Name(i,j))
          call DefParo3v3Hlp1(i,j,'T2',T2Name(i,j))
          call DefParo3v3Hlp1(i,j,'I2',I2Name(i,j))
          call DefParo3v3Hlp1(i,j,'I3',I3Name(i,j))
          call DefParo3v3Hlp1(i,j,'X1',Tmp1Name(i,j))
          call DefParo3v3Hlp1(i,j,'X2',Tmp2Name(i,j))
        end do
        end do
c
c6.2    def L1Name,I1Name
c
        do i=1,MaxGrp
          call DefParo3v3Hlp2 (i,'L1vc',L1Name(i))
          call DefParo3v3Hlp2 (i,'I1in',I1Name(i))
        end do
c
        return
        end
c
c       ---------------------
c
        subroutine DefParo3v3Hlp1 (i,j,Schem,Nomen)
c
c       help routine to DefParo2v4, producing names of Disc files
c       ex: Schem='XY', i=1, j=3  ->  Nomen='XY0103'
c        N.B. suspendovana rutina
c
        implicit none
        integer i,j
        character*2 Schem
        character*6 Nomen
c
c       help variables
        character*1 Chr(1:6)
        character*2 digit(1:64)
        character*2 ichr,jchr
        character*2 baza
        character*6 meno
c
        equivalence (Chr(1),meno)
        equivalence (Chr(1),baza)
        equivalence (Chr(3),ichr)
        equivalence (Chr(5),jchr)
c
c
c        quite a porno this piece
        digit(1)='01'
        digit(2)='02'
        digit(3)='03'
        digit(4)='04'
        digit(5)='05'
        digit(6)='06'
        digit(7)='07'
        digit(8)='08'
        digit(9)='09'
        digit(10)='10'
        digit(11)='11'
        digit(12)='12'
        digit(13)='13'
        digit(14)='14'
        digit(15)='15'
        digit(16)='16'
        digit(17)='17'
        digit(18)='18'
        digit(19)='19'
        digit(20)='20'
        digit(21)='21'
        digit(22)='22'
        digit(23)='23'
        digit(24)='24'
        digit(25)='25'
        digit(26)='26'
        digit(27)='27'
        digit(28)='28'
        digit(29)='29'
        digit(30)='30'
        digit(31)='31'
        digit(32)='32'
        digit(33)='33'
        digit(34)='34'
        digit(35)='35'
        digit(36)='36'
        digit(37)='37'
        digit(38)='38'
        digit(39)='39'
        digit(40)='40'
        digit(41)='41'
        digit(42)='42'
        digit(43)='43'
        digit(44)='44'
        digit(45)='45'
        digit(46)='46'
        digit(47)='47'
        digit(48)='48'
        digit(49)='49'
        digit(50)='50'
        digit(51)='51'
        digit(52)='52'
        digit(53)='53'
        digit(54)='54'
        digit(55)='55'
        digit(56)='56'
        digit(57)='57'
        digit(58)='58'
        digit(59)='59'
        digit(60)='60'
        digit(61)='61'
        digit(62)='62'
        digit(63)='63'
        digit(64)='64'
c
c
        baza=Schem
        ichr=digit(i)
        jchr=digit(j)
        Nomen=meno
c
        return
        end
c
c       ---------------------
c
        subroutine DefParo3v3Hlp2 (i,Schem,Nomen)
c
c       help routine to DefParo2v4, producing names of Disc files
c       ex: Schem='XYZQ', i=1  ->  Nomen='XYZQ01'
c
        implicit none
        integer i
        character*4 Schem
        character*6 Nomen
c
c       help variables
        character*1 Chr(1:6)
        character*2 digit(1:64)
        character*2 ichr
        character*4 baza
        character*6 meno
c
        equivalence (Chr(1),meno)
        equivalence (Chr(1),baza)
        equivalence (Chr(5),ichr)
c
c
        digit(1)='01'
        digit(2)='02'
        digit(3)='03'
        digit(4)='04'
        digit(5)='05'
        digit(6)='06'
        digit(7)='07'
        digit(8)='08'
        digit(9)='09'
        digit(10)='10'
        digit(11)='11'
        digit(12)='12'
        digit(13)='13'
        digit(14)='14'
        digit(15)='15'
        digit(16)='16'
        digit(17)='17'
        digit(18)='18'
        digit(19)='19'
        digit(20)='20'
        digit(21)='21'
        digit(22)='22'
        digit(23)='23'
        digit(24)='24'
        digit(25)='25'
        digit(26)='26'
        digit(27)='27'
        digit(28)='28'
        digit(29)='29'
        digit(30)='30'
        digit(31)='31'
        digit(32)='32'
        digit(33)='33'
        digit(34)='34'
        digit(35)='35'
        digit(36)='36'
        digit(37)='37'
        digit(38)='38'
        digit(39)='39'
        digit(40)='40'
        digit(41)='41'
        digit(42)='42'
        digit(43)='43'
        digit(44)='44'
        digit(45)='45'
        digit(46)='46'
        digit(47)='47'
        digit(48)='48'
        digit(49)='49'
        digit(50)='50'
        digit(51)='51'
        digit(52)='52'
        digit(53)='53'
        digit(54)='54'
        digit(55)='55'
        digit(56)='56'
        digit(57)='57'
        digit(58)='58'
        digit(59)='59'
        digit(60)='60'
        digit(61)='61'
        digit(62)='62'
        digit(63)='63'
        digit(64)='64'
c
        baza=Schem
        ichr=digit(i)
        Nomen=meno
c
        return
        end
c
c       ------------------------------------
c
        subroutine DistMemo3v3jk (NvGrp,maxdim,
     c        PossV1,PossV2,PossV3,PossV4,
     c        PossH1,PossH2,PossH3,PossH4,PossH5,
     c        PossK,PossQ,
     c        PossT)

c
c       This routine do:
c       define initial possitions of H,V,XY
c       described in o3v3jk routine
c
c
c       I/O parameter description:
c       NvGrp    - # of groups in a,b,be,ga set (I)
c       maxdim   - maximal dimension of V'
c       Possx    - initial possitinos of arrays (O-all)
c       PossT    - initial and last possition (I/O)
c
c        requirements for o3v3jk:
c        H1 - max {v'o}
c        H2 - max {v'o}
c        H3 - max {v'v'}
c        H4 - max {v'o}
c        H5 - max {v'o}
c        V1 - max {v'ooo, v'v'oo, o2oo}
c        V2 - max {oooo, v'v'oo, v'ooo}
c        V3 - max {v'v'oo}
c        V4 - max {v'ooo}
c        PX - max {v'v'oo}
c        QY - max {v'v'oo}
c
        implicit none
#include "chcc1.fh"
c
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4,PossH5
        integer PossK,PossQ
        integer PossT
c
c       help variables
        integer length
c
c
c1      Q,K (used also as X,Y)
c
        length=no*no*maxdim*maxdim
c
        PossQ=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Q  ',PossQ,length
        end if
        PossK=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM K  ',PossK,length
        end if
c
c2.1    V1 file - max {v'ooo, v'v'oo, o2oo}
c
        length=no*no*maxdim*maxdim
        if (no*maxdim*nc.gt.length) then
          length=no*maxdim*nc
        end if
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
        if (no*no*no*(no+1)/2.gt.length) then
          length=no*no*no*(no+1)/2
        end if
c
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
c
c2.2    V2 files - max {oooo, v'v'oo, v'ooo}
c
        length=no*no*maxdim*maxdim
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
        if (no*no*no*no.gt.length) then
          length=no*no*no*no
        end if
c
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
c
c2.3    V3 file - max {v'v'oo}
c
        length=no*no*maxdim*maxdim
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
c
c2.4    V4 file - max {v'ooo}
c
        length=no*no*no*maxdim
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
c
c3      H1,2 files
c
        length=no*maxdim
c
        PossH1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H1 ',PossH1,length
        end if
        PossH2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H2 ',PossH2,length
        end if
c
c3.2    H3 file
c
        length=maxdim*maxdim
        if (no*maxdim.gt.length) then
          length=no*maxdim
        end if
c
        PossH3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H3 ',PossH3,length
        end if
c
c3.3        H4,H5 file
c
        length=no*maxdim
c
        PossH4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H4 ',PossH4,length
        end if
        PossH5=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H5 ',PossH5,length
        end if
c
c
        if (printkey.ge.10) then
        write (6,*) 'PossT ',PossT
        end if
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(NvGrp)
        end
c
c       ------------------------------------
c
        subroutine DistMemo3v3t2 (NvGrp,maxdim,
     c        PossV1,PossV2,PossV3,PossV4,
     c        PossH1,PossH2,PossH3,PossH4,
     c        PossK,PossQ,
     c        PossT)

c
c       This routine do:
c       define initial possitions of H,V, QK
c
c
c       I/O parameter description:
c       NvGrp    - # of groups in a,b,be,ga set (I)
c       maxdim   - maximal dimension of V'
c       Possx    - initial possitinos of arrays (O-all)
c       PossT    - initial and last possition (I/O)
c
c        requirements for o3v3t2:
c        H1 - max {v'o}
c        H2 - max {v'o}
c        H3 - max {v'v',ooo}
c        H4 - max {v'o}
c        V1 - max {v'ov'o, vv'}
c        V2 - max {v'ov'o}
c        V3 - max {v'ov'o, o2oo}
c        V4 - max {v'ov'o}
c        PX - max {v'ov'o}
c        QY - max {v'ov'o}
c
        implicit none
#include "chcc1.fh"
c
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4
        integer PossK,PossQ
        integer PossT
c
c       help variables
        integer length
c
c
c1      Q,K (used also as X,Y)
c
        length=no*no*maxdim*maxdim
c
        PossQ=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Q  ',PossQ,length
        end if
        PossK=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM K  ',PossK,length
        end if
c
c2.1    V1 file - max {v'ov'o, vv'}
c
        length=no*no*maxdim*maxdim
        if (nv*maxdim.gt.length) then
          length=maxdim*nv
        end if
c
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
c
c2.2    V2-4 files
c        V2 - max {v'ov'o}
c        V3 - max {v'ov'o, o2oo}
c        V4 - max {v'ov'o}
c
        length=no*no*maxdim*maxdim
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
c
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
c
        length=no*no*maxdim*maxdim
        if (no*(no+1)*no*no/2.gt.length) then
          length=no*no*no*(no+1)/2
        end if
c
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
c
        length=no*no*maxdim*maxdim
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
c
c3      H1,2 files
c
        length=no*maxdim
c
        PossH1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H1 ',PossH1,length
        end if
        PossH2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H2 ',PossH2,length
        end if
c
c3.2    H3 file
c
        length=maxdim*maxdim
        if (no*maxdim.gt.length) then
          length=no*maxdim
        end if
        if (no*no*no.gt.length) then
          length=no*no*no
        end if
c
        PossH3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H3 ',PossH3,length
        end if
c
c3.3        H4 file
c
        length=no*maxdim
c
        PossH4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H4 ',PossH4,length
        end if
c
c
        if (printkey.ge.10) then
        write (6,*) 'PossT ',PossT
        end if
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(NvGrp)
        end
c
c       ------------------------------------
c
        subroutine DistMemo3v3chol (NvGrp,maxdim,
     c        PossV1,PossV2,PossV3,PossV4,
     c        PossH1,PossH2,PossH3,PossH4,
     c        PossM1,PossM2,PossM3,PossM4,PossM5,
     c        PossK,PossQ,
     c        PossT)

c
c       This routine do:
c       define initial possitions of T,L,M and W arrays,
c       used in routine o3v3chol
c
c
c       I/O parameter description:
c       NvGrp    - # of groups in a,b,be,ga set (I)
c       maxdim   - maximal dimension of V'
c       Possx    - initial possitinos of arrays (O-all)
c       PossT    - initial and last possition (I/O)
c
c        requirements of o3v3chol step
c        K  - max{v'v'oo}
c        Q  - max{v'v'oo}
c        V1 - max{v'oo2, mv'v', v'v'oo, mv'o, moo, vv}
c        V2 - max{v'ooo, mv'v', v'v'oo, vo}
c        V3 - max{mv'o, v'v'oo, v'ooo}
c        V4 - max{v'ooo}
c        H1 - max{v'o}
c        H2 - max{v'o}
c        H3 - max{v'o}
c        H4 - max{v'o}
c        M1 - max{moo}
c        M2 - max{mv'o}
c        M3 - max{mv'o}
c        M4 - max{moo}
c        M5 - max{mv'o}
c
        implicit none
#include "chcc1.fh"
c
        integer NvGrp,maxdim
        integer PossV1,PossV2,PossV3,PossV4
        integer PossH1,PossH2,PossH3,PossH4
        integer PossM1,PossM2,PossM3,PossM4,PossM5
        integer PossK,PossQ
        integer PossT
c
c       help variables
        integer length
c
c
c1      Q,K (used also as X,Y)
c
        length=no*no*maxdim*maxdim
c
        PossQ=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Q  ',PossQ,length
        end if
        PossK=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM K  ',PossK,length
        end if
c
c2.1    V1 file max{v'oo2, mv'v', v'v'oo, mv'o, moo, vv}
        length=maxdim*no*no*(no+1)/2
        if (maxdim*maxdim*nc.gt.length) then
          length=maxdim*maxdim*nc
        end if
        if (no*no*maxdim*maxdim.gt.length) then
          length=no*no*maxdim*maxdim
        end if
        if (no*nc*maxdim.gt.length) then
          length=no*nc*maxdim
        end if
        if (no*nc*no.gt.length) then
          length=no*no*nc
        end if
        if (nv*nv.gt.length) then
          length=nv*nv
        end if
c
        PossV1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V1 ',PossV1,length
        end if
c
c2.2    V2 file - max{v'ooo, mv'v', v'v'oo, vo}
c
        length=no*no*maxdim*maxdim
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
        if (nc*maxdim*maxdim.gt.length) then
          length=nc*maxdim*maxdim
        end if
        if (no*nv.gt.length) then
          length=no*nv
        end if
c
        PossV2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V2 ',PossV2,length
        end if
c
c2.2    V3 files - max{mv'o, v'v'oo, v'ooo}
c
        length=no*no*maxdim*maxdim
        if (no*no*no*maxdim.gt.length) then
          length=no*no*no*maxdim
        end if
        if (no*nc*maxdim.gt.length) then
          length=nc*no*maxdim
        end if
c
        PossV3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V3 ',PossV3,length
        end if
c
c2.4    V4 files - max{v'ooo}
c
        length=no*no*no*maxdim
c
        PossV4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM V4 ',PossV4,length
        end if
c
c3      H1-4 files - max{v'o}
c
        length=no*maxdim
c
        PossH1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H1 ',PossH1,length
        end if
        PossH2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H2 ',PossH2,length
        end if
        PossH3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H3 ',PossH3,length
        end if
        PossH4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM H4 ',PossH4,length
        end if
c
c4      M1-5 files
c        M1 - max{moo}
c        M2 - max{mv'o}
c        M3 - max{mv'o}
c        M4 - max{moo}
c        M5 - max{mv'o}
c
        length=no*no*nc
        PossM1=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M1 ',PossM1,length
        end if
        length=no*nc*maxdim
        PossM2=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M2 ',PossM2,length
        end if
        length=no*nc*maxdim
        PossM3=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M3 ',PossM3,length
        end if
        length=no*nc*no
        PossM4=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M4 ',PossM4,length
        end if
        length=no*nc*maxdim
        PossM5=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM M5 ',PossM5,length
        end if
c
c
        if (printkey.ge.10) then
        write (6,*) 'PossT ',PossT
        end if
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(NvGrp)
        end
c
c       --------------------------------
c
        subroutine MkV_Hvv2 (Va,V,dima,dimb,no)
c
c        this routine do:
c        Va(b',i,j,a') <- 2(a'i|b'j)-(a'j|b'i)  V(b'I|a'J)
c
        implicit none
        integer dima,dimb,no
        real*8 Va(1:dimb,1:no,1:no,1:dima)
        real*8 V(1:dimb,1:no,1:dima,1:no)
c
c        help variables
        integer i,j,a,b
c
c
        do a=1,dima
        do j=1,no
        do i=1,no
        do b=1,dimb
          Va(b,i,j,a)=2.0d0*V(b,j,a,i)-V(b,i,a,j)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine AdH_Hvv2 (H,Hvv,dima,dimbe,adda,addbe,nv)
c
c        this routine do:
c        Hvv(a,be) <<- - H(be',a')
c
        implicit none
        integer dima,dimbe,nv,adda,addbe
        real*8 H(1:dimbe,1:dima)
        real*8 Hvv(1:nv,1:nv)
c
c        help variables
        integer a,be
c
c
        do be=1,dimbe
          do a=1,dima
            Hvv(adda+a,addbe+be)=Hvv(adda+a,addbe+be)-H(be,a)
          end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine ExH_X2 (Gvv,H,dima,dimbe,nv,adda,addbe)
c
c        this routine do:
c        H(a',be') <- Gvv(be,a)
c
        implicit none
        integer dima,dimbe,nv,adda,addbe
        real*8 H(1:dima,1:dimbe)
        real*8 Gvv(1:nv,1:nv)
c
c        help variables
        integer a,be,bev
c
        do be=1,dimbe
        bev=addbe+be
          do a=1,dima
            H(a,be)=Gvv(bev,adda+a)
          end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine MkT_T17 (Ta,T,dimb,dimbe,no)
c
c       this routine do:
c       Ta(i,b',be',u) <- 2 T(be',b',u,i) - T(be',b',i,u)
c
c       N.B. Qvajt odflaknute
c
        implicit none
        integer dimb,dimbe,no
        real*8 T(1:dimbe,1:dimb,1:no,1:no)
        real*8 Ta(1:no,1:dimb,1:dimbe,1:no)
c
c       help variables
        integer i,u,be,b
c
        do u=1,no
        do i=1,no
        do b=1,dimb
        do be=1,dimbe
          Ta(i,b,be,u)=2.0d0*T(be,b,u,i)-T(be,b,i,u)
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
        subroutine MkV_Hvo2 (V,V2,dimbe,dima,no)
c
c       this routine do:
c       Make AntiSymetric integrals
c        V2(a',i,be',j) <- [2 V(be',j|a'i) - V(be',i|a'j)]
c
        implicit none
        integer dimbe,dima,no
        real*8 V2(1:dima,1:no,1:dimbe,1:no)
        real*8 V(1:dimbe,1:no,1:dima,1:no)
c
c       help variables
        integer a,be,i,j
c
        do j=1,no
        do be=1,dimbe
        do i=1,no
        do a=1,dima
          V2(a,i,be,j)=2.0d0*V(be,j,a,i)-V(be,i,a,j)
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
        subroutine MkV_Hoo2 (V2,V,dima,dimb,no)
c
c       this routine do:
c       Make AntiSymetric integrals
c        V2(i,a',b',j) <- 2 V(b'i|a'j) - V(b'j|a'i)
c
        implicit none
        integer dimb,dima,no
        real*8 V2(1:no,1:dima,1:dimb,1:no)
        real*8 V(1:dimb,1:no,1:dima,1:no)
c
c       help variables
        integer a,b,i,j
c
        do j=1,no
        do b=1,dimb
        do a=1,dima
        do i=1,no
          V2(i,a,b,j)=2.0d0*V(b,j,a,i)-V(b,i,a,j)
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
        subroutine MkV_Goo3 (V,V2,dima,no)
c
c       this routine do:
c       Make AntiSymetric integrals
c        V2(a',j,i,u) <- 2 V(a,j|iu) - (a,i|ju)
c
        implicit none
        integer dima,no
        real*8 V2(1:dima,1:no,1:no,1:no)
        real*8 V(1:dima,1:no,1:no*(no+1)/2)
c
c       help variables
        integer a,i,j,u,iu,ju
c
c
        do u=1,no
c
          do i=1,no
            if (i.gt.u) then
            iu=(i-1)*i/2+u
            else
            iu=(u-1)*u/2+i
            end if
c
            do j=1,no
c
              if (j.gt.u) then
              ju=(j-1)*j/2+u
              else
              ju=(u-1)*u/2+j
              end if
c
              do a=1,dima
c
                V2(a,j,i,u)=2.0d0*V(a,j,iu)-V(a,i,ju)
c
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
        subroutine AdV_A23 (V1,A,dimij,no)
c
c       this routine do:
c       A(ij,u,v) <<- V1(j,iu,v) + V1(i,jv,u)
c
c        Velmi odflaknute, da sa to urobit podstatne lepsie, ale
c        o4 proces osrat fok
c
        implicit none
        integer dimij,no
        real*8 A(1:dimij,1:no,1:no)
        real*8 V1(1:no,1:dimij,1:no)
c
c       help variables
        integer i,j,ij,u,v,iu,jv
c
        do v=1,no
        do u=1,no
c
          ij=0
c
          do i=1,no
          if (i.ge.u) then
          iu=i*(i-1)/2+u
          else
          iu=u*(u-1)/2+i
          end if
c
          do j=1,i
          ij=ij+1
          if (j.ge.v) then
          jv=j*(j-1)/2+v
          else
          jv=v*(v-1)/2+j
          end if
c
            A(ij,u,v)=A(ij,u,v)+V1(j,iu,v)+V1(i,jv,u)
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
c       --------------------------------
c
        subroutine MkV_A1 (Ve,V1,dimo2,no)
c
c       this routine do:
c       Ve(ij,u,v) <<- V1(iu|jv)
c
        implicit none
        integer dimo2,no
        real*8 Ve(1:dimo2,1:no,1:no)
        real*8 V1(1:dimo2,1:dimo2)
c
c       help variables
        integer i,j,ij,u,v,iu,jv
c
c
        do v=1,no
        do u=1,no
c
          ij=0
          do i=1,no
          do j=1,i
            ij=ij+1
c
            if (i.gt.u) then
            iu=(i-1)*i/2+u
            else
            iu=(u-1)*u/2+i
            end if
c
            if (j.gt.v) then
            jv=(j-1)*j/2+v
            else
            jv=(v-1)*v/2+j
            end if
c
            Ve(ij,u,v)=Ve(ij,u,v)+V1(iu,jv)
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
c       --------------------------------
c
        subroutine MkV_A4 (Vp,V,dimb,dima,no,dimij)
c
c       this routine do:
c       Vp(a,b,ij) <- (ai|bj) from V(b,j,a,i)
c
        implicit none
        integer dima,dimb,no,dimij
        real*8 Vp(1:dima,1:dimb,1:dimij)
        real*8 V(1:dimb,1:no,1:dima,1:no)
c
c       help variables
        integer i,j,ij,a,b
c
        ij=0
        do i=1,no
        do j=1,i
          ij=ij+1
          do b=1,dimb
          do a=1,dima
            Vp(a,b,ij)=V(b,j,a,i)
          end do
          end do
        end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine ExV_X41 (Vp,V,dimab,no)
c
c       this routine do:
c       Vp(a_b,ij) <- V(a_b,i,j) for i>=j
c
        implicit none
        integer dimab,no
        real*8 Vp(1:dimab,1:no*(no+1)/2)
        real*8 V(1:dimab,1:no,1:no)
c
c       help variables
        integer i,j,ij,ab
c
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
          do ab=1,dimab
            Vp(ab,ij)=V(ab,i,j)
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
        subroutine ExV_X43 (Vp,V,dimab,no)
c
c       this routine do:
c       Vp(a_b,ij) <- V(a_b,j,i) for i>=j
c
        implicit none
        integer dimab,no
        real*8 Vp(1:dimab,1:no*(no+1)/2)
        real*8 V(1:dimab,1:no,1:no)
c
c       help variables
        integer i,j,ij,ab
c
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
          do ab=1,dimab
            Vp(ab,ij)=V(ab,j,i)
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
        subroutine ExV_X42 (Vp,V,dimab,no)
c
c       this routine do:
c       Vp(a,b,i) <- V(a,b,i,i)
c
        implicit none
        integer dimab,no
        real*8 Vp(1:dimab,1:no)
        real*8 V(1:dimab,1:no,1:no)
c
c       help variables
        integer i,ab
c
        do i=1,no
          do ab=1,dimab
            Vp(ab,i)=V(ab,i,i)
          end do
        end do
c
c
        return
        end
c
c       --------------------------------
c
        subroutine ExA_X4 (A,Ap,no)
c
c       this routine do:
c       Ap(i,u,v) <- A(ii,u,v)
c
        implicit none
        integer no
        real*8 Ap(1:no,1:no*no)
        real*8 A(1:no*(no+1)/2,1:no*no)
c
c       help variables
        integer uv,i,ii
c
        do uv=1,no*no
        do i=1,no
          ii=i*(i+1)/2
            Ap(i,uv)=A(ii,uv)
          end do
        end do
c
c
        return
        end
c
c       --------------------------------
c
        subroutine MkV_T18 (Va,V,dima,no)
c
c       this routine do:
c        Va(a,j,i,u) = - [2V(ai|ju)-V(aj|iu)]
c        from V(o_a,P,Q,a)
c
c        N.B. Kvajto odflaknute, ozaj ze hnus, treba sa zamysliet
c        ci je nutva permutacia s a-ckom na konci - bod QK2.2 @@
c
c
        implicit none
        integer dima,no
        real*8 Va(1:dima,1:no,1:no,1:no)
        real*8 V(1:no,1:no,1:no,1:dima)
c
c       help variables
        integer a,j,i,u
c
        do u=1,no
          do i=1,no
            do j=1,no
              do a=1,dima
                Va(a,j,i,u)=V(j,i,u,a)-2.0d0*V(i,j,u,a)
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
        subroutine AdT_T17 (T1,Q,dimat,dimaq,no,adda,f)
c
c       this routine do:
c       T1(a,i) <<- f . Q(a',i)
c
        implicit none
        integer dimat,dimaq,no,adda
        real*8 T1(1:dimat,1:no)
        real*8 Q(1:dimaq,no)
        real*8 f
c
c       help variables
        integer i,a
c
        do i=1,no
          do a=1,dimaq
            t1(adda+a,i)=t1(adda+a,i)+f*q(a,i)
          end do
        end do
c
        return
        end
c
c       --------------------------------
c
        subroutine MkT_T15 (Tp,T2,T11,T12,dimbe,dima,no)
c
c       this routine do:
c       Tp(be',u,a',i) <- 2 t2(be,a,u,i)
c                       - t2(be,a,i,u)+t12(a,u).t11(be,i)/2
c        N.B. mozno sa to da este vylepsit
c
        implicit none
        integer dimbe,dima,no
        real*8 Tp(1:dimbe,1:no,1:dima,1:no)
        real*8 T2(1:dimbe,1:dima,1:no,1:no)
        real*8 T11(1:dimbe,1:no)
        real*8 T12(1:dima,1:no)
c
c       help variables
        integer be,a,u,i
        real*8 c1

        do i=1,no
          do a=1,dima
            do u=1,no
            c1=T12(a,u)
              do be=1,dimbe
                 Tp(be,u,a,i)=2.0d0*T2(be,a,u,i)
     c                       -T2(be,a,i,u)+c1*T11(be,i)
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
        subroutine DfH_Hvv1 (Hvv,Fvv,nv,dimbe,addbe)
c
c        this routine do:
c        Hvv(a,be') <- Fvv(a,be)
c
        implicit none
        integer dimbe,nv,addbe
        real*8 Hvv(1:nv,1:dimbe)
        real*8 Fvv(1:nv,1:nv)
c
c        help variables
        integer a,be,bev
c
        do be=1,dimbe
          bev=be+addbe
          do a=1,nv
            Hvv(a,be)=Fvv(a,bev)
          end do
        end do
c
c
        return
        end
c
c       --------------------------------
c
        subroutine ExH_T13 (V,Hvv,dimbe,addbe,nv)
c
c        this routine do:
c       Extract V1(a,be') <- Hvv(a,be)
c        Hvv(a,be') <- Fvv(a,be)
c
        implicit none
        integer dimbe,nv,addbe
        real*8 Hvv(1:nv,1:nv)
        real*8 V(1:nv,1:dimbe)
c
c        help variables
        integer a,be,bev
c
        do be=1,dimbe
          bev=be+addbe
          do a=1,nv
            V(a,be)=Hvv(a,bev)
          end do
        end do
c
c
        return
        end
c
c       --------------------------------
c
        subroutine Ext_Aex (Aex,VV,no)
c
c        this routine def:
c        VV(i,u,v,j) <- Aex(ij,u,v)
c        for Aex(i,j,u,v) = Aex(j,i,v,u), Aex stored only for i>=j
c
        implicit none
        integer no
        real*8 VV(1:no,1:no,1:no,1:no)
        real*8 Aex(1:no*(no+1)/2,1:no,1:no)
c
c        help variables
        integer i,j,ij,u,v
c
c
        do u=1,no
        do v=1,no
          ij=0
          do i=1,no
          do j=1,i
          ij=ij+1
c
            VV(i,u,v,j)=Aex(ij,u,v)
            VV(j,v,u,i)=Aex(ij,u,v)
c
          end do
          end do
        end do
        end do
c
c
        return
        end
