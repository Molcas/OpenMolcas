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
c       This package contains following files:
c
c       DefParo2v4
c       DistMemo2v4
c       GetTau
c          GetTauHlp1
c          GetTauHlp2
c       GetChV
c         getChVHlp1
c         getChVHlp2
c         getChVHlp3
c         getChVHlp4
c       VanishT2n
c       ExtractM
c       MakeT2p
c       MakeT2m
c         makeT2pHlp1
c         makeT2pHlp2
c         makeT2pHlp3
c         makeT2ptHlp1
c         makeT2ptHlp2
c         makeT2ptHlp3
c       MakeWw
c         makeWwHlp1
c         makeWwHlp2
c         makeWwHlp3
c         makeWwHlp4
c       MakeT2pd
c         makeT2pdHlp
c       MakeWwd
c         makeWwdHlp1
c         makeWwdHlp2
c        Calc_addSG
c        Calc_addG
c
c        Vint approach routines
c        Mk_T1t
c        ReaW3
c          MkNameV3
c          ReaW3hlp1
c          ReaW3hlp2
c          ReaW3hlp3
c        ReaW4
c          DefW4abcd
c          DefW4abdc
c          DefW4bacd
c          DefW4badc
c          DefW4cdab
c          DefW4dcab
c          DefW4cdba
c          DefW4dcba
c          MkNameV4
c
c        Xo2v4ctl (kokotske meno)
c        InsReqo2v4
c          InsReaW3
c          InsReaW4
c
c
c       ------------------------------------
c
        subroutine DefParo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                         mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
c
c       This routine do:
c       define parameters in o2v4.fh using NaGrp,NbeGrp,NaSGrp,NbeSgrp
c
c       I/O parameter description:
c       NaGrp    - # of groups in a set (I)
c       NbeGrp   - # of groups in be set (I)
c       NaSGrp   - # of subgroups in each (a)' group (I)
c       NbeSGrp  - # of subgroups in each (be)' group (I)
c       mdGrpa   - # maximal dimension od (a)' group (O)
c       mdGrpbe  - # maximal dimension od (be)' group (O)
c       mdSGrpa  - # maximal dimension od (a)" subgroup (O)
c       mdSGrpbe - # maximal dimension od (be)" subgroup (O)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
c
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
c
c       help variables
c
        real*8 rdim
        integer i,j,ij
c
c
c1      define parameters of Groups of a(b) set
c
        rdim=1.0d0*nv/(1.0d0*NaGrp)
c
        Grpalow(1)=1
        Grpaup(1)=NaSGrp
c
        do i=1,NaGrp
c
           if (i.eq.1) then
             UpGrpa(i)=int(rdim*i)
             LowGrpa(i)=1
           else if (i.eq.NaGrp) then
             UpGrpa(i)=nv
             LowGrpa(i)=UpGrpa(i-1)+1
           else
             UpGrpa(i)=int(rdim*i)
             LowGrpa(i)=UpGrpa(i-1)+1
           end if
c
           DimGrpa(i)=(UpGrpa(i)-LowGrpa(i))+1
c
           if (i.gt.1) then
             Grpalow(i)=Grpalow(i-1)+NaSGrp
             Grpaup(i)=Grpaup(i-1)+NaSGrp
           end if
c
        end do
c
c
c2      define parameters of Groups of be(ga) set
c
        rdim=1.0d0*nv/(1.0d0*NbeGrp)
c
        Grpbelow(1)=1
        Grpbeup(1)=NbeSGrp
c
        do i=1,NbeGrp
c
           if (i.eq.1) then
             UpGrpbe(i)=int(rdim*i)
             LowGrpbe(i)=1
           else if (i.eq.NbeGrp) then
             UpGrpbe(i)=nv
             LowGrpbe(i)=UpGrpbe(i-1)+1
           else
             UpGrpbe(i)=int(rdim*i)
             LowGrpbe(i)=UpGrpbe(i-1)+1
           end if
c
           DimGrpbe(i)=(UpGrpbe(i)-LowGrpbe(i))+1
c
           if (i.gt.1) then
             Grpbelow(i)=Grpbelow(i-1)+NbeSGrp
             Grpbeup(i)=Grpbeup(i-1)+NbeSGrp
           end if
c
        end do
c
c
c3      define parameters of SubGroups of a(b) set
c
        ij=0
        do j=1,NaGrp

        rdim=1.0d0*DimGrpa(j)/(1.0d0*NaSGrp)
c
        do i=1,NaSGrp
        ij=ij+1
c
           if (i.eq.1) then
             UpSGrpa(ij)=(LowGrpa(j)-1)+int(rdim*i)
             LowSGrpa(ij)=LowGrpa(j)
           else if (i.eq.NaSGrp) then
             UpSGrpa(ij)=UpGrpa(j)
             LowSGrpa(ij)=UpSGrpa(ij-1)+1
           else
             UpSGrpa(ij)=(LowGrpa(j)-1)+int(rdim*i)
             LowSGrpa(ij)=UpSGrpa(ij-1)+1
           end if
c
          DimSGrpa(ij)=(UpSGrpa(ij)-LowSGrpa(ij))+1
c
        end do
        end do
c
c
c4      define parameters of SubGroups of be(ga) set
c
        ij=0
        do j=1,NbeGrp

        rdim=1.0d0*DimGrpbe(j)/(1.0d0*NbeSGrp)
c
        do i=1,NbeSGrp
        ij=ij+1
c
           if (i.eq.1) then
             UpSGrpbe(ij)=(LowGrpbe(j)-1)+int(rdim*i)
             LowSGrpbe(ij)=LowGrpbe(j)
           else if (i.eq.NbeSGrp) then
             UpSGrpbe(ij)=UpGrpbe(j)
             LowSGrpbe(ij)=UpSGrpbe(ij-1)+1
           else
             UpSGrpbe(ij)=(LowGrpbe(j)-1)+int(rdim*i)
             LowSGrpbe(ij)=UpSGrpbe(ij-1)+1
           end if
c
          DimSGrpbe(ij)=(UpSGrpbe(ij)-LowSGrpbe(ij))+1
c
        end do
        end do
c
c
c5      find maximal dimensions mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
c
        mdGrpa=DimGrpa(1)
        do i=1,NaGrp
          if (DimGrpa(i).gt.mdGrpa) then
          mdGrpa=DimGrpa(i)
          end if
        end do
c
        mdGrpbe=DimGrpbe(1)
        do i=1,NbeGrp
          if (DimGrpbe(i).gt.mdGrpbe) then
          mdGrpbe=DimGrpbe(i)
          end if
        end do
c
        mdSGrpa=DimSGrpa(1)
        do ij=1,NaGrp*NaSGrp
          if (DimSGrpa(ij).gt.mdSGrpa) then
          mdSGrpa=DimSGrpa(ij)
          end if
        end do
c
        mdSGrpbe=DimSGrpbe(1)
        do ij=1,NbeGrp*NbeSGrp
          if (DimSGrpbe(ij).gt.mdSGrpbe) then
          mdSGrpbe=DimSGrpbe(ij)
          end if
        end do
c
c
c6      def L2Names
c
        do i=1,NaGrp
        do j=1,NbeGrp
          call DefParo3v3Hlp1(i,j,'L2',L2Name(i,j))
        end do
        end do
c
c7      def Tmp3Names
c
        do i=1,MaxSGrp
        do j=1,MaxSGrp
          call DefParo3v3Hlp1(i,j,'X3',Tmp3Name(i,j))
        end do
        end do
c
c
        return
        end
c
c       ------------------------------------
c
        subroutine DistMemo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                 mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,
     c                 PossTau,PossT2n1,PossT2n2,PossT2w,
     c                 PossL11,PossL12,
     c                 PossL21,PossL22,PossL23,PossL24,PossL2W,
     c                 PossH1,PossH2,
     c                 PossM1,PossM2,PossW1,PossW2,PossWw,PossWx,
     c                 PossT,NL2)
c
c       This routine do:
c       define initial possitions of T,L,M and W arrays,
c       described in o2v4ctl routine
c
c
c       I/O parameter description:
c       NaGrp    - # of groups in a set (I)
c       NbeGrp   - # of groups in be set (I)
c       NaSGrp   - # of subgroups in each (a)' group (I)
c       NbeSGrp  - # of subgroups in each (be)' group (I)
c       mdGrpa   - # maximal dimension od (a)' group (I)
c       mdGrpbe  - # maximal dimension od (be)' group (I)
c       mdSGrpa  - # maximal dimension od (a)" subgroup (I)
c       mdSGrpbe - # maximal dimension od (be)" subgroup (I)
c       PossX    - initial possitinos of arrays (O-all)
c       PossT    - next free possition (O)
c       NL2      - # of L2 vectors (O)
c
        implicit none
#include "chcc1.fh"
c
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
        integer PossTau,PossT2n1,PossT2n2,PossT2w
        integer PossL11,PossL12
        integer PossH1,PossH2
        integer PossL21,PossL22,PossL23,PossL24,PossL2W
        integer PossM1,PossM2,PossW1,PossW2,PossWw,PossWx
        integer PossT,NL2
c
c       help variables
        integer length,lenab,lenbega
c
c1      Tau
c
        PossTau=possT
        if (NaGrp.eq.1) then
          length=no*no*nv*(nv+1)/2
        else
          length=no*no*mdGrpa*mdGrpa
        end if
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,99) 'DM Tau',PossTau
        end if
c
c2      T2n files
c
        possT2n1=PossT
        if ((NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
          length=no*(no+1)*nv*(nv+1)/4
        else
          length=no*(no+1)*mdSGrpbe*mdSGrpbe/2
        end if
        PossT=PossT+length
c
        possT2n2=PossT
        if ((NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
          length=no*(no-1)*nv*(nv-1)/4
        else
          length=no*(no-1)*mdSGrpbe*mdSGrpbe/2
        end if
        PossT=PossT+length
c
        possT2w=PossT
        if ((NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
          length=no*(no+1)*nv*(nv+1)/4
        else
          length=no*(no+1)*mdSGrpbe*mdSGrpbe/2
        end if
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,99) 'DM T2 ',PossT2n1,PossT2n2,PossT2w
        end if
c
c       L1 files
c
        if (intkey.eq.1) then
c        integral based
          length=0
        else
c        choleski based
          length=nc*mdGrpa*no
        end if
c
        if (NaGrp.eq.1) then
          PossL11=PossT
          PossL12=PossT
          PossT=PossT+length
        else
          PossL11=PossT
          PossT=PossT+length
          PossL12=PossT
          PossT=PossT+length
        end if
        if (printkey.ge.10) then
        write (6,99) 'DM L1 ',PossL11,PossL12
        end if
c
c       L2 files
c
        if (intkey.eq.1) then
c        integral based
          length=0
        else
c        choleski based
          length=nc*mdGrpa*mdGrpbe
        end if
c
        if (NaGrp.eq.1) then
          if (NbeGrp.eq.1) then
c         All L2 are identical
            PossL21=PossT
            PossL22=PossT
            PossL23=PossT
            PossL24=PossT
            PossT=PossT+length
            NL2=1
          else
c         L21=L23 and L22=L24
            PossL21=PossT
            PossL23=PossT
            PossT=PossT+length
            PossL22=PossT
            PossL24=PossT
            PossT=PossT+length
            NL2=2
          end if
        else
          if (NbeGrp.eq.1) then
c         L21=L22 and L23=L24
            PossL21=PossT
            PossL22=PossT
            PossT=PossT+length
            PossL23=PossT
            PossL24=PossT
            PossT=PossT+length
            NL2=2
          else
c         all L2 files are different
            PossL21=PossT
            PossT=PossT+length
            PossL22=PossT
            PossT=PossT+length
            PossL23=PossT
            PossT=PossT+length
            PossL24=PossT
            PossT=PossT+length
            NL2=4
          end if
        end if
c
        PossL2W=PossT
        if (intkey.eq.1) then
c        integral based
          length=0
        else
c        choleski based
        if ((NaGrp.eq.1).and.(NbeGrp.eq.1)) then
          length=nc*nv*(nv+1)/2
        else
          length=nc*mdGrpa*mdGrpbe
        end if
        if ((no*mdGrpbe).gt.length) then
          length=no*mdGrpbe
        end if
        if (nc*no*mdGrpa.gt.length) then
          length=nc*no*mdGrpa
        end if
        end if
c
        PossT=PossT+length
c
        if (printkey.ge.10) then
        write (6,99) 'DM L2 ',PossL21,PossL22,PossL23,PossL24,
     c                       PossL2W
        end if
c
c       H files
c
        length=no*mdSGrpbe
        if (intkey.eq.1) then
          PossH1=PossT
          PossT=PossT+length
          PossH2=PossT
          PossT=PossT+length
        else
          PossH1=PossT
          PossH2=PossT
        end if
        if (printkey.ge.10) then
        write (6,99) 'DM H  ',PossH1,PossH2
        end if
c
c       M files
c
        length=nc*mdSGrpa*mdSGrpbe
c
        if ((NaGrp.eq.1).and.(NaSGrp.eq.1).and.
     c      (NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
c       M1=M2
          PossM1=PossT
          PossM2=PossT
          PossT=PossT+length
        else
c       M files are different
          PossM1=PossT
          PossT=PossT+length
          PossM2=PossT
          PossT=PossT+length
        end if
        if (printkey.ge.10) then
        write (6,99) 'DM M  ',PossM1,PossM2
        end if
c
c       W1,W2 files
c
        if (NaGrp*NaSGrp.eq.1) then
c       W1 and W2 are identical
          length=mdSGrpa*mdSGrpbe*mdSGrpa*mdSGrpbe
          if ((no.gt.mdSGrpbe).and.(intkey.eq.1)) then
          length=mdSGrpa*mdSGrpbe*mdSGrpa*no
          end if
          PossW1=PossT
          PossW2=PossT
          PossT=PossT+length
        else
c       W1 and W2 are different
          length=mdSGrpa*mdSGrpbe*mdSGrpa*mdSGrpbe
          if ((no.gt.mdSGrpbe).and.(intkey.eq.1)) then
          length=mdSGrpa*mdSGrpbe*mdSGrpa*no
          end if
          PossW1=PossT
          PossT=PossT+length
          PossW2=PossT
          PossT=PossT+length
        end if
c
c       Ww file
c
        PossWw=PossT
c
        if ((NaGrp.eq.1).and.(NaSGrp.eq.1)) then
          lenab=nv*(nv+1)/2
        else
          lenab=mdSGrpa*mdSGrpa
        end if
        if ((NbeGrp.eq.1).and.(NbeSGrp.eq.1)) then
          lenbega=nv*(nv+1)/2
        else
          lenbega=mdSGrpbe*mdSGrpbe
        end if
c
        length=lenab*lenbega
        if (intkey.eq.1) then
          length=mdSGrpa*mdSGrpbe*mdSGrpa*mdSGrpbe
          if (no.gt.mdSGrpbe) then
          length=mdSGrpa*mdSGrpbe*mdSGrpa*no
          end if
        end if
        PossT=PossT+length
c
c        Wx file
c
        PossWx=PossT
        length=mdSGrpa*mdSGrpbe*mdSGrpa*no
        if (intkey.eq.0) then
          length=0
        end if
        PossT=PossT+length
c
        if (printkey.ge.10) then
        write (6,99) 'DM W  ',PossW1,PossW2,PossWw,PossWx
c
        write (6,99) 'PossT ',PossT
        end if
c
99        format (a7,10(i10,1x))
c
        return
        end
c
c       ------------------------------------
c
        subroutine GetTau (Tau,T1,
     c                     aGrp,bGrp,dima,dimb,adda,addb,lunTau)
c
c       this routine do:
c       1) read the block of T2((a,b)',ij) from T2Name(aGrp,bGrp)
c        2) Make Tau ((a,b)',ij) in T2((a,b)',ij) array
c
c       I/O parameter description:
c       Tau      - array for Tau((a,b),ij) (O)
c       aGrp     - group of a index
c       bGrp     - group of b index
c       dima     - dimension group of a' index
c       dimb     - dimension group of b' index
c       adda     - shift of a' in full virtual space
c       addb     - shift of b' in full virtual space
c       lunTau   - Lun of opened file, where Tau is stored
c
c
        implicit none
#include "chcc1.fh"
#include "chcc_files.fh"
c
        real*8 Tau(1)
        real*8 T1(1)
        integer aGrp,bGrp,dima,dimb,adda,addb,lunTau
c
c       help variables
        integer length
c
c
c1      def legth
c
        if (aGrp.eq.bGrp) then
c       groups of a and b are equal, reading for a'>=b'
          length=no*no*Dima*(Dima+1)/2
        else
c       aGrp>bGrp, reading for a',b' in given groups
          length=no*no*Dima*Dimb
        end if
c
c
c2      read block of T2 amplitudes
c
        call GetX (Tau(1),length,LunTau,T2Name(aGrp,bGrp),1,1)
c
c
c3        make Tau
c
        if (aGrp.ne.bGrp) then
           call GetTauHlp1 (Tau(1),T1(1),dima,dimb,adda,addb,no,nv)
        else
           call GetTauHlp2 (Tau(1),T1(1),dima,adda,no,nv)
        end if


        return
        end
c
c       ------------------------------------
c
        subroutine GetTauHlp1 (Tau,T1,dima,dimb,adda,addb,no,nv)
c
c        Make Tau for aGrp.ne.bGrp
c
        implicit none
        integer dima,dimb,adda,addb,no,nv
        real*8 Tau(1:dima,1:dimb,1:no,1:no)
        real*8 T1(1:nv,1:no)
c
c       help variables
        integer i,j,a,b
        real*8 c
c
c
        do j=1,no
          do i=1,no
            do b=1,dimb
              c=t1(addb+b,j)
              do a=1,dima
                Tau(a,b,i,j)=Tau(a,b,i,j)+c*t1(adda+a,i)
              end do
            end do
          end do
        end do
c
c
        return
        end
c
c       ------------------------------------
c
         subroutine GetTauHlp2 (Tau,T1,dima,adda,no,nv)
c
c        Make Tau for aGrp.eq.bGrp
c
        implicit none
         integer dima,adda,no,nv
        real*8 Tau(1:dima*(dima+1)/2,1:no,1:no)
        real*8 T1(1:nv,1:no)
c
c       help variables
        integer i,j,a,b,ab
        real*8 c
c
c
        do j=1,no
          do i=1,no
            ab=0
            do a=1,dima
              c=t1(adda+a,i)
              do b=1,a
                ab=ab+1
                Tau(ab,i,j)=Tau(ab,i,j)+c*t1(adda+b,j)
              end do
            end do
          end do
        end do
c
c
        return
        end
c
c       ------------------------------------
c
        subroutine GetChV (wrk,wrksize,
     c                     aGrp,bGrp,beGrp,gaGrp,NL2,L2Status,
     c                     pL21,pL22,pL23,pL24,pL2W,
     c                     PossL11,PossL12,LunAux)
c
c       this rotune do:
c       read L2 files from disc file
c*      L21 (m,a',be')
c       L22 (m,a',ga')
c       L23 (m,b',be')
c       L24 (m,b',ga')
c        and completed it by the term
c        L2x (m,c',de') <<- L2(m,c',i) . T1(T)(de',i)
c
c       Routine define pL2x, which are values, which show, which index
c       i (correspondig to L2Status(i,1-3) correspond to ginev L2x.
c       Array L2Status (neaningful only for i=1,NL2) stores informations
c       about really loaded L2's in memory. If some of the blocks were
c       already defined in previous spep(s), they are not readed repetidly.
c
c       N.B.1 Rutina sa da spravit aj trivialne, ze sa pokazde vsetky
c       L21-4 nacitaju, toto riesenie je menej prehladne, ale setri IO
c       N.B.2 Efektivitu ovplyvnuje hlavne, ako vybera rutine getChVHlp2
c       kam sa ma nove L2 nacitat, najma ak ma vybrat, ktory s
c       existujucich rekordov premaze
c
c       I/O parameter description:
c       aGrp, bGrp   - groups of a,b (I)
c       beGrp, GaGrp - groups of be,ga (I)
c       NL2          - real number of L2 declared in memory(1,2or4) (I)
c       L2Status    - L2 status matrix (I/O)
c                     L2Status(i,1) - cGrp
c                     L2Status(i,2) - deGrp
c                     L2Status(i,3) - Poss
c       pL2x        - possition of L2x in L2Status (i.e. index i
c                     in L2Status(i,1-3) for given L2x (O)
c        PossL11     - Possition of L1(m,a',i) (I)
c        PossL12     - Possition of L1(m,b',i) (I)
c       LunAux      - lun for auxiliary reading (I)
c
c
        implicit none
#include "wrk.fh"
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,bGrp,beGrp,gaGrp,NL2
        integer pL21,pL22,pL23,pL24,pL2W,PossL11,PossL12,LunAux
        integer L2Status(1:4,1:3)
c
c       help variables
        integer Used(1:4)
        integer WhatNeedToRead(1:4,1:2)
        integer HowMany
        integer Which(1:4)
        integer Kde(1:4)
        integer i,kery,kam,yes,ToDo
        integer cGrp,deGrp,dimc,dimde,addde
c
c
c1      Define
c         a) what need to be read (WhatNeedToRead)
c         b) How many different L2 need to be defined (HowMany)
c         c) How L21-4 is mapped to the unique ones, that need to
c            be read (Which)
c
        if (aGrp.eq.bGrp) then
          if (beGrp.eq.gaGrp) then
c         case L21=L22=L23=L24
            HowMany=1
            WhatNeedToRead(1,1)=aGrp
            WhatNeedToRead(1,2)=beGrp
            Which(1)=1
            Which(2)=1
            Which(3)=1
            Which(4)=1
          else
c         case L21=L23, L22=L24
            HowMany=2
            WhatNeedToRead(1,1)=aGrp
            WhatNeedToRead(1,2)=beGrp
            WhatNeedToRead(2,1)=aGrp
            WhatNeedToRead(2,2)=gaGrp
            Which(1)=1
            Which(2)=2
            Which(3)=1
            Which(4)=2
          end if
        else
          if (beGrp.eq.gaGrp) then
c         case L21=L22,L23=L24
            HowMany=2
            WhatNeedToRead(1,1)=aGrp
            WhatNeedToRead(1,2)=beGrp
            WhatNeedToRead(2,1)=bGrp
            WhatNeedToRead(2,2)=beGrp
            Which(1)=1
            Which(2)=1
            Which(3)=2
            Which(4)=2
          else
c         all L2s are different
            HowMany=4
            WhatNeedToRead(1,1)=aGrp
            WhatNeedToRead(1,2)=beGrp
            WhatNeedToRead(2,1)=aGrp
            WhatNeedToRead(2,2)=gaGrp
            WhatNeedToRead(3,1)=bGrp
            WhatNeedToRead(3,2)=beGrp
            WhatNeedToRead(4,1)=bGrp
            WhatNeedToRead(4,2)=gaGrp
            Which(1)=1
            Which(2)=2
            Which(3)=3
            Which(4)=4
          end if
        end if
c
c
c2      Load neccesarry L2 from disk
c
c2.1    SetUp
        do i=1,NL2
          Kde(i)=0
          Used(i)=0
        end do
c
1       ToDo=HowMany
c
c2.2    Test, which ones of those, that need to be read are
c       already loaded and where (Kde) and find, how many of
c       L2 need to be readed really (ToDo)
        do i=1,HowMany
          call getChVHlp1
     c    (WhatNeedToRead(i,1),WhatNeedToRead(i,2),yes,NL2,L2Status)
          Kde(i)=yes
          if (yes.ne.0) then
            used(yes)=1
            ToDo=ToDo-1
          end if
        end do
c
c2.3    if there is atleast one L2 to read, read only one
c       and run 2.2 section again
        if (ToDo.gt.0) then
c
c         zistit kere L2 treba nacitat (kery)
          kery=0
          do i=1,HowMany
            if (Kde(i).eq.0) then
              kery=i
            end if
          end do
c
c         found, where to read it (kam)
          call getChVHlp2 (L2Status,NL2,used,kam)
c
c         read L2 (kery) into (kam) (+ expand, f needed)
          L2Status(kam,1)=WhatNeedToRead(kery,1)
          L2Status(kam,2)=WhatNeedToRead(kery,2)
          cGrp=L2Status(kam,1)
          deGrp=L2Status(kam,2)
          call getChVHlp3 (wrk(L2Status(kam,3)),wrk(pL2W),
     c                     L2Status(kam,1),L2Status(kam,2),LunAux)
c
c          Create L2W(i,de') <- T1(de,i)
          dimde=DimGrpbe(deGrp)
          addde=0
          do i=1,deGrp-1
            addde=addde+DimGrpbe(i)
          end do
          call getChVHlp4 (wrk(pL2W),wrk(PossT1o),dimde,addde)
c
c          upgrade L2(m,c',de') <<- - L1(m,c',i). L2W(i,de')
          dimc=DimGrpa(cGrp)
          if (cGrp.eq.aGrp) then
            call mc0c2a3b (nc*dimc,no,no,dimde,nc*dimc,dimde,
     c                     nc*dimc,no,dimde,
     c                     wrk(PossL11),wrk(pL2W),wrk(L2Status(kam,3)))
          else if (cGrp.eq.bGrp) then
            call mc0c2a3b (nc*dimc,no,no,dimde,nc*dimc,dimde,
     c                     nc*dimc,no,dimde,
     c                     wrk(PossL12),wrk(pL2W),wrk(L2Status(kam,3)))
          else
            write (6,*) ' Nieje dobre, c nieje ani a ani b :-( Ch. K.'
            call Abend()
          end if
c
          used(kam)=1
c
          goto 1
        end if
c
c
c3      At this point all neccesarry L2 are read - def pL2x
c
        i=Which(1)
        pL21=Kde(i)
        i=Which(2)
        pL22=Kde(i)
        i=Which(3)
        pL23=Kde(i)
        i=Which(4)
        pL24=Kde(i)
c
c
        return
        end
c
c       ---------------------
c
        subroutine GetChVHlp1 (cGrp,deGrp,yes,NL2,L2Status)
c
c       this routine do:
c       check, if Choleski vector block cGrp,deGrp is actually situated
c       in the memory as one of L2x
c
c       descrition of parameters:
c       cGrp, deGrp - required Group of c, delta in L2(m,c',de') (I)
c       yes    -  0 - not loaded  (O)
c                 x - loaded, value of L2x (x=1-4)
c       NL2         - Number of L2 arrays really reserved in memory (I)
c       L2Status    - L2 status matrix (I)
c                     L2Status(i,1) - cGrp
c                     L2Status(i,2) - deGrp
c                     L2Status(i,3) - Poss
c
        implicit none
        integer cGrp,deGrp,yes,NL2
        integer L2Status(1:4,1:3)
c
c       help variables
        integer i
c
        yes=0
        do i=1,NL2
        if ((cGrp.eq.L2Status(i,1)).and.(deGrp.eq.L2Status(i,2))) then
        yes=i
        end if
        end do
c
        return
        end
c
c       ---------------------
c
        subroutine GetChVHlp2 (L2Status,NL2,used,kam)
c
c       this routine do:
c       Define, which one of unused L2 arrays allocated in memory
c       will be used for placing newly needed L2x
c       Criteria of present approach:
c       1) First those possitions that are allocated, but never used yet
c       2) if no 1), take first unused in this step, according L2Status
c          N.B. 2) druhy krok moze byt eventuelne vylepseny, zatial
c               takto odflaknute
c
c       description ov variables:
c       L2Status  - L2 Status array (I)
c       NL2       - number of L2 arrays reserver in memory (I)
c       used      - array indicating, which L2x are already used
c                   in this step (those cannot be touched) (I)
c       kam       - index of possition, where new L2 can be placed (O)
c
        implicit none
        integer L2Status(1:4,1:3)
        integer used(1:4)
        integer NL2,kam
c
c       help variables
        integer i
c
c
c1      search, if there are never used possitions
c
        do i=1,NL2
          if (L2Status(i,1).eq.0) then
          kam=i
          return
          end if
        end do
c
c2      find first unused in this step
c
        do i=1,NL2
          if (used(i).eq.0) then
          kam=i
          return
          end if
        end do
c
c       Jaj nieje dobre ak sme sa dostali az sem
        write (6,*) ' Sorry fish getChVHlp2 '
        call Abend()
c
        return
        end
c
c       ---------------------
c
        subroutine GetChVHlp3 (L2,Tmp,cGrp,deGrp,LunAux)
c
c       this routine do:
c       read L2(m,c'de') into L2 from disk file
c       structure of files: for each c',de' one file with
c       name L2Name(cGrp,deGrp)
c       @ citanie zatial odflaknute
c
c       parameter description:
c       L2     - Array for L2 (O)
c       Tmp    - Temporary array of L2 size, used for mapping, if needed
c       xGrp   - Groups of c, delta (I)
c       LunAux - lun for auxiliary reading (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
c
        real*8 L2(1)
        real*8 Tmp(1)
        integer cGrp,deGrp,LunAux
        character*6 LunName
c
c       help variables
        integer length
c
c        nacitanie (+expanzia, ak treba)
        if (cGrp.gt.deGrp) then
          LunName=L2Name(cGrp,deGrp)
          length=nc*DimGrpa(cGrp)*DimGrpbe(deGrp)
          call GetX (L2(1),length,LunAux,LunName,1,1)
        else if (cGrp.eq.deGrp) then
          LunName=L2Name(cGrp,deGrp)
          length=nc*DimGrpa(cGrp)*(DimGrpbe(deGrp)+1)/2
          call GetX (Tmp(1),length,LunAux,LunName,1,1)
          length=DimGrpa(cGrp)*(DimGrpbe(deGrp)+1)/2
          call Exp1 (Tmp(1),L2(1),nc,length,DimGrpa(cGrp))
        else
          LunName=L2Name(deGrp,cGrp)
          length=nc*DimGrpa(cGrp)*DimGrpbe(deGrp)
          call GetX (Tmp(1),length,LunAux,LunName,1,1)
          call Map3_132 (Tmp(1),L2(1),nc,DimGrpbe(deGrp),DimGrpa(cGrp))
        end if
c
        return
        end
c
c       ---------------------
c
        subroutine GetChVHlp4 (H,T1,dima,adda)
c
c       this routine do:
c       Extract H(i,a') <- T1(a,i) for given aGrp
c
c       parameter description:
c       H       - Output file (O)
c       T1      - T1 amplitudes (I)
c       dima    - dimension of given Group (I)
c       adda    - shift of a' in full a set (I)
c
c       N.B. Kvajt odflaknute, je to ExtT1 len transponovane
c
        implicit none
#include "chcc1.fh"
        integer dima,adda
        real*8 T1(1:nv,1:no)
        real*8 H(1:no,1:dima)
c
c       help variables
        integer a,aa,i
c
        do a=1,dima
        aa=adda+a
          do i=1,no
            H(i,a)=T1(aa,i)
          end do
        end do
c
        return
        end
c
c       ------------------------------------
c
        subroutine VanishT2n (T2n1,T2n2,beSGrp,gaSGrp)
c
c       this routine do:
c       vanish space for T2n = T2n(-(+)) (i>(>=)j,(be>(>=)ga)")
c
c       parameter description:
c       T2nx    - arrays for T2+- (O)
c       xSGrp   - SubGroups of be,ga (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 T2n1(1)
        real*8 T2n2(1)
        integer beSGrp,gaSGrp
c
c       help variables
        integer length1,length2
c
c1      calc legths
        if (beSGrp.eq.gaSGrp) then
          length1=no*(no+1)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/4
          length2=no*(no-1)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)-1)/4
        else
          length1=no*(no+1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)/2
          length2=no*(no-1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)/2
        end if
c
c2      vanish
        call mv0zero (length1,length1,T2n1(1))
        call mv0zero (length2,length2,T2n2(1))
c
        return
        end
c
c       ------------------------------------
c
        subroutine GetT2n (T2n1,T2n2,beSGrp,gaSGrp,LunAux)
c
c       this routine do:
c       Read T2n = T2n(-(+)) (i>(>=)j,(be>(>=)ga)")
c        from Tmp3Name(be",ga")
c
c       parameter description:
c       T2nx    - arrays for T2+- (O)
c       xSGrp   - SubGroups of be,ga (I)
c        LunAux  - Lun (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
c
        real*8 T2n1(1)
        real*8 T2n2(1)
        integer beSGrp,gaSGrp,LunAux
c
c       help variables
        integer length1,length2
        character*6 LunName
c
c1      calc legths
        if (beSGrp.eq.gaSGrp) then
          length1=no*(no+1)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/4
        else
          length1=no*(no+1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)/2
        end if
c
        if (beSGrp.eq.gaSGrp) then
          length2=no*(no-1)*DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)-1)/4
        else
          length2=no*(no-1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)/2
        end if
c
c2      get T2n1, T2n2
        LunName=Tmp3Name(beSGrp,gaSGrp)
        call GetX (T2n1(1),length1,LunAux,LunName,1,0)
        call GetX (T2n2(1),length2,LunAux,LunName,0,1)
c
        return
        end
c
c       ------------------------------------
c
        subroutine ExtractM (M,L2,cGrp,deGrp,cSGrp,deSGrp)
c
c       this routine do:
c       extract M(m,c",de") from L2(m,c',de')
c
c       parameter description:
c       M     - M file (O)
c       L     - L2 file (O)
c       xGrp  - c,delta Group (I)
c       xSGrp - c,delta Group (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 M(*),L2(*)
        integer cGrp,deGrp,cSGrp,deSGrp
c
c       help variables
        integer lenMCpp
        integer possM0,possL20,incL20
        integer i,k,depp
c
c
c1      Initial settings
c
c1.1    def length of mc" block (=increment for possM0)
        lenMCpp=nc*DimSGrpa(cSGrp)
c
c1.2    def possM0
        possM0=1
c
c1.3    def possL20
        PossL20=1
        k=0
        if (deSGrp.gt.Grpbelow(deGrp)) then
          do i=Grpbelow(deGrp),deSGrp-1
            k=k+DimSGrpbe(i)
          end do
        end if
        possL20=PossL20+nc*DimGrpa(cGrp)*k
c
        k=0
        if (cSGrp.gt.Grpalow(cGrp)) then
          do i=Grpalow(cGrp),cSGrp-1
            k=k+DimSGrpa(i)
          end do
        end if
        possL20=PossL20+nc*k
c
c1.4    def increment for possL20
        incL20=dimGrpa(cGrp)*nc
c
c2      cycle over de"
        do depp=1,DimSgrpbe(deSGrp)
c
c2.1      copy block of #mc" size
          call mv0u (lenMCpp,L2(possL20),1,M(possM0),1)
c
c2.2      upgrade possM0 and possL20 for next use
          possM0=possM0+lenMCpp
          possL20=possL20+incL20
c
        end do
c
c
        return
        end
c
c       ------------------------------------
c
        subroutine MakeT2p (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,keyT)
c
c       this routine do:
c       Make T2p(i>=j,(a>b)" )  = Tau(i,j,(a>=b)")+Tau(j,i,a>=b)")
c                           from  Tau((a>=b)',i,j)
c        or Transposed (T(ab",ij)
c
c       parameter description:
c       T2p    - T2+ array (O)
c       Tau    - Tau array (I)
c       xGrp   - Group of a,b (I)
c       xSGrp  - SubGroup of a,b (I)
c        keyT   - 0 - make T(ij,ab")
c                1 - make T(ab",ij)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 T2p(1)
        real*8 Tau(1)
        integer aGrp,bGrp,aSGrp,bSGrp,keyT
c
c       help variables
        integer dimi,dimij,dimap,dimbp,dimapp,dimbpp,dimabp,dimabpp
c
        dimi=no
        dimij=no*(no+1)/2
c
        dimap=DimGrpa(aGrp)
        dimbp=DimGrpa(bGrp)
        if (aGrp.eq.bGrp) then
          dimabp=dimap*(dimap+1)/2
        else
          dimabp=dimap*dimbp
        end if
c
        dimapp=DimSGrpa(aSGrp)
        dimbpp=DimSGrpa(bSGrp)
        if (aSGrp.eq.bSGrp) then
          dimabpp=dimapp*(dimapp-1)/2
        else
          dimabpp=dimapp*dimbpp
        end if
c
        if (keyT.eq.0) then
c        T+(ij,ab") case
c
        if (aGrp.eq.bGrp) then
          if (aSGrp.eq.bSGrp) then
            call makeT2pHlp1 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                        dimi,dimij,dimapp,dimabpp,dimap,dimabp)
          else
            call makeT2pHlp2 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                        dimi,dimij,dimapp,dimbpp,dimapp,dimabp)
          end if
        else
          call makeT2pHlp3 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                      dimi,dimij,dimapp,dimbpp,dimap,dimbp)
        end if
c
        else
c        T+(ab",ij) case
c
        if (aGrp.eq.bGrp) then
          if (aSGrp.eq.bSGrp) then
            call makeT2ptHlp1 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                        dimi,dimij,dimapp,dimabpp,dimap,dimabp)
          else
            call makeT2ptHlp2 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                        dimi,dimij,dimapp,dimbpp,dimapp,dimabp)
          end if
        else
          call makeT2ptHlp3 (T2p(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,0,
     c                      dimi,dimij,dimapp,dimbpp,dimap,dimbp)
        end if
c
        end if

c
        return
        end
c
c       ------------------------------------
c
        subroutine MakeT2m (T2m,Tau,aGrp,bGrp,aSGrp,bSGrp,keyT)
c
c       this routine do:
c       Make T2p(i>,(a>b)")  = Tau(i,j,(a>=b)")-Tau(j,i,a>=b)")
c                     from     Tau((a>=b)',i,j)
c        or Transposed (T(ab",ij)
c
c       parameter description:
c       T2m    - T2- array (O)
c       Tau    - Tau array (I)
c       xGrp   - Group of a,b (I)
c       xSGrp  - SubGroup of a,b (I)
c        keyT   - 0 - make T(ij,ab")
c                1 - make T(ab",ij)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 T2m(1)
        real*8 Tau(1)
        integer aGrp,bGrp,aSGrp,bSGrp,keyT
c
c       help variables
        integer dimi,dimij,dimap,dimbp,dimapp,dimbpp,dimabp,dimabpp
c
        dimi=no
        dimij=no*(no-1)/2
c
        dimap=DimGrpa(aGrp)
        dimbp=DimGrpa(bGrp)
        if (aGrp.eq.bGrp) then
          dimabp=dimap*(dimap+1)/2
        else
          dimabp=dimap*dimbp
        end if
c
        dimapp=DimSGrpa(aSGrp)
        dimbpp=DimSGrpa(bSGrp)
        if (aSGrp.eq.bSGrp) then
          dimabpp=dimapp*(dimapp-1)/2
        else
          dimabpp=dimapp*dimbpp
        end if
c
        if (keyT.eq.0) then
c        T-(ij,ab") case
c
        if (aGrp.eq.bGrp) then
          if (aSGrp.eq.bSGrp) then
            call makeT2pHlp1 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,
     c                        dimi,dimij,dimapp,dimabpp,dimap,dimabp)
          else
            call makeT2pHlp2 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,
     c                        dimi,dimij,dimapp,dimbpp,dimapp,dimabp)
          end if
        else
          call makeT2pHlp3 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,
     c                      dimi,dimij,dimapp,dimbpp,dimap,dimbp)
        end if
c
        else
c        T-(ab",ij) case
c
        if (aGrp.eq.bGrp) then
          if (aSGrp.eq.bSGrp) then
            call makeT2ptHlp1 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,
     c                        dimi,dimij,dimapp,dimabpp,dimap,dimabp)
          else
            call makeT2ptHlp2 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,
     c                        dimi,dimij,dimapp,dimbpp,dimapp,dimabp)
          end if
        else
          call makeT2ptHlp3 (T2m(1),Tau(1),aGrp,bGrp,aSGrp,bSGrp,1,
     c                      dimi,dimij,dimapp,dimbpp,dimap,dimbp)
        end if
c
        end if
c
        return
        end
c
c       ---------------------
c
        subroutine makeT2pHlp1 (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,key,
     c                          dimi,dimij,dimapp,dimabpp,dimap,dimabp)
c
c       this routine do:
c       define T2(+-)(ij,(ab)") = Tau((ab)',i,j) +- Tau((ab)',j,i)
c       for the case: aGrp=bGrp and aSGrp=bSGrp
c
c       parameter description:
c       T2p     - array for T2+ (O)
c       Tau     - array for Tau (I)
c       xGrp    - Groups of a',b' (I)
c       xSGrp   - SubGroups of a",b" (I)
c       key     - 0 - T2+; 1 = T2- will be produced
c       dimx    - Dimension of i,(i>=j),a",(a>b)",a',(a>=b)' (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,bGrp,aSGrp,bSGrp,key
        integer dimi,dimij,dimapp,dimabpp,dimap,dimabp
        real*8 T2p(1:dimij,1:dimabpp)
        real*8 Tau(1:dimabp,1:dimi,1:dimi)
c
c
c       help variables
        integer i,j,ij,app,bpp,abpp,ap,abp,appAdd,bppAdd
c
c
c1      def appAdd,bppAdd
c
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
c
        bppAdd=0
        if (bSGrp.ne.Grpalow(bGrp)) then
          do i=Grpalow(bGrp),bSGrp-1
          bppAdd=bppAdd+DimSGrpa(i)
          end do
        end if
c
c
c2      define T2(+,-)
c
        if (key.eq.0) then
c
c2.1    define T2+
        abpp=0
        do app=2,dimapp
        ap=appAdd+app
        abp=ap*(ap-1)/2+bppAdd
        do bpp=1,app-1
        abp=abp+1
        abpp=abpp+1
          ij=0
          do i=1,dimi
          do j=1,i
            ij=ij+1
            T2p(ij,abpp)=Tau(abp,i,j)+Tau(abp,j,i)
          end do
          end do
        end do
        end do
c
        else
c
c2.2    define T2-
        abpp=0
        do app=2,dimapp
        ap=appAdd+app
        abp=ap*(ap-1)/2+bppAdd
        do bpp=1,app-1
        abp=abp+1
        abpp=abpp+1
          ij=0
          do i=2,dimi
          do j=1,i-1
            ij=ij+1
            T2p(ij,abpp)=Tau(abp,i,j)-Tau(abp,j,i)
c           T2p(ij,abpp)=T2c(ap,bpp,i,j)-T2c(ap,bpp,j,i)
          end do
          end do
        end do
        end do
c
        end if
c
        call mv0sv (dimij*dimabpp,dimij*dimabpp,
     c              T2p(1,1),0.5d0)
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_integer(dimap)
        end
c
c       ---------------------
c
        subroutine makeT2pHlp2 (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,key,
     c                          dimi,dimij,dimapp,dimbpp,dimap,dimabp)
c
c       this routine do:
c       define T2(+-)(ij,(ab)") = Tau((ab)',i,j) +- Tau((ab)',j,i)
c       for the case: aGrp=bGrp but aSGrp>bSGrp
c
c       parameter description:
c       T2p     - array for T2+ (O)
c       Tau     - array for Tau (I)
c       xGrp    - Groups of a',b' (I)
c       xSGrp   - SubGroups of a",b" (I)
c       key     - 0 - T2+; 1 = T2- will be produced
c       dimx    - Dimension of i,(i>=j),a",b",a',(a>=b)' (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,bGrp,aSGrp,bSGrp,key
        integer dimi,dimij,dimapp,dimbpp,dimap,dimabp
        real*8 T2p(1:dimij,1:dimapp,1:dimbpp)
        real*8 Tau(1:dimabp,1:dimi,1:dimi)
c
c
c       help variables
        integer i,j,ij,app,bpp,ap,abp,appAdd,bppAdd
c
c
c1      def appAdd,bppAdd
c
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
c
        bppAdd=0
        if (bSGrp.ne.Grpalow(bGrp)) then
          do i=Grpalow(bGrp),bSGrp-1
          bppAdd=bppAdd+DimSGrpa(i)
          end do
        end if
c
c
c2      define T2(+,-)
c
        if (key.eq.0) then
c
c2.1    define T2+
        do app=1,dimapp
        ap=appAdd+app
        abp=ap*(ap-1)/2+bppAdd
        do bpp=1,dimbpp
        abp=abp+1
          ij=0
          do i=1,dimi
          do j=1,i
            ij=ij+1
            T2p(ij,app,bpp)=Tau(abp,i,j)+Tau(abp,j,i)
          end do
          end do
        end do
        end do
c
        else
c
c2.2    define T2-
        do app=1,dimapp
        ap=appAdd+app
        abp=ap*(ap-1)/2+bppAdd
        do bpp=1,dimbpp
        abp=abp+1
          ij=0
          do i=2,dimi
          do j=1,i-1
            ij=ij+1
            T2p(ij,app,bpp)=Tau(abp,i,j)-Tau(abp,j,i)
          end do
          end do
        end do
        end do
c
        end if
c
        call mv0sv (dimij*dimapp*dimbpp,dimij*dimapp*dimbpp,
     c              T2p(1,1,1),0.5d0)
c
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_integer(dimap)
        end
c
c       ---------------------
c
        subroutine makeT2pHlp3 (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,key,
     c                          dimi,dimij,dimapp,dimbpp,dimap,dimbp)
c
c       this routine do:
c       define T2(+-)(ij,(ab)") = Tau((ab)',i,j) +- Tau((ab)',j,i)
c       for the case: aGrp>bGrp (=> aSGrp>bSGrp)
c
c       parameter description:
c       T2p     - array for T2+ (O)
c       Tau     - array for Tau (I)
c       xGrp    - Groups of a',b' (I)
c       xSGrp   - SubGroups of a",b" (I)
c       key     - 0 - T2+; 1 = T2- will be produced
c       dimx    - Dimension of i,(i>=j),a",b",a',b' (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,bGrp,aSGrp,bSGrp,key
        integer dimi,dimij,dimapp,dimbpp,dimap,dimbp
        real*8 T2p(1:dimij,1:dimapp,1:dimbpp)
        real*8 Tau(1:dimap,1:dimbp,1:dimi,1:dimi)
c
c       help variables
        integer i,j,ij,app,bpp,ap,bp,appAdd,bppAdd
c
c
c1      def appAdd,bppAdd
c
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
c
        bppAdd=0
        if (bSGrp.ne.Grpalow(bGrp)) then
          do i=Grpalow(bGrp),bSGrp-1
          bppAdd=bppAdd+DimSGrpa(i)
          end do
        end if
c
c
c2      define T2(+,-)
c
        if (key.eq.0) then
c
c2.1    define T2+
        do bpp=1,dimbpp
        bp=bppAdd+bpp
        do app=1,dimapp
        ap=appAdd+app
          ij=0
          do i=1,dimi
          do j=1,i
            ij=ij+1
            T2p(ij,app,bpp)=Tau(ap,bp,i,j)+Tau(ap,bp,j,i)
          end do
          end do
        end do
        end do
c
        else
c
c2.2    define T2-
        do bpp=1,dimbpp
        bp=bppAdd+bpp
        do app=1,dimapp
        ap=appAdd+app
          ij=0
          do i=2,dimi
          do j=1,i-1
            ij=ij+1
            T2p(ij,app,bpp)=Tau(ap,bp,i,j)-Tau(ap,bp,j,i)
          end do
          end do
        end do
        end do
c
        end if
c
        call mv0sv (dimij*dimapp*dimbpp,dimij*dimapp*dimbpp,
     c              T2p(1,1,1),0.5d0)
c
c
        return
        end
c
c       ---------------------
c
        subroutine makeT2ptHlp1 (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,key,
     c                          dimi,dimij,dimapp,dimabpp,dimap,dimabp)
c
c       this routine do:
c       define T2(+-)((ab)",ij) = Tau((ab)',i,j) +- Tau((ab)',j,i)
c       for the case: aGrp=bGrp and aSGrp=bSGrp
c
c       parameter description:
c       T2p     - array for T2+ (O)
c       Tau     - array for Tau (I)
c       xGrp    - Groups of a',b' (I)
c       xSGrp   - SubGroups of a",b" (I)
c       key     - 0 - T2+; 1 = T2- will be produced
c       dimx    - Dimension of i,(i>=j),a",(a>(or>=)b)",a',(a>=b)' (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,bGrp,aSGrp,bSGrp,key
        integer dimi,dimij,dimapp,dimabpp,dimap,dimabp
        real*8 T2p(1:dimabpp,1:dimij)
        real*8 Tau(1:dimabp,1:dimi,1:dimi)
c
c
c       help variables
        integer i,j,ij,app,bpp,abpp,ap,abp,appAdd,bppAdd
c
c
c1      def appAdd,bppAdd
c
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
c
        bppAdd=0
        if (bSGrp.ne.Grpalow(bGrp)) then
          do i=Grpalow(bGrp),bSGrp-1
          bppAdd=bppAdd+DimSGrpa(i)
          end do
        end if
c
c
c2      define T2(+,-)
c
        if (key.eq.0) then
c
c2.1    define T2+
        ij=0
        do i=1,dimi
        do j=1,i
        ij=ij+1
          abpp=0
          do app=2,dimapp
          ap=appAdd+app
          abp=ap*(ap-1)/2+bppAdd
          do bpp=1,app-1
          abp=abp+1
          abpp=abpp+1
            T2p(abpp,ij)=Tau(abp,i,j)+Tau(abp,j,i)
          end do
          end do
        end do
        end do
c
        else
c
c2.2    define T2-
        ij=0
        do i=2,dimi
        do j=1,i-1
        ij=ij+1
          abpp=0
          do app=2,dimapp
          ap=appAdd+app
          abp=ap*(ap-1)/2+bppAdd
          do bpp=1,app-1
          abp=abp+1
          abpp=abpp+1
            T2p(abpp,ij)=Tau(abp,i,j)-Tau(abp,j,i)
          end do
          end do
        end do
        end do
c
        end if
c
        call mv0sv (dimij*dimabpp,dimij*dimabpp,
     c              T2p(1,1),0.5d0)
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_integer(dimap)
        end
c
c       ---------------------
c
        subroutine makeT2ptHlp2 (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,key,
     c                          dimi,dimij,dimapp,dimbpp,dimap,dimabp)
c
c       this routine do:
c       define T2(+-)((ab)",ij) = Tau((ab)',i,j) +- Tau((ab)',j,i)
c       for the case: aGrp=bGrp but aSGrp>bSGrp
c
c       parameter description:
c       T2p     - array for T2+ (O)
c       Tau     - array for Tau (I)
c       xGrp    - Groups of a',b' (I)
c       xSGrp   - SubGroups of a",b" (I)
c       key     - 0 - T2+; 1 = T2- will be produced
c       dimx    - Dimension of i,(i>=j),a",b",a',(a>=b)' (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,bGrp,aSGrp,bSGrp,key
        integer dimi,dimij,dimapp,dimbpp,dimap,dimabp
        real*8 T2p(1:dimapp,1:dimbpp,1:dimij)
        real*8 Tau(1:dimabp,1:dimi,1:dimi)
c
c
c       help variables
        integer i,j,ij,app,bpp,ap,abp,appAdd,bppAdd
c
c
c1      def appAdd,bppAdd
c
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
c
        bppAdd=0
        if (bSGrp.ne.Grpalow(bGrp)) then
          do i=Grpalow(bGrp),bSGrp-1
          bppAdd=bppAdd+DimSGrpa(i)
          end do
        end if
c
c
c2      define T2(+,-)
c
        if (key.eq.0) then
c
c2.1    define T2+
        ij=0
        do i=1,dimi
        do j=1,i
        ij=ij+1
          do app=1,dimapp
          ap=appAdd+app
          abp=ap*(ap-1)/2+bppAdd
          do bpp=1,dimbpp
          abp=abp+1
            T2p(app,bpp,ij)=Tau(abp,i,j)+Tau(abp,j,i)
          end do
          end do
        end do
        end do
c
        else
c
c2.2    define T2-
        ij=0
        do i=2,dimi
        do j=1,i-1
        ij=ij+1
          do app=1,dimapp
          ap=appAdd+app
          abp=ap*(ap-1)/2+bppAdd
          do bpp=1,dimbpp
          abp=abp+1
            T2p(app,bpp,ij)=Tau(abp,i,j)-Tau(abp,j,i)
          end do
          end do
        end do
        end do
c
        end if
c
        call mv0sv (dimij*dimapp*dimbpp,dimij*dimapp*dimbpp,
     c              T2p(1,1,1),0.5d0)
c
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_integer(dimap)
        end
c
c       ---------------------
c
        subroutine makeT2ptHlp3 (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,key,
     c                          dimi,dimij,dimapp,dimbpp,dimap,dimbp)
c
c       this routine do:
c       define T2(+-)((ab)",ij) = Tau((ab)',i,j) +- Tau((ab)',j,i)
c       for the case: aGrp>bGrp (=> aSGrp>bSGrp)
c
c       parameter description:
c       T2p     - array for T2+ (O)
c       Tau     - array for Tau (I)
c       xGrp    - Groups of a',b' (I)
c       xSGrp   - SubGroups of a",b" (I)
c       key     - 0 - T2+; 1 = T2- will be produced
c       dimx    - Dimension of i,(i>=j),a",b",a',b' (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,bGrp,aSGrp,bSGrp,key
        integer dimi,dimij,dimapp,dimbpp,dimap,dimbp
        real*8 T2p(1:dimapp,1:dimbpp,1:dimij)
        real*8 Tau(1:dimap,1:dimbp,1:dimi,1:dimi)
c
c       help variables
        integer i,j,ij,app,bpp,ap,bp,appAdd,bppAdd
c
c
c1      def appAdd,bppAdd
c
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
c
        bppAdd=0
        if (bSGrp.ne.Grpalow(bGrp)) then
          do i=Grpalow(bGrp),bSGrp-1
          bppAdd=bppAdd+DimSGrpa(i)
          end do
        end if
c
c
c2      define T2(+,-)
c
        if (key.eq.0) then
c
c2.1    define T2+
        ij=0
        do i=1,dimi
        do j=1,i
        ij=ij+1
          do bpp=1,dimbpp
          bp=bppAdd+bpp
          do app=1,dimapp
          ap=appAdd+app
            T2p(app,bpp,ij)=Tau(ap,bp,i,j)+Tau(ap,bp,j,i)
          end do
          end do
        end do
        end do
c
        else
c
c2.2    define T2-
        ij=0
        do i=2,dimi
        do j=1,i-1
        ij=ij+1
          do bpp=1,dimbpp
          bp=bppAdd+bpp
          do app=1,dimapp
          ap=appAdd+app
            T2p(app,bpp,ij)=Tau(ap,bp,i,j)-Tau(ap,bp,j,i)
          end do
          end do
        end do
        end do
c
        end if
c
        call mv0sv (dimij*dimapp*dimbpp,dimij*dimapp*dimbpp,
     c              T2p(1,1,1),0.5d0)
c
c
        return
        end
c
c       ------------------------------------
c
        subroutine MakeWw (Ww,W1,W2,aSGrp,bSGrp,beSGrp,gaSGrp,key)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c                                     and W2(b",be",a",ga")
c        where index ab" is always a">b" (not a"=b")
c        and   index bega" is be">=ga" for W+ and be">ga" for W-
c       for all cases aSGrp >= bSGrp
c       N.B. for a"=b" W2 is not used (although it is passed through
c            the header (because of uniform shape) but no matter what
c            is there)
c       N.B. rutinky MakeWwHlpx niesu prilis vymakane @
c
c       parameter description:
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       W2     - array for W1(b",be",a",ga") (I)
c       xSGrp  - SubGroup of a",b",be",ga" (I)
c       key    - 1 - calc Ww+, 2 - calc Ww- (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 Ww(1)
        real*8 W1(1)
        real*8 W2(1)
        integer aSGrp,bSGrp,beSGrp,gaSGrp,key
c
c       help variables
        integer dimbe,dimga,dimbega,dima,dimb,dimab
c
c
c1      def dimensions
c
        dimbe=DimSGrpbe(beSGrp)
        dimga=DimSGrpbe(gaSGrp)
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpa(bSGrp)
c
        if (beSGrp.eq.gaSGrp) then
          if (key.eq.1) then
            dimbega=dimbe*(dimbe+1)/2
          else
            dimbega=dimbe*(dimbe-1)/2
          end if
        else
          dimbega=dimbe*dimga
        end if
c
        if (aSGrp.eq.bSGrp) then
          dimab=dima*(dima-1)/2
        else
          dimab=dima*dimb
        end if
c
c
c2      Make Ww matrix
c
        if (aSGrp.eq.bSGrp) then
          if (beSGrp.eq.gaSGrp) then
            call MakeWwHlp1 (Ww(1),W1(1),
     c                       dima,dimb,dimab,dimbe,dimga,dimbega,key)
          else
            call MakeWwHlp2 (Ww(1),W1(1),
     c                       dima,dimb,dimab,dimbe,dimga,key)
          end if
        else
          if (beSGrp.eq.gaSGrp) then
            call MakeWwHlp3 (Ww(1),W1(1),W2(1),
     c                       dima,dimb,dimbe,dimga,dimbega,key)
          else
            call MakeWwHlp4 (Ww(1),W1(1),W2(1),
     c                       dima,dimb,dimbe,dimga,key)
          end if
        end if
c
c
        return
        end
c
c       ---------------------
c
        subroutine MakeWwHlp1 (Ww,W1,
     c                       dima,dimb,dimab,dimbe,dimga,dimbega,key)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c       for the case a"=b" , be"=ga"
c       N.B. algoritmus nieje prilis vymakany
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       dimx   - dimension of a",b",ab",be",ga",bega" (I)
c       key    - 1 - calc Ww+, 2 - calc Ww- (I)
c
        implicit none
        integer dima,dimb,dimab,dimbe,dimga,dimbega,key
        real*8 Ww(1:dimab,1:dimbega)
        real*8 W1(1:dima,1:dimbe,1:dimb,1:dimga)
c
c       help variables
        integer a,b,ab,be,ga,bega
c
c
        if (key.eq.1) then
          bega=0
          do be=1,dimbe
          do ga=1,be
          bega=bega+1
          ab=0
          do a=2,dima
          do b=1,a-1
          ab=ab+1
            Ww(ab,bega)=W1(a,be,b,ga)+W1(b,be,a,ga)
          end do
          end do
          end do
          end do
        else
          bega=0
          do be=2,dimbe
          do ga=1,be-1
          bega=bega+1
          ab=0
          do a=2,dima
          do b=1,a-1
          ab=ab+1
            Ww(ab,bega)=W1(a,be,b,ga)-W1(b,be,a,ga)
          end do
          end do
          end do
          end do
        end if
c
c        Cely clen ma Faktor 2, tu teda nevydelim 2
c        call mv0sv (dimab*dimbega,dimab*dimbega,
c    c              Ww(1,1),0.5d0)
c
        return
        end
c
c       ---------------------
c
        subroutine MakeWwHlp2 (Ww,W1,
     c                       dima,dimb,dimab,dimbe,dimga,key)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c       for the case a"=b" , be".ne.ga"
c       N.B. algoritmus nieje prilis vymakany
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       dimx   - dimension of a",b",ab",be",ga" (I)
c       key    - 1 - calc Ww+, 2 - calc Ww- (I)
c
        implicit none
        integer dima,dimb,dimab,dimbe,dimga,key
        real*8 Ww(1:dimab,1:dimbe,1:dimga)
        real*8 W1(1:dima,1:dimbe,1:dimb,1:dimga)
c
c       help variables
        integer a,b,ab,be,ga
c
c
        if (key.eq.1) then
          do ga=1,dimga
          do be=1,dimbe
          ab=0
          do a=2,dima
          do b=1,a-1
          ab=ab+1
            Ww(ab,be,ga)=W1(a,be,b,ga)+W1(b,be,a,ga)
          end do
          end do
          end do
          end do
        else
          do ga=1,dimga
          do be=1,dimbe
          ab=0
          do a=2,dima
          do b=1,a-1
          ab=ab+1
            Ww(ab,be,ga)=W1(a,be,b,ga)-W1(b,be,a,ga)
          end do
          end do
          end do
          end do
        end if
c
c        Cely clen ma Faktor 2, tu teda nevydelim 2
c        call mv0sv (dimab*dimbe*dimga,dimab*dimbe*dimga,
c    c              Ww(1,1,1),0.5d0)
c
        return
        end
c
c       ---------------------
c
        subroutine MakeWwHlp3 (Ww,W1,W2,
     c                       dima,dimb,dimbe,dimga,dimbega,key)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c                                     and W2(b",be",a",ga")
c       for the case a".ne.b" , be"=ga"
c       N.B. algoritmus nieje prilis vymakany, prva cast sa da
c            urobit pomodou dcopy
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       W2     - array for W1(b",be",a",ga") (I)
c       dimx   - dimension of a",b",be",ga",be>(>=)ga" (I)
c       key    - 1 - calc Ww+, 2 - calc Ww- (I)
c
        implicit none
        integer dima,dimb,dimbe,dimga,dimbega,key
        real*8 Ww(1:dima,1:dimb,1:dimbega)
        real*8 W1(1:dima,1:dimbe,1:dimb,1:dimga)
        real*8 W2(1:dimb,1:dimbe,1:dima,1:dimga)
c
c       help variables
        integer a,b,be,ga,bega
c
c
        if (key.eq.1) then
          bega=0
          do be=1,dimbe
          do ga=1,be
          bega=bega+1
          do b=1,dimb
          do a=1,dima
            Ww(a,b,bega)=W1(a,be,b,ga)+W2(b,be,a,ga)
          end do
          end do
          end do
          end do
        else
          bega=0
          do be=2,dimbe
          do ga=1,be-1
          bega=bega+1
          do b=1,dimb
          do a=1,dima
            Ww(a,b,bega)=W1(a,be,b,ga)-W2(b,be,a,ga)
          end do
          end do
          end do
          end do
        end if
c
c        Cely clen ma Faktor 2, tu teda nevydelim 2
c        call mv0sv (dima*dimb*dimbega,dima*dimb*dimbega,
c    c              Ww(1,1,1),0.5d0)
c
        return
        end
c
c       ---------------------
c
        subroutine MakeWwHlp4 (Ww,W1,W2,
     c                       dima,dimb,dimbe,dimga,key)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c                                     and W2(b",be",a",ga")
c       for the case a".ne.b" , be".ne.ga"
c       N.B. algoritmus nieje prilis vymakany, prva cast sa da
c            urobit pomodou dcopy
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       W2     - array for W1(b",be",a",ga") (I)
c       dimx   - dimension of a",b",be",ga" (I)
c       key    - 1 - calc Ww+, 2 - calc Ww- (I)
c
        implicit none
        integer dima,dimb,dimbe,dimga,key
        real*8 Ww(1:dima,1:dimb,1:dimbe,1:dimga)
        real*8 W1(1:dima,1:dimbe,1:dimb,1:dimga)
        real*8 W2(1:dimb,1:dimbe,1:dima,1:dimga)
c
c       help variables
        integer a,b,be,ga
c
c
        if (key.eq.1) then
          do ga=1,dimga
          do be=1,dimbe
          do b=1,dimb
          do a=1,dima
            Ww(a,b,be,ga)=W1(a,be,b,ga)+W2(b,be,a,ga)
          end do
          end do
          end do
          end do
        else
          do ga=1,dimga
          do be=1,dimbe
          do b=1,dimb
          do a=1,dima
            Ww(a,b,be,ga)=W1(a,be,b,ga)-W2(b,be,a,ga)
          end do
          end do
          end do
          end do
        end if
c
c        Cely clen ma Faktor 2, tu teda nevydelim 2
c        call mv0sv (dima*dimb*dimbe*dimga,dima*dimb*dimbe*dimga,
c    c              Ww(1,1,1,1),0.5d0)
c
        return
        end
c
c       ------------------------------------
c
        subroutine MakeT2pd (T2p,Tau,aGrp,aSGrp)
c
c       this routine do:
c       Make T2p((a>b)",i>=j )  = Tau(i,j,(a>=b)")+Tau(j,i,a>=b)")
c                           from  Tau((a>=b)',i,j)
c
c       parameter description:
c       T2p    - T2+ array (O)
c       Tau    - Tau array (I)
c       xGrp   - Group of a (I)
c       xSGrp  - SubGroup of a (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 T2p(1)
        real*8 Tau(1)
        integer aGrp,aSGrp
c
c       help variables
        integer dimi,dimij,dimap,dimapp,dimabp
c
        dimi=no
        dimij=no*(no+1)/2
c
        dimap=DimGrpa(aGrp)
        dimabp=dimap*(dimap+1)/2
c
        dimapp=DimSGrpa(aSGrp)
c
c
        call makeT2pdHlp (T2p(1),Tau(1),aGrp,aSGrp,
     c                   dimi,dimij,dimapp,dimap,dimabp)
c
        return
        end
c
c       ---------------------
c
        subroutine makeT2pdHlp (T2p,Tau,aGrp,aSGrp,
     c                          dimi,dimij,dimapp,dimap,dimabp)
c
c       this routine do:
c       define T2(+)((aa)",ij) = Tau((ab)',i,j) + Tau((ab)',j,i)
c       for the case: aGrp=bGrp and aSGrp=bSGrp
c
c       parameter description:
c       T2p     - array for T2+ (O)
c       Tau     - array for Tau (I)
c       xGrp    - Groups of a',b' (I)
c       xSGrp   - SubGroups of a",b" (I)
c       dimx    - Dimension of i,(i>=j),a",a',(a>=b)' (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,aSGrp
        integer dimi,dimij,dimapp,dimap,dimabp
        real*8 T2p(1:dimapp,1:dimij)
        real*8 Tau(1:dimabp,1:dimi,1:dimi)
c
c
c       help variables
        integer i,j,ij,app,ap,abp,appAdd
c
c
c1      def appAdd
c
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
c
c2        define T2+(aa",ij)
c
        ij=0
        do i=1,dimi
        do j=1,i
        ij=ij+1
          do app=1,dimapp
          ap=appAdd+app
          abp=ap*(ap+1)/2
            T2p(app,ij)=Tau(abp,i,j)+Tau(abp,j,i)
          end do
        end do
        end do
c
        call mv0sv (dimij*dimapp,dimij*dimapp,
     c              T2p(1,1),0.5d0)
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_integer(dimap)
        end
c
c       ------------------------------------
c
        subroutine MakeWwd (Ww,W1,aSGrp,beSGrp,gaSGrp)
c
c       this routine do:
c       Make  Ww(+)((aa)",(bega)") from W1(a",be",a",ga")
c        and   index bega" is be">=ga"
c       for all cases aSGrp >= bSGrp
c       N.B. rutinky MakeWwHlpx niesu prilis vymakane @
c
c       parameter description:
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",a",ga") (I)
c       xSGrp  - SubGroup of a",be",ga" (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 Ww(1)
        real*8 W1(1)
        integer aSGrp,beSGrp,gaSGrp
c
c       help variables
        integer dimbe,dimga,dimbega,dima
c
c
c1      def dimensions
c
        dimbe=DimSGrpbe(beSGrp)
        dimga=DimSGrpbe(gaSGrp)
        dima=DimSGrpa(aSGrp)
c
        if (beSGrp.eq.gaSGrp) then
          dimbega=dimbe*(dimbe+1)/2
        else
          dimbega=dimbe*dimga
        end if
c
c
c2      Make Ww matrix
c
        if (beSGrp.eq.gaSGrp) then
          call MakeWwdHlp1 (Ww(1),W1(1),
     c                   dima,dimbe,dimbega)
        else
          call MakeWwdHlp2 (Ww(1),W1(1),
     c                   dima,dimbe,dimga)
        end if

c
c
        return
        end
c
c       ---------------------
c
        subroutine MakeWwdHlp1 (Ww,W1,
     c                          dima,dimbe,dimbega)
c
c       this routine do:
c       Make  Ww(+)((aa)",(bega)") from W1(a",be",b",ga")
c       for the case beSGrp=gaSGrp
c       N.B. algoritmus nieje prilis vymakany
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       dimx   - dimension of a",ga",bega" (I)
c
        implicit none
        integer dima,dimbe,dimbega
        real*8 Ww(1:dima,1:dimbega)
        real*8 W1(1:dima,1:dimbe,1:dima,1:dimbe)
c
c       help variables
        integer a,be,ga,bega
c
c
          bega=0
          do be=1,dimbe
          do ga=1,be
          bega=bega+1
          do a=1,dima
            Ww(a,bega)=W1(a,be,a,ga)
          end do
          end do
          end do
c
        return
        end
c
c       ---------------------
c
        subroutine MakeWwdHlp2 (Ww,W1,
     c                       dima,dimbe,dimga)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c       for the case a"=b" , beSGrp.ne.gaSGrp
c       N.B. algoritmus nieje prilis vymakany
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       dimx   - dimension of a",be",ga" (I)
c
        implicit none
        integer dima,dimbe,dimga
        real*8 Ww(1:dima,1:dimbe,1:dimga)
        real*8 W1(1:dima,1:dimbe,1:dima,1:dimga)
c
c       help variables
        integer a,be,ga
c
c
CVpV 2014 Fix for Intel Compiler v14.*
#ifdef __INTEL_COMPILER
          do a=1,dima
          do ga=1,dimga
cDEC$ VECTOR UNALIGNED
          do be=1,dimbe
            Ww(a,be,ga)=W1(a,be,a,ga)
          end do
          end do
          end do
#else
          do ga=1,dimga
          do be=1,dimbe
          do a=1,dima
            Ww(a,be,ga)=W1(a,be,a,ga)
          end do
          end do
          end do
#endif
c
        return
        end
c
c        -------------------------
c
        subroutine Calc_addSG (aSGrp,adda)
c
c        calc add constant from SubGroup
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
        integer aSGrp,adda
c
c        help var
        integer i
c
        adda=0
        do i=1,aSgrp-1
          adda=adda+DimSGrpa(i)
        end do
c
        return
        end
c
c        -------------------------
c
        subroutine Calc_addG (aGrp,adda)
c
c        calc add constant from group
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
        integer aGrp,adda
c
c        help var
        integer i
c
        adda=0
        do i=1,agrp-1
          adda=adda+DimGrpa(i)
        end do
c
        return
        end
c
c        -------------------------
c        -------------------------
c
        subroutine  Mk_T1t (T1,H,dimbepp,no,nv,addbepp)
c
c        this routine do:
c        H(i,be") <- T1o(be,i)
c
        implicit none
        integer dimbepp,no,nv,addbepp
        real*8 T1(1:nv,1:no)
        real*8 H(1:no,1:dimbepp)
c
c        help variables
        integer i,bepp
c
        do i=1,no
        do bepp=1,dimbepp
          H(i,bepp)=T1(addbepp+bepp,i)
        end do
        end do
c
        return
        end
c
c        -------------------------
c
        subroutine ReaW3 (Ww,Wx,aSGrp,beSGrp,bSGrp,LunAux)
c
c        this routine do:
c        define Ww(a",be",b",i) <- (a",be"|b",i)
c
c        integrals (a",be"|b",i) are stored in files V3xxyyzz for xx>=yy
c
c        Storing of integrals in V3files
c        do i=1,no
c          record of V3(a"be",b",_i): a">=be",b" for given i
c        end do
c
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aSGrp,beSGrp,bSGrp,LunAux
        real*8 Ww(1)
        real*8 Wx(1)
c
c        help variables
c
        integer dima,dimb,dimbe
        character*8 LunName

c0        def dimensions
c
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpa(bSGrp)
        dimbe=DimSGrpbe(beSGrp)
c
c
        if (aSGrp.gt.beSGrp) then
c1        case aSGrp>beSGrp, integrals (a",be"|b",i) to be read
c1.1          make Name
          call MkNameV3 (aSGrp,beSGrp,bSGrp,'W3',LunName)
c1.2          read (a",be"|b",i) from V3:(a",be"|b",i)
          call ReaW3hlp1 (Ww,Wx,dima,dimbe,dimb,no,LunName,LunAux)
c
        else if (aSGrp.eq.beSGrp) then
c2        case aSGrp=beSGrp, integrals (a">=be"|b",i)to be read and expand
c2.1          make Name
          call MkNameV3 (aSGrp,beSGrp,bSGrp,'W3',LunName)
c2.2          read (a",be"|b",i) from V3:(a">=be"|b",i)
          call ReaW3hlp2 (Ww,Wx,dima,dimb,no,LunName,LunAux)
c
        else
c3        case aSGrp<beSGrp, integrals (be",a"|b",i) to be read and mapped
c3.1          make Name
          call MkNameV3 (beSGrp,aSGrp,bSGrp,'W3',LunName)
c3.2          read (a",be"|b",i) from V3:(be",a"|b",i)
          call ReaW3hlp3 (Ww,Wx,dima,dimbe,dimb,no,LunName,LunAux)

        end if
c
        return
        end
c
c        -------------------
c
        subroutine MkNameV3 (i,j,k,Schem,Nomen)
c
c       help routine to ReaW3, producing name of V3file
c       ex: Schem='XY', i=1, j=3, k=5 ->  Nomen='XY010305'
c
        implicit none
        integer i,j,k
        character*2 Schem
        character*8 Nomen
c
c       help variables
        character*1 Chr(1:8)
        character*2 digit(1:64)
        character*2 ichr,jchr,kchr
        character*2 baza
        character*8 meno
c
        equivalence (Chr(1),meno)
        equivalence (Chr(1),baza)
        equivalence (Chr(3),ichr)
        equivalence (Chr(5),jchr)
        equivalence (Chr(7),kchr)
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
        kchr=digit(k)
        Nomen=meno
c
        return
        end
c
c        -------------------
c
        subroutine ReaW3hlp1 (Ww,Wx,dima,dimbe,dimb,no,LunName,LunAux)
c
c        this routine do:
c        reconstruct  Ww(a",be',b,i)  for aSGrp>beSGrp
c        from (a",be'|b,i) records in V3 file LunName
c
        implicit none
        integer dima,dimbe,dimb,no,LunAux
        character*8 LunName
        real*8 Ww(1)
        real*8 Wx(1)
c
c        help variables
c
        integer length

*       open (unit=LunAux,file=LunName,form='unformatted')
        Call Molcas_BinaryOpen_Vanilla(LunAux,LunName)
c
        length=dima*dimbe*dimb*no
c
c        read block (a",be'|b,_i)
        call rea_chcc (LunAux,length,Ww(1))
cmp        call mv0zero (length,length,Ww(1))
c
        close (LunAux)
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_real_array(Wx)
        end
c
c        -------------------
c
        subroutine ReaW3hlp2 (Ww,Wx,dima,dimb,no,LunName,LunAux)
c
c        this routine do:
c        reconstruct  Ww(a",be',b,i)  for aSGrp=beSGrp
c        from (a">=be"|b,i) records in V3 file LunName
c
        implicit none
        integer dima,dimb,no,LunAux
        character*8 LunName
        real*8 Ww(1:dima,1:dima,1:dimb,1:no)
        real*8 Wx(*)

c        help variables
c
        integer i,a,be,abebi,b,length

c        read block (a">=be"|b"_i)
        length=(no*dima*(dima+1)*dimb)/2
*       open (unit=LunAux,file=LunName,form='unformatted')
        Call Molcas_BinaryOpen_Vanilla(LunAux,LunName)
        call rea_chcc (LunAux,length,Wx(1))
cmp        call mv0zero (length,length,Wx(1))
        close (LunAux)
c
c          Expand and Set Ww(a",be",b",i) <- Wx(a">=be"|b",i)
        abebi=0
c
        do i=1,no
          do b=1,dimb
          do a=1,dima
          do be=1,a
            abebi=abebi+1
            Ww(a,be,b,i)=Wx(abebi)
            Ww(be,a,b,i)=Wx(abebi)
          end do
          end do
          end do
        end do
c
c
        return
        end
c
c        -------------------
c
        subroutine ReaW3hlp3 (Ww,Wx,dima,dimbe,dimb,no,LunName,LunAux)
c
c        this routine do:
c        reconstruct  Ww(a",be',b,i)  for aSGrp<beSGrp
c        from (be",a"|b,i) records in V3 file LunName
c
        implicit none
        integer dima,dimbe,dimb,no,LunAux
        character*8 LunName
        real*8 Ww(1:dima,1:dimbe,1:dimb,1:no)
        real*8 Wx(*)

c        help variables
c
        integer i,a,be,b,beabi,length
c
c        read block (be",a"|b"_i)
c
        length=dima*dimbe*dimb*no
*       open (unit=LunAux,file=LunName,form='unformatted')
        Call Molcas_BinaryOpen_Vanilla(LunAux,LunName)
        call rea_chcc (LunAux,length,Wx(1))
cmp        call mv0zero (length,length,Wx(1))
        close (LunAux)
c
c          Map Ww(a",be",b",i) <- Wx(be",a"|b"_i)
c
        beabi=0
        do i=1,no
          do b=1,dimb
          do a=1,dima
          do be=1,dimbe
          beabi=beabi+1
            Ww(a,be,b,i)=Wx(beabi)
          end do
          end do
          end do
        end do
c
c
        return
        end
c
c        -------------------
c
        subroutine ReaW4 (W,Wa,aSGrp,bSGrp,cSGrp,dSGrp,LunAux)
c
c        this routine do:
c        Add W(a",b",c",d") <<- (a",b"|c",d") for any a",b",c",d"
c        via reading records W4 (p">=q")>=(r">=s")
c
c        Wa is an auxiliary file
c
c        Structure of records:
c        1) only records for pSGrp>=qSGrp & rSGrp>=sSGrp & pqSGrp>=rsSGrp
c           exists
c        2) in individual records p">=q",r">=s" are stored
c           however p">=q")>=(r">=a") reduction is not yet used
c
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aSGrp,bSGrp,cSGrp,dSGrp,LunAux
        real*8 W(1)
        real*8 Wa(1)
c
c        help variables
        integer abPerm,cdPerm,abcdPerm
        integer abSGrp,cdSGrp,abcdSGrp
        integer abLen,cdLen,abcdLen
        integer dima,dimb,dimc,dimd
        integer pSGrp,qSGrp,rSGrp,sSGrp,dimp,dimq,dimr,dims
        character*10 LunName
        integer i
c
c
c        -------- part I - define basic parameters
c
c1        Def dima,dimb,dimc,dimc
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpbe(bSGrp)
        dimc=DimSGrpa(cSGrp)
        dimd=DimSGrpbe(dSGrp)
c
c
c        In steps 2.1 - 2.3 also dimp-dims and pSGrp-sSGrp
c2.1        Def abPerm, abSGrp
        if (aSGrp.ge.bSGrp) then
          abPerm=0
          abSGrp=(aSGrp*(aSGrp-1)/2)+bSGrp
          pSGrp=aSGrp
          qSGrp=bSGrp
          dimp=dima
          dimq=dimb
        else
          abPerm=1
          abSGrp=(bSGrp*(bSGrp-1)/2)+aSGrp
          pSGrp=bSGrp
          qSGrp=aSGrp
          dimp=dimb
          dimq=dima
        end if
c
c2.2        Def cdPerm, cdSGrp
        if (cSGrp.ge.dSGrp) then
          cdPerm=0
          cdSGrp=(cSGrp*(cSGrp-1)/2)+dSGrp
          rSGrp=cSGrp
          sSGrp=dSGrp
          dimr=dimc
          dims=dimd
        else
          cdPerm=1
          cdSGrp=(dSGrp*(dSGrp-1)/2)+cSGrp
          rSGrp=dSGrp
          sSGrp=cSGrp
          dimr=dimd
          dims=dimc
        end if
c
c2.3        Def abcdPerm, abcdSGrp
        if (abSGrp.ge.cdSGrp) then
          abcdPerm=0
          abcdSGrp=(abSGrp*(abSGrp-1)/2)+cdSGrp
        else
          abcdPerm=1
          abcdSGrp=(cdSGrp*(cdSGrp-1)/2)+abSGrp
          i=pSGrp
          pSGrp=rSGrp
          rSGrp=i
          i=qSGrp
          qSGrp=sSGrp
          sSGrp=i
          i=dimp
          dimp=dimr
          dimr=i
          i=dimq
          dimq=dims
          dims=i
        end if
c
c
c3.1        Def abLen
        if (aSGrp.eq.bSGrp) then
          abLen=dima*(dima+1)/2
        else
          abLen=dima*dimb
        end if
c
c3.2        Def cdLen
        if (cSGrp.eq.dSGrp) then
          cdLen=dimc*(dimc+1)/2
        else
          cdLen=dimc*dimd
        end if
c
c3.3        Def abcdLen
        abcdLen=abLen*cdLen
c
c
c        -------- part II - read proper integrals (pq|rs) from disc
c
c4        Def proper LunName
        call MkNameV4 (pSGrp,qSGrp,rSGrp,sSGrp,'W4',LunName)
c
c
c5        Read Wa(pqrs) - velka udalost
C       open (unit=LunAux,file=LunName,form='unformatted')
        Call Molcas_BinaryOpen_Vanilla(LunAux,LunName)
        call rea_chcc (LunAux,abcdLen,Wa(1))
        close (LunAux)
c
c
c        -------- part III - Upgrade W(a",b",c",d") from Wx(pqrs)
c        symbolics used: [xy] means xy or yx
c
        if (abcdPerm.eq.0) then
c          case ([ab]|[cd])
c
          if ((abPerm.eq.0).and.(cdPerm.eq.0)) then
c            subcase (ab|cd)
            call DefW4abcd (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else if ((abPerm.eq.0).and.(cdPerm.eq.1)) then
c            subcase (ab|dc)
            call DefW4abdc (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else if ((abPerm.eq.1).and.(cdPerm.eq.0)) then
c            subcase (ba|cd)
            call DefW4bacd (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else
c            subcase (ba|cd)
            call DefW4badc (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          end if
c
        else
c        case ([cd]|[ab])
c
          if ((abPerm.eq.0).and.(cdPerm.eq.0)) then
c            subcase (cd|ab)
            call DefW4cdab (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else if ((abPerm.eq.0).and.(cdPerm.eq.1)) then
c            subcase (cd|ab)
            call DefW4dcab (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else if ((abPerm.eq.1).and.(cdPerm.eq.0)) then
c            subcase (cd|ba)
            call DefW4cdba (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else
c            subcase (dc|ba)
            call DefW4dcba (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          end if
c
        end if
c
c
        return
        end
c
c        -------------------
c
        subroutine DefW4abcd (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (ab|cd)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:abLen,1:cdLen)
c
c        help variables
        integer a,b,c,d,ab,cd
c
        if ((aSGrp.eq.bSGrp).and.(cSGrp.eq.dSGrp)) then
c        case (a=b|c=d)
          do c=2,dimc
          cd=c*(c-1)/2
          do d=1,c-1
          cd=cd+1
            do a=2,dima
            ab=a*(a-1)/2
            do b=1,a-1
            ab=ab+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ab,cd)
              W(b,a,c,d)=W(b,a,c,d)+Wx(ab,cd)
              W(a,b,d,c)=W(a,b,d,c)+Wx(ab,cd)
              W(b,a,d,c)=W(b,a,d,c)+Wx(ab,cd)
            end do
            end do
            do a=1,dima
            ab=a*(a+1)/2
              W(a,a,c,d)=W(a,a,c,d)+Wx(ab,cd)
              W(a,a,d,c)=W(a,a,d,c)+Wx(ab,cd)
            end do
          end do
          end do
c
          do c=1,dimc
          cd=c*(c+1)/2
            do a=2,dima
            ab=a*(a-1)/2
            do b=1,a-1
            ab=ab+1
              W(a,b,c,c)=W(a,b,c,c)+Wx(ab,cd)
              W(b,a,c,c)=W(b,a,c,c)+Wx(ab,cd)
            end do
            end do
            do a=1,dima
            ab=a*(a+1)/2
              W(a,a,c,c)=W(a,a,c,c)+Wx(ab,cd)
            end do
          end do
c
        else if ((aSGrp.eq.bSGrp).and.(cSGrp.ne.dSGrp)) then
c        case (a=b|c,d)
          cd=0
          do d=1,dimd
          do c=1,dimc
          cd=cd+1
            do a=2,dima
            ab=a*(a-1)/2
            do b=1,a-1
            ab=ab+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ab,cd)
              W(b,a,c,d)=W(b,a,c,d)+Wx(ab,cd)
            end do
            end do
            do a=1,dima
            ab=a*(a+1)/2
              W(a,a,c,d)=W(a,a,c,d)+Wx(ab,cd)
            end do
          end do
          end do
c
        else if ((aSGrp.ne.bSGrp).and.(cSGrp.eq.dSGrp)) then
c        case (a,b|c=d)
          do c=2,dimc
          do d=1,c-1
          cd=c*(c-1)/2+d
            ab=0
            do b=1,dimb
            do a=1,dima
            ab=ab+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ab,cd)
              W(a,b,d,c)=W(a,b,d,c)+Wx(ab,cd)
            end do
            end do
          end do
          end do
c
          do c=1,dimc
          cd=c*(c+1)/2
            ab=0
            do b=1,dimb
            do a=1,dima
            ab=ab+1
              W(a,b,c,c)=W(a,b,c,c)+Wx(ab,cd)
            end do
            end do
          end do
c
        else
c        case (a,b|c,d)
          cd=0
          do d=1,dimd
          do c=1,dimc
          cd=cd+1
            ab=0
            do b=1,dimb
            do a=1,dima
            ab=ab+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ab,cd)
            end do
            end do
          end do
          end do
c
        end if
c
c
        return
        end
c
c        ---------------------
c
        subroutine DefW4abdc (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (ab|dc)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:abLen,1:cdLen)
c
c        help variables
        integer a,b,c,d,ab,dc
c
        if (aSGrp.eq.bSGrp) then
c        case (a=b|d,c)
          dc=0
          do c=1,dimc
          do d=1,dimd
          dc=dc+1
            do a=2,dima
            ab=a*(a-1)/2
            do b=1,a-1
            ab=ab+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ab,dc)
              W(b,a,c,d)=W(b,a,c,d)+Wx(ab,dc)
            end do
            end do
            do a=1,dima
            ab=a*(a+1)/2
              W(a,a,c,d)=W(a,a,c,d)+Wx(ab,dc)
            end do
          end do
          end do
c
        else
c        case (a,b|d,c)
          dc=0
          do c=1,dimc
          do d=1,dimd
          dc=dc+1
            ab=0
            do b=1,dimb
            do a=1,dima
            ab=ab+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ab,dc)
            end do
            end do
          end do
          end do
c
        end if
c
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(cSGrp)
        call Unused_integer(dSGrp)
      end if
        end
c
c        ---------------------
c
        subroutine DefW4bacd (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (ba|cd)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:abLen,1:cdLen)
c
c        help variables
        integer a,b,c,d,ba,cd
c
        if (cSGrp.eq.dSGrp) then
c        case (b,a|c=d)
          do c=2,dimc
          cd=c*(c-1)/2
          do d=1,c-1
          cd=cd+1
            ba=0
            do a=1,dima
            do b=1,dimb
            ba=ba+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ba,cd)
              W(a,b,d,c)=W(a,b,d,c)+Wx(ba,cd)
            end do
            end do
          end do
          end do
c
          do c=1,dimc
          cd=c*(c+1)/2
            ba=0
            do a=1,dima
            do b=1,dimb
            ba=ba+1
              W(a,b,c,c)=W(a,b,c,c)+Wx(ba,cd)
            end do
            end do
          end do
c
        else
c        case (b,a|c,d)
          cd=0
          do d=1,dimd
          do c=1,dimc
          cd=cd+1
            ba=0
            do a=1,dima
            do b=1,dimb
            ba=ba+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ba,cd)
            end do
            end do
          end do
          end do
c
        end if
c
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(aSGrp)
        call Unused_integer(bSGrp)
      end if
        end
c
c        ---------------------
c
        subroutine DefW4badc (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (ba|dc)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:abLen,1:cdLen)
c
c        help variables
        integer a,b,c,d,ba,dc
c
c        case (b,a|c,d)
          dc=0
          do c=1,dimc
          do d=1,dimd
          dc=dc+1
            ba=0
            do a=1,dima
            do b=1,dimb
            ba=ba+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(ba,dc)
            end do
            end do
          end do
          end do
c
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(aSGrp)
        call Unused_integer(bSGrp)
        call Unused_integer(cSGrp)
        call Unused_integer(dSGrp)
      end if
        end
c
c        ---------------------
c
        subroutine DefW4cdab (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (cd|ab)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:cdLen,1:abLen)
c
c        help variables
        integer a,b,c,d,ab,cd
c
        if ((aSGrp.eq.bSGrp).and.(cSGrp.eq.dSGrp)) then
c        case (c=d|a=b)
          do a=2,dima
          ab=a*(a-1)/2
          do b=1,a-1
          ab=ab+1
            do c=2,dimc
            cd=c*(c-1)/2
            do d=1,c-1
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ab)
              W(b,a,c,d)=W(b,a,c,d)+Wx(cd,ab)
              W(a,b,d,c)=W(a,b,d,c)+Wx(cd,ab)
              W(b,a,d,c)=W(b,a,d,c)+Wx(cd,ab)
            end do
            end do
            do c=1,dimc
            cd=c*(c+1)/2
              W(a,b,c,c)=W(a,b,c,c)+Wx(cd,ab)
              W(b,a,c,c)=W(b,a,c,c)+Wx(cd,ab)
            end do
          end do
          end do
c
          do a=1,dima
          ab=a*(a+1)/2
            do c=2,dimc
            cd=c*(c-1)/2
            do d=1,c-1
            cd=cd+1
              W(a,a,c,d)=W(a,a,c,d)+Wx(cd,ab)
              W(a,a,d,c)=W(a,a,d,c)+Wx(cd,ab)
            end do
            end do
            do c=1,dimc
            cd=c*(c+1)/2
              W(a,a,c,c)=W(a,a,c,c)+Wx(cd,ab)
            end do
          end do
c
        else if ((aSGrp.eq.bSGrp).and.(cSGrp.ne.dSGrp)) then
c        case (c,d|a=b)
          do a=2,dima
          ab=a*(a-1)/2
          do b=1,a-1
          ab=ab+1
            cd=0
            do d=1,dimd
            do c=1,dimc
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ab)
              W(b,a,c,d)=W(b,a,c,d)+Wx(cd,ab)
            end do
            end do
          end do
          end do
c
          do a=1,dima
          ab=a*(a+1)/2
            cd=0
            do d=1,dimd
            do c=1,dimc
            cd=cd+1
              W(a,a,c,d)=W(a,a,c,d)+Wx(cd,ab)
            end do
            end do
          end do
c
        else if ((aSGrp.ne.bSGrp).and.(cSGrp.eq.dSGrp)) then
c        case (c=d|a,b)
          ab=0
          do b=1,dimb
          do a=1,dima
          ab=ab+1
            do c=2,dimc
            cd=c*(c-1)/2
            do d=1,c-1
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ab)
              W(a,b,d,c)=W(a,b,d,c)+Wx(cd,ab)
            end do
            end do
            do c=1,dimc
            cd=c*(c+1)/2
              W(a,b,c,c)=W(a,b,c,c)+Wx(cd,ab)
            end do
          end do
          end do
c
        else
c        case (c,d|a,b)
          ab=0
          do b=1,dimb
          do a=1,dima
          ab=ab+1
            cd=0
            do d=1,dimd
            do c=1,dimc
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ab)
            end do
            end do
          end do
          end do
c
        end if
c
c
        return
        end
c
c        ---------------------
c
        subroutine DefW4dcab (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (dc|ab)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:cdLen,1:abLen)
c
c        help variables
        integer a,b,c,d,ab,dc
c
        if (aSGrp.eq.bSGrp) then
c        case (d,c|a=b)
          do a=2,dima
          ab=a*(a-1)/2
          do b=1,a-1
          ab=ab+1
            dc=0
            do c=1,dimc
            do d=1,dimd
            dc=dc+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(dc,ab)
              W(b,a,c,d)=W(b,a,c,d)+Wx(dc,ab)
            end do
            end do
          end do
          end do
c
          do a=1,dima
          ab=a*(a+1)/2
            dc=0
            do c=1,dimc
            do d=1,dimd
            dc=dc+1
              W(a,a,c,d)=W(a,a,c,d)+Wx(dc,ab)
            end do
            end do
          end do
c
        else
c        case (d,c|a,b)
          ab=0
          do b=1,dimb
          do a=1,dima
          ab=ab+1
            dc=0
            do c=1,dimc
            do d=1,dimd
            dc=dc+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(dc,ab)
            end do
            end do
          end do
          end do
c
        end if
c
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(cSGrp)
        call Unused_integer(dSGrp)
      end if
        end
c
c        -------------------
c
        subroutine DefW4cdba (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (cd|ba)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:cdLen,1:abLen)
c
c        help variables
        integer a,b,c,d,ba,cd
c
        if (cSGrp.eq.dSGrp) then
c        case (c=d|b,a)
          ba=0
          do a=1,dima
          do b=1,dimb
          ba=ba+1
            do c=2,dimc
            cd=c*(c-1)/2
            do d=1,c-1
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ba)
              W(a,b,d,c)=W(a,b,d,c)+Wx(cd,ba)
            end do
            end do
            do c=1,dimc
            cd=c*(c+1)/2
              W(a,b,c,c)=W(a,b,c,c)+Wx(cd,ba)
            end do
          end do
          end do
c
        else
c        case (c,d|b,a)
          ba=0
          do a=1,dima
          do b=1,dimb
          ba=ba+1
            cd=0
            do d=1,dimd
            do c=1,dimc
            cd=cd+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(cd,ba)
            end do
            end do
          end do
          end do
c
        end if
c
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(aSGrp)
        call Unused_integer(bSGrp)
      end if
        end
c
c        -------------------
c
        subroutine DefW4dcba (W,Wx,dima,dimb,dimc,dimd,abLen,cdLen,
     c                        aSGrp,bSGrp,cSGrp,dSGrp)
c
c        define W(a,b,c,d) from (dc|ba)
c
        implicit none
        integer dima,dimb,dimc,dimd,abLen,cdLen,aSGrp,bSGrp,cSGrp,dSGrp
        real*8 W(1:dima,1:dimb,1:dimc,1:dimd)
        real*8 Wx(1:cdLen,1:abLen)
c
c        help variables
        integer a,b,c,d,ba,dc
c
c        case (c,d|b,a)
          ba=0
          do a=1,dima
          do b=1,dimb
          ba=ba+1
            dc=0
            do c=1,dimc
            do d=1,dimd
            dc=dc+1
              W(a,b,c,d)=W(a,b,c,d)+Wx(dc,ba)
            end do
            end do
          end do
          end do
c
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(aSGrp)
        call Unused_integer(bSGrp)
        call Unused_integer(cSGrp)
        call Unused_integer(dSGrp)
      end if
        end
c
c        -------------------
c
        subroutine MkNameV4 (i,j,k,l,Schem,Nomen)
c
c       help routine to ReaW4, producing name of V4file
c       ex: Schem='XY', i=1, j=3, k=5, l=07 ->  Nomen='XY01030507'
c
        implicit none
        integer i,j,k,l
        character*2 Schem
        character*10 Nomen
c
c       help variables
        character*1 Chr(1:10)
        character*2 digit(1:64)
        character*2 ichr,jchr,kchr,lchr
        character*2 baza
        character*10 meno
c
        equivalence (Chr(1),meno)
        equivalence (Chr(1),baza)
        equivalence (Chr(3),ichr)
        equivalence (Chr(5),jchr)
        equivalence (Chr(7),kchr)
        equivalence (Chr(9),lchr)
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
        kchr=digit(k)
        lchr=digit(l)
        Nomen=meno
c
        return
        end
c
c       ------------------------------------
c
c
        subroutine Xo2v4ctl (NvGrp,NvSGrp,LunAux)
c
c!      drajver procesu na testovanie ktory W3/W4 file
c        treba na ktorom node
c        N.B. upraveny drajver o2v4 procesu
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "chcc_parcc.fh"
#include "para_info.fh"
c
        integer NvGrp,NvSGrp,LunAux
c
c        help variables
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
c##
        integer aGrp,bGrp,proc,i,j
        integer nJobs,addJobs,actJobs
c
c
c1      Inicializacia premennych (Predbezna)
c
        NaGrp=NvGrp
        NbeGrp=NvGrp
        NaSGrp=NvSGrp
        nbeSGrp=NvSGrp
c
c
c2      define all groups and ssungroup parameters
c
        call DefParo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
c
c3        distribute work among nodes (def ABID)
c
        do proc=0,(nProcs-1)
        do aGrp=1,NaGrp
        do bGrp=1,NaGrp
          ABID(proc,aGrp,bGrp)=0
        end do
        end do
        end do
c
c
        if ((nProcs.eq.1).or.(NaGrp.eq.1)) then
c3.1        nP = <1>
c        all jobs to one node
c
          do aGrp=1,NaGrp
          do bGrp=1,NaGrp
            ABID(0,aGrp,bGrp)=1
          end do
          end do
c
c
        else
c3.2        nP = <2, inf)
c        N.B. cele zatial trochu odflaknute, lebo sa neberie do uvahy,
c        ze pre half joby je lepsie ked nadvazuju na plne joby s
c        rovnakymi indexami
c
c3.2.1        Full Jobs
c
          i=(NaGrp*(NaGrp-1))/2
          nJobs=int(i/nProcs)
          addJobs=mod(i,nProcs)
c
c          first nodes: 0-addJobs
c             will have int(N'(N'-1)/2 /nProc)+1 Full Jobs
c          rest nodes: addJobs-nProcs-1
c               will have int(N'(N'-1)/2 /nProc)
c
          proc=0
          actJobs=nJobs
          do aGrp=2,NaGrp
          do bGrp=1,aGrp-1
            ABID(proc,aGrp,bGrp)=1
            actJobs=actJobs-1
            if (actJobs.eq.-1) then
              proc=proc+1
              actJobs=nJobs
            else if (actJobs.eq.0) then
              if (addJobs.gt.0) then
                addJobs=addJobs-1
              else
                proc=proc+1
                actJobs=nJobs
              end if
            end if
          end do
          end do
c
c
c3.2.2        Half Jobs
c
c        distribution of Half Jobs
c          - first distribution - addjobs-nProcs-1
c          - second distribution - addjobs-nProcs-1
c          - other distributions - 0 - nProcs-1
c
          addJobs=mod(i,nProcs)
          proc=addJobs
          j=0
c
          do aGrp=1,NaGrp
            ABID(proc,aGrp,aGrp)=1
            proc=proc+1
            if (proc.eq.nProcs) then
              if (j.eq.0) then
c              first distribution - addjobs-nProcs-1
                proc=addJobs
                j=1
              else if (j.eq.1) then
c              second distribution - addjobs-nProcs-1
                proc=addJobs
                j=2
              else
c              other distributions - 0 - nProcs-1
                proc=0
              end if
            end if
          end do
c
        end if
c
c
c        Printing ABID
        do proc=0,nProcs-1
        if (printkey.ge.10) then
        write (6,*)   ' For myRank = ',proc
        end if
          do aGrp=1,NaGrp
          do bGrp=1,aGrp
          if (ABID(proc,aGrp,bGrp).eq.1) then
          if (printkey.ge.10) then
          write (6,*) '    aGrp,bGrp ',aGrp,bGrp
          end if
          end if
          end do
          end do
        end do
c
c
c
c4      A ideme na to
c
        call InsReqo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c             mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
c
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_integer(LunAux)
        end
c
c       -------------------
c
        subroutine InsReqo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
c
c       this routine do:
c        Inspect W3 and W4 files requirements of o2v4 procedure
c        on this node. It checks which of the W3 and W4 files
c        are used on this node
c
        implicit none
#include "chcc1.fh"
#include "chcc_parcc.fh"
#include "para_info.fh"
#include "o2v4.fh"
c
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
c
c        help variables
c
        integer aGrp,bGrp,gaGrp,beGrp,aSGrp,bSGrp,gaSGrp,beSGrp
        integer dima,dimb,adda,addb
        integer dimbe,dimga,addbe,addga,addbepp,addgapp
        integer bSGrpUp,gaSGrpUp
        integer LenW3,LenW4
        integer i,j,NSGrp
c
c
c*        Initial Set InqW3,InqW4 = F
        NSGrp=NaGrp*NaSGrp
        do i=1,(NSGrp*(NSGrp+1))/2
          do j=1,NSGrp
          InqW3(i,j)=.False.
          end do
          do j=1,(NSGrp*(NSGrp+1))/2
          InqW4(i,j)=.False.
          end do
        end do
        LenW3=0
        LenW4=0
c
c
c*      cycle over all groups of (a>=b)
        adda=0
        do aGrp=1,naGrp
        dima=DimGrpa(aGrp)
c
c##        test, if on this node atleast one combination with this
c        aGrp is scheduled. Skip if no
        i=0
        do j=1,NaGrp
          i=i+ABID(myRank,aGrp,j)
        end do
        if (i.eq.0) goto 12
c
c
        addb=0
        do bGrp=1,aGrp
        dimb=DimGrpa(bGrp)
c
c##     test, if this a'b' combination is planed to be run on this node
        if (ABID(myRank,aGrp,bGrp).eq.0) goto 11
c
c
cx        cycle over all groups of (be>=ga)
          addbe=0
          do beGrp=1,nbeGrp
          dimbe=DimGrpbe(beGrp)
          addga=0
          do gaGrp=1,beGrp
          dimga=DimGrpbe(gaGrp)
c
c
cxx         cycle over all subgroups of (be>=ga)'
            addbepp=addbe
            do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
c
            if (beGrp.eq.gaGrp) then
              gaSGrpUp=beSGrp
            else
              gaSGrpUp=GrpbeUp(gaGrp)
            end if
c
            addgapp=addga
            do gaSGrp=GrpbeLow(gaGrp),gaSGrpUp
c
c
cxxx          cycle over all subgroups of (a>=b)'
              do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
              if (aGrp.eq.bGrp) then
                bSGrpUp=aSGrp
              else
                bSGrpUp=GrpaUp(bGrp)
              end if
              do bSGrp=GrpaLow(bGrp),bSGrpUp
c
c
c*              clasical reading or generating of (VV|VV) integrals
c
cxxxx                term W1(a",be",b",ga") <- -(b"ga"|a"i) . T1(i,be")
cxxxx1          Get Ww(b",ga",a",i) (Wx used as Aux)
                call InsReaW3 (bSGrp,gaSGrp,aSGrp,LenW3)
c
cxxxx           Upgrade W1(a",be",b",ga") <<- (a",be"|b",ga")
                call InsReaW4 (aSGrp,beSGrp,bSGrp,gaSGrp,LenW4)
c
c
                if (aSGrp.ne.bSGrp) then
cxxxx                  term W2(b",be",a",ga") <- -(a"ga"|b"i) . T1(i,ga")
cxxxx1            Get Ww(a",ga",b",i) (Wx used as Aux)
                  call InsReaW3 (aSGrp,gaSGrp,bSGrp,LenW3)
c
cxxxx                  term W2(b",be",a",ga") <<- -(b"be"|a"i) . T1(i,ga")
cxxxx1            Get Ww(b",be",a",i) (Wx used as Aux)
                  call InsReaW3 (bSGrp,beSGrp,aSGrp,LenW3)
c
cxxxx             Add W2(b",be",a",ga") <<- (b",be"Ia",ga")
c                  Upgrade:W2, destroy:Ww
                  call InsReaW4 (bSGrp,beSGrp,aSGrp,gaSGrp,LenW4)
                end if
c
c
c
cxxx          end cycle over all subgroups of (a>=b)'
              end do
              end do
c
cxx         end cycle over all subgroups of (be>=ga)'
            addgapp=addgapp+DimSGrpbe(gaSGrp)
            end do
            addbepp=addbepp+DimSGrpbe(beSGrp)
            end do
c
cx        end cycle over all groups of (be>=ga)
          addga=addga+dimga
          end do
          addbe=addbe+dimbe
          end do
c
c*      end cycle over all groups of (a>=b)
11        addb=addb+dimb
        end do
12        adda=adda+dima
        end do
c
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(NbeSgrp)
        call Unused_integer(mdGrpa)
        call Unused_integer(mdGrpbe)
        call Unused_integer(mdSGrpa)
        call Unused_integer(mdSGrpbe)
      end if
        end
c
c        -------------------------
c
        subroutine InsReaW3 (aSGrp,beSGrp,bSGrp,length)
c
c        this routine do:
c        Check which W3 file corresponds to given combination
c        of indexes,
c        - increase the parameter, checking the overal number of
c        required integrals if these block is conted first time
c        - and set corresponding InqW3 to T
c
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "chcc_parcc.fh"
c
        integer aSGrp,beSGrp,bSGrp,length
c
c        help variables
c
        integer dima,dimb,dimbe
        integer abeSGrp
c
c
c0        def dimensions
c
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpa(bSGrp)
        dimbe=DimSGrpbe(beSGrp)
c
        if (aSGrp.gt.beSGrp) then
c1        case aSGrp>beSGrp, integrals (a",be"|b",i) to be read
          abeSGrp=(aSGrp*(aSGrp-1)/2)+beSGrp
          if (InqW3(abeSGrp,bSGrp).eqv..False.) then
            InqW3(abeSGrp,bSGrp)=.True.
            length=length+dima*dimbe*dimb*no
          end if
c
        else if (aSGrp.eq.beSGrp) then
c2        case aSGrp=beSGrp, integrals (a">=be"|b",i)to be read and expand
          abeSGrp=(aSGrp*(aSGrp-1)/2)+beSGrp
          if (InqW3(abeSGrp,bSGrp).eqv..False.) then
            InqW3(abeSGrp,bSGrp)=.True.
            length=length+(no*dima*(dima+1)*dimb)/2
          end if
c
        else
c3        case aSGrp<beSGrp, integrals (be",a"|b",i) to be read and mapped
          abeSGrp=(beSGrp*(beSGrp-1)/2)+aSGrp
          if (InqW3(abeSGrp,bSGrp).eqv..False.) then
            InqW3(abeSGrp,bSGrp)=.True.
            length=length+dima*dimbe*dimb*no
          end if
        end if
c
        return
        end
c
c        -------------------
c
        subroutine InsReaW4 (aSGrp,bSGrp,cSGrp,dSGrp,length)
c
c        this routine do:
c        Check which W4 file corresponds to given combination
c        of indexes,
c        - increase the parameter, checking the overal number of
c        required integrals if these block is conted first time
c        - and set corresponding InqW4 to T
c
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "chcc_parcc.fh"
c
        integer aSGrp,bSGrp,cSGrp,dSGrp,length
c
c        help variables
        integer abPerm,cdPerm,abcdPerm
        integer abSGrp,cdSGrp,abcdSGrp
        integer abLen,cdLen,abcdLen
        integer dima,dimb,dimc,dimd
        integer pSGrp,qSGrp,rSGrp,sSGrp,dimp,dimq,dimr,dims
        integer pqSGrp,rsSGrp
        integer i
c
c
c        -------- part I - define basic parameters
c
c1        Def dima,dimb,dimc,dimc
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpbe(bSGrp)
        dimc=DimSGrpa(cSGrp)
        dimd=DimSGrpbe(dSGrp)
c
c
c        In steps 2.1 - 2.3 also dimp-dims and pSGrp-sSGrp
c2.1        Def abPerm, abSGrp
        if (aSGrp.ge.bSGrp) then
          abPerm=0
          abSGrp=(aSGrp*(aSGrp-1)/2)+bSGrp
          pSGrp=aSGrp
          qSGrp=bSGrp
          dimp=dima
          dimq=dimb
        else
          abPerm=1
          abSGrp=(bSGrp*(bSGrp-1)/2)+aSGrp
          pSGrp=bSGrp
          qSGrp=aSGrp
          dimp=dimb
          dimq=dima
        end if
c
c2.2        Def cdPerm, cdSGrp
        if (cSGrp.ge.dSGrp) then
          cdPerm=0
          cdSGrp=(cSGrp*(cSGrp-1)/2)+dSGrp
          rSGrp=cSGrp
          sSGrp=dSGrp
          dimr=dimc
          dims=dimd
        else
          cdPerm=1
          cdSGrp=(dSGrp*(dSGrp-1)/2)+cSGrp
          rSGrp=dSGrp
          sSGrp=cSGrp
          dimr=dimd
          dims=dimc
        end if
c
c2.3        Def abcdPerm, abcdSGrp
        if (abSGrp.ge.cdSGrp) then
          abcdPerm=0
          abcdSGrp=(abSGrp*(abSGrp-1)/2)+cdSGrp
        else
          abcdPerm=1
          abcdSGrp=(cdSGrp*(cdSGrp-1)/2)+abSGrp
          i=pSGrp
          pSGrp=rSGrp
          rSGrp=i
          i=qSGrp
          qSGrp=sSGrp
          sSGrp=i
          i=dimp
          dimp=dimr
          dimr=i
          i=dimq
          dimq=dims
          dims=i
        end if
c
c
c3.1        Def abLen
        if (aSGrp.eq.bSGrp) then
          abLen=dima*(dima+1)/2
        else
          abLen=dima*dimb
        end if
c
c3.2        Def cdLen
        if (cSGrp.eq.dSGrp) then
          cdLen=dimc*(dimc+1)/2
        else
          cdLen=dimc*dimd
        end if
c
c3.3        Def abcdLen
        abcdLen=abLen*cdLen
c
c
c        -------- part II - read proper integrals (pq|rs) from disc
c
        pqSGrp=(pSGrp*(pSGrp-1)/2)+qSGrp
        rsSGrp=(rSGrp*(rSGrp-1)/2)+sSGrp
c
        if (InqW4(pqSGrp,rsSGrp).eqv..False.) then
          InqW4(pqSGrp,rsSGrp)=.True.
          length=length+abcdLen
        end if
c
c
        return
        end
