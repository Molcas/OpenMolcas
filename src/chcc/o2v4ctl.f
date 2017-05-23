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
        subroutine o2v4ctl (wrk,wrksize,NvGrp,NvSGrp,LunAux)
c
c!      drajver o2v4 procesu
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "wrk.fh"
#include "chcc_parcc.fh"
#include "para_info.fh"
#include "chcc_casy.fh"
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
c@c
c        return
c@c
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
        if (printkey.ge.10) then
        write (6,*) NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
        end if
cmp
        Call CWTime(TCpu,TWall)
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
cmp
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
        call o2v4 (wrk(1),wrksize,
     c             NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c             mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,
     c             LunAux)

c
        return
        end
c
c       -------------------
c
        subroutine o2v4 (wrk,wrksize,
     c                   NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,
     c                   LunAux)
c
c       this routine do:
c
c       t(u,v,bega) <- sum(a,b) [ b(a,b,be,ga) . Tau(u,v,a,b) ]
c
c       b(a,b,be,ga) =             (a,be|b,ga)
c                    - sum (i)   [ (a,be|b,i) . t(i,ga) ]
c                    - sum (i)   [ (a,i|b,ga) . t(i,be) ]
c                    + sum (i,j) [ (a,i|b,j) . t(i,be) . t(j,ga) ]
c
c       Tau(u,v,ab) = t(u,v,ab) + t(i,a) . t(j,b) : stored for: i,j,a>=b
c
c
cA      List of main arrays:
c
c*      tau(i,j,(a,b)') - PossTau
c
c*      t2n1(ij,(bega)")- PossT2n1
c       t2n2(ij,(bega)")- PossT2n2
c       t2w(ij,(ab)")   - PossT2w - work file for T(+-)
c
c*      L21 (m,a',be')  - PossL21
c       L22 (m,a',ga')  - PossL22
c       L23 (m,b',be')  - PossL23
c       L24 (m,b',ga')  - PossL24
c        L2W (m,c',de')  - PossL2W - Used for Tmp file in GetCHV
c       some (or all) of L2s might be identical
c
c*      M1(m,a",be")  - PossM1, used also for M(m,b'',be'')
c       M2(m,b",ga")  - PossM2, used also for M(m,a'',ga'')
c       Both Ms might be identical (mozno + M3 ako manipulacne)
c
c*      W1(a",be",b",ga") - PossW1 - b(a",be"|b",ga")
c       W2(b",be",a",ga") - PossW2 - b(b",be"|a",ga")
c       Ww((ab)",(bega)") - PossWw - work file for b(+-)
c
c
cB      list of additional arrays
c
cC      list of used pre- and suffixes and assignements
c
c       _Grp         - Group
c       _SGrp        - SubGroup
c       Poss_        - Possition
c       PsAc_        - Actual Possition
c       (...)'       - Group, Block
c       (...)"       - SubGroup, SubBlock
c
cD      list of most important variables
c
        implicit none
#include "wrk.fh"
#include "chcc1.fh"
#include "chcc_parcc.fh"
#include "para_info.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp,LunAux

        integer PossTau,PossT2n1,PossT2n2,PossT2w
        integer PossL11,PossL12
        integer PossL21,PossL22,PossL23,PossL24,PossL2W
        integer PossH1,PossH2
        integer PossM1,PossM2,PossW1,PossW2,PossWw,PossWx
        integer PsAcL21,PsAcL22,PsAcL23,PsAcL24
c       integer PsAcM1,PsAcM2
        integer pL21,pL22,pL23,pL24
        integer L2Status(1:4,1:3)
        integer PossT,dim1,dim2,dim3,lent2n1,lent2n2
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,NL2
c
        integer aGrp,bGrp,gaGrp,beGrp,aSGrp,bSGrp,gaSGrp,beSGrp
        integer dima,dimb,adda,addb
        integer dimbe,dimga,addbe,addga,addbepp,addgapp
        integer bSGrpUp,gaSGrpUp
        integer i,j
        integer choleskikey
        integer FirstT2n(1:maxSGrp,1:maxSGrp)
        character*6 LunName
c
        if (intkey.eq.1) then
          choleskikey=0
        else
          choleskikey=1
        end if
c
cx      distribute memory
c
        PossT=PossFree
        call DistMemo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     &                 mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,
     &                 PossTau,PossT2n1,PossT2n2,PossT2w,
     &                 PossL11,PossL12,
     &                 PossL21,PossL22,PossL23,PossL24,PossL2W,
     &                 PossH1,PossH2,
     &                 PossM1,PossM2,PossW1,PossW2,PossWw,PossWx,
     &                 PossT,NL2)
        if (printkey.ge.10) then
        write (6,*) ' Last Value :',PossT,wrksize
        end if
        if (PossT.gt.wrksize) then
cmp!          write (6,*) ' Nieje dobre - o2v4, Dr. Ch. Kokotopuloss',
          write (6,*) ' Not Enough memory in o2v4 step!'//
     & 'Increase large and/or small segmentation ',
     &    (1.0d0*PossT)/(1.0d0*wrksize)
          call abend
        end if
c
c@@
c        call Calc_Bc
c@@
c
cx      initialize L2Status
c
        do i=1,NL2
        do j=1,3
          L2Status(i,j)=0
        end do
        end do
c
        if (NL2.eq.1) then
          L2Status(1,3)=possL21
        else if (NL2.eq.2) then
          if (NbeGrp.eq.1) then
            L2Status(1,3)=possL21
            L2Status(2,3)=possL23
          else
            L2Status(1,3)=possL21
            L2Status(2,3)=possL22
          end if
        else
          L2Status(1,3)=possL21
          L2Status(2,3)=possL22
          L2Status(3,3)=possL23
          L2Status(4,3)=possL24
        end if
c
c
cx        initialize FirstT2n and Generate Tmp3Name s
c##        initialize T2o2v4yes
c
c
        do beGrp=1,nbeGrp
        do gaGrp=1,beGrp
c
          do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
          if (beGrp.eq.gaGrp) then
            gaSGrpUp=beSGrp
          else
            gaSGrpUp=GrpbeUp(gaGrp)
          end if
          do gaSGrp=GrpbeLow(gaGrp),gaSGrpUp
c
          FirstT2n(beSGrp,gaSGrp)=1
c
          T2o2v4yes(beSGrp,gaSGrp)=0
c
          end do
          end do
c
        end do
        end do
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
        if (choleskikey.eq.1) then
cx          read L2W(m,i,a') <- L1(m,i,a')
          LunName=L1Name(aGrp)
          dim1=nc*dima*no
          call GetX (wrk(PossL2W),dim1,LunAux,LunName,1,1)
cx          Map L11(m,a',i) <- L2W(m,i,a')
          call Map3_132(wrk(PossL2W),wrk(PossL11),nc,no,dima)
        end if
c
        addb=0
        do bGrp=1,aGrp
        dimb=DimGrpa(bGrp)
c
c##     test, if this a'b' combination is planed to be run on this node
        if (ABID(myRank,aGrp,bGrp).eq.0) goto 11
        if (printkey.ge.10) then
        write (6,*) ' O2V4 MyRank, aGrp, bGrp',myRank, aGrp,bGrp
        end if
c
c
          if (choleskikey.eq.1) then
cx            read L12(m,b',i) <- L1(m,b',i)
              LunName=L1Name(bGrp)
            dim1=nc*dimb*no
            if (NaGrp.gt.1) then
              if (aGrp.eq.bGrp) then
                call mv0u (dim1,wrk(PossL11),1,wrk(PossL12),1)
              else
                call GetX (wrk(PossL2W),dim1,LunAux,LunName,1,1)
                call Map3_132(wrk(PossL2W),wrk(PossL12),nc,no,dimb)
              end if
            end if
          end if
c
c
cx        read the block of  amplitudes T2((ab)',ij) for given aGrp,bGrp
c          and make Tau from them
          call getTau (wrk(PossTau),wrk(PossT1o),
     c                 aGrp,bGrp,dima,dimb,adda,addb,LunAux)
c
cx        cycle over all groups of (be>=ga)
          addbe=0
          do beGrp=1,nbeGrp
          dimbe=DimGrpbe(beGrp)
          addga=0
          do gaGrp=1,beGrp
          dimga=DimGrpbe(gaGrp)
c
            if (choleskikey.eq.1) then
cxx         read Choleski vectors
c*          L21 (m,a',be')
c           L22 (m,a',ga')
c           L23 (m,b',be')
c           L24 (m,b',ga')
            call getChV (wrk,wrksize,
     c                     aGrp,bGrp,beGrp,gaGrp,NL2,L2Status,
     c                     pL21,pL22,pL23,pL24,PossL2W,
     c                     PossL11,PossL12,LunAux)
            PsAcL21=L2Status(pL21,3)
            PsAcL22=L2Status(pL22,3)
            PsAcL23=L2Status(pL23,3)
            PsAcL24=L2Status(pL24,3)
            end if
c
cxx         cycle over all subgroups of (be>=ga)'
            addbepp=addbe
            do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)
c
            if (intkey.eq.1) then
cxx              get H1(i,be") <- t1o(be,i)
              call Mk_T1t (wrk(PossT1o),wrk(PossH1),
     c        DimSGrpbe(beSGrp),no,nv,addbepp)
            end if
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
              if (intkey.eq.1) then
cxxx                get H2(i,ga") <- t1o(ga,i)
                call Mk_T1t (wrk(PossT1o),wrk(PossH2),
     c          DimSGrpbe(gaSGrp),no,nv,addgapp)
              end if
c
cxxx          vanish arrays for new (final) amplitudes, if this is a
c               first use, or read from Tmp3Name(be",ga") actual stage
c             T2n1((bega)",ij)
c             T2n2((bega)",ij)
              if (FirstT2n(beSGrp,gaSgrp).eq.1) then
                call VanishT2n (wrk(PossT2n1),wrk(PossT2n2),
     c                          beSGrp,gaSGrp)
                FirstT2n(beSGrp,gaSgrp)=0
              else
                call GetT2n (wrk(PossT2n1),wrk(PossT2n2),
     c                       beSGrp,gaSGrp,LunAux)
              end if
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
        if (printkey.ge.10) then
        write (6,99)aGrp,bGrp,beGrp,gaGrp,aSGrp,bSGrp,beSGrp,gaSGrp
        end if
99      format (8(i3,1x))
c
                if (choleskikey.eq.1) then
c
c*              choleski generation of (VV|VV) integrals
c
cxxxx             Extract M1(m,a",be") from L21(m,a',be')
                  call  ExtractM (wrk(PossM1),wrk(PsAcL21),
     c                            aGrp,beGrp,aSGrp,beSGrp)
cxxxx             Extract M2(m,b",ga") from L24(m,b',ga')
                  if ((aSGrp.eq.bSGrp).and.(beSGrp.eq.gaSGrp)) then
                    if (PossM1.ne.PossM2) then
                      dim1=nc*DimSGrpa(aSGrp)*DimSGrpbe(beSGrp)
                      call mv0u (dim1,wrk(PossM1),1,wrk(PossM2),1)
                    end if
                  else
                    call  ExtractM (wrk(PossM2),wrk(PsAcL24),
     c                              bGrp,gaGrp,bSGrp,gaSGrp)
                  end if
cxxxx             Calc W1(a",be",b",ga")=M1(m,a",be")(T).M2(m,b",ga")
                  dim1=DimSGrpa(aSGrp)*DimSGrpbe(beSGrp)
                  dim2=DimSGrpa(bSGrp)*DimSGrpbe(gaSGrp)
                  call mv0zero (dim1*dim2,dim1*dim2,wrk(possW1))
                  call mc0c1at3b(nc,dim1,nc,dim2,dim1,dim2,dim1,nc,dim2,
     c                          wrk(possM1),wrk(possM2),wrk(possW1))
c
                  if (aSGrp.ne.bSGrp) then
cxxxx             Extract M1(m,b",be") from L23(m,b',be')
                  call  ExtractM (wrk(PossM1),wrk(PsAcL23),
     c                          bGrp,beGrp,bSGrp,beSGrp)
cxxxx             Extract M2(m,a",ga") from L22(m,a',ga')
                  call  ExtractM (wrk(PossM2),wrk(PsAcL22),
     c                          aGrp,gaGrp,aSGrp,gaSGrp)
cxxxx             Calc W2(b",be",a",ga")=M1(m,b",be")(T).M2(m,a",ga")
                  dim1=DimSGrpa(bSGrp)*DimSGrpbe(beSGrp)
                  dim2=DimSGrpa(aSGrp)*DimSGrpbe(gaSGrp)
                  call mv0zero (dim1*dim2,dim1*dim2,wrk(possW2))
                  call mc0c1at3b (nc,dim1,nc,dim2,dim1,dim2,
     c                          dim1,nc,dim2,
     c                          wrk(possM1),wrk(possM2),wrk(possW2))
                  end if
c
                else
c
c*              clasical reading or generating of (VV|VV) integrals
c
cxxxx                term W1(a",be",b",ga") <- -(b"ga"|a"i) . T1(i,be")
c                Def:W1 destroy:Ww
cxxxx1          Get Ww(b",ga",a",i) (Wx used as Aux)
                call ReaW3 (wrk(PossWw),wrk(PossWx),
     c                      bSGrp,gaSGrp,aSGrp,LunAux)
cxxxx2                Map V1(a",ga",b",i) <- Ww(b",ga",a",i)
                call Map4_3214 (wrk(PossWw),wrk(PossW1),
     c                          DimSGrpa(bSGrp),DimSGrpbe(gaSGrp),
     c                          DimSGrpa(aSGrp),no)
cxxxx3          Calc Ww(a",ga",b",be") <- - W1(a",ga",b",i) . H1(i,be")
                dim1=DimSGrpa(bSGrp)*DimSGrpa(aSGrp)*DimSGrpbe(gaSGrp)
                dim2=DimSGrpbe(beSGrp)
                call mv0zero (dim1*dim2,dim1*dim2,wrk(PossWw))
                call mc0c2a3b (dim1,no,no,dim2,dim1,dim2,dim1,no,dim2,
     c                         wrk(PossW1),wrk(PossH1),wrk(PossWw))
cxxxx4                Map W1(a",be",b",ga") <<- Ww(a",ga",b",be")
                call Map4_1432 (wrk(PossWw),wrk(PossW1),
     c                          DimSGrpa(aSGrp),DimSGrpbe(gaSGrp),
     c                          DimSGrpa(bSGrp),DimSGrpbe(beSGrp))
c
cxxxx                term W1(a",be",b",ga") <<- -(a"be"|b"i) . T1(i,ga")
c                Upgrade: W1, destroy: Ww
cxxxx1          Get Ww(a",be",b",i) (Wx used as Aux)
                call ReaW3 (wrk(PossWw),wrk(PossWx),
     c                      aSGrp,beSGrp,bSGrp,LunAux)
cxxxx2          Add W1(a",be",b",ga") <<- - Ww(a",be",b",i) . H2(i,ga")
                dim1=DimSGrpa(aSGrp)*DimSGrpa(bSGrp)*DimSGrpbe(beSGrp)
                dim2=DimSGrpbe(gaSGrp)
                call mc0c2a3b (dim1,no,no,dim2,dim1,dim2,dim1,no,dim2,
     c                         wrk(PossWw),wrk(PossH2),wrk(PossW1))
c
c
cxxxx           Upgrade W1(a",be",b",ga") <<- (a",be"|b",ga")
c                Upgrade:W1 destroy:Ww
                call ReaW4 (wrk(PossW1),wrk(PossWw),
     c                     aSGrp,beSGrp,bSGrp,gaSGrp,LunAux)
c
c
                if (aSGrp.ne.bSGrp) then

cxxxx                  term W2(b",be",a",ga") <- -(a"ga"|b"i) . T1(i,ga")
c                  Def:W2 destroy: Ww
cxxxx1            Get Ww(a",ga",b",i) (Wx used as Aux)
                  call ReaW3 (wrk(PossWw),wrk(PossWx),
     c                        aSGrp,gaSGrp,bSGrp,LunAux)
cxxxx2                  Map W2(b",ga",a",i) <- Ww(a",ga",b",i) @@Nepreskusany
                  call Map4_3214 (wrk(PossWw),wrk(PossW2),
     c                            DimSGrpa(aSGrp),DimSGrpbe(gaSGrp),
     c                            DimSGrpa(bSGrp),no)
cxxxx3            Calc Ww(b",ga",a",be") <- -W2(b",ga",a",i) . H1(i,be")
                  dim1=DimSGrpa(aSGrp)*DimSGrpbe(gaSGrp)*DimSGrpa(bSGrp)
                  dim2=DimSGrpbe(beSGrp)
                  call mv0zero (dim1*dim2,dim1*dim2,wrk(PossWw))
                  call mc0c2a3b (dim1,no,no,dim2,dim1,dim2,dim1,no,dim2,
     c                           wrk(PossW2),wrk(PossH1),wrk(PossWw))
cxxxx4                  Map W2 (b",be",a",ga") <- Ww(b",ga",a",be")
                  call Map4_1432 (wrk(PossWw),wrk(PossW2),
     c                            DimSGrpa(bSGrp),DimSGrpbe(gaSGrp),
     c                            DimSGrpa(aSGrp),DimSGrpbe(beSGrp))
c
cxxxx                  term W2(b",be",a",ga") <<- -(b"be"|a"i) . T1(i,ga")
c                  Upgrade:W2 destroy: 0
cxxxx1            Get Ww(b",be",a",i) (Wx used as Aux)
                  call ReaW3 (wrk(PossWw),wrk(PossWx),
     c                        bSGrp,beSGrp,aSGrp,LunAux)
cxxxx2            Add W2(b",be",a",ga") <<- - Ww(b",be",a",i) .H2(i,ga")
                  dim1=DimSGrpa(aSGrp)*DimSGrpa(bSGrp)*DimSGrpbe(beSGrp)
                  dim2=DimSGrpbe(gaSGrp)
                  call mc0c2a3b (dim1,no,no,dim2,dim1,dim2,dim1,no,dim2,
     c                           wrk(PossWw),wrk(PossH2),wrk(PossW2))
c
cxxxx             Add W2(b",be",a",ga") <<- (b",be"Ia",ga")
c                  Upgrade:W2, destroy:Ww
                  call ReaW4 (wrk(PossW2),wrk(PossWw),
     c                       bSGrp,beSGrp,aSGrp,gaSGrp,LunAux)
                end if
c
                end if
c
c
c
cxxxx           Make t2w(+)((ab)",ij) from Tau((ab)',ij)
                call MakeT2p (wrk(PossT2w),wrk(PossTau),
     c                        aGrp,bGrp,aSGrp,bSGrp,1)
cxxxx           Make  Ww(+)((ab)",(bega)") from W1(a",be",b",ga")
c                                           and W2(b",be",a",ga")
c                     N.B. for a"=b" we have only W1 defined (=W2)
c                          but never mind (solved within MakeWw)
                call MakeWw (wrk(PossWw),wrk(PossW1),wrk(PossW2),
     c                       aSGrp,bSGrp,beSGrp,gaSGrp,1)
cxxxx           Add T2n1((bega)",ij) <<-
cxxxxc              Ww(+)(T)((ab)",(bega)").t2w(+)((ab)",ij)
                dim1=no*(no+1)/2
                if (aSGrp.eq.bSGrp) then
                  dim2=DimSGrpa(aSGrp)*(DimSGrpa(bSGrp)-1)/2
                else
                  dim2=DimSGrpa(aSGrp)*DimSGrpa(bSGrp)
                end if
                if (beSGrp.eq.gaSGrp) then
                  dim3=DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/2
                else
                  dim3=DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
                end if
                call mc0c1at3b (dim2,dim3,dim2,dim1,dim3,dim1,
     c                         dim3,dim2,dim1,
     c                         wrk(possWw),wrk(possT2w),wrk(possT2n1))
c
c
                 if (aSGrp.eq.bSGrp) then
cxxxx           Make t2w(+)((aa)",ij) from Tau((ab)',ij)
                call MakeT2pd (wrk(PossT2w),wrk(PossTau),aGrp,aSGrp)
cxxxx           Make  Ww(+)((aa)",(bega)") from W1(a",be",a",ga")
                call MakeWwd (wrk(PossWw),wrk(PossW1),
     c                       aSGrp,beSGrp,gaSGrp)
cxxxx           Add T2n1((bega)",ij) <<-
cxxxxc              Ww(+)(T)((aa)",(bega)").t2w(+)((aa)",ij)
                dim1=no*(no+1)/2
                dim2=DimSGrpa(aSGrp)
                if (beSGrp.eq.gaSGrp) then
                  dim3=DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/2
                else
                  dim3=DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
                end if
                call mc0c1at3b (dim2,dim3,dim2,dim1,dim3,dim1,
     c                         dim3,dim2,dim1,
     c                         wrk(possWw),wrk(possT2w),wrk(possT2n1))
                end if
c
c
cxxxx           Make t2w(-)((ab)",ij) from Tau((ab)',ij)
                call MakeT2m (wrk(PossT2w),wrk(PossTau),
     c                        aGrp,bGrp,aSGrp,bSGrp,1)
cxxxx           Make  Ww(-)((ab)",(bega)") from W1(a",be",b",ga")
c                                           and W2(b",be",a",ga")
c                          but never mind (solved within MakeWw)
                call MakeWw (wrk(PossWw),wrk(PossW1),wrk(PossW2),
     c                       aSGrp,bSGrp,beSGrp,gaSGrp,2)
c       write (6,*) 'V'
cxxxx           Add T2n2((bega)",ij) <<-
c                   Ww(-)(T)((ab)",(bega)").t2w(-)((ab)",ij)
                dim1=no*(no-1)/2
                if (aSGrp.eq.bSGrp) then
                  dim2=DimSGrpa(aSGrp)*(DimSGrpa(bSGrp)-1)/2
                else
                  dim2=DimSGrpa(aSGrp)*DimSGrpa(bSGrp)
                end if
                if (beSGrp.eq.gaSGrp) then
                  dim3=DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)-1)/2
                else
                  dim3=DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)
                end if
                call mc0c1at3b (dim2,dim3,dim2,dim1,dim3,dim1,
     c                         dim3,dim2,dim1,
     c                         wrk(possWw),wrk(possT2w),wrk(possT2n2))
c
c
cxxx          end cycle over all subgroups of (a>=b)'
              end do
              end do
c
cxxx          Save T2n1(bega)",ij) and T2n2((bega)",ij)
cxxxc              in recreated form, corresponding to standard
cxxxc              storage in T2 amplitudes
              LunName=Tmp3Name(beSGrp,gaSGrp)
              if (beSGrp.eq.gaSGrp) then
                lent2n1=no*(no+1)*
     c                  DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)+1)/4
                lent2n2=no*(no-1)*
     c                  DimSGrpbe(beSGrp)*(DimSGrpbe(gaSGrp)-1)/4
              else
                lent2n1=no*(no+1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)/2
                lent2n2=no*(no-1)*DimSGrpbe(beSGrp)*DimSGrpbe(gaSGrp)/2
              end if
              T2o2v4yes(beSGrp,gaSGrp)=1
              call SaveX (wrk(PossT2n1),lent2n1,LunAux,LunName,1,0)
              call SaveX (wrk(PossT2n2),lent2n2,LunAux,LunName,0,1)
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
        end
c
c        ---------------------
c
        subroutine Chck_W1 (W1,aSGrp,beSGrp,bSGrp,gaSGrp)
c
c        cek W1
        implicit none
#include "chcc1.fh"
        integer aSGrp,beSGrp,bSGrp,gaSGrp
        real*8 W1(1:16*31,1:16*33)
c
        integer a,b,be,ga,bad,ap,bp,bep,gap,ab,bega
        real*8 s
c
        if (aSGrp.eq.2) then
          ap=nv/2
        else
          ap=0
        end if
c
        if (bSGrp.eq.2) then
          bp=nv/2
        else
          bp=0
        end if
c
        if (gaSGrp.eq.2) then
          gap=nv/2
        else
          gap=0
        end if
c
        if (beSGrp.eq.2) then
          bep=nv/2
        else
          bep=0
        end if
c
        bad=0
        bega=0
        do be=1,nv/2
        do ga=1,be
        bega=bega+1
        ab=0
        do a=2,nv/2
        do b=1,a-1
        ab=ab+1
          s=(Q4(ap+a,bep+be,bp+b,gap+ga)
     c      +Q4(ap+a,gap+ga,bp+b,bep+be))/1
          if (abs(W1(ab,bega)-s).gt.1.0d-10) then
          bad=bad+1
c          write (6,99) a,b,be,ga,ab,bega,s,W1(a,be,ga)
c99        format(4(i2,1x),2(i6,1x),2(f15.10))
          end if
        W1(ab,bega)=s
        end do
        end do
        end do
        end do
c
        if (bad.eq.0) then
        write (6,*) ' Chck W OK ', bad
        else
        write (6,*) ' Chck W Bug !!!!!!! ', bad
        end if
c
        return
        end
c
c        ---------------------
c
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
c
c        ---------------------
c
        subroutine Chck_T21od (T21,beSGrp,gaSGrp)
c
c        test T2n+
c
        implicit none
#include "chcc1.fh"
        integer beSGrp,gaSGrp
        real*8 T21(1:32,1:32,1:no*(no-1)/2)
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
          do be=1,nv/2
          do ga=1,nv/2
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
          if (abs(T21(be,ga,uv)-s).gt.1.0d-10) then
            bad=bad+1
c        write (6,99) be,ga,u,v
c99        format (4(i3,1x))
          end if
          T21(be,ga,uv)=s
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
c
c        ---------------------
c
        subroutine Chck_T2p (T21,aSGrp,bSGrp)
c
c        test T2+
c
        implicit none
#include "chcc1.fh"
        real*8 T21(1:16*31,1:no*(no+1)/2)
        integer aSGrp,bSGrp
        integer u,v,a,b,ab,uv,bad,ap,bp
        real*8 s
c
c
        if (aSGrp.eq.2) then
          ap=nv/2
        else
          ap=0
        end if
c
        if (bSGrp.eq.2) then
          bp=nv/2
        else
          bp=0
        end if
c
        bad=0
c
        uv=0
        do u=1,no
        do v=1,u
        uv=uv+1
c
          ab=0
          do a=2,nv/2
          do b=1,a-1
          ab=ab+1
c
            s=(T2c(bp+b,ap+a,v,u)+T2c(bp+b,ap+a,u,v))/2
c
          if (abs(T21(ab,uv)-s).gt.1.0d-10) then
            bad=bad+1
          end if
          T21(ab,uv)=s
c
          end do
          end do
c
        end do
        end do
c
c
        if (bad.eq.0) then
        write (6,*) ' Chck T2+ OK ', bad
        else
        write (6,*) ' Chck T2+ Bug !!!!!!! ', bad
        end if
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_W2p (W2)
c
c        test W2+
c
        implicit none
#include "chcc1.fh"
        real*8 W2(1:nv,1:nv*(nv+1)/2)
        integer a,be,ga,bega,bad
        real*8 s
c
        bad=0
c
          bega=0
          do be=1,nv
          do ga=1,be
          bega=bega+1
c
            do a=1,nv
c
            s=Q4(a,ga,a,be)/2
            if (abs(W2(a,bega)-s).gt.1.0d-10) then
              bad=bad+1
            end if
            W2(a,bega)=s
c
            end do
c
          end do
          end do
c
        write (6,*) ' W2+ chck ',bad
c
        return
        end
c
c        ---------------------
c
        subroutine Calc_Bc
c
c        this routine calc Bc
c        Bc(a,b,be,ga) =        (a,be|b,ga)
c                     - S(i)   (a,be,|b,i).t1(ga,i)
c                     - S(i)   (a,i,|b,ga).t1(be,i)
c                     + S(i,j) (a,i|b,j).t1(be,i).t1(ga,j)
c
        implicit none
#include "chcc1.fh"
c
c        help var
        integer i,a,b,be,ga
c       integer j
        real*8 s
c
        do ga=1,nv
        do be=1,nv
        do b=1,nv
        do a=1,nv
c
c
c1            (a,be|b,ga)
          s=Q4(be,a,ga,b)
c
c2,3          - S(i)   (a,be,|b,i).t1(ga,i)
c         - S(i)   (a,i,|b,ga).t1(be,i)
          do i=1,no
           s=s-Q3(a,be,b,i)*T1c(ga,i)
           s=s-Q3(b,ga,a,i)*T1c(be,i)
          end do

c4        + S(i,j) (a,i|b,j).t1(be,i).t1(ga,j)
c          do j=1,no
c          do i=1,no
c          s=s+Q21(a,i,b,j)*T1c(be,i)*T1c(ga,j)
c          end do
c          end do
c
c
          Bc(a,b,be,ga)=s
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
c        ---------------------
c
        subroutine Chck_T1g (T1,dima,adda)
c
c        this routine test T1g
c
        implicit none
#include "chcc1.fh"
        integer dima,adda
         real*8 T1(1:no,1:dima)
c
c        help var
        integer a,bad,i,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do i=1,no
        do a=adda+1,adda+dima
c
           s=T1c(a,i)
c
           if (abs(T1(i,a-adda)-s).gt.1.0d-10) then
            bad=bad+1
            T1(i,a-adda)=s
          end if
          ntot=ntot+1

        end do
        end do
c
        write (6,*) ' T1g   ',bad,ntot
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_L1 (L1,dima,adda)
c
c        this routine test L1
c
        implicit none
#include "chcc1.fh"
        integer dima,adda
         real*8 L1(1:nc,1:dima,1:no)
c
c        help var
        integer m,a,bad,i,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do i=1,no
        do a=adda+1,adda+dima
        do m=1,nc
c
           s=L1k(m,i,a)
c
           if (abs(L1(m,a-adda,i)-s).gt.1.0d-10) then
            bad=bad+1
            L1(m,a-adda,i)=s
          end if
          ntot=ntot+1

        end do
        end do
        end do
c
        write (6,*) ' L1   ',bad,ntot
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_L2 (L2,dima,dimb,adda,addb)
c
c        this routine test L2
c
        implicit none
#include "chcc1.fh"
        integer dima,dimb,adda,addb
         real*8 L2(1:nc,1:dima,1:dimb)
c
c        help var
        integer m,a,b,bad,i,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do b=addb+1,addb+dimb
        do a=adda+1,adda+dima
        do m=1,nc
c
          s=L2k(m,a,b)
           do i=1,no
           s=s-L1k(m,i,a)*T1c(b,i)
           end do
c
           if (abs(L2(m,a-adda,b-addb)-s).gt.1.0d-10) then
            bad=bad+1
             L2(m,a-adda,b-addb)=s
          end if
          ntot=ntot+1

        end do
        end do
        end do
c
        write (6,*) ' L2   ',bad,ntot
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_B (BB,
     c             dima,dimb,dimbe,dimga,adda,addb,addbe,addga)
c
c        this routine test L2
c
        implicit none
#include "chcc1.fh"
        integer dima,dimb,dimbe,dimga,adda,addb,addbe,addga
         real*8 BB(1:dima,1:dimbe,1:dimb,1:dimga)
c
c        help var
        integer a,b,be,ga,bad,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do ga=addga+1,addga+dimga
        do b=addb+1,addb+dimb
        do be=addbe+1,addbe+dimbe
        do a=adda+1,adda+dima
c
          s=Bc(a,b,be,ga)
c
           if (abs(BB(a-adda,be-addbe,b-addb,ga-addga)-s)
     &              .gt.1.0d-10) then
            bad=bad+1
            BB(a-adda,be-addbe,b-addb,ga-addga)=s
          end if
          ntot=ntot+1

        end do
        end do
        end do
        end do
c
        write (6,*) ' B test ',bad,ntot
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_W4
     c  (W4,dima,dimbe,dimb,dimga,adda,addbe,addb,addga)
c    c  (W4,dima,dimga,dimb,dimbe,adda,addga,addb,addbe)
c
c        this routine test W4
c
        implicit none
#include "chcc1.fh"
        integer dima,dimbe,dimb,dimga,adda,addbe,addb,addga
         real*8 W4(1:dima,1:dimbe,1:dimb,1:dimga)
c        real*8 W4(1:dima,1:dimga,1:dimb,1:dimbe)
c
c        help var
        integer a,b,be,ga,i,bad,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do ga=1,dimga
        do b=1,dimb
        do be=1,dimbe
        do a=1,dima
          s=Q4(a+adda,be+addbe,b+addb,ga+addga)
          s=0.0d0
           do i=1,no
             s=s-Q3(a+adda,be+addbe,b+addb,i)*T1c(ga+addga,i)
             s=s-Q3(b+addb,ga+addga,a+adda,i)*T1c(be+addbe,i)
           end do
           if (abs(W4(a,be,b,ga)-s).gt.1.0d-10) then
c          if (abs(W4(a,ga,b,be)-s).gt.1.0d-10) then
            bad=bad+1
c            W4(a,be,b,ga)=s
          end if
          ntot=ntot+1
        end do
        end do
        end do
        end do
c
        write (6,*) ' W4 test ',bad,ntot
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_W3
     c  (W3,dima,dimbe,dimb,adda,addbe,addb)
c
c        this routine test W3 (a,be|b,i)
c
        implicit none
#include "chcc1.fh"
        integer dima,dimbe,dimb,adda,addbe,addb
         real*8 W3(1:dima,1:dimbe,1:dimb,1:no)
c
c        help var
        integer a,b,be,i,bad,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do i=1,no
        do b=1,dimb
        do be=1,dimbe
        do a=1,dima
          s=Q3(a+adda,be+addbe,b+addb,i)
           if (abs(W3(a,be,b,i)-s).gt.1.0d-10) then
            bad=bad+1
c            W3(a,be,b,i)=s
          end if
          ntot=ntot+1
        end do
        end do
        end do
        end do
c
        write (6,*) ' W3 test ',bad,ntot
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_W31
     c  (W3,dima,dimbe,dimb,adda,addbe,addb)
c
c        this routine test W3 (a,be,b,i) = (be,b|a,i)
c
        implicit none
#include "chcc1.fh"
        integer dima,dimbe,dimb,adda,addbe,addb
         real*8 W3(1:dima,1:dimbe,1:dimb,1:no)
c
c        help var
        integer a,b,be,i,bad,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do i=1,no
        do b=1,dimb
        do be=1,dimbe
        do a=1,dima
          s=Q3(b+addb,be+addbe,a+adda,i)
           if (abs(W3(a,be,b,i)-s).gt.1.0d-10) then
            bad=bad+1
c            W3(a,be,b,i)=s
          end if
          ntot=ntot+1
        end do
        end do
        end do
        end do
c
        write (6,*) ' W31 tst ',bad,ntot
c
        return
        end
c
c        ----------------
c
        subroutine CalcAddpp (aSGrp,addapp)
c
c        this routine calc addapp
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
        integer aSGrp,addapp
c
c        help var
        integer i
c
          addapp=0
        if (aSGrp.gt.1) then
          do i=1,aSGrp-1
            addapp=addapp+DimSGrpa(i)
          end do
        end if
c
        return
        end
c
c        ---------------------
c
        subroutine Chck_H1 (H1,dim,add)
c
c        this routine test H1(i,a") = t1o(a,i)
c
        implicit none
#include "chcc1.fh"
        integer dim,add
         real*8 H1(1:no,1:dim)
c
c        help var
        integer a,i,bad,ntot
        real*8 s
c
        bad=0
        ntot=0
c
        do a=1,dim
        do i=1,no
          s=T1c(a+add,i)
           if (abs(H1(i,a)-s).gt.1.0d-10) then
            bad=bad+1
c            H1(i,a)=s
          end if
          ntot=ntot+1
        end do
        end do
c
        write (6,*) ' H1 test ',bad,ntot
c
        return
        end
