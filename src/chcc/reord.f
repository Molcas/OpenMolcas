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
        subroutine Reord_chcc (wrk,wrksize,
     c             NaGrpR,NaSGrpR,NchBlk,LunAux)
c
c       This routine do:
c       1) Read local CD1 file of Choleski vectors from MC
c       2) Reorder ChV from pq,mloc to mloc,pq
c           Make L0-L2 files (local, dimensioned as ncLoc)
c           (if needed: (JoinLkey=1)
c       3) Make I0-I3 integrals
c          Calc E MP2
c           Make T20 (fist estimation (ai|jb)/Dijab
c           Make L0-L2 files (Global, dimensioned as nc)
c           (if needed: (JoinLkey=2)
c
c        for integral based approach (intkey=1)
c        5.1) W3 (vv|vo) integrals
c        5.2) W4 (vv|vv) integrals
c        5.3) Make L0-L2 files (Global, dimensioned as nc)
c           (if needed: JoinLkey=3)
c
c
c
c1      Structure of files, where selected group of (pq|rs) are
c       stored (V'O|OO) - I1 ; (V'O|V'O) - I2 ; (V'V'|OO) - I3
c
c       (IJ |KL)  I0intg
c
c       (A'I|JK)  I1inxx xx - Group of A'
c
c       (A'I|B'J) I2xxyy xx - Group of A'
c                        yy - Group of B'
c
c       (A'B'|IJ) I3xxyy xx - Group of A'
c                        yy - Group of B'
c
c
c2      Structure of Choleski vector files
c
c       L0(m,IJ)    L0vctr  I>=J
c
c       L1(m,I ,A') L1vcxx xx - Group of A'
c@@        kokot som, ze som to takto urobil, prerobit na L(m,a,i) to treba
c
c       L2(m,A'B')  L2xxyy xx - Group of A', A'>=B'
c                          yy - Group of B'
c
c3      Structure of Amplitude file
c       t2(A'B',I,J)  T2xxyy xx - Group of A'
c                            yy - Group of B'
c
c       Memory requirements:
c        choleski:
c       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm;  oom}
c       V2   - max {ov'ov'; v'v'm, ov'm; oom}
c       V3   - max {ov'm; oom}
c        V4   -      oom
c        integral based:
c       V1   - max {ov'ov'; nbas.nbas.m'; ov'm; v'v'm; oom; V"V"V"V"}
c       V2   - max {ov'ov'; v'v'm, ov'm; oom}
c       V3   - max {ov'm; oom; V'V'M}
c       V4   - oom
c       M1   - V"V"m
c       M2   - max {V"V"M; OV"M)
c
        implicit none
#include "chcc1.fh"
#include "chcc_reord.fh"
#include "o2v4.fh"
#include "chcc_files.fh"
#include "wrk.fh"
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "chcc_parcc.fh"
#endif
c
        integer NaGrpR,NaSGrpR,NChBlk
        integer LunAux

c       help variables
        integer dim1,dimij
        integer PossV1,PossV2,PossV3,PossV4,PossM1,PossM2
        integer PossT,maxdim
#ifdef _MOLCAS_MPP_
        integer j,NSGrp
#endif
        integer i,nbs,PossX
        integer ncLoc
c
        integer DimCh(1:100)
        integer LunChVF,NCh,ChLow,ChUp,idisk
c
        character*6 LunName
c
        character*24 Label
        Logical      Found
        integer      nOrbE
c
        integer dima,dimb,dimc,dimd,dimab,dimcd,dimci
        integer dimapp,dimbpp,dimcpp,dimdpp,dimabpp,dimcdpp
        integer adda,addb,addc,addd,addapp,addbpp,addcpp,adddpp
        integer aGrp,bGrp,cGrp,dGrp,aSGrp,bSGrp,cSGrp,dSGrp
        integer abGrp,cdGrp,abSGrp,cdSGrp
        integer bSGrpUp,dSGrpUp
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
        character*8 LunName8
        character*10 LunName10
        integer Nv4Ints
        integer w3aby,w4aby,w3abcy,w4abcdy
c
        real*8 e2,e2os
c
        integer isfreeunit
        integer idum(1)
c
c       Def parameters
        call DefParReord (NaGrpR,maxdim)
        NaGrp=NaGrpR
        NaSGrp=NaSGrpR
        NbeGrp=NaGrpR
        NbeSGrp=NaSGrpR
        call DefParo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
        if (printkey.ge.10) then
        write (6,*) ' Maxdim',maxdim,mdSGrpa
        end if
c
c        Distribute memory
        PossT=PossFree
        call DistMemReord (NaGrpR,maxdim,mdSGrpa,NchBlk,
     c       PossV1,PossV2,PossV3,PossV4,PossM1,PossM2,
     c       PossT)
        if (printkey.ge.10) then
        write (6,*) ' Last Value :',PossT,wrksize
        end if
        if (PossT.gt.wrksize) then
cmp!          write (6,*) ' Nieje dobre - Reord_chcc, Dr. Ch. Kokotopuloss',
        write (6,*) ' Not Enough memory in Reord_chcc step!',
     & 'Increase large and/or small segmentation ',
     c    (1.0d0*PossT)/(1.0d0*wrksize)
          call abend()
        end if
c
c*      Get Orbital energies
c
c       nOrbE=nfr+no+nv ! wrong size if ndel.ne.0
        Call Get_iArray('nBas',idum,1) ! must read always nBas fr runf
        nOrbE=idum(1)
        Label='OrbE'
        Call qpg_dArray(Label,Found,nOrbE)
        If(.not.Found .or. nOrbE.eq.0) Then
          Call SysAbendMsg('get_orbe','Did not find:',Label)
        End If
        call Get_dArray(Label,wrk(PossOE),nOrbE)
        if (printkey.ge.10) then ! toto som si nie isty
        do i=1,nfr+no+nv
        write (6,*) i,wrk(PossOE+i-1)
        end do
        end if
c
c*      skip frozen OE
        PossOE=PossOE+nfr
c
c*        Make Foo,Fvv,Fov
c        N.B. Ked bude OE(q) ine ako Fqq,  alebo F nediagonalny,
c            bude treba urobit inak
c
        dim1=no*no
        call mv0zero (dim1,dim1,wrk(PossFoo))
        dim1=nv*nv
        call mv0zero (dim1,dim1,wrk(PossFvv))
        dim1=nv*no
        call mv0zero (dim1,dim1,wrk(PossFvo))
c
c*        Escape, if Reord is not needed
c
        if (generkey.eq.0) then
          return
        end if
c
c*      ------- part 1, read  data form _CD file
c               and distributed in the form L(p',q',m') ------
c
c
c*      open _CD file and initialize parameters
c
cmp!<new 21/04/09
        LunChVF = 80
        LunChVF = isfreeunit(LunChVF)
cmp!>
        call DaName_mf_wa (LunChVF,'CD1tmp')
c
#ifdef _MOLCAS_MPP_
        ncLoc=NChLoc(myRank)
#else
        ncLoc=nc
#endif
c
        Nch=1
        ChLow=1
        ChUp=NChBlk
        DimCh(1)=ChUp-ChLow+1
        idisk=1
cmp
        if (printkey.ge.10) then
        write (6,*)
        write (6,*) 'ncLoc ',ncLoc
        write (6,*)
        end if
cmp
c
c*      read the block into V1(p,q,m') from _CD1
1       dim1=(no+nv)*(no+nv)*DimCh(NCh)
        if (printkey.ge.10) then
        write (6,*) 'Read _CD1',DimCh(NCh),dim1,idisk
        end if
        call ddafile (LunChVF,2,wrk(PossV1),dim1,idisk)
c
c
c1.1    Extract L0(ij,m')
c
        dimc=DimCh(NCh)
        nbs=no+nv
        dimij=no*(no+1)/2
        call Ext_L0 (wrk(PossV1),wrk(PossV2),
     c               no,dimij,dimc,nbs)
c
        LunName=L0Name
        dim1=dimij*dimc
        if (ChLow.eq.1) then
          call SaveX (wrk(PossV2),dim1,LunAux,LunName,1,1)
        else
          call SaveX (wrk(PossV2),dim1,LunAux,LunName,3,1)
        end if
c
c
c1.2    Extract L1(i,a',m')
c
        adda=0
        do aGrp=1,NaGrpR
        dima=DimGrpaR(aGrp)
c
          dimc=DimCh(NCh)
          nbs=no+nv
          call Ext_L1 (wrk(PossV1),wrk(PossV2),
     c                 no,dima,dimc,adda+no,nbs)
c
          LunName=L1Name(aGrp)
          dim1=dima*no*dimc
          if (ChLow.eq.1) then
            call SaveX (wrk(PossV2),dim1,LunAux,LunName,1,1)
          else
            call SaveX (wrk(PossV2),dim1,LunAux,LunName,3,1)
          end if
c
        adda=adda+dima
        end do
c
c
c1.3    Extract L2(a'b',m')
c
        adda=0
        do aGrp=1,NaGrpR
        dima=DimGrpaR(aGrp)
c
          addb=0
          do bGrp=1,aGrp
          dimb=DimGrpaR(bGrp)
c
            dimc=DimCh(NCh)
            nbs=no+nv
            if (aGrp.eq.bGrp) then
              dimab=dima*(dima+1)/2
              call Ext_L2s (wrk(PossV1),wrk(PossV2),
     c                      dima,dimab,dimc,adda+no,addb+no,nbs)
            else
              dimab=dima*dimb
              call Ext_L2u (wrk(PossV1),wrk(PossV2),
     c                      dima,dimb,dimc,adda+no,addb+no,nbs)
            end if
c
            LunName=L2Name(aGrp,bGrp)
            dim1=dimab*dimc
            if (ChLow.eq.1) then
              call SaveX (wrk(PossV2),dim1,LunAux,LunName,1,1)
            else
              call SaveX (wrk(PossV2),dim1,LunAux,LunName,3,1)
            end if
c
          addb=addb+dimb
          end do
        adda=adda+dima
        end do
c
c
c*      Upgrade parameters, if needed
        if (ChUp.lt.ncLoc) then
          Nch=Nch+1
          ChLow=ChLow+DimCh(Nch-1)
          if ((ChUp+NChBlk).gt.ncLoc) then
            ChUp=ncLoc
          else
            ChUp=ChUp+NChBlk
          end if
          DimCh(Nch)=ChUp-ChLow+1
          goto 1
        end if
c
c*      close _CD file
c
        call DaClos (LunChVF)
c
c
c
c*      ------- part 2, read data in the form L(p',q',m')
c               and redistributed in the form L(m,p'q')   --------
c
c
c*
c2.1    Reorder L0(ij,m') to L0(ml,ij)
c
          LunName=L0Name
          dimij=no*(no+1)/2
c
c2.1.1  get L0(ij,ml) into V1(ij,ml)
c
        if (NCh.eq.1) then
c       all in one
          dim1=dimij*ncLoc
          call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
        else if (NCh.eq.2) then
c       all in two
          PossX=PossV1
          dim1=dimij*DimCh(1)
          call GetX (wrk(PossX),dim1,LunAux,LunName,1,0)
          PossX=PossV1+dim1
          dim1=dimij*DimCh(2)
          call GetX (wrk(PossX),dim1,LunAux,LunName,0,1)
c
        else
c       more than 2 records
          PossX=PossV1
          dim1=dimij*DimCh(1)
          call GetX (wrk(PossX),dim1,LunAux,LunName,1,0)
          do i=2,NCh-1
            PossX=PossX+dim1
            dim1=dimij*DimCh(i)
            call GetX (wrk(PossX),dim1,LunAux,LunName,0,0)
          end do
          PossX=PossX+dim1
          dim1=dimij*DimCh(NCh)
          call GetX (wrk(PossX),dim1,LunAux,LunName,0,1)
c
        end if
c
c2.1.2  map V1(ij,ml) -> V2(ml,ij)
c
        call Map2_21 (wrk(PossV1),wrk(PossV2),dimij,ncLoc)
c
c2.1.3  write L0 = V2(ml,ij) back into file
c
        dim1=dimij*ncLoc
        call SaveX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
c
c*
c2.2    Reorder L1(i,a',m') to L1(ml,i,a')
c
        do aGrp=1,NaGrpR
        dima=DimGrpaR(aGrp)
        LunName=L1Name(aGrp)
c
c2.2.1  get L1(i,a',ml) into V1(i,a',ml)
c
        if (NCh.eq.1) then
c       all in one
          dim1=no*dima*ncLoc
          call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
        else if (NCh.eq.2) then
c       all in two
          PossX=PossV1
          dim1=no*dima*DimCh(1)
          call GetX (wrk(PossX),dim1,LunAux,LunName,1,0)
          PossX=PossV1+dim1
          dim1=no*dima*DimCh(2)
          call GetX (wrk(PossX),dim1,LunAux,LunName,0,1)
c
        else
c       more than 2 records
          PossX=PossV1
          dim1=no*dima*DimCh(1)
          call GetX (wrk(PossX),dim1,LunAux,LunName,1,0)
          do i=2,NCh-1
            PossX=PossX+dim1
            dim1=no*dima*DimCh(i)
            call GetX (wrk(PossX),dim1,LunAux,LunName,0,0)
          end do
          PossX=PossX+dim1
          dim1=no*dima*DimCh(NCh)
          call GetX (wrk(PossX),dim1,LunAux,LunName,0,1)
c
        end if
c
c2.2.2  map V1(i,a',ml) -> V2(ml,i,a')
c
        call Map2_21 (wrk(PossV1),wrk(PossV2),no*dima,ncLoc)
c
c2.2.3  write L1 = V2(ml,i,a') back into file
c
        dim1=no*dima*ncLoc
        call SaveX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
        end do
c
c
c*
c2.3    Reorder L2(a'b',m') to L2(ml,a'b')
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
c2.3.1  get L2(a'b',ml) into V1(a'b',ml)
c
        if (NCh.eq.1) then
c       all in one
          dim1=dimab*ncLoc
          call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
        else if (NCh.eq.2) then
c       all in two
          PossX=PossV1
          dim1=dimab*DimCh(1)
          call GetX (wrk(PossX),dim1,LunAux,LunName,1,0)
          PossX=PossV1+dim1
          dim1=dimab*DimCh(2)
          call GetX (wrk(PossX),dim1,LunAux,LunName,0,1)
c
        else
c       more than 2 records
          PossX=PossV1
          dim1=dimab*DimCh(1)
          call GetX (wrk(PossX),dim1,LunAux,LunName,1,0)
          do i=2,NCh-1
            PossX=PossX+dim1
            dim1=dimab*DimCh(i)
            call GetX (wrk(PossX),dim1,LunAux,LunName,0,0)
          end do
          PossX=PossX+dim1
          dim1=dimab*DimCh(NCh)
          call GetX (wrk(PossX),dim1,LunAux,LunName,0,1)
c
        end if
c
c2.3.2  map V1(a'b',ml) -> V2(ml,a'b')
c
        call Map2_21 (wrk(PossV1),wrk(PossV2),dimab,ncLoc)
c
c2.3.3  write L2 = V2(ml,a'b') back into file
c
        dim1=dimab*ncLoc
        call SaveX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
        end do
        end do
c
c
c2.4    reconstructing L0-L2 (Global, m=nc) ---
c           from Local L0-L2 (ml=NChLoc(myRank))
c        if needed (JoinLKey=1)
c
        if (JoinLkey.eq.1) then
          if (ncLoc.lt.nc) then
            call JoinLvec (wrk,wrksize,
     c                     PossV1,PossV2,NaGrpR,LunAux)
            ncLoc=nc
          end if
        end if
c
c
c
c*      ------- part 3, produce I integrals from Choleski vectors
c
c3.1    Generate I0(ij,kl)
c
c3.124.1 read V4(ml,ij) <- L0(ml,ij)
        LunName=L0name
        dimij=no*(no+1)/2
        call GetX (wrk(PossV4),ncLoc*dimij,LunAux,LunName,1,1)
c
c3.1.2        V1(ij,kl) = V4(T)(ml,ij) . V4(ml,kl)
        dim1=dimij*dimij
        call mv0zero(dim1,dim1,wrk(PossV1))
        call mc0c1at3b (ncLoc,dimij,ncLoc,dimij,dimij,dimij,
     c                  dimij,ncLoc,dimij,
     c                  wrk(PossV4),wrk(PossV4),wrk(PossV1))
c
#ifdef _MOLCAS_MPP_
c##        Synchronizacny bod
c3.1.3        Allreduce V1
        if (ncLoc.lt.nc) then
          dim1=dimij*dimij
          call gadgop (wrk(PossV1),dim1,'+')
        end if
#endif
c
c3.1.4  write V1(ij,kl)
        LunName=I0name
        dim1=dimij*dimij
        call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c
c3.23   generate I1(a'i,jk), I2(a'i,b'j) :-)
c
c
        e2=0.0d0
        e2os=0.0d0
c
        addb=0
        do bGrp=1,NaGrpR
        dimb=DimGrpaR(bGrp)
c
c3.23.2   read V2(ml,i,b') <- L1(ml,i,b')
          LunName=L1name(bGrp)
          dim1=ncLoc*no*dimb
          call GetX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
c3.23.3   Map V3(ml,b',i) <- V2(ml,i,b')
          call Map3_132 (wrk(PossV2),wrk(PossV3),ncLoc,no,dimb)
c
c3.2.4    calc V1(b',i,kl) <<- V3(T)(ml,b',i) . V4(ml,kl)
cBug      dim1=no*dima*dimij
          dim1=no*dimb*dimij
          call mv0zero(dim1,dim1,wrk(PossV1))
          call mc0c1at3b (ncLoc,dimb*no,ncLoc,dimij,dimb*no,dimij,
     c                    dimb*no,ncLoc,dimij,
     c                    wrk(PossV3),wrk(PossV4),wrk(PossV1))
c
#ifdef _MOLCAS_MPP_
c##          Synchronizacny bod
c3.2.5          Allreduce V1
          if (ncLoc.lt.nc) then
            dim1=no*dimb*dimij
            call gadgop (wrk(PossV1),dim1,'+')
          end if
#endif
c
c3.2.6    write V1(b',i,kl)
          LunName=I1name(bGrp)
          dim1=no*dimb*dimij
          call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c
          adda=0
          do aGrp=1,NaGrpR
          dima=DimGrpaR(aGrp)
          if (printkey.ge.10) then
          write (6,*) aGrp,bGrp,dima,dimb
          end if
c
c3.3.4      read V1(ml,i,a)
            LunName=L1name(aGrp)
            dim1=ncLoc*no*dima
            call GetX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c3.3.5      Map V2(ml,a',i) <- V1(ml,i,a')
            call Map3_132 (wrk(PossV1),wrk(PossV2),ncLoc,no,dima)
c
c3.3.6      V1(a',i,b',j) <<- V2(T)(ml,a',i) . V3(ml,b',j)
            dim1=no*no*dima*dimb
            call mv0zero(dim1,dim1,wrk(PossV1))
            call mc0c1at3b (ncLoc,dima*no,ncLoc,dimb*no,dima*no,dimb*no,
     c                      dima*no,ncLoc,dimb*no,
     c                      wrk(PossV2),wrk(PossV3),wrk(PossV1))
c
#ifdef _MOLCAS_MPP_
c##            Synchronizacny bod
c3.3.7            Allreduce V1
            if (ncLoc.lt.nc) then
              dim1=no*no*dima*dimb
              call gadgop (wrk(PossV1),dim1,'+')
            end if
#endif
c
c3.3.8      write V1(a',i,b',j)
            LunName=I2name(aGrp,bGrp)
            dim1=no*dima*no*dimb
            call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
c3.3.x      cierny vypocet E2
            call CVE2 (wrk(PossV1),wrk(PossOE),
     c                 dima,dimb,adda+no,addb+no,no,e2,e2os)
c
c
c3.3.x            Make V2(ab',i,j) <- V1(a',i,b',j)/Dijab
            if (aGrp.eq.bGrp) then
               call MkT20p (wrk(PossV2),wrk(PossV1),wrk(PossOE),
     c                     dima,adda+no,no)
            else
               call MkT20u (wrk(PossV2),wrk(PossV1),wrk(PossOE),
     c                     dima,dimb,adda+no,addb+no,no)
            end if
c
c3.3.x            Save V2 - T2(0) into T2Name
            LunName=T2Name(aGrp,bGrp)
            if (aGrp.eq.bGrp) then
              dim1=dima*(dima+1)*no*no/2
            else
              dim1=dima*dimb*no*no
            end if
            call SaveX (wrk(PossV2),dim1,LunAux,LunName,1,1)
c
          adda=adda+dima
          end do
c
        addb=addb+dimb
        end do
c
        write (6,*)
        write (6,92) ' E2 MP2    :',e2
        write (6,92) ' E2 ss     :',e2-e2os
        write (6,92) ' E2 os     :',e2os
        write (6,*)
92      format (a12,1x,f15.12)

c
c
c3.4    generate I3(a'b'|ij)
c
        adda=0
        do aGrp=1,NaGrpR
        dima=DimGrpaR(aGrp)
c
          addb=0
          do bGrp=1,aGrp
          dimb=DimGrpaR(bGrp)
c
c3.4.2      read V2(ml,a'b')
            LunName=L2name(aGrp,bGrp)
            if (aGrp.eq.bGrp) then
              dimab=dima*(dima+1)/2
            else
              dimab=dima*dimb
            end if
            call GetX (wrk(PossV2),ncLoc*dimab,LunAux,LunName,1,1)
c
c3.4.3      V1(a'b',ij) <<- V2(T)(ml,a'b') . V4(ml,ij)
            call mv0zero(dimab*dimij,dimab*dimij,wrk(PossV1))
            call mc0c1at3b (ncLoc,dimab,ncLoc,dimij,dimab,dimij,
     c                      dimab,ncLoc,dimij,
     c                      wrk(PossV2),wrk(PossV4),wrk(PossV1))
c
#ifdef _MOLCAS_MPP_
c##            Synchronizacny bod
c3.4.4            Allreduce V1
            if (ncLoc.lt.nc) then
              dim1=dimij*dimab
              call gadgop (wrk(PossV1),dim1,'+')
            end if
#endif
c
c3.4.5      write V1(a',i,b',j)
            LunName=I3name(aGrp,bGrp)
            dim1=dimij*dimab
            call SaveX (wrk(PossV1),dim1,LunAux,LunName,1,1)
c
          adda=adda+dima
          end do
c
        addb=addb+dimb
        end do
c
c
c3.5    reconstructing L0-L2 (Global, m=nc) ---
c           from Local L0-L2 (ml=NChLoc(myRank))
c        if needed (JoinLKey=2)
c
        if (JoinLkey.eq.2) then
          if (ncLoc.lt.nc) then
            call JoinLvec (wrk,wrksize,
     c                     PossV1,PossV2,NaGrpR,LunAux)
            ncLoc=nc
          end if
        end if
c
c
c        in W3 and W4 files are not generated, finish
        if (intkey.eq.0) then
          return
        end if
c
c
c*      ------- part 5, produce (vv|vo) and (vv|vv) integrals
c                        for integral based approach
c                        all terminology in o2v4 language
c
c        Extra cost: Mult. read - Na*(Na+1)/2 . L2(m,cd)
c                              - Na*(Na+1)/2 . L1(m,ci)
c
c        N.B.
c        - Ak je L2 uz spojene (ncLoc=nc) a W34DistKey=1
c          potom sa obskakuju ani nepisu tie kombinacie, ktore na
c          tomto node netreba (vhodne pre velke nProcs aj pri rychlej
c          sieti, pre pomalu siet aj pre mensie nProcs)
c        - Ak L2 este nieje spojene (ncLoc<nc) a W34DistKey=1
c          potom sa nic neobskakuje, vsetky W3 a W4 sa pocitaju
c          a allreducuju, ale sa nepisu na disk ak ich na tomto
c          node netreba, cim sa setri miesto
c          (vhodne pre mensie nProcs a rychlu siet)
c        - Ak L2 este nieje spojene (ncLoc<nc) a W34DistKey=0
c          potom sa nic neobskakuje, vsetky W3 a W4 sa pocitaju
c          a allreducuju, a pisu sa na disk vsade,
c          cim sa nesetri miesto, ale vsade su vsetky W3,4
c         (zbytocny pripad, pokial sa nechceme hrat s load-ballance
c           a nepotrebujeme mat vsade vsetko)
c        - Paralelny pripad (ncLoc=nc) a W34DistKey=0
c          bude pocitat iba kokot
c        - Osetrene na paralelny Nprocs=1, aj na skalarny case
c

#ifdef _MOLCAS_MPP_
c
        if (W34DistKey.eq.1) then
c*          case: Ditributed W34 files
c*.1          make a map, which W3 and W4 files are needed on this node
c          i.e. def InqW3, InqW4
          call Xo2v4ctl (NaGrp,NaSGrp,LunAux)
c
        else
c*          case: All W34 files on each node
c*.1          set InqW3,InqW4 - True
          NSGrp=NaGrp*NaSGrp
          do i=1,(NSGrp*(NSGrp+1))/2
            do j=1,NSGrp
            InqW3(i,j)=.True.
            end do
            do j=1,(NSGrp*(NSGrp+1))/2
            InqW4(i,j)=.True.
            end do
          end do
c
        end if
c
#endif
c
        Nv4Ints=0
c
c        cycle over a'>=b' groups
        adda=0
        do aGrp=1,NaGrp
        dima=DimGrpa(aGrp)
        addb=0
        do bGrp=1,aGrp
        dimb=DimGrpa(bGrp)
        abGrp=aGrp*(aGrp-1)/2+bGrp
c
        if (printkey.ge.10) then
        write (6,*) ' W3 + W4 ',aGrp,bGrp
        end if
c
c        test, if atleast one file of W3/W4 integrals
c        need to be calculated on this node for given a',b'
c        skip, if there is nothing to calculate for this a'b'
c        on this node and L is joined (ncLoc=nc)
        call DefW34y (aGrp,bGrp,w3aby,w4aby,NaGrp)
        if ((w3aby+w4aby).eq.0) then
          if (ncLoc.eq.nc) goto 41
        end if
c
c         read V2(ml,a'b') = L2(ml,a'b')
          LunName=L2name(aGrp,bGrp)
          if (aGrp.eq.bGrp) then
            dimab=dima*(dima+1)/2
          else
            dimab=dima*dimb
             end if
          call GetX (wrk(PossV2),ncLoc*dimab,LunAux,LunName,1,1)
c
c
c*          ---- Part 5.1 - Generation of (ab|ci) integrals ----
c               (Svincaja morda)
c
c
c          skip, if there is no W3 file to calculate for this a'b'
c          on this node and L is joined (ncLoc=nc)
          if ((w3aby).eq.0) then
            if (ncLoc.eq.nc) goto 42
          end if
c
c          cycle over c' groups
          addc=0
          do cGrp=1,NaGrp
          dimc=DimGrpa(cGrp)
c
c            test, if atleast one file of W3 integrals
c            need to be calculated on this node for given a',b',c'
c            and skip, if there is nothing to calculate for this a',b',c'
c              on this node and L is joined (ncLoc=nc)
            call DefW3y (aGrp,bGrp,cGrp,w3abcy)
            if ((w3abcy).eq.0) then
              if (ncLoc.eq.nc) goto 43
            end if


c            Get V3(ml,c',i)
c            Read V1(ml,i,c') <- L1(ml,i,c)
            dimci=dimc*no
            LunName=L1Name(cGrp)
            call GetX (wrk(PossV1),ncLoc*dimci,LunAux,LunName,1,1)
c            Map V3(ml,c',i) <- V1(ml,i,c')
            call Map3_132 (wrk(PossV1),wrk(PossV3),ncLoc,no,dimc)
c
c
c            cycle over a">=b" subgroups
            addapp=0
            do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
            dimapp=DimSGrpa(aSGrp)
            if (aGrp.eq.bGrp) then
              bSGrpUp=aSGrp
            else
              bSGrpUp=GrpaUp(bGrp)
            end if
            addbpp=0
            do bSGrp=GrpaLow(bGrp),bSGrpUp
            dimbpp=DimSGrpa(bSGrp)
            abSGrp=aSGrp*(aSGrp-1)/2+bSGrp
c
c
c              Extract M1(ml,a"b") <- V2(ml,a'b')
              if (aSGrp.eq.bSGrp) then
                dimabpp=dimapp*(dimapp+1)/2
              else
                dimabpp=dimapp*dimbpp
              end if
              call Ext_W4 (wrk(PossV2),wrk(PossM1),
     c                     ncLoc,dima,dimb,dimab,dimapp,dimbpp,dimabpp,
     c                     addapp,addbpp,aGrp,bGrp,aSGrp,bSGrp)
c
c
c              cycle over c" subgroups
              addcpp=0
              do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
              dimcpp=DimSGrpa(cSGrp)
c
#ifdef _MOLCAS_MPP_
c                skip, if W3(a"b",c"i) is not needed on this node
c               and L is joined (ncLoc=nc)
                if (InqW3(abSGrp,cSGrp).eqv..False.) then
                  if (ncLoc.eq.nc) goto 44
                end if
#endif
c
c                Extract M2(ml,c",i) <- V3(ml,c',i)
                call Ext_W3 (wrk(PossV3),wrk(PossM2),
     c                       ncLoc,no,dimc,dimcpp,addcpp)
c
c                Calc V1(a"b",c"i) <- M1(T)(ml,a"b") . M2(ml,c",i)
                dim1=dimcpp*no
                call mv0zero (dimabpp*dim1,dimabpp*dim1,wrk(PossV1))
                call mc0c1at3b (ncLoc,dimabpp,ncLoc,dim1,dimabpp,dim1,
     c                          dimabpp,ncLoc,dim1,
     c                          wrk(PossM1),wrk(PossM2),wrk(PossV1))
c
#ifdef _MOLCAS_MPP_
c##                Synchronizacny bod
c                     Allreduce V1
                if (ncLoc.lt.nc) then
                  dim1=dimabpp*dimcpp*no
                  call gadgop (wrk(PossV1),dim1,'+')
                end if
c                skip, if W3(a"b",c"i) is not needed on this node
                if (InqW3(abSGrp,cSGrp).eqv..False.) then
                  goto 44
                end if
#endif
c
c                Create Proper LunName8
                call MkNameV3 (aSGrp,bSGrp,cSGrp,'W3',LunName8)
c
c                Write integral block to proper file
*                open (unit=LunAux,file=LunName8,form='unformatted')
                 Call MOLCAS_BinaryOpen_Vanilla(LunAux,LunName8)
                call wri_chcc (LunAux,dimabpp*dimcpp*no,wrk(PossV1))
                close (LunAux)
c
c
c              end cycle over c" subgroups
#ifdef _MOLCAS_MPP_
44              continue
#endif
                addcpp=addcpp+dimcpp
              end do
c
c
c            end cycle over a">=b" subgroups
            addbpp=addbpp+dimbpp
            end do
            addapp=addapp+dimapp
            end do
c
c
c          end cycle over c' groups
43          addc=addc+dimc
          end do
c
c
c*          ---- Part 5.2 - Generation of (ab|cd) integrals ----
c               (Kobyljacaja sraka)
c
c
c          skip, if there is no W4 file to calculate for this a'b'
c          on this node and L is joined (ncLoc=nc)
42        if ((w4aby).eq.0) then
            if (ncLoc.eq.nc) goto 41
          end if
c
c          cycle over c'>=d' groups
          addc=0
          do cGrp=1,NaGrp
          dimc=DimGrpa(cGrp)
          addd=0
          do dGrp=1,cGrp
          dimd=DimGrpa(dGrp)
          cdGrp=cGrp*(cGrp-1)/2+dGrp
             if (cGrp.gt.aGrp) goto 49
c
c            test, if atleast one file of W4 integrals
c            need to be calculated on this node for given a',b',c',d'
c            and skip, if there is nothing to calculate for this
c                  a',b',c',d' on this node and L is joined (ncLoc=nc)
            call DefW4y (aGrp,bGrp,cGrp,dGrp,w4abcdy)
            if ((w4abcdy).eq.0) then
              if (ncLoc.eq.nc) goto 49
            end if
c
c           read V3(ml,c'd') = L2(ml,c'd')
            if (abGrp.eq.cdGrp) then
              dimcd=dimab
              call mv0u (ncLoc*dimcd,wrk(PossV2),1,wrk(PossV3),1)
            else
              LunName=L2name(cGrp,dGrp)
              if (cGrp.eq.dGrp) then
                dimcd=dimc*(dimc+1)/2
              else
                dimcd=dimc*dimd
                 end if
              call GetX (wrk(PossV3),ncLoc*dimcd,LunAux,LunName,1,1)
            end if

c
c            cycle over a">=b" subgroups
            addapp=0
            do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
            dimapp=DimSGrpa(aSGrp)
            if (aGrp.eq.bGrp) then
              bSGrpUp=aSGrp
            else
              bSGrpUp=GrpaUp(bGrp)
            end if
            addbpp=0
            do bSGrp=GrpaLow(bGrp),bSGrpUp
            dimbpp=DimSGrpa(bSGrp)
            abSGrp=aSGrp*(aSGrp-1)/2+bSGrp
c
c              test, if atleast one file of W4 integrals
c              need to be calculated on this node for given a",b",c',d'
c              and skip, if there is nothing to calculate for this
c                    a",b",c',d' on this node and L is joined (ncLoc=nc)
              call DefW4y2 (aSGrp,bSGrp,cGrp,dGrp,w4abcdy)
              if ((w4abcdy).eq.0) then
                if (ncLoc.eq.nc) goto 45
              end if
c
c              Extract M1(ml,a"b") <- V2(ml,a'b')
              if (aSGrp.eq.bSGrp) then
                dimabpp=dimapp*(dimapp+1)/2
              else
                dimabpp=dimapp*dimbpp
              end if
              call Ext_W4 (wrk(PossV2),wrk(PossM1),
     c                     ncLoc,dima,dimb,dimab,dimapp,dimbpp,dimabpp,
     c                     addapp,addbpp,aGrp,bGrp,aSGrp,bSGrp)
c
c
c              cycle over c">=d" subgroups
              addcpp=0
              do cSGrp=GrpaLow(cGrp),GrpaUp(cGrp)
              dimcpp=DimSGrpa(cSGrp)
              if (cGrp.eq.dGrp) then
                dSGrpUp=cSGrp
              else
                dSGrpUp=GrpaUp(dGrp)
              end if
              adddpp=0
              do dSGrp=GrpaLow(dGrp),dSGrpUp
              dimdpp=DimSGrpa(dSGrp)
              cdSGrp=cSGrp*(cSGrp-1)/2+dSGrp
              if (cdSGrp.gt.abSGrp) goto 48
c
#ifdef _MOLCAS_MPP_
c                skip, if W4(a"b",c"d") is not needed on this node
c                and L is joined (ncLoc=nc)
                if (InqW4(abSGrp,cdSGrp).eqv..False.)  then
                  if (ncLoc.eq.nc) goto 48
                end if
#endif
c
c                Extract M2(ml,a"b") <- V3(ml,a'b')
                if (cSGrp.eq.dSGrp) then
                  dimcdpp=dimcpp*(dimcpp+1)/2
                else
                  dimcdpp=dimcpp*dimdpp
                end if
                call Ext_W4 (wrk(PossV3),wrk(PossM2),
     c                      ncLoc,dimc,dimd,dimcd,dimcpp,dimdpp,dimcdpp,
     c                       addcpp,adddpp,cGrp,dGrp,cSGrp,dSGrp)
c
c                Calc V1(a"b",c"d") <- M1(T)(ml,a"b") . M2(ml,c"d")
                dim1=dimabpp*dimcdpp
                call mv0zero (dim1,dim1,wrk(PossV1))
                call mc0c1at3b (ncLoc,dimabpp,ncLoc,dimcdpp,
     c                          dimabpp,dimcdpp,dimabpp,ncLoc,dimcdpp,
     c                          wrk(PossM1),wrk(PossM2),wrk(PossV1))
c
#ifdef _MOLCAS_MPP_
c##                Synchronizacny bod
c                     Allreduce V1
                if (ncLoc.lt.nc) then
                  dim1=dimabpp*dimcdpp
                  call gadgop (wrk(PossV1),dim1,'+')
                end if
c                skip, if W4(a"b",c"d") not needed on this node
                if (InqW4(abSGrp,cdSGrp).eqv..False.)  then
                  goto 48
                end if
#endif
c
c               Def proper LunName10
                call MkNameV4 (aSGrp,bSGrp,cSGrp,dSGrp,'W4',LunName10)
c
c                Write integral block to proper file
*                open (unit=LunAux,file=LunName10,form='unformatted')
                 Call MOLCAS_BinaryOpen_Vanilla(LunAux,LunName10)
                call wri_chcc (LunAux,dimabpp*dimcdpp,wrk(PossV1))
                close (LunAux)
                Nv4Ints=Nv4Ints+dimabpp*dimcdpp
c
c
c              end cycle over c">=d" subgroups
48              continue
              adddpp=adddpp+dimdpp
              end do
              addcpp=addcpp+dimcpp
              end do
c
c            end cycle over a">=b" subgroups
45          addbpp=addbpp+dimbpp
            end do
            addapp=addapp+dimapp
            end do

c          end cycle over c'>=d' groups
49          continue
          addd=addd+dimd
          end do
          addc=addc+dimc
          end do
c
c        end cycle over a'>=b' groups
41        addb=addb+dimb
        end do
        adda=adda+dima
        end do
c
        write (6,*) ' V4 ints compressing factor: ',
     c              1.0d0*(nv*nv*nv*nv)/Nv4Ints

c
c5.3    reconstructing L0-L2 (Global, m=nc) ---
c        if needed (JoinLKey=3)
c
        if (JoinLkey.eq.3) then
          if (ncLoc.lt.nc) then
            call JoinLvec (wrk,wrksize,
     c                     PossV1,PossV2,NaGrpR,LunAux)
            ncLoc=nc
          end if
        end if
c
c
        return
        end
