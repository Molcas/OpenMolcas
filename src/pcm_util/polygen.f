************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Christian Pomelli                                      *
************************************************************************
      Subroutine PolyGen(MaxT,Ptype,Pflag,TsAre,TsNum,
     +          XEN,YEN,ZEN,REN,TsEff,CV,JTR)
      implicit real*8 (d,g)
c
c Polygen: a program to generate spherical polyhedra with triangular
c faces by C.S. Pomelli (cris@ibm550.icqem.cnr.it)
c An equilater division algorithm is used.
      parameter (Icosa = 1)
      parameter (NVI=12, NTI=20, NEI=30)
      parameter (NVP=32, NTP=60, NEP=90)
      parameter (NVT= 4, NTT=4 , NET= 6)
      integer i,j,k,l,m,n,ii,jj,NV,NT0,NE0,NF
      integer NTPT,NVPT,Ptype,TsNum,TsEff,Pflag
      integer NFI,NFP,NFT,NDI,NDP,NDT
      integer trvo(60,3),treo(60,3),edo(90,2)
      integer JTR(3,*)
*                                                                      *
************************************************************************
*                                                                      *
*
*     Change trnew to be allocated dynamically.
*
C     integer oldtr(100,100),ednew(90,100),trnew(60,100,100) ! old code
*
#include "WrkSpc.fh"
      integer oldtr(100,100),ednew(90,100)
*                                                                      *
************************************************************************
*                                                                      *
      integer tettrv(3,4), tettre(3,4), teted(2,6)
      integer icotrv(3,20),icotre(3,20),icoed(2,30)
      integer pentrv(3,60),pentre(3,60),pened(2,90)
      Real*8 tetve(3,4),icove(3,12),penve(3,32),dnorm,TsAre,pi,v1(3),
     $  v2(3),v3(3),alpha,beta,costheta,theta,cos1,cos2,sintheta,XEN,
     $  YEN,ZEN,REN,CV(3,*)
c     data pi/3.141592653589793d0/
c Data for original polyhedra
c 1) Tetrahedron
      data teted/1,2,2,3,1,3,1,4,3,4,2,4/
      data tettre/1,2,3,1,6,4,3,5,4,2,5,6/
      data tettrv/1,2,3,1,2,4,1,3,4,2,3,4/
      data tetve/
     &  -.577350269D+00, .577350269D+00, .577350269D+00,
     &  .577350269D+00,-.577350269D+00, .577350269D+00,
     &  -.577350269D+00,-.577350269D+00,-.577350269D+00,
     &  .577350269D+00, .577350269D+00,-.577350269D+00/
c 2) Icosaedron
      data icoed/
     &  1, 2, 2, 3, 1, 3, 2, 6, 1, 6, 3, 4, 1, 4, 4, 5, 1, 5, 5, 6, 3,
     &  7, 2, 7, 6,11, 2,11, 7,11, 4, 8, 3, 8, 7, 8, 5, 9, 4, 9, 8, 9,
     &  6,10, 5,10, 9,10,10,11, 8,12, 7,12,11,12, 9,12,10,12/
      data icotre/
     &  1, 2, 3, 1, 4, 5, 3, 6, 7, 7, 8, 9, 9,10, 5, 2,11,12,
     &  4,13,14,12,15,14, 6,16,17,11,18,17, 8,19,20,16,21,20,10,22,23,
     &  19,24,23,22,25,13,18,26,27,15,28,27,21,29,26,24,30,29,25,28,30/
      data icotrv/
     &  1 , 2, 3, 1, 2, 6, 1, 3, 4, 1, 4, 5, 1, 5, 6, 2, 3, 7,
     &  2, 6,11, 2, 7,11, 3, 4, 8, 3, 7, 8, 4, 5, 9, 4, 8, 9, 5, 6,10,
     &  5, 9,10, 6,10,11, 7, 8,12, 7,11,12, 8, 9,12, 9,10,12,10,11,12/
      data icove/
     &  .000000000D+00,.000000000D+00, .100000000D+01,.276393202D+00,
     &  .850650808D+00, .447213595D+00,-.723606798D+00, .525731112D+00,
     &  .447213595D+00,-.723606798D+00,-.525731112D+00, .447213595D+00,
     & .276393202D+00,-.850650808D+00, .447213595D+00,.894427191D+00,-
     & .219071479D-15, .447213595D+00,-.276393202D+00, .850650808D+00,-
     & .447213595D+00,-.894427191D+00, .109535740D-15,-.447213595D+00,-
     & .276393202D+00,-.850650808D+00,-.447213595D+00,.723606798D+00,-
     & .525731112D+00,-.447213595D+00,.723606798D+00, .525731112D+00,-
     & .447213595D+00,.000000000D+00, .000000000D+00,-.100000000D+01/
c 3) Pentakisdodecahedron
      data pened/
     &  1 , 2, 2, 3, 1, 3, 2, 6, 1, 6, 3, 4, 1, 4, 4, 5, 1, 5, 5, 6, 3,
     &  8,2, 8, 6, 7, 2, 7, 7,12, 2,12, 8,12, 4, 9, 3, 9, 8,13, 3,13, 9,
     &  13, 5,10, 4,10, 9,14, 4,14,10,14, 6,11, 5,11,10,15, 5,15,11,15,
     &  7,16, 6,16,11,16,12,17, 7,17,16,17,12,18, 8,18,13,18,13,19, 9,
     &  19,14,19,14,20,10,20,15,20,15,21,11,21,16,21,17,22,12,22,18,22,
     &  18,23,13,23,19,23,19,24,14,24,20,24,20,25,15,25,21,25,17,26,16,
     &  26,21,26,22,27,17,27,26,27,22,28,18,28,23,28,23,29,19,29,24,29,
     &  24,30,20,30,25,30,25,31,21,31,26,31,27,28,28,29,29,30,30,31,27,
     &  31,28,32,27,32,31,32,29,32,30,32/
      data pentre/
     &  1, 2, 3, 1, 4, 5, 3, 6, 7, 7, 8, 9, 9,10, 5, 2,11,12, 4,13,14,
     &  14,15,16,12,17,16, 6,18,19,11,20,21,19,22,21, 8,23,24,18,25,26,
     &  24,27,26,10,28,29,23,30,31,29,32,31,13,33,34,28,35,34,15,36,37,
     &  33,38,37,17,39,40,20,41,40,22,42,43,25,44,43,27,45,46,30,47,46,
     &  32,48,49,35,50,49,36,51,52,39,53,52,41,54,55,42,56,55,44,57,58,
     &  45,59,58,47,60,61,48,62,61,38,63,64,50,65,64,51,66,67,63,68,67,
     &  53,69,70,54,71,70,56,72,73,57,74,73,59,75,76,60,77,76,62,78,79,
     &  65,80,79,66,81,69,71,82,72,74,83,75,77,84,78,68,85,80,81,86,87,
     &  85,88,87,82,89,86,83,90,89,84,88,90/
      data pentrv/
     &  1, 2, 3, 1, 2, 6, 1, 3, 4, 1, 4, 5, 1, 5, 6, 2, 3, 8, 2, 6, 7,
     &  2, 7,12, 2, 8,12, 3, 4, 9, 3, 8,13, 3, 9,13, 4, 5,10, 4, 9,14,
     &  4,10,14, 5, 6,11, 5,10,15, 5,11,15, 6, 7,16, 6,11,16, 7,12,17,
     &  7,16,17, 8,12,18, 8,13,18, 9,13,19, 9,14,19,10,14,20,10,15,20,
     &  11,15,21,11,16,21,12,17,22,12,18,22,13,18,23,13,19,23,14,19,24,
     &  14,20,24,15,20,25,15,21,25,16,17,26,16,21,26,17,22,27,17,26,27,
     &  18,22,28,18,23,28,19,23,29,19,24,29,20,24,30,20,25,30,21,25,31,
     &  21,26,31,22,27,28,23,28,29,24,29,30,25,30,31,26,27,31,27,28,32,
     &  27,31,32,28,29,32,29,30,32,30,31,32/
      data penve/
     &  .000000000D+00, .000000000D+00, .100000000D+01, .491123473D+00,
     &  .356822090D+00, .794654472D+00,-.187592474D+00, .577350269D+00,
     &  .794654472D+00,-.607061998D+00, .540159669D-09, .794654472D+00,-
     &  .187592475D+00,-.577350269D+00, .794654472D+00, .491123473D+00,-
     &  .356822091D+00, .794654472D+00, .894427191D+00, .000000000D+00,
     &  .447213595D+00, .276393203D+00, .850650808D+00, .447213595D+00,-
     &  .723606797D+00, .525731113D+00, .447213595D+00,-.723606799D+00,-
     &  .525731111D+00, .447213595D+00, .276393201D+00,-.850650809D+00,
     &  .447213595D+00, .794654472D+00, .577350269D+00, .187592474D+00,-
     &  .303530999D+00, .934172359D+00, .187592474D+00,-.982246946D+00,
     &  .873996703D-09, .187592474D+00,-.303531000D+00,-.934172359D+00,
     &  .187592474D+00, .794654471D+00,-.577350271D+00, .187592474D+00,
     &  .982246946D+00, .000000000D+00,-.187592474D+00, .303530999D+00,
     &  .934172359D+00,-.187592474D+00,-.794654472D+00, .577350270D+00,-
     &  .187592474D+00,-.794654473D+00,-.577350268D+00,-.187592474D+00,
     &  .303530997D+00,-.934172359D+00,-.187592474D+00, .723606798D+00,
     &  .525731112D+00,-.447213596D+00,-.276393202D+00, .850650808D+00,-
     &  .447213596D+00,-.894427191D+00, .795855278D-09,-.447213596D+00,-
     &  .276393203D+00,-.850650808D+00,-.447213596D+00, .723606797D+00,-
     &  .525731113D+00,-.447213596D+00, .607061998D+00, .000000000D+00,-
     &  .794654472D+00, .187592474D+00, .577350269D+00,-.794654472D+00,-
     &  .491123473D+00, .356822090D+00,-.794654472D+00,-.491123474D+00,-
     &  .356822089D+00,-.794654472D+00, .187592473D+00,-.577350269D+00,-
     &  .794654472D+00,.000000000D+00, .000000000D+00,-.100000000D+01/
c
c set pi
c
      PI = DBLE(4)*ATan(DBLE(1))
      NT0= 0 ! dummy initialize
      NV = 0 ! dummy initialize
      NE0= 0 ! dummy initialize
c
c  A seconda delle due opzioni di funzionamento (Area ottimale o
c  numero di tessere ottimale ) vengono stabilite il tipo di poliedro
c  (Ptype) e la  frequenza di divisione (NF)
c
c  Pflag = 0 Numero di tessere: si genera il poliedro disponibile
c  con il numero di tessere piu' prossimo a TsNum
c
c  Pflag  =1  Area: si genera il poliedro con il numero di tessere tali
c  da avere un area il piu' simile possibile
      if (Pflag.eq.1) then
        TsNum = INT (4.0d0 * pi * REN**2 / TsAre + 0.5 )
      endif
C
c  -si controlla che TsNum <= MxTs
C
      if (TsNum.gt.MaxT) then
        write(6,'(''Too many tesserae in PolyGen'')')
        Call Abend()
      endif
C
c  -calcola il valore di NF ottimale per ciascuna famiglia di poliedri.
C
      NFI = INT( SQRT(TsNum / ( NTI * 1.0d0 )) + 0.5d0)
      NFP = INT( SQRT(TsNum / ( NTP * 1.0d0 )) + 0.5d0)
      NFT = INT( SQRT(TsNum / ( NTT * 1.0d0 )) + 0.5d0)
C
c  -trova lo scostamento del numero di tessere in input
C
      NDI = IABS( TsNum - NTI * NFI**2 )
      NDP = IABS( TsNum - NTP * NFP**2 )
      NDT = IABS( TsNum - NTT * NFT**2 )
C
c  -sceglie tipo e frequenza di divisione
c  a parita' di differenza in numero di tessere sceglie quello con NF
c  minore
C
      if     ((NDP.le.NDI).and.(NDP.le.NDT)) then
        Ptype = 2
        NF = NFP
        TsEff = NTP * NFP**2
      else
        if (NDI.le.NDT) then
          Ptype = Icosa
          NF = NFI
          TsEff = NTI * NFI**2
        else
          Ptype = 3
          NF = NFT
          TsEff = NTT * NFT**2
        endif
      endif
C
c  -carica gli opportuni valori iniziali
C
      if (Ptype.eq.1) then
C
c  icosaedro
C
        NT0=NTI
        NV=NVI
        NE0=NEI
        do 2000 i=1, NT0
          do 2010 k=1,3
            trvo(i,k)=icotrv(k,i)
            treo(i,k)=icotre(k,i)
 2010 continue
 2000 continue
        do 2020 i=1,NV
          do 2030 k=1,3
            CV(k,i)=icove(k,i)
 2030 continue
 2020 continue
        do 2040 i=1,NE0
          do 2050 k=1,2
            edo(i,k)=icoed(k,i)
 2050 continue
 2040 continue
      elseif (Ptype.eq.2) then
C
c  -pentakisdodecaedro
C
        NT0=NTP
        NV=NVP
        NE0=NEP
        do 2060 i=1, NT0
          do 2070 k=1,3
            trvo(i,k)=pentrv(k,i)
            treo(i,k)=pentre(k,i)
 2070 continue
 2060 continue
        do 2080 i=1,NV
          do 2090 k=1,3
            CV(k,i)=penve(k,i)
 2090 continue
 2080 continue
        do 2100 i=1,NE0
          do 2110 k=1,2
            edo(i,k)=pened(k,i)
 2110 continue
 2100 continue
      elseif (Ptype.eq.3) then
C
c  tetraedro
C
        NT0=NTT
        NV=NVT
        NE0=NET
        do 2120 i=1, NT0
          do 2130 k=1,3
            trvo(i,k)=tettrv(k,i)
            treo(i,k)=tettre(k,i)
 2130 continue
 2120 continue
        do 2140 i=1,NV
          do 2150 k=1,3
            CV(k,i)=tetve(k,i)
 2150 continue
 2140 continue
        do 2160 i=1,NE0
          do 2170 k=1,2
            edo(i,k)=teted(k,i)
 2170 continue
 2160 continue
      endif
      NVPT=NV+1
      If(NVPT.gt.1000) then
        Write(6,*)'NVPT out of range in polygen',nvpt
        Call Abend()
      EndIf
C
c  -nuovi vertici posti lungo i vecchi spigoli
c  -le regole di calcolo derivano da  calcoli di algebra
c   lineare e trigonometria
c  -i vertici vengono memorizzati in ordine progressivo ed riferiti
c  tramite ednew al vertice di appartenenza
C
      do j=1,NE0
        do k=1,3
          v1(k)=CV(k,edo(j,1))
          v2(k)=CV(k,edo(j,2))
        end do
        costheta=
     &    ( v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)) /
     &    (sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2) *
     &     sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2))
        theta=acos(costheta)
        sintheta=sin(theta)
        do l=1,NF-1
          m=NF-l
          cos1=cos(theta*l/NF)
          cos2=cos(theta*(NF-l)/NF)
          alpha=(cos1-costheta*cos2)/sintheta**2
          beta=(cos2-costheta*cos1)/sintheta**2
          do k=1,3
            v3(k)=alpha*v1(k)+beta*v2(k)
          end do
          dnorm = sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)
          do k=1,3
            v3(k) = v3(k) /dnorm
          end do
          do k=1,3
            CV(k,NVPT)=v3(k)
          end do
          ednew(j,l+1)=NVPT
          NVPT=NVPT+1
          If(NVPT.gt.1000) then
            Write(6,*)'NVPT out of range in polygen',nvpt
            Call Abend()
          EndIf
      end do
      end do
C
c  -nuovi vertici non posti lungo i vecchi spigoli
c  -a partire dai vertici in ednew secondo regole analoghe alle
c  precedenti, vengono memorizzati a seconda del triangolo in trnew
c  trnew(triangolo,fila,n ordine)=nvertice
C
*
*     Allocate TrNew. Note that we change the order from NT0*NF*NF
*     to NF*NF*NT0, (ixx,jxx,kxx) to (kxx,jxx,ixx)
*
      nTrNew=NT0*NF*NF
      Call GetMem('TrNew','Allo','Real',ipTrNew,nTrNew)
*
      do 2240 j=1,NT0
        ii=treo(j,1)
        jj=treo(j,3)
        do 2250 l=3,NF
          do 2260 k=1,3
            v1(k)=CV(k,ednew(ii,l))
            v2(k)=CV(k,ednew(jj,l))
 2260 continue
        costheta=
     &    ( v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)) /
     &    (sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2) *
     &     sqrt(v2(1)**2 + v2(2)**2 + v2(3)**2))
          theta=acos(costheta)
          sintheta=sin(theta)
          do 2270 n=1,l-2
            cos1=cos(theta*n/(l-1))
            cos2=cos(theta*(l-1-n)/(l-1))
            alpha=(cos1-costheta*cos2)/sintheta**2
            beta=(cos2-costheta*cos1)/sintheta**2
            do k=1,3
              v3(k)=alpha*v1(k)+beta*v2(k)
            end do
            dnorm = sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)
            do k=1,3
              v3(k) = v3(k) /dnorm
            end do
            do k=1,3
              CV(k,NVPT)=v3(k)
            end do
*
C           trnew(j,l,n+1)=NVPT ! Old code
*
            ixx=j
            jxx=l
            kxx=n+1
            kj  = (jxx-1)*NF    + kxx
            kji = (ixx-1)*NF**2 + kj
            Work(ipTrNew+kji-1)=NVPT
*
            NVPT=NVPT+1
            If(NVPT.gt.1000) then
              Write(6,*)'NVPT out of range in polygen',nvpt
              Call Abend()
            EndIf
 2270 continue
 2250 continue
 2240 continue
      NV=NVPT-1
C
c  -ora per ogni triangolo originario vengono posti nella matrice oldtr
c  i numeri di ordine dei vertici originali,creati lungo gli spigoli e
c  no dei vecchi triangoli secondo lo schema:
c
c          11                    La disposizione e' secondo
c          | \                   la posizione geometrica del
c          |  \                  vertice.
c          21--22
c          | \  |\               i vertici (1,1)-(NF+1,1)-(NF+1,NF+1) sono
c          |  \ | \              quelli originari.
c          31--32--33
c          | \  |  | \           i nuovi triangoli sono:
c          |  \ |  |  \          (i,j)-(i+1,j)-(i+1,j+1)
c          41--42--43--44        i=1,...,NF j=1,...,i
c
c                                (i,j)--(i,j+1),(i+1,j+1)
c                                i=2,...,NF j=1,...,i-1
c
      NTPT=1
      do 2310 n=1,NT0
C
c  -1 vecchi spigoli
C
        oldtr(1,1)=trvo(n,1)
        oldtr(NF+1,1)=trvo(n,2)
        oldtr(NF+1,NF+1)=trvo(n,3)
C
c  -2 nuovi vertici lungo i vecchi spigoli
C
        do 2320 l=2,NF
          oldtr(l,1)=ednew(treo(n,1),l)
          oldtr(NF+1,l)=ednew(treo(n,2),l)
          oldtr(l,l)=ednew(treo(n,3),l)
 2320 continue
C
c  -3 nuovi vertici non lungo i vecchi spigoli
C
        do l=3,NF
          do m=2,l-1
*
C           oldtr(l,m)=trnew(n,l,m) ! old code
*
            ixx=n
            jxx=l
            kxx=m
            kj  = (jxx-1)*NF    + kxx
            kji = (ixx-1)*NF**2 + kj
            oldtr(l,m)=int(Work(ipTrNew-1+kji))
*
      end do
      end do
C
c  -ora si creano i nuovi triangoli
C
        Do 2350 i=1,NF
          Do 2360 j=1,i
            JTR(1,NTPT)=oldtr(i,j)
            JTR(2,NTPT)=oldtr(i+1,j)
            JTR(3,NTPT)=oldtr(i+1,j+1)
            NTPT=NTPT+1
 2360 continue
 2350 continue
        Do i=2,NF
          Do j=1,i-1
            JTR(1,NTPT)=oldtr(i,j)
            JTR(2,NTPT)=oldtr(i,j+1)
            JTR(3,NTPT)=oldtr(i+1,j+1)
            NTPT=NTPT+1
          End Do
        End Do
 2310 continue
      Call GetMem('TrNew','Free','Real',ipTrNew,nTrNew)
C
c scrittura della JVT
C
      do 2390 i=1,NV
         CV(1,i) = CV(1,i) * REN + XEN
         CV(2,i) = CV(2,i) * REN + YEN
         CV(3,i) = CV(3,i) * REN + ZEN
 2390 continue
      return
      end
