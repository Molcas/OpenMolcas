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
c
       subroutine sumovera (wrk,wrksize,
     & lunt2o1,lunt2o2,lunt2o3,
     & lunw3aaaa,lunw3baab,lunw3bbaa,lunw3bbbb,lunw3abba,lunw3aabb)
c
c     This routine do:
c     I) Prepair phase:
c
c     II) sum over A - alfa
c     FI3, WIII1, WIII2, T15, T16, T27, T2Ex
c
c     FI3
cF13.1FI(a,e)aa <- -sum(m,f-aa) [ <ma||ef>aaaa . T1o(f,m)aa ] *
cF13.2FI(a,e)aa <- -sum(m,f-bb) [ <ma||ef>baab . T1o(f,m)bb ] *
cF13.3FI(a,e)bb <-  sum(m,f-bb) [ <ma||fe>bbbb . T1o(f,m)bb ] +
cf13.4FI(a,e)bb <-  sum(m,f-aa) [ <ma||fe>abab . T1o(f,m)aa ] +
c
c     WIII1
cW31.1WIII(m,e,b,j)aaaa <= <mb||ej>aaaa *
cW31.2WIII(m,e,b,j)bbbb <= <mb||ej>bbbb +
cW31.3WIII(m,e,b,j)aabb <= <mb||ej>abab +
cW31.4WIII(m,e,b,j)abba <= <mb||ej>abba +
cW31.5WIII(m,e,b,j)baab <= <mb||ej>baab *
cW31.6WIII(m,e,b,j)bbaa <= <mb||ej>baba *
c
c     WIII2
cW32.1WIII(m,e,b,j)aaaa <- + sum(f-a) [  <mb||ef>aaaa . T1o(f,j)aa ] *
cW32.2WIII(m,e,b,j)bbbb <- _ sum(f-b) [  <mb||ef>bbbb . T1o(f,j)bb ] +
cW32.3WIII(m,e,b,j)aabb <- + sum(f-b) [  <mb||ef>abab . T1o(f,j)bb ] +
cW32.4WIII(m,e,b,j)abba <- - sum(f-a) [  <mb||fe>abab . T1o(f,j)aa ] +
cW32.5WIII(m,e,b,j)baab <-  sum(f-b) [  <mb||ef>baab . T1o(f,j)bb ] *
cW32.6WIII(m,e,b,j)bbaa <- - sum(f-a) [  <mb||fe>baab . T1o(f,j)aa ] *
c
c     T15
cT15.1T1n(a,i)aa <- sum(n,f-aa) [ <na||fi>aaaa . T1o(f,n)aa ] *
cT15.2T1n(a,i)aa <- sum(n,f-bb) [ <na||fi>baba . T1o(f,n)bb ] *
cT15.3T1n(a,i)bb <- sum(n,f-bb) [ <na||fi>bbbb . T1o(f,n)bb ] +
cT15.4T1n(a,i)bb <- sum(n,f-aa) [ <na||fi>abab . T1o(f,n)aa ] +
c
c     T16
cT16.1T1n(a,i)aa <- - sum(m,e>f-aaa) [ T2o(ef,i,m)aaaa  . <ma||ef>aaaa ] *
cT16.2T1n(a,i)aa <- - sum(m,e,f-bab) [ T2o(e,f,i,m)abab . <ma||ef>baab ] *
cT16.3T1n(a,i)bb <- - sum(m,e>f-bbb) [ T2o(ef,i,m)bbbb  . <ma||ef>bbbb ] +
cT16.4T1n(a,i)bb <- + sum(m,e,f-aab) [ T2o(e,f,m,i)abab . <ma||ef>abab ] +
c
c     T2Ex
cT2E.12  R(a,m,ij)aaaa   <= - sum(e>f-aa) [ <ma||ef>aaaa . Tau(ef,ij)aaaa ] *
cT2E.1T2n(ab,ij)aaaa   <- - sum(m-a)    [ T1o(b,m)aa   .  R(a,m,ij)aaaa ] *
cT2E.2T2n(ab,ij)aaaa         <-   sum(m-a)    [ T1o(a,m)aa   .  R(b,m,ij)aaaa ] *
cT2E.34         R(a,m,ij)bbbb   <= - sum(e>f-bb) [ <ma||ef>bbbb . Tau(ef,ij)bbbb ] +
cT2E.3T2n(ab,ij)bbbb   <- - sum(m-b)    [ T1o(b,m)bb   .  R(a,m,ij)bbbb ] +
cT2E.4T2n(ab,ij)bbbb         <-   sum(m-b)    [ T1o(a,m)bb   .  R(b,m,ij)bbbb ] +
cT2E.5R1(a,m,i,j)abab <= - sum(e,f-ab) [ <ma||ef>baab . Tau(e,f,i,j)abab ] *
cT2E.5T2n(a,b,i,j)abab <- - sum(m-b)    [ T1o(m,b)bb . R1(a,m,i,j)baab ] *
cT2E.6R2(a,m,i,j)baab <= - sum(e,f-ab) [ <ma||ef>abab . Tau(e,f,i,j)abab ] +
cT2E.6T2n(a,b,i,j)abab <-   sum(m-a)    [ T1o(m,a)aa . R2(a,m,i,j)abab ] +
c
c     T27
cT27.1G(b,m,i,j)aaaa  <= - sum(e-a)  [ <mb||ej>aaaa . T1o(e,i)aa ]  *
cT27.1Q(b,m,ij)aaaa   <= G(b,m,i,j)aaaa - G(b,m,j,i)aaaa            *
cT27.1R(c,d,ij)aaaa   <= sum(m-a)    [ T1o(c,m)aa   . Q(d,m,ij)aa ] *
cT27.1T2n(ab,ij)aaaa   <- R(a,b,ij)aaaa - R(b,a,ij)aaaa              *
cT27.2G(b,m,i,j)bbbb  <= - sum(e-b)  [ <mb||ej>bbbb . T1o(e,i)bb ]  +
cT27.2Q(b,m,ij)bbbb   <= G(b,m,i,j)bbbb - G(b,m,j,i)bbbb            +
cT27.2R(c,d,ij)bbbb   <= sum(m-b)    [ T1o(c,m)bb   . Q(d,m,ij)bb ] +
cT27.2T2n(ab,ij)bbbb   <- R(a,b,ij)bbbb - R(b,a,ij)bbbb              +
cT27.3G1(b,m,i,j)abab <= + sum(e-a)  [ <mb||ej>baab . T1o(e,i)aa ]  *
cT27.3<- - sum(e-b)  [ <mb||ei>baba . T1o(e,j)bb ]  *
cT27.3T2n(a,b,i,j)abab <- sum(m-b)    [ T1o(b,m)bb   . G1(a,m,i,j)abab ] *
cT27.4G2(b,m,i,j)baab <= - sum(e-a)  [ <mb||ej>abab . T1o(e,i)aa ] +
cT27.4<- + sum(e-b)  [ <mb||ei>abba . T1o(e,j)bb ] +
cT27.4T2n(a,b,i,j)abab <- sum(m-a)    [ T1o(a,m)aa   . G2(b,m,i,j)baab ] +
c
c
c     General status
c
c     Temporary files lunt2o1,lunt2o2,lunt2o3
c     with T2o4, T2o4b and T2o2b must be opened
c
c     Integrals are strored as follows:
c
c     1) <ma||ef> for a - alfa
c     lun = n1aalfa
c     integrals : <ma||ef>aaaa, <ma||ef>baab
c
c     do syma=1,nsym
c     ! one record with mapd,mapi of <ma||ef>aaaa - always
c     ! one record with mapd,mapi of <ma||ef>baab - always
c     do a=1,nva(syma)
c     ! one record with <ma||ef>aaaa if any
c     ! one record with <ma||ef>baab if any
c     end do
c     end do
c
c
c     2) <ma||ef> for a - beta
c     lun = n1abeta
c     integrals : <ma||ef>bbbb, <ma||ef>abab
c
c     do syma=1,nsym
c     ! one record with mapd,mapi of <ma||ef>bbbb - always
c     ! one record with mapd,mapi of <ma||ef>abab - always
c     do a=1,nvb(syma)
c     ! one record with <ma||ef>bbbb if any
c     ! one record with <ma||ef>abab if any
c     end do
c     end do
c
c
c     3) <ma||ej> for a - alfa
c     lun = n2aalfa
c     integrals : <ma||ej>aaaa, <ma||ej>baab, <ma||ej>baba
c
c     do syma=1,nsym
c     ! one record with mapd,mapi of <ma||ej>aaaa - always
c     ! one record with mapd,mapi of <ma||ej>baab - always
c     ! one record with mapd,mapi of <ma||ej>baba - always
c     do a=1,nva(syma)
c     ! one record with <ma||ej>aaaa if any
c     ! one record with <ma||ej>baab if any
c     ! one record with <ma||ej>baba if any
c     end do
c     end do
c
c
c     4) <ma||ej> for a - beta
c     lun = n2abeta
c     integrals : <ma||ej>bbbb, <ma||ej>abab, <ma||ej>abba
c
c     do syma=1,nsym
c     ! one record with mapd,mapi of <ma||ej>bbbb - always
c     ! one record with mapd,mapi of <ma||ej>abba - always
c     ! one record with mapd,mapi of <ma||ej>abab - always
c     do a=1,nva(syma)
c     ! one record with <ma||ej>bbbb if any
c     ! one record with <ma||ej>abba if any
c     ! one record with <ma||ej>abab if any
c     end do
c     end do
c
c     N.B. use and destroy : all
c     N.B. # of reads      : 4
c
c
c        Paralell status
c
c        Alternatives:
c
c        A) 1 proc - idaaaa,idbbbb,idaabb,idabba,idbaab,idbbaa = 1
c
c        B) 2 proc - idaaaa,idbaab,idbbaa = 1
c                   idbbbb,idaabb,idabba = 2
c
c        C) 4 proc - idaaaa        = 1
c                   idbbaa,idbaab = 2
c                   idbbbb        = 3
c                   idabba,idaabb = 4
c
c        D) 6 proc - idaaaa = 1
c                   idbaab = 2
c                   idbbaa = 3
c                   idbbbb = 4
c                   idaabb = 5
c                   idabba = 6
c
c        Pilot nodes (in the case of multiplicity contribution to Tn
c       is evaluated only on pivot node)
c                      : idbaab for alpha, idaabb for beta
c                 i.e. : A-1; B-1,2; C-2,4; D-2,5
c
c        In case D) contributions F13 are calculated separately. In such
c        a case on idaaaa, and idbbbb nodes contr. F13.1 and F13.3 resp.
c        are set (nod add) directly to corresp. F1. On pilot nodes,
c       cont. F13.2 and F13.4 are standardly added to F1's
c
        implicit none
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "paralell.fh"
#include "wrk.fh"
c
       integer lunt2o1,lunt2o2,lunt2o3
       integer lunw3aaaa,lunw3baab,lunw3bbaa,lunw3bbbb,lunw3abba,
     & lunw3aabb
c
c     help variables
c
       integer n1aalfa,n1abeta,n2aalfa,n2abeta
       integer m1lenght,m2lenght,h1lenght,h2lenght,h3lenght
       integer syma,a
       integer rc,posst
       integer ssm3,ssm4,ssh4
c
c     paralell parameters
       integer yesa,yesb
c
c
cA0   paralell
c
c     I.par.1  - escape, if this node is not reserved for sumoverb_a
c
      if ((myRank.eq.idaaaa).or.
     & (myRank.eq.idbaab).or.(myRank.eq.idbbaa)) then
        yesa=1
      else
        yesa=0
      end if
c
      if ((myRank.eq.idbbbb).or.
     & (myRank.eq.idaabb).or.(myRank.eq.idabba)) then
        yesb=1
      else
        yesb=0
      end if
c
        if ((yesa.eq.0).and.(yesb.eq.0)) then
        return
        end if
c
c     ****************   A - prepair phase part 1 -------------------
c
cpar
      if (yesa.eq.1) then
c
cA.1  read T2o4a from opened lunt2o1 file
c     V1(ef,ij)    <- T2o(ef,ij)aaaa
       call filemanager (2,lunt2o1,rc)
       call getmediate (wrk,wrksize,
     & lunt2o1,possv10,mapdv1,mapiv1,rc)
c
cA.2  V2(ij,ef) <- V1(ef,ij)
       call map (wrk,wrksize,
     & 4,3,4,1,2,
     & mapdv1,mapiv1,1,mapdv2,mapiv2,possv20,posst,rc)
c
cA.3  mktau V1(ef,ij) <- V1(ef,ij)
       call mktau (wrk,wrksize,
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,1.0d0,
     &             rc)
c
cA.4  expand V3(i,j,ef) <- V2(ij,ef)
       call expand (wrk,wrksize,
     & 4,5,mapdv2,mapiv2,1,possv30,mapdv3,mapiv3,rc)
c
cA.5
        if ((idaaaa.ne.idbaab).and.(myRank.eq.idaaaa)) then
        call set0 (wrk,wrksize,
     &             mapdf11,mapif11)
        end if
cparend
      end if
c
cpar
      if ((yesa.eq.1).or.(yesb.eq.1)) then
c
cA.6  read T2o2b from opened lunt2o3 file
c     V2(e,f,i,j)  <- T2o(e,f,i,j)abab
       call filemanager (2,lunt2o3,rc)
       call getmediate (wrk,wrksize,
     & lunt2o3,possv20,mapdv2,mapiv2,rc)
c
cA.6  V4(i,j,e,f) <- V2(e,f,i,j)
       call map (wrk,wrksize,
     & 4,3,4,1,2,
     & mapdv2,mapiv2,1,mapdv4,mapiv4,possv40,posst,rc)
c
cA.7  mktau V2(e,f,i,j) <- V2(e,f,i,j)
       call mktau (wrk,wrksize,
     & mapdv2,mapiv2,mapdt11,mapit11,mapdt12,mapit12,1.0d0,
     &             rc)
cparend
      end if
c
c
c     now there is:
c     V1(ef,ij)   = Tau(ef,ij)aaaa
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(i,j,ef)  = T2o(ef,ij)aaaa
c     V4(i,j,e,f) = T2o(e,f,i,j)abab
c     free: M1-M4, H1-H4
c
c
c
c     ****************   B - sum over a alpha ---------------------
c
cpar
       if (yesa.eq.1) then
c
cB.1  open n1aalfa,n2aalfa
       n1aalfa=11
       n2aalfa=13
       call filemanager (4,n1aalfa,rc)
       call filemanager (4,n2aalfa,rc)
c
cB.2  open temp files lunw3aaaa, lunw3baab, lunw3bbaa
cpar
       if (myRank.eq.idaaaa) then
       call filemanager (1,lunw3aaaa,rc)
       end if
c
       if (myRank.eq.idbaab) then
       call filemanager (1,lunw3baab,rc)
       end if
c
       if (myRank.eq.idbbaa) then
       call filemanager (1,lunw3bbaa,rc)
       end if
c
c
       do 3000 syma=1,nsym
c
        if (fullprint.ge.3) then
        write (6,*) ' SymA alfa ',syma
        end if
c
c     storing of mediates:
c
c     V1(ef,ij)   = Tau(ef,ij)aaaa
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(i,j,ef)  = T2o(ef,ij)aaaa
c     V4(i,j,e,f) = T2o(e,f,i,j)abab
c     M1(m,ef)  - <ma||ef>aaaa
c     M2(m,e,f) - <ma||ef>baab
c     H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
c     H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
c     H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa
c     free: M3,M4,H4
c
cB.3  get mapd and mapi for M1,M2
       call getmap (n1aalfa,possm10,m1lenght,mapdm1,mapim1,rc)
       call getmap (n1aalfa,possm20,m2lenght,mapdm2,mapim2,rc)
c
cB.4  get mapd and mapi for H1 - H3
       call getmap (n2aalfa,possh10,h1lenght,mapdh1,mapih1,rc)
       call getmap (n2aalfa,possh20,h2lenght,mapdh2,mapih2,rc)
       call getmap (n2aalfa,possh30,h3lenght,mapdh3,mapih3,rc)
c
cB.5  write mapd and mapi of W3aaaa(H1), W3baab(H2) and W3bbaa(H3)
cpar
       if (myRank.eq.idaaaa) then
       call wrtmap (lunw3aaaa,mapdh1,mapih1,rc)
       end if
c
       if (myRank.eq.idbaab) then
       call wrtmap (lunw3baab,mapdh2,mapih2,rc)
       end if
c
       if (myRank.eq.idbbaa) then
       call wrtmap (lunw3bbaa,mapdh3,mapih3,rc)
       end if
c
c
cB.6  skip cycle over a if lenght of all files is zero
       if ((m1lenght+m2lenght+h1lenght+h2lenght+h3lenght).eq.0) goto
     & 3000
c
       do 2500 a=1,nva(syma)
        if (fullprint.ge.3) then
        write (6,*) ' A alfa ',a
        end if
c
c
       if (h1lenght.gt.0) then
cW31.1.1    read H1(m,e,j) = <ma||ej>aaaa
cpar
       if (myRank.eq.idaaaa) then
       call rea (n2aalfa,h1lenght,wrk(possh10))
       else
       call reajalovy (n2aalfa,h1lenght,wrk(possh10))
       end if
c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)aaaa
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(i,j,ef)  = T2o(ef,ij)aaaa
c     V4(i,j,e,f) = T2o(e,f,i,j)abab
c     M1(m,ef)  - <ma||ef>aaaa (free, but reserved map's)
c     M2(m,e,f) - <ma||ef>baab (free, but reserved map's)
c     H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
c     H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab (free, but reserved map's)
c     H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa (free, but reserved map's)
c     free: M3,M4,H4
c
cpar
        if (myRank.eq.idaaaa) then
c
c     ------- cont to T15
cT15.1T1n(a,i)aa <- sum(n,f-aa) [ <na||fi>aaaa . T1o(f,n)aa ]
c
cT15.1.1    H4(i,f,n) <- H1(n,f,i)
       call map (wrk,wrksize,
     & 3,3,2,1,0,
     & mapdh1,mapih1,syma,mapdh4,mapih4,possh40,posst,rc)
c
cT15.1.2    M3(i) <- H4(i,f,n) . T1o(f,n)aa
       call mult (wrk,wrksize,
     & 3,2,1,2,
     & mapdh4,mapih4,syma,mapdt11,mapit11,1,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT15.1.3    T1n(a,i)aa <- 1.0d0 . M3(i)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,1.0d0,
     & mapdm3,syma,mapdt13,mapit13,1,rc)
c
c     ------- cont to T27
cT27.1G(b,m,i,j)aaaa  <= - sum(e-a)  [ <mb||ej>aaaa . T1o(e,i)aa ]
cT27.1Q(b,m,ij)aaaa   <= G(b,m,i,j)aaaa - G(b,m,j,i)aaaa
cT27.1R(c,d,ij)aaaa   <= sum(m-a)    [ T1o(c,m)aa   . Q(d,m,ij)aa ]
cT27.1T2n(ab,ij)aaaa   <- R(a,b,ij)aaaa - R(b,a,ij)aaaa
c
cT27.1.1    M3(m,j,e) <- H1(m,e,j)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdh1,mapih1,syma,mapdm3,mapim3,possm30,posst,rc)
c
cT27.1.2    H4(m,j,i) <- M3(m,j,e) . T1o(e,i)aa
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm3,mapim3,syma,mapdt11,mapit11,1,mapdh4,mapih4,ssh4,
     & possh40,rc)
c
cT27.1.3    M3(m,ji) <- H4(m,j,i) - H4(m,i,j)
       call fack (wrk,wrksize,
     & 3,2,mapdh4,syma,mapih4,mapdm3,mapim3,possm30,rc)
c
cT27.1.4    M4(c,ji) <- T1o(c,m)aa . M3(m,ji)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt11,mapit11,1,mapdm3,mapim3,syma,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cT27.1.5    T2n(ab,ij) <- - (M4(a,ji)-M4(b,ji) = M4(b,ij)-M4(a,ij)
c     since -1 was skipped in first step (T27.1)
       call add (wrk,wrksize,
     & 3,4,1,2,a,0,syma,1,1.0d0,
     & mapdm4,syma,mapdt21,mapit21,1,rc)
c
cparend
       end if
c
       end if
c
c
c
c
       if (h2lenght.gt.0) then
cW31.5read H2(m,e,j) = <ma||ej>baab
cpar
        if (myRank.eq.idbaab) then
       call rea (n2aalfa,h2lenght,wrk(possh20))
        else
       call reajalovy (n2aalfa,h2lenght,wrk(possh20))
        end if

c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)aaaa
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(i,j,ef)  = T2o(ef,ij)aaaa
c     V4(i,j,e,f) = T2o(e,f,i,j)abab
c     M1(m,ef)  - <ma||ef>aaaa (free, but reserved map's)
c     M2(m,e,f) - <ma||ef>baab (free, but reserved map's)
c     H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
c     H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
c     H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa (free, but reserved map's)
c     free: M3,M4,H4
c
cpar
        if (myRank.eq.idbaab) then
c
c     ------- cont to T27 (part)
cT27.3p            G1(b,m,i,j)abab <= + sum(e-a)  [ <mb||ej>baab . T1o(e,i)aa ]
c
cT27.3.1    M4(m,j,e) <- H2(m,e,j)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdh2,mapih2,syma,mapdm4,mapim4,possm40,posst,rc)
c
cT27.3.2    M3(m,j,i) <- M4(m,j,e) . T1o(e,i)
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm4,mapim4,syma,mapdt11,mapit11,1,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT27.3.3    M4(m,i,j) <- M3(m,j,i)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdm3,mapim3,syma,mapdm4,mapim4,possm40,posst,rc)
c
cpar
        if (idbaab.ne.idbbaa) then
c        if idbaab and idbbaa are different, we need to add contribution
c       from this part of G1 to T2n, since they will be not present on
c       idbbaa
c
cT27.3  T2n(a,b,i,j)abab <- sum(m-b)    [ T1o(b,m)bb   . G1(a,m,i,j)abab ]
c
cT27.3.7    M3(b,i,j) <- T1o(b,m)bb . M4(m,i,j)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt12,mapit12,1,mapdm4,mapim4,syma,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT27.3.8    T2n(a,b,i,j) <- 1.0d0 . M3(b,i,j)
       call add (wrk,wrksize,
     & 3,4,1,1,a,0,syma,1,1.0d0,
     & mapdm3,syma,mapdt23,mapit23,1,rc)
cparend
        end if
c
cparend
       end if
c
       end if
c
c
c
c
       if (h3lenght.gt.0) then
cW31.6read H3(m,e,j) = <ma||ej>baba
cpar
         if (myRank.eq.idbbaa) then
       call rea (n2aalfa,h3lenght,wrk(possh30))
        else
       call reajalovy (n2aalfa,h3lenght,wrk(possh30))
        end if
c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)aaaa
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(i,j,ef)  = T2o(ef,ij)aaaa
c     V4(i,j,e,f) = T2o(e,f,i,j)abab
c     M1(m,ef)  - <ma||ef>aaaa (free, but reserved map's)
c     M2(m,e,f) - <ma||ef>baab (free, but reserved map's)
c     M4(m,j,i) - part of G1 intermediat from previos step (only if idbaab=idbbaa)
c     H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
c     H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
c     H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa
c     free: M3,H4
c
cpar
        if (myRank.eq.idbbaa) then
c
c     ------- cont to T27 (continue)
cT27.3G1(b,m,i,j)abab <- - sum(e-b)  [ <mb||ei>baba . T1o(e,j)bb ]
cT27.3T2n(a,b,i,j)abab <- sum(m-b)    [ T1o(b,m)bb   . G1(a,m,i,j)abab ]
c     part of G1 is in M4(m,i,j) from part T27.3.3
c
cT27.3.4    M3(m,i,e) <- H3(m,e,i)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdh3,mapih3,syma,mapdm3,mapim3,possm30,posst,rc)
c
cT27.3.5    H4(m,i,j) <- M3(m,i,e) . T1o(e,j)bb
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm3,mapim3,syma,mapdt12,mapit12,1,mapdh4,mapih4,ssh4,
     & possh40,rc)
c
cpar
        if (idbaab.eq.idbbaa) then
c        if idbaab=idbbaa add H4 to M4 from prev step and realize
c        contribution from both parts of G1 in one step
c
cT27.3.6    (G1) M4(m,i,j) <- -H4(m,i,j)
       call add (wrk,wrksize,
     & 3,3,0,0,0,0,1,1,-1.0d0,
     & mapdh4,syma,mapdm4,mapim4,syma,rc)
c
cT27.3.7    M3(b,i,j) <- T1o(b,m)bb . M4(m,i,j)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt12,mapit12,1,mapdm4,mapim4,syma,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT27.3.8    T2n(a,b,i,j) <- 1.0d0 . M3(b,i,j)
       call add (wrk,wrksize,
     & 3,4,1,1,a,0,syma,1,1.0d0,
     & mapdm3,syma,mapdt23,mapit23,1,rc)
c
        else
c        if idbaab.ne.idbbaa we need to calc only contribution from
c        2.nd part of G1
c
cT27.3.7    M3(b,i,j) <- T1o(b,m)bb . H4(m,i,j)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt12,mapit12,1,mapdh4,mapih4,syma,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT27.3.8    T2n(a,b,i,j) <- -1.0d0 . M3(b,i,j)
c           (minus sign, since we skip sign in previous step)
       call add (wrk,wrksize,
     & 3,4,1,1,a,0,syma,1,-1.0d0,
     & mapdm3,syma,mapdt23,mapit23,1,rc)
cparend
        end if
c
c     ------- cont to T15
cT15.2T1n(a,i)aa <- sum(n,f-bb) [ <na||fi>baba . T1o(f,n)bb ]
c
cT15.2.1    H4(i,f,n) <- H3(n,f,i)
       call map (wrk,wrksize,
     & 3,3,2,1,0,
     & mapdh3,mapih3,syma,mapdh4,mapih4,possh40,posst,rc)
c
cT15.2.2    M3(i) <- H4(i,f,n) . T1o(f,n)bb
       call mult (wrk,wrksize,
     & 3,2,1,2,
     & mapdh4,mapih4,syma,mapdt12,mapit12,1,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT15.2.3    T1n(a,i)aa <- 1.0d0 . M3(i)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,1.0d0,
     & mapdm3,syma,mapdt13,mapit13,1,rc)
c
cparend
       end if
c
       end if
c
c
c
c
       if (m1lenght.gt.0) then
c
cB.7  read M1(m,ef) = <ma||ef>aaaa
cpar
        if (myRank.eq.idaaaa) then
       call rea (n1aalfa,m1lenght,wrk(possm10))
        else
       call reajalovy (n1aalfa,m1lenght,wrk(possm10))
        end if
c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)aaaa
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(i,j,ef)  = T2o(ef,ij)aaaa
c     V4(i,j,e,f) = T2o(e,f,i,j)abab
c     M1(m,ef)  - <ma||ef>aaaa
c     M2(m,e,f) - <ma||ef>baab (free, but reserved map's)
c     H1(m,e,j) - <ma||ej>aaaa = W3(m,e,a,j)aaaa
c     H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
c     H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa
c     free: M3,M4,H4
c
cpar
        if (myRank.eq.idaaaa) then
c
c     ------- cont to T16
cT16.1T1n(a,i)aa <- - sum(m,e>f-aaa) [ T2o(ef,i,m)aaaa  . <ma||ef>aaaa ]
c
cT16.1.1    M3(i) = V3(i,m,ef). M1(m,ef)
       call mult (wrk,wrksize,
     & 4,3,1,3,
     & mapdv3,mapiv3,1,mapdm1,mapim1,syma,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT16.1.2    t1n(a,i)aa <- -1.0d0 . M3(i)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,-1.0d0,
     & mapdm3,syma,mapdt13,mapit13,1,rc)
c
c     ------- cont to T2Ex
cT2E.12             R(a,m,ij)aaaa   <= - sum(e>f-aa) [ <ma||ef>aaaa . Tau(ef,ij)aaaa ]
cT2E.1T2n(ab,ij)aaaa   <- - sum(m-a)    [ T1o(b,m)aa   .  R(a,m,ij)aaaa ]
cT2E.2T2n(ab,ij)aaaa   <-   sum(m-a)    [ T1o(a,m)aa   .  R(b,m,ij)aaaa ]
c
cT2E.12.1   M3(m,ij) = M1(m,ef) . V1(ef,ij)
       call mult (wrk,wrksize,
     & 3,4,3,2,
     & mapdm1,mapim1,syma,mapdv1,mapiv1,1,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT2E.1,2.1  M4(b,ij) = T1o(b,m)aa . M3(m,ij)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt11,mapit11,1,mapdm3,mapim3,ssm3,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cT2E.1,2.2    T2n(ab,ij)aaaa <- 1.0d0 . M4(b,ij) (@@ pozor na factor @@)
       call add (wrk,wrksize,
     & 3,4,1,1,a,0,syma,1,1.0d0,
     & mapdm4,syma,mapdt21,mapit21,1,rc)
c
c     ------- cont to W32
cW32.1WIII(m,e,b,j)aaaa <- + sum(f-a) [  <mb||ef>aaaa . T1o(f,j)aa ]
c
cW32.1.1    M3(m,e,f) <- M1(m,ef)
       call expand (wrk,wrksize,
     & 3,2,mapdm1,mapim1,syma,possm30,mapdm3,mapim3,rc)
c
cW32.1.2    M4(m,e,j) <- M3(m,e,f) . T1o(f,j)aa
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm3,mapim3,syma,mapdt11,mapit11,1,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cW32.1.3    H1(W3aaaa)(m,e,j) <- +1.0d0 . M4(m,e,j)
       call add (wrk,wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0,
     & mapdm4,syma,mapdh1,mapih1,syma,rc)
c
cW3.1 write W3aaaa(m,e,j) to  lunw3aaaa file if size is not zero
       if (h1lenght.gt.0) then
       call wri (lunw3aaaa,h1lenght,wrk(possh10))
       end if
c
c     ------- cont to FI3
cF13.1FI(a,e)aa <- -sum(m,f-aa) [ <ma||ef>aaaa . T1o(f,m)aa ]
c
cF13.1.1    M4(e,f,m) <- M3(m,e,f)  (M3 from previous part is not destroyed)
       call map (wrk,wrksize,
     & 3,3,1,2,0,
     & mapdm3,mapim3,syma,mapdm4,mapim4,possm40,posst,rc)
c
cF13.1.2    H4(e) <- M4(e,f,m) . T1o(f,m)aa
       call mult (wrk,wrksize,
     & 3,2,1,2,
     & mapdm4,mapim4,syma,mapdt11,mapit11,1,mapdh4,mapih4,ssh4,
     & possh40,rc)
c
cF13.1.3    F1(a,e)aa <- -1.0d0 . H4(e)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,-1.0d0,
     & mapdh4,syma,mapdf11,mapif11,1,rc)
c
cparend
       end if
c
       end if
c
c
c
c
       if (m2lenght.gt.0) then
cB.8  read M2(m,e,f) = <ma||e,f>baab
cpar
        if ((myRank.eq.idbaab).or.(myRank.eq.idbbaa)) then
       call rea (n1aalfa,m2lenght,wrk(possm20))
        else
       call reajalovy (n1aalfa,m2lenght,wrk(possm20))
        end if
c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)aaaa
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(i,j,ef)  = T2o(ef,ij)aaaa
c     V4(i,j,e,f) = T2o(e,f,i,j)abab
c     M1(m,ef)  - <ma||ef>aaaa (free, but reserved map's)
c     M2(m,e,f) - <ma||ef>baab
c     H1(m,e,j) - <ma||ej>aaaa  (free, but reserved map's)
c     H2(m,e,j) - <ma||ej>baab = W3(m,e,a,j)baab
c     H3(m,e,j) - <ma||ej>baba = W3(m,e,a,j)bbaa
c     free: M3,M4,H4
c
c
cpar
        if (myRank.eq.idbaab) then
c
c     ------- cont to T16
cT16.2T1n(a,i)aa <- - sum(m,e,f-bab) [ T2o(e,f,i,m)abab . <ma||ef>baab ]
c
cT16.2.1    M3(i) = V4(i,m,e,f) . M2(m,e,f)
       call mult (wrk,wrksize,
     & 4,3,1,3,
     & mapdv4,mapiv4,1,mapdm2,mapim2,syma,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT16.2.2    t1n(a,i)aa <- -1.0d0 . M3(i)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,-1.0d0,
     & mapdm3,syma,mapdt13,mapit13,1,rc)
c
c     ------- cont to T2Ex
cT2E.5R1(a,m,i,j)abab <= - sum(e,f-ab) [ <ma||ef>baab . Tau(e,f,i,j)abab ]
cT2E.5T2n(a,b,i,j)abab <- - sum(m-b)    [ T1o(m,b)bb . R1(a,m,i,j)baab ]
c
cT2E.5.1    M4(m,i,j) <- M2(m,e,f) . V2(e,f,i,j)
       call mult (wrk,wrksize,
     & 3,4,3,2,
     & mapdm2,mapim2,syma,mapdv2,mapiv2,1,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cT2E.5.2    M3(b,i,j) <- T1o(b,m)bb . M4(m,i,j)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt12,mapit12,1,mapdm4,mapim4,syma,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT2E.5.3    T2n(a,b,i,j)abab <- 1.0d0 M3(b,i,j)
       call add (wrk,wrksize,
     & 3,4,1,1,a,0,syma,1,1.0d0,
     & mapdm3,ssm3,mapdt23,mapit23,1,rc)
cparend
        end if
c
c
c     ------- cont to W32
cW32.5WIII(m,e,b,j)baab <-  sum(f-b) [  <mb||ef>baab . T1o(f,j)bb ]
cW32.6WIII(m,e,b,j)bbaa <- - sum(f-a) [  <mb||fe>baab . T1o(f,j)aa ]
c
cpar
        if (myRank.eq.idbaab) then
cW32.5.1    M4(m,e,j) = M2(m,e,f) . T1o(f,j)bb
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm2,mapim2,syma,mapdt12,mapit12,1,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cW32.5.2    H2(W3baab)(m,e,j) <- 1.0d0 M4(m,e,j)
       call add (wrk,wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0,
     & mapdm4,syma,mapdh2,mapih2,syma,rc)
c
cW3.5 write W3baab(m,e,j) to  lunw3baab file if size is not zero
       if (h2lenght.gt.0) then
       call wri (lunw3baab,h2lenght,wrk(possh20))
       end if
cendpar
        end if
c
c
cpar
        if (myRank.eq.idbbaa) then
cW32.6.1    M4(m,e,f) <- M2(m,f,e)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdm2,mapim2,syma,mapdm4,mapim4,possm40,posst,rc)
c
cW32.6.2    M3(m,e,j) <- M4(m,e,f) . T1o(f,j)aa
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm4,mapim4,syma,mapdt11,mapit11,1,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cW32.6.3    H3(W3bbaa)(m,e,j) <- -1.0d0 M3(m,e,j)
       call add (wrk,wrksize,
     & 3,3,0,0,0,0,1,1,-1.0d0,
     & mapdm3,syma,mapdh3,mapih3,syma,rc)
c
cW3.6 write W3bbaa(m,e,j) to  lunw3bbaa file if size is not zero
       if (h3lenght.gt.0) then
       call wri (lunw3bbaa,h3lenght,wrk(possh30))
       end if
cparend
        end if
c
c     ------- cont to FI3
cpar
        if (myRank.eq.idbaab) then
cF13.2FI(a,e)aa <- -sum(m,f-bb) [ <ma||ef>baab . T1o(f,m)bb ]
c
cF13.2.1    M3(e,f,m) <- M2(m,e,f)
       call map (wrk,wrksize,
     & 3,3,1,2,0,
     & mapdm2,mapim2,syma,mapdm3,mapim3,possm30,posst,rc)
c
cF13.2.2    M4(e)     <- M3(e,f,m) . T1o(f,m)bb
       call mult (wrk,wrksize,
     & 3,2,1,2,
     & mapdm3,mapim3,syma,mapdt12,mapit12,1,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cF13.2.3    F1(a,e)aa <- -1.0d0 M4(e)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,-1.0d0,
     & mapdm4,syma,mapdf11,mapif11,1,rc)
       end if
cparend
        end if
c
 2500   continue
c
 3000   continue
c
cB.9  close n1aalfa,n2aalfa
       call filemanager (3,n1aalfa,rc)
       call filemanager (3,n2aalfa,rc)
cparend
       end if
c
c
c     ****************   A - prepair phase part 2 -------------------
c
c
cpar
       if (yesb.eq.1) then
cA.10 V3(j,i,e,f) <- V4(i,j,e,f)
       call map (wrk,wrksize,
     & 4,2,1,3,4,
     & mapdv4,mapiv4,1,mapdv3,mapiv3,possv30,posst,rc)
c
cA.11 read T2o4b from opened lunt2o1 file
c     V4(ef,ij)    <- T2o(ef,ij)aaaa
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv40,mapdv4,mapiv4,rc)
c
cA.12 V1(ij,ef) <- V4(ef,ij)
       call map (wrk,wrksize,
     & 4,3,4,1,2,
     & mapdv4,mapiv4,1,mapdv1,mapiv1,possv10,posst,rc)
c
cA.13 expand V4(i,j,ef) <- V1(ij,ef)
       call expand (wrk,wrksize,
     & 4,5,mapdv1,mapiv1,1,possv40,mapdv4,mapiv4,rc)
c
cA.14 read T2o4b from opened lunt2o1 file
c     V1(ef,ij)    <- T2o(ef,ij)bbbb
       call filemanager (2,lunt2o2,rc)
       call getmediate (wrk,wrksize,
     & lunt2o2,possv10,mapdv1,mapiv1,rc)
c
cA.15 mktau V1(ef,ij) <- V1(ef,ij)
       call mktau (wrk,wrksize,
     & mapdv1,mapiv1,mapdt11,mapit11,mapdt12,mapit12,1.0d0,
     &             rc)
cA.16
        if ((idbbbb.ne.idaabb).and.(myRank.eq.idbbbb)) then
       call set0 (wrk,wrksize,
     &             mapdf12,mapif12)
        end if
c
c
c     now there is:
c     V1(ef,ij)   = Tau(ef,ij)bbbb
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(j,i,e,f) = T2o(e,f,i,j)abab
c     V4(i,j,ef)  = T2o(ef,ij)bbbb
c     free: M1-M4, H1-H4
c
c
c     ****************   C - sum over a beta ---------------------
c
cC.1  open n1abeta,n2abeta
       n1abeta=12
       n2abeta=14
       call filemanager (4,n1abeta,rc)
       call filemanager (4,n2abeta,rc)
c
cC.2  open temp files lunw3bbbb, lunw3abba, lunw3aabb
cpar
       if (myRank.eq.idbbbb) then
       call filemanager (1,lunw3bbbb,rc)
       end if
c
       if (myRank.eq.idabba) then
       call filemanager (1,lunw3abba,rc)
       end if
c
       if (myRank.eq.idaabb) then
       call filemanager (1,lunw3aabb,rc)
       end if
c
c
       do 6000 syma=1,nsym
c
        if (fullprint.ge.3) then
        write (6,*) ' SymA beta ',syma
        end if
c
c     storing of mediates:
c
c     V1(ef,ij)   = Tau(ef,ij)bbbb
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(j,i,e,f) = T2o(e,f,i,j)abab
c     V4(i,j,ef)  = T2o(ef,ij)bbbb
c     M1(m,ef)  - <ma||ef>bbbb
c     M2(m,e,f) - <ma||ef>abab
c     H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
c     H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
c     H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb
c     free: M3,M4,H4
c
cC.3  get mapd and mapi for M1,M2
       call getmap (n1abeta,possm10,m1lenght,mapdm1,mapim1,rc)
       call getmap (n1abeta,possm20,m2lenght,mapdm2,mapim2,rc)
c
cC.4  get mapd and mapi for H1 - H3
       call getmap (n2abeta,possh10,h1lenght,mapdh1,mapih1,rc)
       call getmap (n2abeta,possh20,h2lenght,mapdh2,mapih2,rc)
       call getmap (n2abeta,possh30,h3lenght,mapdh3,mapih3,rc)
c
cC.5  write mapd and mapi of W3bbbb(H1), W3abba(H2) and W3aabb(H3)
cpar
        if (myRank.eq.idbbbb) then
       call wrtmap (lunw3bbbb,mapdh1,mapih1,rc)
        end if
c
        if (myRank.eq.idabba) then
       call wrtmap (lunw3abba,mapdh2,mapih2,rc)
        end if
c
        if (myRank.eq.idaabb) then
       call wrtmap (lunw3aabb,mapdh3,mapih3,rc)
        end if
cparend
c
c
cC.6  skip cycle over a if lenght of all files is zero
       if ((m1lenght+m2lenght+h1lenght+h2lenght+h3lenght).eq.0) goto
     & 6000
c
       do 5500 a=1,nvb(syma)
        if (fullprint.ge.3) then
        write (6,*) ' A beta ',a
        end if
c
c
       if (h1lenght.gt.0) then
cW31.2.1    read H1(m,e,j) = <ma||ej>bbbb
cpar
        if (myRank.eq.idbbbb) then
       call rea (n2abeta,h1lenght,wrk(possh10))
        else
       call reajalovy (n2abeta,h1lenght,wrk(possh10))
        end if
c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)bbbb
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(j,i,e,f) = T2o(e,f,i,j)abab
c     V4(i,j,ef)  = T2o(ef,ij)bbbb
c     M1(m,ef)  - <ma||ef>bbbb (free, but reserved map's)
c     M2(m,e,f) - <ma||ef>abab (free, but reserved map's)
c     H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
c     H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba (free, but reserved map's)
c     H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb (free, but reserved map's)
c     free: M3,M4,H4
c
cpar
        if (myRank.eq.idbbbb) then
c
c     ------- cont to T15
cT15.3T1n(a,i)bb <- sum(n,f-bb) [ <na||fi>bbbb . T1o(f,n)bb ]
c
cT15.3.1    H4(i,f,n) <- H1(n,f,i)
       call map (wrk,wrksize,
     & 3,3,2,1,0,
     & mapdh1,mapih1,syma,mapdh4,mapih4,possh40,posst,rc)
c
cT15.3.2    M3(i) <- H4(i,f,n) . T1o(f,n)bb
       call mult (wrk,wrksize,
     & 3,2,1,2,
     & mapdh4,mapih4,syma,mapdt12,mapit12,1,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT15.3.3    T1n(a,i)aa <- 1.0d0 . M3(i)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,1.0d0,
     & mapdm3,syma,mapdt14,mapit14,1,rc)
c
c     ------- cont to T27
cT27.2G(b,m,i,j)bbbb  <= - sum(e-b)  [ <mb||ej>bbbb . T1o(e,i)bb ]
cT27.2Q(b,m,ij)bbbb   <= G(b,m,i,j)bbbb - G(b,m,j,i)bbbb
cT27.2R(c,d,ij)bbbb   <= sum(m-b)    [ T1o(c,m)bb   . Q(d,m,ij)bb ]
cT27.2T2n(ab,ij)bbbb   <- R(a,b,ij)bbbb - R(b,a,ij)bbbb
c
cT27.2.1    M3(m,j,e) <- H1(m,e,j)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdh1,mapih1,syma,mapdm3,mapim3,possm30,posst,rc)
c
cT27.2.2    H4(m,j,i) <- M3(m,j,e) . T1o(e,i)bb
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm3,mapim3,syma,mapdt12,mapit12,1,mapdh4,mapih4,ssh4,
     & possh40,rc)
c
cT27.2.3    M3(m,ji) <- H4(m,j,i) - H4(m,i,j)
       call fack (wrk,wrksize,
     & 3,2,mapdh4,syma,mapih4,mapdm3,mapim3,possm30,rc)
c
cT27.2.4    M4(c,ji) <- T1o(c,m)bb . M3(m,ji)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt12,mapit12,1,mapdm3,mapim3,syma,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cT27.2.5    T2n(ab,ij) <- - (M4(a,ji)-M4(b,ji) = M4(b,ij)-M4(a,ij)
c     since -1 was skipped in first step (T27.1)
       call add (wrk,wrksize,
     & 3,4,1,2,a,0,syma,1,1.0d0,
     & mapdm4,syma,mapdt22,mapit22,1,rc)
cparend
       end if
c
       end if
c
c
c
c
       if (h2lenght.gt.0) then
cW31.4read H2(m,e,j) = <ma||ej>abba
cpar
        if (myRank.eq.idabba) then
       call rea (n2abeta,h2lenght,wrk(possh20))
        else
       call reajalovy (n2abeta,h2lenght,wrk(possh20))
        end if
c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)bbbb
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(j,i,e,f) = T2o(e,f,i,j)abab
c     V4(i,j,ef)  = T2o(ef,ij)bbbb
c     M1(m,ef)  - <ma||ef>bbbb (free, but reserved map's)
c     M2(m,e,f) - <ma||ef>abab (free, but reserved map's)
c     H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
c     H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
c     H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb (free, but reserved map's)
c     free: M3,M4,H4
c
cpar
        if (myRank.eq.idabba) then
c
c     ------- cont to T27 (part)
cT27.4p      G2(b,m,i,j)baab <- + sum(e-b)  [ <mb||ei>abba . T1o(e,j)bb ] +
c
cT27.4.1    M3(m,i,e) <- H2(m,e,i)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdh2,mapih2,syma,mapdm3,mapim3,possm30,posst,rc)
c
cT27.4.2    M4(m,i,j) <- M3(m,i,e) . T1o(e,j)bb
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm3,mapim3,syma,mapdt12,mapit12,1,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cpar
        if (idaabb.ne.idabba) then
c        if idaabb and idabba are different, we need to add contribution
c       from this part of G2 to T2n, since they will be not present on
c       idaabb
c
cT27.4        T2n(a,b,i,j)abab <- sum(m-a)    [ T1o(a,m)aa   . G2(b,m,i,j)baab ] +
c
cT27.4.7    H4(a,i,j) <- T1o(a,m)aa . M4(m,i,j)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt11,mapit11,1,mapdm4,mapim4,syma,mapdh4,mapih4,ssh4,
     & possh40,rc)
c
cT27.4.8    T2n(c,a,i,j) <- 1.0d0 . H4(c,i,j)
       call add (wrk,wrksize,
     & 3,4,1,2,a,0,syma,1,1.0d0,
     & mapdh4,syma,mapdt23,mapit23,1,rc)
cparend
        end if
c
cparend
       end if
c
       end if
c
c
c
c
       if (h3lenght.gt.0) then
cW31.3read H3(m,e,j) = <ma||ej>abab
cpar
        if (myRank.eq.idaabb) then
       call rea (n2abeta,h3lenght,wrk(possh30))
        else
       call reajalovy (n2abeta,h3lenght,wrk(possh30))
        end if
c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)bbbb
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(j,i,e,f) = T2o(e,f,i,j)abab
c     V4(i,j,ef)  = T2o(ef,ij)bbbb
c     M1(m,ef)  - <ma||ef>bbbb (free, but reserved map's)
c     M2(m,e,f) - <ma||ef>abab (free, but reserved map's)
c     M4(m,j,i) - part of G2 intermediat from previos step
c     H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
c     H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
c     H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb
c     free: M3,H4
c
cpar
        if (myRank.eq.idaabb) then
c
c     ------- cont to T27 (continue)
cT27.4G2(b,m,i,j)baab <= - sum(e-a)  [ <mb||ej>abab . T1o(e,i)aa ] +
cT27.4T2n(a,b,i,j)abab <- sum(m-a)    [ T1o(a,m)aa   . G2(b,m,i,j)baab ] +
c     N.B. part of G2(m,i,j) is in M4 from part T27.4.2
c
cT27.4.3    M3(m,j,e) <- H3(m,e,j)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdh3,mapih3,syma,mapdm3,mapim3,possm30,posst,rc)
c
cT27.4.4    H4(m,j,i) <- M3(m,j,e) . T1o(e,i)aa
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm3,mapim3,syma,mapdt11,mapit11,1,mapdh4,mapih4,ssh4,
     & possh40,rc)
c
cT27.4.5    M3(m,i,j) <- H4(m,j,i)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdh4,mapih4,syma,mapdm3,mapim3,possm30,posst,rc)
c
cpar
        if (idaabb.eq.idabba) then
c        if idaabb=idabba add H4 to M4 from prev step and realize
c        contribution from both parts of G2 in one step

cT27.4.6    (G2) M4(m,i,j) <- - M3(m,i,j)
       call add (wrk,wrksize,
     & 3,3,0,0,0,0,1,1,-1.0d0,
     & mapdm3,syma,mapdm4,mapim4,syma,rc)
c
cT27.4.7    H4(a,i,j) <- T1o(a,m)aa . M4(m,i,j)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt11,mapit11,1,mapdm4,mapim4,syma,mapdh4,mapih4,ssh4,
     & possh40,rc)
c
cT27.4.8    T2n(c,a,i,j) <- 1.0d0 . H4(c,i,j)
       call add (wrk,wrksize,
     & 3,4,1,2,a,0,syma,1,1.0d0,
     & mapdh4,syma,mapdt23,mapit23,1,rc)
c
        else
c        if idabba.ne.idaabb we need to calc only contribution from
c        2.nd part of G2
c
cT27.4.7    H4(a,i,j) <- T1o(a,m)aa . M3(m,i,j)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt11,mapit11,1,mapdm3,mapim3,syma,mapdh4,mapih4,ssh4,
     & possh40,rc)
c
cT27.4.8    T2n(c,a,i,j) <- 1.0d0 . -H4(c,i,j)
c           (minus sign, since we skip sign in previous step)
       call add (wrk,wrksize,
     & 3,4,1,2,a,0,syma,1,-1.0d0,
     & mapdh4,syma,mapdt23,mapit23,1,rc)
cparend
        end if
c
c
c     ------- cont to T15
cT15.4T1n(a,i)bb <- sum(n,f-aa) [ <na||fi>abab . T1o(f,n)aa ]
c
cT15.4.1    H4(i,f,n) <- H3(n,f,i)
       call map (wrk,wrksize,
     & 3,3,2,1,0,
     & mapdh3,mapih3,syma,mapdh4,mapih4,possh40,posst,rc)
c
cT15.4.2    M3(i) <- H4(i,f,n) . T1o(f,n)aa
       call mult (wrk,wrksize,
     & 3,2,1,2,
     & mapdh4,mapih4,syma,mapdt11,mapit11,1,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT15.4.3    T1n(a,i)bb <- 1.0d0 . M3(i)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,1.0d0,
     & mapdm3,syma,mapdt14,mapit14,1,rc)
cparend
       end if
c
       end if
c
c
c
c
       if (m1lenght.gt.0) then
c
cC.7  read M1(m,ef) = <ma||ef>bbbb
cpar
        if (myRank.eq.idbbbb) then
       call rea (n1abeta,m1lenght,wrk(possm10))
        else
       call reajalovy (n1abeta,m1lenght,wrk(possm10))
        end if
c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)bbbb
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(j,i,e,f) = T2o(e,f,i,j)abab
c     V4(i,j,ef)  = T2o(ef,ij)bbbb
c     M1(m,ef)  - <ma||ef>bbbb
c     M2(m,e,f) - <ma||ef>abab (free, but reserved map's)
c     H1(m,e,j) - <ma||ej>bbbb = W3(m,e,a,j)bbbb
c     H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
c     H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb
c     free: M3,M4,H4
c
cpar
        if (myRank.eq.idbbbb) then
c
c     ------- cont to T16
cT16.3T1n(a,i)bb <- - sum(m,e>f-bbb) [ T2o(ef,i,m)bbbb  . <ma||ef>bbbb ]
c
cT16.3.1*   M3(i) = V4(i,m,ef). M1(m,ef)
       call mult (wrk,wrksize,
     & 4,3,1,3,
     & mapdv4,mapiv4,1,mapdm1,mapim1,syma,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT16.3.2    t1n(a,i)bb <- -1.0d0 . M3(i)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,-1.0d0,
     & mapdm3,syma,mapdt14,mapit14,1,rc)
c
c     ------- cont to T2Ex
cT2E.34             R(a,m,ij)bbbb   <= - sum(e>f-bb) [ <ma||ef>bbbb . Tau(ef,ij)bbbb ] +
cT2E.3T2n(ab,ij)bbbb   <- - sum(m-b)    [ T1o(b,m)bb   .  R(a,m,ij)bbbb ] +
cT2E.4T2n(ab,ij)bbbb   <-   sum(m-b)    [ T1o(a,m)bb   .  R(b,m,ij)bbbb ] +
c
cT2E.34.1   M3(m,ij) = M1(m,ef) . V1(ef,ij)
       call mult (wrk,wrksize,
     & 3,4,3,2,
     & mapdm1,mapim1,syma,mapdv1,mapiv1,1,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT2E.3,4.1  M4(b,ij) = T1o(b,m)bb . M3(m,ij)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt12,mapit12,1,mapdm3,mapim3,ssm3,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cT2E.3,4.2    T2n(ab,ij)bbbb <- 1.0d0 . M4(b,ij) (@@ pozor na factor @@)
       call add (wrk,wrksize,
     & 3,4,1,1,a,0,syma,1,1.0d0,
     & mapdm4,ssm4,mapdt22,mapit22,1,rc)
c
c     ------- cont to W32
cW32.2WIII(m,e,b,j)bbbb <- + sum(f-b) [  <mb||ef>bbbb . T1o(f,j)bb ]
c
cW32.2.1    M3(m,e,f) <- M1(m,ef)
       call expand (wrk,wrksize,
     & 3,2,mapdm1,mapim1,syma,possm30,mapdm3,mapim3,rc)
c
cW32.2.2    M4(m,e,j) <- M3(m,e,f) . T1o(f,j)bb
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm3,mapim3,syma,mapdt12,mapit12,1,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cW32.2.3    H1(W3bbbb)(m,e,j) <- +1.0d0 . M4(m,e,j)
       call add (wrk,wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0,
     & mapdm4,syma,mapdh1,mapih1,syma,rc)
c
cW3.2 write W3bbbb(m,e,j) to  lunw3bbbb file if size is not zero
       if (h1lenght.gt.0) then
       call wri (lunw3bbbb,h1lenght,wrk(possh10))
       end if
c
c     ------- cont to FI3
cF13.3FI(a,e)bb <-  sum(m,f-bb) [ <ma||fe>bbbb . T1o(f,m)bb ]
c
cF13.3.1    M4(e,f,m) <- M3(m,f,e)  (M3 from previous part is not destroyed)
       call map (wrk,wrksize,
     & 3,3,2,1,0,
     & mapdm3,mapim3,syma,mapdm4,mapim4,possm40,posst,rc)
c
cF13.3.2    H4(e) <- M4(e,f,m) . T1o(f,m)bb
       call mult (wrk,wrksize,
     & 3,2,1,2,
     & mapdm4,mapim4,syma,mapdt12,mapit12,1,mapdh4,mapih4,ssh4,
     & possh40,rc)
c
cF13.3.3    F1(a,e)bb <- 1.0d0 . H4(e)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,1.0d0,
     & mapdh4,syma,mapdf12,mapif12,1,rc)
cparend
       end if
c
       end if
c
c
c
c
       if (m2lenght.gt.0) then
cB.8  read M2(m,e,f) = <ma||e,f>abab
        if ((myRank.eq.idaabb).or.(myRank.eq.idabba)) then
       call rea (n1abeta,m2lenght,wrk(possm20))
        else
       call reajalovy (n1abeta,m2lenght,wrk(possm20))
        end if
c
c     Status:
c     V1(ef,ij)   = Tau(ef,ij)bbbb
c     V2(e,f,i,j) = Tau(e,f,i,j)abab
c     V3(j,i,e,f) = T2o(e,f,i,j)abab
c     V4(i,j,ef)  = T2o(ef,ij)bbbb
c     M1(m,ef)  - <ma||ef>bbbb (free, but reserved map's)
c     M2(m,e,f) - <ma||ef>abab
c     H1(m,e,j) - <ma||ej>bbbb  (free, but reserved map's)
c     H2(m,e,j) - <ma||ej>abba = W3(m,e,a,j)abba
c     H3(m,e,j) - <ma||ej>abab = W3(m,e,a,j)aabb
c     free: M3,M4,H4
c
cpar
        if (myRank.eq.idaabb) then
c
c     ------- cont to T16
cT16.4T1n(a,i)bb <- + sum(m,e,f-aab) [ T2o(e,f,m,i)abab . <ma||ef>abab ]
c
cT16.4.1    M3(i) = V3(i,m,e,f) . M2(m,e,f)
       call mult (wrk,wrksize,
     & 4,3,1,3,
     & mapdv3,mapiv3,1,mapdm2,mapim2,syma,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT16.4.2    t1n(a,i)bb <- +1.0d0 . M3(i)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,1.0d0,
     & mapdm3,syma,mapdt14,mapit14,1,rc)
c
c     ------- cont to T2Ex
cT2E.6R2(a,m,i,j)baab <= - sum(e,f-ab) [ <ma||ef>abab . Tau(e,f,i,j)abab ]
cT2E.6T2n(a,b,i,j)abab <-   sum(m-a)    [ T1o(m,a)aa . R2(b,m,i,j)abab ]
c
cT2E.6.1    M4(m,i,j) <- M2(m,e,f) . V2(e,f,i,j)
       call mult (wrk,wrksize,
     & 3,4,3,2,
     & mapdm2,mapim2,syma,mapdv2,mapiv2,1,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cT2E.6.2    M3(a,i,j) <- T1o(a,m)aa . M4(m,i,j)
       call mult (wrk,wrksize,
     & 2,3,3,1,
     & mapdt11,mapit11,1,mapdm4,mapim4,syma,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cT2E.6.3    T2n(a,b,i,j)abab <- -1.0d0 M3(a,i,j)
       call add (wrk,wrksize,
     & 3,4,1,2,a,0,syma,1,-1.0d0,
     & mapdm3,ssm3,mapdt23,mapit23,1,rc)
cparend
       end if
c
c     ------- cont to W32
cW32.3WIII(m,e,b,j)aabb <- + sum(f-b) [  <mb||ef>abab . T1o(f,j)bb ]
cW32.4WIII(m,e,b,j)abba <- - sum(f-a) [  <mb||fe>abab . T1o(f,j)aa ]
c
cpar
        if (myRank.eq.idaabb) then
cW32.3.1    M4(m,e,j) = M2(m,e,f) . T1o(f,j)bb
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm2,mapim2,syma,mapdt12,mapit12,1,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cW32.3.2    H3(W3aabb)(m,e,j) <- 1.0d0 M4(m,e,j)
       call add (wrk,wrksize,
     & 3,3,0,0,0,0,1,1,1.0d0,
     & mapdm4,syma,mapdh3,mapih3,syma,rc)
c
cW3.3 write W3aabb(m,e,j) to  lunw3baab file if size is not zero
       if (h3lenght.gt.0) then
       call wri (lunw3aabb,h3lenght,wrk(possh30))
       end if
cparend
        end if
c
c
cpar
        if (myRank.eq.idabba) then
cW32.4.1    M4(m,e,f) <- M2(m,f,e)
       call map (wrk,wrksize,
     & 3,1,3,2,0,
     & mapdm2,mapim2,syma,mapdm4,mapim4,possm40,posst,rc)
c
cW32.4.2    M3(m,e,j) <- M4(m,e,f) . T1o(f,j)aa
       call mult (wrk,wrksize,
     & 3,2,3,1,
     & mapdm4,mapim4,syma,mapdt11,mapit11,1,mapdm3,mapim3,ssm3,
     & possm30,rc)
c
cW32.4.3    H2(W3abba)(m,e,j) <- -1.0d0 M3(m,e,j)
       call add (wrk,wrksize,
     & 3,3,0,0,0,0,1,1,-1.0d0,
     & mapdm3,syma,mapdh2,mapih2,syma,rc)
c
cW3.4 write W3abba(m,e,j) to  lunw3baba file if size is not zero
       if (h2lenght.gt.0) then
       call wri (lunw3abba,h2lenght,wrk(possh20))
       end if
cparend
        end if
c
c     ------- cont to FI3
cpar
        if (myRank.eq.idaabb) then
cF13.4FI(a,e)bb <-  sum(m,f-aa) [ <ma||fe>abab . T1o(f,m)aa ]
c
cF13.4.1    M3(e,f,m) <- M2(m,f,e)
       call map (wrk,wrksize,
     & 3,3,2,1,0,
     & mapdm2,mapim2,syma,mapdm3,mapim3,possm30,posst,rc)
c
cF13.4.2    M4(e)     <- M1(e,f,m) . T1o(f,m)aa
       call mult (wrk,wrksize,
     & 3,2,1,2,
     & mapdm3,mapim3,syma,mapdt11,mapit11,1,mapdm4,mapim4,ssm4,
     & possm40,rc)
c
cF13.4.3    F1(a,e)bb <- -1.0d0 M4(e)
       call add (wrk,wrksize,
     & 1,2,1,1,a,0,syma,1,1.0d0,
     & mapdm4,syma,mapdf12,mapif12,1,rc)
       end if
cparend
        end if
c
 5500   continue
c
 6000   continue
c
cC.9  close n1abeta,n2abeta
       call filemanager (3,n1abeta,rc)
       call filemanager (3,n2abeta,rc)
cparend
       end if
c
       return
       end
