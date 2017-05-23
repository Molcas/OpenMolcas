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
* Copyright (C) 2006, Pavel Neogrady                                   *
************************************************************************
c
c       this package contains:
c       sumoerab
c         multstack
c         grc0stack
c         extstack
c           extstackhlp1
c
c     -----------------------------
c
      subroutine sumoverab (wrk,wrksize,
     & lunt2o1, lunt2o2, lunt2o3,nabstack,possabstack,
     & niter)
c
c     this routine realize sumation over ab
c     and calculate following contributions:
c
c     T25
c     T2n(ab,ij)aaaa   <- sum(e>f-aa) [ <ab||ef>aaaa . Tau(ef,ij)aaaa ]
c     T2n(ab,ij)bbbb   <- sum(e>f-bb) [ <ab||ef>bbbb . Tau(ef,ij)bbbb ]
c     T2n(a,b,i,j)abab <- sum(e,f-ab) [ <ab||ef>abab . Tau(ef,ij)abab ]
c
c     T28
c      Q(ab,k,l)aaaa   <= sum(e-a)    [ <ab||ke>aaaa . T1o(e,l)aa ]
c     T2n(ab,ij)aaaa   <- Q(ab,i,j)aaaa - Q(ab,j,i)aaaa
c      Q(ab,k,l)bbbb   <= sum(e-b)    [ <ab||ke>bbbb . T1o(e,l)bb ]
c     T2n(ab,ij)bbbb   <- Q(ab,i,j)bbbb - Q(ab,j,i)bbbb
c     T2n(a,b,i,j)abab <- - sum(e-a)  [ <ab||je>abba . T1o(e,i)aa ]
c        <- + sum(e-b)  [ <ab||ie>abab . T1o(e,j)bb ]
c
c     General Status:

c     1) Tau ampitudes are stored and properly mapped in V1-V3
c     V1(ij,ef)   = Tau(ef,ij)aaaa
c     V2(ij,ef)   = Tau(ef,ij)bbbb
c     V3(i,j,e,f) = Tau(e,f,i,j)abab
c
c     2) Integrals are stored in file nab
c
c     syma=1,nsym
c       symb=1,syma
c       ! one record mapdn,mapin - in all cases
c         a=1,nvb(syma)
c         if (syma.eq.symb) then
c         limb=a
c         else
c         limb=nvb(symb)
c         end if
c           b=1,limb
c             ! one record N _a_b(p,q) - if any
c           end b
c         end a
c       end sumb
c     end syma
c
c     N.B. nab musi byt dakde deklarovane (daky common)
c
c     N.B. II. in paralell.fh there are:
c     nprocab - number of nodeschoosed for 'sumoverab' part
c     idab    - array of id's of these nodes
c
        implicit none
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
#include "paralell.fh"
      INTEGER lunt2o1,lunt2o2,lunt2o3
      integer nabstack,possabstack
      integer niter
c
c     help variables
c
      INTEGER key,aeqb,yes
      INTEGER a,b,syma,symb,ssh3,ssn,ssh1
      INTEGER nlenght,posst,rc
      INTEGER lunab,nsa,nsb,limb
      integer nab,nadd
      integer idtot(1:maxproc)
c
      integer iab,nabnow,left,bstart,bstop,todo,possab,bb
      integer nhelp
c
      todo=0
      possab=possabstack
c
c
c     I.paralell
c
c     I.par.1  - escape, if this node is not reserved for sumoverab
      yes=0
      do a=1,nprocab
        if (idab(a).eq.myRank) then
c       if ((idab(a).eq.myRank).and.(ideffab(a).gt.0.0d0)) then
        yes=a
        end if
      end do
c
      if (yes.eq.0) then
      return
      end if
c
c      I.par.2 - cal overal number of ab records
      nab=0
      DO  syma=1,nsym
        DO  symb=1,syma
        IF (syma.eq.symb) THEN
        nab=nab+(nvb(syma)*(nvb(syma)+1))/2
        ELSE
        nab=nab+nvb(syma)*nvb(symb)
        END IF
      end do
      end do
c
c       I.par.3 - distribute and calc add parameter
      call sumabdistt (nab,idtot)
        nab=idtot(yes)
c        escape, if there is nothing to do on this node
        if (nab.eq.0) return
        nadd=0
        do b=1,yes-1
        nadd=nadd+idtot(b)
        end do
c
c      I.1  open nab file
c
      lunab=10
      CALL filemanager (4, lunab, rc)
c
c      I.2  prepair
c     V1(ij,ef)   = Tau(ef,ij)aaaa
c     V2(ij,ef)   = Tau(ef,ij)bbbb
c     V3(i,j,e,f) = Tau(e,f,i,j)abab
c
      CALL filemanager (2, lunt2o1, rc)
      CALL getmediate (wrk,wrksize,
     &lunt2o1, possv40, mapdv4, mapiv4, rc)
      CALL mktau (wrk,wrksize,
     &mapdv4, mapiv4, mapdt11, mapit11, mapdt12, mapit12,
     11.0d0, rc)
      CALL map (wrk,wrksize,
     &4, 3, 4, 1, 2, mapdv4, mapiv4, 1, mapdv1, mapiv1,
     1possv10, posst, rc)
c
      CALL filemanager (2, lunt2o2, rc)
      CALL getmediate (wrk,wrksize,
     &lunt2o2, possv40, mapdv4, mapiv4, rc)
      CALL mktau (wrk,wrksize,
     &mapdv4, mapiv4, mapdt11, mapit11, mapdt12, mapit12,
     11.0d0, rc)
      CALL map (wrk,wrksize,
     &4, 3, 4, 1, 2, mapdv4, mapiv4, 1, mapdv2, mapiv2,
     1possv20, posst, rc)
c
      CALL filemanager (2, lunt2o3, rc)
      CALL getmediate (wrk,wrksize,
     &lunt2o3, possv40, mapdv4, mapiv4, rc)
      CALL mktau (wrk,wrksize,
     &mapdv4, mapiv4, mapdt11, mapit11, mapdt12, mapit12,
     11.0d0, rc)
      CALL map (wrk,wrksize,
     &4, 3, 4, 1, 2, mapdv4, mapiv4, 1, mapdv3, mapiv3,
     1possv30, posst, rc)
c
c      II   sum over ab
      DO 20 syma=1,nsym
        DO 20 symb=1,syma
        if (fullprint.ge.2) then
        write (6,*) ' SymA, SymB ',syma,symb
        end if
c
c      II.1   read mapdn,mapin
        CALL getmap (lunab, possn0, nlenght, mapdn, mapin, rc)
c
c      II.3     skip sum over a,b if N is empty
        IF (nlenght.eq.0) GO TO 20
c
c      II.4   skip sum over a,b if # of ab is zero
        IF ((nvb(syma)*nvb(symb)).eq.0) GO TO 20
c
c      II.5   def symmetry of N
        ssn=mmul(syma,symb)
c
c      II.5.1   def mapd and mapi for R1 _a,_b(j,e)aaaa = <ab||je>  (pos
        CALL grc0 (2, 0, 1, 3, 0, 0, ssn, possm10, posst, mapdr1,
     1   mapir1)
c      II.5.2   def mapd and mapi for R2 _a,_b(j,e)bbbb = <ab||je>  (pos
        CALL grc0 (2, 0, 2, 4, 0, 0, ssn, possm20, posst, mapdr2,
     1   mapir2)
c      II.5.3   def mapd and mapi for R3 _a,_b(j,e)abba = <ab||je>  (pos
        CALL grc0 (2, 0, 2, 3, 0, 0, ssn, possm30, posst, mapdr3,
     1   mapir3)
c      II.5.4   def mapd and mapi for R4 _b,_a(j,e)abba = <ba||je>  (pos
        CALL grc0 (2, 0, 2, 3, 0, 0, ssn, possm40, posst, mapdr4,
     1   mapir4)
c      II.5.5   def mapd and mapi for R5 _a,_b(j,e)abab = <ab||je>  (pos
        CALL grc0 (2, 0, 1, 4, 0, 0, ssn, possh10, posst, mapdr5,
     1   mapir5)
c      II.5.6   def mapd and mapi for R6 _b,_a(j,e)abab = <ba||je>  (pos
        CALL grc0 (2, 0, 1, 4, 0, 0, ssn, possh20, posst, mapdr6,
     1   mapir6)
c
c      II.6.1   def mapd and mapi for M1 _a,_b(ef)aaaa = <ab||ef>
        CALL grc0 (2, 1, 3, 3, 0, 0, ssn, possm10, posst, mapdm1,
     1   mapim1)
c      II.6.2   def mapd and mapi for M2 _a,_b(ef)bbbb = <ab||ef>
        CALL grc0 (2, 1, 4, 4, 0, 0, ssn, possm20, posst, mapdm2,
     1   mapim2)
c      II.6.3   def mapd and mapi for M3 _a,_b(e,f)abab = <ab||ef>
        CALL grc0 (2, 0, 3, 4, 0, 0, ssn, possm30, posst, mapdm3,
     1   mapim3)
c      II.6.4   def mapd and mapi for M4 _b,_a(e,f)abab = <ba||fe>
        CALL grc0 (2, 0, 3, 4, 0, 0, ssn, possm40, posst, mapdm4,
     1   mapim4)
c
c      II.7   def # of S orbitals in syma and symb
        nsa=nvb(syma)-nva(syma)
        nsb=nvb(symb)-nva(symb)
c
        DO 10 a=1,nvb(syma)
c
        if (fullprint.ge.3) then
        write (6,*) ' A', a
        end if
c
c      II.7     def limitations for b
          IF (syma.eq.symb) THEN
            limb=a
          ELSE
            limb=nvb(symb)
          END IF
c
c@@
cpar    for nodes, selected for sumab process:
c       1) only blank read (instead of _a_b(e,f) form lunab file for records,
c       pre-standing to the portion, that is going to be done on this node
c       2) prepair setup (def bstart,limb), if needed
c
        if (nadd.ge.limb) then
c       all cycle over b need to be skipped
          do b=1,limb
c           blank read and nothing more
            call reajalovy (lunab, nlenght, wrk(possn0))
            nadd=nadd-1
          end do
c       no setup needed
          goto 10
c
        else if (nadd.gt.0) then
c       setup first
          bstart=nadd
          if (nab.lt.(limb-bstart)) then
            limb=bstart+nab
          end if
c       part of the cycle over b need to be skipped
          do b=1,nadd
c           blank read and nothing more
            call reajalovy (lunab, nlenght, wrk(possn0))
          end do
          nadd=0
c
        else
c       nabadd musi byt 0, inak nieje dobre
c       no part of the cycle over b need to be skipped
c       only setup
          bstart=0
          if (nab.lt.limb) then
            limb=nab
          end if
        end if
c
c
c*      define conditions for reading to stack
c       bstart - lower limit of b-sumation
c       bstop  - upper limit of b-sumation
c       nabnow - actual size af stack (if=1, no stacking)
c
c*.1    initialization
cStare  bstart=0
c
c*.1    set unstacking conditions
4       bstart=bstart+1
        bstop=bstart
        nabnow=1
c
c*.2    without stacking if: 1)a=b; 2)a or b \in S; 3) nabnow.le.3
c       (N.B. it was found inconvenient to execute stacking if
c             size of stack is less then 4,or even more)
c
c       ad 1)
        if ((bstart.eq.a).and.(syma.eq.symb)) goto 6
c       ad 2)
        if ((bstart.le.nsb).or.(a.le.nsa)) goto 6
c       ad 3)
        left=limb-bstart+1
        if (syma.eq.symb) left=left-1
        if (left.le.3) goto 6
c
c*.3    set stacking conditions for real stacking
        if (left.le.nabstack) then
          nabnow=left
        else
          nabnow=nabstack
        end if
        bstop=bstart+nabnow-1
        todo=1
        possab=possabstack
c
c*.     reading to N or to stack
6       if (nabnow.eq.1) then
          CALL rea (lunab, nlenght, wrk(possn0))
        else
          do iab=1,nabnow
          CALL rea (lunab, nlenght, wrk(possabstack+(iab-1)*nlenght))
          if ((possabstack+iab*nlenght).ge.possm10) then
          end if
          end do
c         qq=wrk(possabstack)
        end if
c
c
          do 8 b=bstart,bstop
c
c      III.3       read N
cStare    CALL rea (lunab, nlenght, wrk(possn0))
c
c      III.4       def aeqb
          IF ((syma.eq.symb).and.(a.eq.b)) THEN
            aeqb=1
          ELSE
            aeqb=0
          END IF
c
c      III.5       def key
          IF ((a.gt.nsa).and.(b.gt.nsb)) THEN
            key=4
          ELSE IF ((a.gt.nsa).and.(b.le.nsb)) THEN
            key=2
          ELSE IF ((a.le.nsa).and.(b.gt.nsb)) THEN
            key=3
          ELSE
            key=1
          END IF
c
c      def yes for T25.1
          IF (syma.eq.symb) THEN
            IF (a.eq.b) THEN
              yes=0
            ELSE IF (key.eq.4) THEN
              yes=1
            ELSE
              yes=0
            END IF
          ELSE
            IF (key.eq.4) THEN
              yes=1
            ELSE
              yes=0
            END IF
          END IF
c
c
c           T25
c           T2n(ab,ij)aaaa   <- sum(e>f-aa) [ <ab||ef>aaaa . Tau(ef,ij)aaaa ]
c           T2n(ab,ij)bbbb   <- sum(e>f-bb) [ <ab||ef>bbbb . Tau(ef,ij)bbbb ]
c           T2n(a,b,i,j)abab <- sum(e,f-ab) [ <ab||ef>abab . Tau(ef,ij)abab ]
c
        if (nabnow.eq.1) then
c       AB stacking is not active
c
c      III.6       get
c           -> M1 _a,_b(ef)aaaa  = <ab||ef>
c           -> M2 _a,_b(ef)bbbb  = <ab||ef>
c           -> M3 _a,_b(e,f)abab = <ab||ef>
c           -> M4 _b,_a(e,f)abab = <ba||fe>
c           free: H1-H4 (V4 reserved for stacking)
          CALL unpackab1 (wrk,wrksize,
     &     mapdn, mapin, mapdm1, mapim1, mapdm2, mapim2,
     1     mapdm3, mapim3, mapdm4, mapim4, ssn, key, aeqb)
c
c
          IF (yes.eq.1) THEN
c      T25.1.1      H1 (ij) = V1(ij,ef) . M1 (ef)
            CALL mult (wrk,wrksize,
     &       4, 2, 2, 2, mapdv1, mapiv1, 1, mapdm1, mapim1,
     1       ssn, mapdh1, mapih1,ssh1,possh10, rc)
c      T25.1.2      T2n(ab,ij)aaaa <- 1.0d0 . H1(ij)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a-nsa, b-nsb, syma, symb, 1.0d0,
     1       mapdh1, ssn, mapdt21, mapit21, 1, rc)
          END IF
c
          IF (aeqb.eq.0) THEN
c      T25.2.1      H1 (ij) = V2(ij,ef) . M2 (ef)
            CALL mult (wrk,wrksize,
     &       4, 2, 2, 2, mapdv2, mapiv2, 1, mapdm2, mapim2,
     1       ssn, mapdh1, mapih1,ssh1,possh10, rc)
c      T25.2.2      T2n(ab,ij)bbbb <- 1.0d0 . H1(ij)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a, b, syma, symb, 1.0d0, mapdh1, ssn,
     1       mapdt22, mapit22, 1, rc)
          END IF
c
          IF ((key.eq.2).or.(key.eq.4)) THEN
c      T25.3.1      H1 (i,j) = V3(i,j,e,f) . M3 (e,f)
            CALL mult (wrk,wrksize,
     &       4, 2, 2, 2, mapdv3, mapiv3, 1, mapdm3, mapim3,
     1       ssn, mapdh1, mapih1,ssh1,possh10, rc)
c      T25.3.2      T2n(a,b,i,j)abab <- 1.0d0 . H1(i,j)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a-nsa, b, syma, symb, 1.0d0, mapdh1,
     1       ssn, mapdt23, mapit23, 1, rc)
          END IF
c
          IF (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) THEN
c      T25.3.3      H1 (i,j) = V3(i,j,e,f) . M4 (e,f)
            CALL mult (wrk,wrksize,
     &       4, 2, 2, 2, mapdv3, mapiv3, 1, mapdm4, mapim4,
     1       ssn, mapdh1, mapih1,ssh1,possh10, rc)
c
c      T25.3.4      T2n(a,b,i,j)abab <- 1.0d0 . H1(i,j)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, b-nsb, a, symb, syma, 1.0d0, mapdh1,
     1       ssn, mapdt23, mapit23, 1, rc)
          END IF
c
c
        else
c       AB stacking is active
c
          if (todo.eq.1) then
c         once it need to be done
          todo=0
c
c      III.6       get
c           -> M1 _a,_b(ef)aaaa  = <ab||ef>
c           -> M2 _a,_b(ef)bbbb  = <ab||ef>
c           -> M3 _a,_b(e,f)abab = <ab||ef>
c           -> M4 _b,_a(e,f)abab = <ba||fe>
c           free: H1-H4 (V4 reserved for stacking)
c
c
c       def mapd and mapi for V4 _a,(ef,Bp)
       call grc0stack (nabnow,1,3,3,3,0,ssn,possv40,posst,mapdv4,mapiv4)
c           -> V4 _a,_b(ef)aaaa  = <ab||ef>
       call unpackab3 (wrk,wrksize,
     &      mapdn,mapin,mapdv4,mapiv4,ssn,nabnow,possabstack,nlenght,1)
c
c      T25.1.1      H2 (ij,Bp) = V1(ij,ef) . V4 (ef,Bp)
        call multstack (wrk,wrksize,
     &                  mapdv1,mapdv4,mapdh2,mapiv1,mapiv4,mapih2,
     &                  1,ssn,possh20,nabnow)

c       def mapd and mapi for H1 _a,_b(ij)aa
        call grc0 (2,1,1,1,0,0,ssn,possh10,posst,mapdh1,mapih1)
        do bb=bstart,bstop
c      T25.1.2      ext H1(ij) <- H2(ij,_Bb) (only data transfer)
        call extstack (wrk,wrksize,
     &                 mapdh1,mapdh2,(bb-bstart+1),nabnow)
c      T25.1.3      T2n(ab,ij)aaaa <- 1.0d0 . H1(ij)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a-nsa, bb-nsb, syma, symb, 1.0d0,
     1       mapdh1, ssn, mapdt21, mapit21, 1, rc)
        end do
c
c
c
c       def mapd and mapi for V4 _a,(ef,Bp)
       call grc0stack (nabnow,1,4,4,4,0,ssn,possv40,posst,mapdv4,mapiv4)
c           -> V4 _a,_b(ef)bbbb  = <ab||ef>
       call unpackab3 (wrk,wrksize,
     &      mapdn,mapin,mapdv4,mapiv4,ssn,nabnow,possabstack,nlenght,2)
c
c      T25.2.1      H2 (ij,Bp) = V2(ij,ef) . V4 (ef,Bp)
        call multstack (wrk,wrksize,
     &                  mapdv2,mapdv4,mapdh2,mapiv2,mapiv4,mapih2,
     &                  1,ssn,possh20,nabnow)
c
c       def mapd and mapi for H1 _a,_b(ij)bb
        call grc0 (2,1,2,2,0,0,ssn,possh10,posst,mapdh1,mapih1)
        do bb=bstart,bstop
c      T25.2.2      ext H1(ij) <- H2(ij,Bb)
        call extstack (wrk,wrksize,
     &                 mapdh1,mapdh2,(bb-bstart+1),nabnow)
c      T25.2.3      T2n(ab,ij)bbbb <- 1.0d0 . H1(ij)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a, bb, syma, symb, 1.0d0, mapdh1, ssn,
     1       mapdt22, mapit22, 1, rc)
        end do
c
c
c
c       def mapd and mapi for V4 _a,(ef,Bp)
       call grc0stack (nabnow,0,3,4,4,0,ssn,possv40,posst,mapdv4,mapiv4)
c           -> V4 _a,_b(e,f)abab = <ab||ef>
       call unpackab3 (wrk,wrksize,
     &      mapdn,mapin,mapdv4,mapiv4,ssn,nabnow,possabstack,nlenght,3)
c
c      T25.3.1      H2 (i,j,Bp) = V3(i,j,e,f) . V4 (e,f,Bp)
        call multstack (wrk,wrksize,
     &                  mapdv3,mapdv4,mapdh2,mapiv3,mapiv4,mapih2,
     &                  1,ssn,possh20,nabnow)
c
c       def mapd and mapi for H1 _a,_b(ij)ab
        call grc0 (2,0,1,2,0,0,ssn,possh10,posst,mapdh1,mapih1)
        do bb=bstart,bstop
c      T25.3.2      ext H1(i,j) <- H2(ij,Bb)
        call extstack (wrk,wrksize,
     &                 mapdh1,mapdh2,(bb-bstart+1),nabnow)
c      T25.3.3      T2n(a,b,i,j)abab <- 1.0d0 . H1(i,j)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a-nsa, bb, syma, symb, 1.0d0, mapdh1,
     1       ssn, mapdt23, mapit23, 1, rc)
        end do
c
c
c
c       def mapd and mapi for V4 _a,(ef,Bp)
       call grc0stack (nabnow,0,3,4,3,0,ssn,possv40,posst,mapdv4,mapiv4)
c           -> V4 _b,_a(e,f)abab = <ba||fe>
       call unpackab3 (wrk,wrksize,
     &      mapdn,mapin,mapdv4,mapiv4,ssn,nabnow,possabstack,nlenght,4)
c
c      T25.3.4      H2 (i,j,Bp) = V3(i,j,e,f) . V4 (e,f,Bp)
        call multstack (wrk,wrksize,
     &                  mapdv3,mapdv4,mapdh2,mapiv3,mapiv4,mapih2,
     &                  1,ssn,possh20,nabnow)
c
c       def mapd and mapi for H1 _a,_b(ij)ab
        call grc0 (2,0,1,2,0,0,ssn,possh10,posst,mapdh1,mapih1)
        do bb=bstart,bstop
c      T25.3.5      ext H1(i,j) <- H2(i,j,Bb)
        call extstack (wrk,wrksize,
     &                 mapdh1,mapdh2,(bb-bstart+1),nabnow)
c      T25.3.6      T2n(a,b,i,j)abab <- 1.0d0 . H1(i,j)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, bb-nsb, a, symb, syma, 1.0d0, mapdh1,
     1       ssn, mapdt23, mapit23, 1, rc)
        end do
c
          end if
        end if
c
c     T28
c      Q(ab,k,l)aaaa   <= sum(e-a)    [ <ab||ke>aaaa . T1o(e,l)aa ]
c     T2n(ab,ij)aaaa   <- Q(ab,i,j)aaaa - Q(ab,j,i)aaaa
c      Q(ab,k,l)bbbb   <= sum(e-b)    [ <ab||ke>bbbb . T1o(e,l)bb ]
c     T2n(ab,ij)bbbb   <- Q(ab,i,j)bbbb - Q(ab,j,i)bbbb
c     T2n(a,b,i,j)abab <- - sum(e-a)  [ <ab||je>abba . T1o(e,i)aa ]
c        <- + sum(e-b)  [ <ab||ie>abab . T1o(e,j)bb ]
c
c
c       put appropriate data from AB_stack to N mediate, if needed
        if (nabnow.gt.1) then
c       transfer data
          do nhelp=0,nlenght-1
          wrk(possn0+nhelp)=wrk(possab+nhelp)
          end do
c       upgrade address
          possab=possab+nlenght
        end if
c
c
c      III.7       get
c           -> M1 _a,_b(j,e)aaaa = <ab||je> (map is R1)
c           -> M2 _a,_b(j,e)bbbb = <ab||je> (map is R2)
c           -> M3 _a,_b(j,e)abba = <ab||je> (map is R3)
c           -> M4 _b,_a(j,e)abba = <ba||je> (map is R4)
c           -> H1 _a,_b(j,e)abab = <ab||je> (map is R5)
c           -> H2 _b,_a(j,e)abab = <ba||je> (map is R6)
c           free: H3,H4 (V4 reserved for stacking)
c
       CALL unpackab2 (wrk,wrksize,
     & mapdn,mapin,mapdr1,mapir1,mapdr2,mapir2,mapdr3,mapir3,
     & mapdr4,mapir4,mapdr5,mapir5,mapdr6,mapir6,ssn,key,aeqb)
c
c
c     T28
c      Q(ab,k,l)aaaa   <= sum(e-a)    [ <ab||ke>aaaa . T1o(e,l)aa ]
c     T2n(ab,ij)aaaa   <- Q(ab,i,j)aaaa - Q(ab,j,i)aaaa
c      Q(ab,k,l)bbbb   <= sum(e-b)    [ <ab||ke>bbbb . T1o(e,l)bb ]
c     T2n(ab,ij)bbbb   <- Q(ab,i,j)bbbb - Q(ab,j,i)bbbb
c     T2n(a,b,i,j)abab <- - sum(e-a)  [ <ab||je>abba . T1o(e,i)aa ]
c        <- + sum(e-b)  [ <ab||ie>abab . T1o(e,j)bb ]
c
c
c      def yes for T28.1
          IF (syma.eq.symb) THEN
            IF (a.eq.b) THEN
              yes=0
            ELSE IF (key.eq.4) THEN
              yes=1
            ELSE
              yes=0
            END IF
          ELSE
            IF (key.eq.4) THEN
              yes=1
            ELSE
              yes=0
            END IF
          END IF
c
c
          IF (yes.eq.1) THEN
c      T28.1.1      H3(k,l) = M1(k,e) . T1o(e,l)aa
            CALL mult (wrk,wrksize,
     &       2, 2, 2, 1, mapdr1, mapir1, ssn, mapdt11,
     1       mapit11, 1, mapdh3, mapih3,ssh3,possh30, rc)
c      T28.1.2      H4(kl)  = H3(k,l)-H3(l,k)
            CALL fack (wrk,wrksize,
     &       2, 1, mapdh3, ssh3, mapih3, mapdh4, mapih4,
     1       possh40, rc)
c      T28.1.3      T2n(ab,ij)aaaa <- 1.0d0 H4(ij)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a-nsa, b-nsb, syma, symb, 1.0d0,
     1       mapdh4, ssn, mapdt21, mapit21, 1, rc)
          END IF
c
          IF (aeqb.eq.0) THEN
c      T28.2.1      H3(k,l) = M2(k,e) . T1o(e,l)bb
            CALL mult (wrk,wrksize,
     &       2, 2, 2, 1, mapdr2, mapir2, ssn, mapdt12,
     1       mapit12, 1, mapdh3, mapih3,ssh3,possh30, rc)
c      T28.2.2      H4(kl)  = H3(k,l)-H3(l,k)
            CALL fack (wrk,wrksize,
     &       2, 1, mapdh3, ssh3, mapih3, mapdh4, mapih4,
     1       possh40, rc)
c      T28.2.3      T2n(ab,ij)bbbb <- 1.0d0 H4(ij)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a, b, syma, symb, 1.0d0, mapdh4, ssn,
     1       mapdt22, mapit22, 1, rc)
          END IF
c
          IF ((key.eq.2).or.(key.eq.4)) THEN
c      T28.3.1      H3(j,i)  = M3(j,e) . T1o(e,i)aa
            CALL mult (wrk,wrksize,
     &       2, 2, 2, 1, mapdr3, mapir3, ssn, mapdt11,
     1       mapit11, 1, mapdh3, mapih3,ssh3,possh30, rc)
c      T28.3.2      H4(i,j)  = H3(j,i)
            CALL map (wrk,wrksize,
     &       2, 2, 1, 0, 0, mapdh3, mapih3, ssn, mapdh4,
     1       mapih4, possh40, posst, rc)
c      T28.3.3      T2n(a,b,i,j)abab <- -1.0d0 H4(i,j)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a-nsa, b, syma, symb, -1.0d0, mapdh4,
     1       ssn, mapdt23, mapit23, 1, rc)
          END IF
c
          IF (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) THEN
c      T28.3.4      H3(j,i)  = M4(j,e) . T1o(e,i)aa
            CALL mult (wrk,wrksize,
     &       2, 2, 2, 1, mapdr4, mapir4, ssn, mapdt11,
     1       mapit11, 1, mapdh3, mapih3,ssh3,possh30, rc)
c      T28.3.5      H4(i,j)  = H3(j,i)
            CALL map (wrk,wrksize,
     &       2, 2, 1, 0, 0, mapdh3, mapih3, ssn, mapdh4,
     1       mapih4, possh40, posst, rc)
c      T28.3.6      T2n(a,b,i,j)abab <- -1.0d0 H4(i,j)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, b-nsb, a, symb, syma, -1.0d0, mapdh4,
     1       ssn, mapdt23, mapit23, 1, rc)
          END IF
c
          IF ((key.eq.2).or.(key.eq.4)) THEN
c      T28.3.7      H3(i,j)  = H1(i,e) . T1o(e,j)bb
            CALL mult (wrk,wrksize,
     &       2, 2, 2, 1, mapdr5, mapir5, ssn, mapdt12,
     1       mapit12, 1, mapdh3, mapih3,ssh3,possh30, rc)
c      T28.3.8      T2n(a,b,i,j)abab <- 1.0d0 H3(i,j)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, a-nsa, b, syma, symb, 1.0d0, mapdh3,
     1       ssn, mapdt23, mapit23, 1, rc)
          END IF
c
          IF (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) THEN
c      T28.3.9      H3(j,i)  = H2(i,e) . T1o(e,j)bb
            CALL mult (wrk,wrksize,
     &       2, 2, 2, 1, mapdr6, mapir6, ssn, mapdt12,
     1       mapit12, 1, mapdh3, mapih3,ssh3,possh30, rc)
c      T28.3.10     T2n(a,b,i,j)abab <- 1.0d0 H3(i,j)
            CALL add (wrk,wrksize,
     &       2, 4, 2, 5, b-nsb, a, symb, syma, 1.0d0, mapdh3,
     1       ssn, mapdt23, mapit23, 1, rc)
          END IF
c
c
  8     CONTINUE
c
cpar
        nab=nab-(bstop-bstart+1)
        if (nab.eq.0) then
        goto 30
        end if
c
        if (bstop.lt.limb) then
          bstart=bstop
          goto 4
        end if
 10     CONTINUE
c
 20   CONTINUE
c
c      IV   close lunab
30    CALL filemanager (3, lunab, rc)
c
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(niter)
      END
c
c     -----------------------------
c
       subroutine multstack (wrk,wrksize,
     &                       mapda,mapdb,mapdc,mapia,mapib,mapic,
     &                       ssa,ssb,possc0,bsize)
c
c        This is a special routine for multiplying of the:
c        C(ij,Bp) = A(ij,cd) . B(cd,Bp)
c        where Bp is a limited (partial) sumation over
c        b index (#b - bsize), namely those, hich are stacked
c
c        This routine is used only in stacking in sumoverab
c        process and is a modification of grc42y routine.
c        Type of index Bp is registered as for standard b, but
c        all lengths of blocks are calculated with
c        bsize, instead of dimm(typb,symb) and the symmetry
c       of b is ignored, B and C are treated as 2index
c
c        P.N. 17.02.06
c
c
#include "ccsd1.fh"
#include "wrk.fh"
c
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer mapdc(0:512,1:6)
c
       integer mapia(1:8,1:8,1:8)
       integer mapib(1:8,1:8,1:8)
       integer mapic(1:8,1:8,1:8)
c
       integer mvec(1:4096,1:7)
       integer possc0,bsize
       integer ssa,ssb
c
c     help variables
c
       integer nhelp1,nhelp2,nhelp3,nhelp4
       integer nhelp21,nhelp22,nhelp41,nhelp42
       integer ntest1,ntest2
       integer sa1,sa2,sa3,sa4,sb1,sb2,sa34,sa134
CLD    integer sa1,sa2,sa3,sa4,sb1,sb2,sb3,sa12,sa34,sa134,sb12
CLD    integer nsyma2
       integer ia,ib,ix,iy
CLD    integer ia,ib,ic,ix,iy
       integer possct
c
c1*
c
c     sctructure A(pq,rs)*B(rs,t)=C(pq,t)
c     sctructure A(pq,rs)*B(rs)=YC(pq)
c
c1.1  define limitations -  p>q,r,s must be tested - ntest1
c     p,q,r>s must be tested - ntest2
c
       if ((mapda(0,6).eq.1).or.(mapda(0,6).eq.4)) then
       ntest1=1
       else
       ntest1=0
       end if
c
       if ((mapda(0,6).eq.3).or.(mapda(0,6).eq.4)) then
       ntest2=1
       else
       ntest2=0
       end if
c
c1.0  prepare mapdc,mapic
c
       call grc0stack (bsize,ntest1,mapda(0,1),mapda(0,2),mapdb(0,3),
     &                 0,mmul(ssa,ssb),possc0,possct,mapdc,mapic)

c
c
c1.2  def symm states and test the limitations
c
       ix=1
       do 100 sb1=1,nsym
       sa3=sb1
c
       sb2=mmul(ssb,sb1)
       sa4=sb2
       sa34=mmul(sa3,sa4)
       if ((ntest2.eq.1).and.(sb1.lt.sb2)) then
c     Meggie out
       goto 100
       end if
c
       do 50 sa1=1,nsym
       sa134=mmul(sa1,sa34)
c
       sa2=mmul(ssa,sa134)
       if ((ntest1.eq.1).and.(sa1.lt.sa2)) then
c     Meggie out
       goto 50
       end if
c
c1.3  def mvec,mapdc and mapdi
c
       ia=mapia(sa1,sa2,sa3)
       ib=mapib(sb1,1,1)
       iy=mapic(sa1,1,1)
c
c     yes/no
       if ((mapda(ia,2).gt.0).and.(mapdb(ib,2).gt.0)) then
       nhelp1=1
       else
       goto 50
       end if
c
c     rowA
       nhelp21=dimm(mapda(0,1),sa1)
       nhelp22=dimm(mapda(0,2),sa2)
       if ((ntest1.eq.1).and.(sa1.eq.sa2)) then
       nhelp2=nhelp21*(nhelp21-1)/2
       else
       nhelp2=nhelp21*nhelp22
       end if
c
c     sum
       nhelp41=dimm(mapda(0,3),sa3)
       nhelp42=dimm(mapda(0,4),sa4)
       if ((ntest2.eq.1).and.(sa3.eq.sa4)) then
       nhelp4=nhelp41*(nhelp41-1)/2
       else
       nhelp4=nhelp41*nhelp42
       end if
c
c     colBp
       nhelp3=bsize
c
       mvec(ix,1)=nhelp1
       mvec(ix,2)=mapda(ia,1)
       mvec(ix,3)=mapdb(ib,1)
       mvec(ix,4)=mapdc(iy,1)
       mvec(ix,5)=nhelp2
       mvec(ix,6)=nhelp4
       mvec(ix,7)=nhelp3
c
       ix=ix+1
c
 50     continue
 100    continue
       ix=ix-1
c
c
c*        multiplying
c
        call multc0 (wrk,wrksize,
     &               mvec,ix,mapdc,1)
c
       return
       end
c
c     -----------------------------
c
       subroutine grc0stack (bsize,typ,typp,typq,typr,typs,stot,
     & poss0,posst,mapd,mapi)
c
c             This routine defines mapd and mapi for specific
c        3 index intermediat A(pq,Bp), needed when stacking
c        (About Bp, see notes in multstack)
c        This routine is a modification of grc0 routine
c
c        P.N. 17.02.06
c     !N.B. (this routine cannot run with +OP2)
c
       integer bsize,typ,typp,typq,typr,typs,stot,poss0,posst
c
#include "ccsd1.fh"
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer sp,sq
CLD    integer sp,sq,sr,ss,spq,spqr
CLD    integer nsymq,nsymr
       integer poss,i,nhelp1,nhelp2,nhelp3
CLD    integer poss,i,nhelp1,nhelp2,nhelp3,nhelp4

c     To get rid of compiler warning
      poss=0
      i=0
c
c     vanishing mapi files
c
       do nhelp1=1,nsym
       do nhelp2=1,nsym
       do nhelp3=1,nsym
       mapi(nhelp3,nhelp2,nhelp1)=0
       end do
       end do
       end do
c
c     matrix A(p,q) or specifilally A(i,j,Bp)
c
       i=1
       poss=poss0
c
       do 100 sp=1,nsym
c
       sq=mmul(stot,sp)
       if ((typ.eq.1).and.(sp.lt.sq)) then
c     Meggie out
       goto 100
       end if
c
       nhelp1=dimm(typp,sp)
       nhelp2=dimm(typq,sq)
c
c     def mapi
       mapi(sp,1,1)=i
c
c     def possition
       mapd(i,1)=poss
c
c     def lenght
       if ((typ.eq.1).and.(sp.eq.sq)) then
       mapd(i,2)=bsize*nhelp1*(nhelp1-1)/2
       else
       mapd(i,2)=bsize*nhelp1*nhelp2
       end if
c
c     def sym p,q
       mapd(i,3)=sp
       mapd(i,4)=sq
       mapd(i,5)=0
       mapd(i,6)=0
c
       poss=poss+mapd(i,2)
       i=i+1
c
 100    continue
c
c
       posst=poss
c
c     definition of other coll
c
       mapd(0,1)=typp
       mapd(0,2)=typq
       mapd(0,3)=typr
       mapd(0,4)=typs
       mapd(0,5)=i-1
       mapd(0,6)=typ
c
       return
       end
c
c     -----------------------------
c
        subroutine extstack (wrk,wrksize,
     &                 mapda,mapdb,b,dimb)
c
c       This routine do:
c       A(ij) <- B(ij,_b) for given b
c
c       A special routine used only in sumoverab for stacking
c       case.
c
c       Yet it is assumed,that blocks in A and B are in the same
c       order. To je pomerne odflaknuty predpoklad, a moze to
c       byt bugous
c
#include "wrk.fh"
c
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer b,dimb
c
c     help variables
c
       integer ii,dimij,possa,possb
c
c
        do ii=1,mapda(0,5)
          dimij=mapda(ii,2)
          possa=mapda(ii,1)
          possb=mapdb(ii,1)
          call extstackhlp1 (wrk(possa),wrk(possb),dimij,dimb,b)
        end do
c
        return
        end
c
c       -----------
c
        subroutine extstackhlp1 (a,b,dimij,dimb,bb)
        integer dimij,dimb,bb,ij
        real*8 a(1:dimij),b(1:dimij,1:dimb)
        do ij=1,dimij
        a(ij)=b(ij,bb)
        end do
        return
        end
c
c     -----------------------------
c
