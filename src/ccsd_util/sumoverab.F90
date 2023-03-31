!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2006, Pavel Neogrady                                   *
!***********************************************************************
      subroutine sumoverab (wrk,wrksize,                                &
     & lunt2o1, lunt2o2, lunt2o3,nabstack,possabstack,                  &
     & niter)
!
!     this routine realize sumation over ab
!     and calculate following contributions:
!
!     T25
!     T2n(ab,ij)aaaa   <- sum(e>f-aa) [ <ab||ef>aaaa . Tau(ef,ij)aaaa ]
!     T2n(ab,ij)bbbb   <- sum(e>f-bb) [ <ab||ef>bbbb . Tau(ef,ij)bbbb ]
!     T2n(a,b,i,j)abab <- sum(e,f-ab) [ <ab||ef>abab . Tau(ef,ij)abab ]
!
!     T28
!      Q(ab,k,l)aaaa   <= sum(e-a)    [ <ab||ke>aaaa . T1o(e,l)aa ]
!     T2n(ab,ij)aaaa   <- Q(ab,i,j)aaaa - Q(ab,j,i)aaaa
!      Q(ab,k,l)bbbb   <= sum(e-b)    [ <ab||ke>bbbb . T1o(e,l)bb ]
!     T2n(ab,ij)bbbb   <- Q(ab,i,j)bbbb - Q(ab,j,i)bbbb
!     T2n(a,b,i,j)abab <- - sum(e-a)  [ <ab||je>abba . T1o(e,i)aa ]
!        <- + sum(e-b)  [ <ab||ie>abab . T1o(e,j)bb ]
!
!     General Status:

!     1) Tau ampitudes are stored and properly mapped in V1-V3
!     V1(ij,ef)   = Tau(ef,ij)aaaa
!     V2(ij,ef)   = Tau(ef,ij)bbbb
!     V3(i,j,e,f) = Tau(e,f,i,j)abab
!
!     2) Integrals are stored in file nab
!
!     syma=1,nsym
!       symb=1,syma
!       ! one record mapdn,mapin - in all cases
!         a=1,nvb(syma)
!         if (syma.eq.symb) then
!         limb=a
!         else
!         limb=nvb(symb)
!         end if
!           b=1,limb
!             ! one record N _a_b(p,q) - if any
!           end b
!         end a
!       end sumb
!     end syma
!
!     N.B. nab musi byt dakde deklarovane (daky common)
!
!     N.B. II. in parallel.fh there are:
!     nprocab - number of nodeschoosed for 'sumoverab' part
!     idab    - array of id's of these nodes
!
       use Para_Info, only: MyRank
        implicit none
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "wrk.fh"
#include "parallel.fh"
      INTEGER lunt2o1,lunt2o2,lunt2o3
      integer nabstack,possabstack
      integer niter
!
!     help variables
!
      INTEGER key,aeqb,yes
      INTEGER a,b,syma,symb,ssh3,ssn,ssh1
      INTEGER nlength,posst,rc
      INTEGER lunab,nsa,nsb,limb
      integer nab,nadd
      integer idtot(1:maxproc)
!
      integer iab,nabnow,left,bstart,bstop,todo,possab,bb
      integer nhelp
!
      todo=0
      possab=possabstack
!
!
!     I.parallel
!
!     I.par.1  - escape, if this node is not reserved for sumoverab
      yes=0
      do a=1,nprocab
        if (idab(a).eq.myRank) then
!       if ((idab(a).eq.myRank).and.(ideffab(a).gt.0.0d0)) then
        yes=a
        end if
      end do
!
      if (yes.eq.0) then
      return
      end if
!
!      I.par.2 - cal overal number of ab records
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
!
!       I.par.3 - distribute and calc add parameter
      call sumabdistt (nab,idtot)
        nab=idtot(yes)
!        escape, if there is nothing to do on this node
        if (nab.eq.0) return
        nadd=0
        do b=1,yes-1
        nadd=nadd+idtot(b)
        end do
!
!      I.1  open nab file
!
      lunab=10
      CALL filemanager (4, lunab, rc)
!
!      I.2  prepair
!     V1(ij,ef)   = Tau(ef,ij)aaaa
!     V2(ij,ef)   = Tau(ef,ij)bbbb
!     V3(i,j,e,f) = Tau(e,f,i,j)abab
!
      CALL filemanager (2, lunt2o1, rc)
      CALL getmediate (wrk,wrksize,                                     &
     &lunt2o1, possv40, mapdv4, mapiv4, rc)
      CALL mktau (wrk,wrksize,                                          &
     &mapdv4, mapiv4, mapdt11, mapit11, mapdt12, mapit12,               &
     &1.0d0, rc)
      CALL map (wrk,wrksize,                                            &
     &4, 3, 4, 1, 2, mapdv4, mapiv4, 1, mapdv1, mapiv1,                 &
     &possv10, posst, rc)
!
      CALL filemanager (2, lunt2o2, rc)
      CALL getmediate (wrk,wrksize,                                     &
     &lunt2o2, possv40, mapdv4, mapiv4, rc)
      CALL mktau (wrk,wrksize,                                          &
     &mapdv4, mapiv4, mapdt11, mapit11, mapdt12, mapit12,               &
     &1.0d0, rc)
      CALL map (wrk,wrksize,                                            &
     &4, 3, 4, 1, 2, mapdv4, mapiv4, 1, mapdv2, mapiv2,                 &
     &possv20, posst, rc)
!
      CALL filemanager (2, lunt2o3, rc)
      CALL getmediate (wrk,wrksize,                                     &
     &lunt2o3, possv40, mapdv4, mapiv4, rc)
      CALL mktau (wrk,wrksize,                                          &
     &mapdv4, mapiv4, mapdt11, mapit11, mapdt12, mapit12,               &
     &1.0d0, rc)
      CALL map (wrk,wrksize,                                            &
     &4, 3, 4, 1, 2, mapdv4, mapiv4, 1, mapdv3, mapiv3,                 &
     &possv30, posst, rc)
!
!      II   sum over ab
      DO 20 syma=1,nsym
        DO 21 symb=1,syma
        if (fullprint.ge.2) then
        write (6,*) ' SymA, SymB ',syma,symb
        end if
!
!      II.1   read mapdn,mapin
        CALL getmap (lunab, possn0, nlength, mapdn, mapin, rc)
!
!      II.3     skip sum over a,b if N is empty
        IF (nlength.eq.0) GO TO 21
!
!      II.4   skip sum over a,b if # of ab is zero
        IF ((nvb(syma)*nvb(symb)).eq.0) GO TO 21
!
!      II.5   def symmetry of N
        ssn=mmul(syma,symb)
!
!      II.5.1   def mapd and mapi for R1 _a,_b(j,e)aaaa = <ab||je>  (pos
        CALL grc0 (2, 0, 1, 3, 0, 0, ssn, possm10, posst, mapdr1,       &
     &   mapir1)
!      II.5.2   def mapd and mapi for R2 _a,_b(j,e)bbbb = <ab||je>  (pos
        CALL grc0 (2, 0, 2, 4, 0, 0, ssn, possm20, posst, mapdr2,       &
     &   mapir2)
!      II.5.3   def mapd and mapi for R3 _a,_b(j,e)abba = <ab||je>  (pos
        CALL grc0 (2, 0, 2, 3, 0, 0, ssn, possm30, posst, mapdr3,       &
     &   mapir3)
!      II.5.4   def mapd and mapi for R4 _b,_a(j,e)abba = <ba||je>  (pos
        CALL grc0 (2, 0, 2, 3, 0, 0, ssn, possm40, posst, mapdr4,       &
     &   mapir4)
!      II.5.5   def mapd and mapi for R5 _a,_b(j,e)abab = <ab||je>  (pos
        CALL grc0 (2, 0, 1, 4, 0, 0, ssn, possh10, posst, mapdr5,       &
     &   mapir5)
!      II.5.6   def mapd and mapi for R6 _b,_a(j,e)abab = <ba||je>  (pos
        CALL grc0 (2, 0, 1, 4, 0, 0, ssn, possh20, posst, mapdr6,       &
     &   mapir6)
!
!      II.6.1   def mapd and mapi for M1 _a,_b(ef)aaaa = <ab||ef>
        CALL grc0 (2, 1, 3, 3, 0, 0, ssn, possm10, posst, mapdm1,       &
     &   mapim1)
!      II.6.2   def mapd and mapi for M2 _a,_b(ef)bbbb = <ab||ef>
        CALL grc0 (2, 1, 4, 4, 0, 0, ssn, possm20, posst, mapdm2,       &
     &   mapim2)
!      II.6.3   def mapd and mapi for M3 _a,_b(e,f)abab = <ab||ef>
        CALL grc0 (2, 0, 3, 4, 0, 0, ssn, possm30, posst, mapdm3,       &
     &   mapim3)
!      II.6.4   def mapd and mapi for M4 _b,_a(e,f)abab = <ba||fe>
        CALL grc0 (2, 0, 3, 4, 0, 0, ssn, possm40, posst, mapdm4,       &
     &   mapim4)
!
!      II.7   def # of S orbitals in syma and symb
        nsa=nvb(syma)-nva(syma)
        nsb=nvb(symb)-nva(symb)
!
        DO 10 a=1,nvb(syma)
!
        if (fullprint.ge.3) then
        write (6,*) ' A', a
        end if
!
!      II.7     def limitations for b
          IF (syma.eq.symb) THEN
            limb=a
          ELSE
            limb=nvb(symb)
          END IF
!
!@@
!par    for nodes, selected for sumab process:
!       1) only blank read (instead of _a_b(e,f) form lunab file for records,
!       pre-standing to the portion, that is going to be done on this node
!       2) prepair setup (def bstart,limb), if needed
!
        if (nadd.ge.limb) then
!       all cycle over b need to be skipped
          do b=1,limb
!           blank read and nothing more
            call reajalovy (lunab, nlength, wrk(possn0))
            nadd=nadd-1
          end do
!       no setup needed
          goto 10
!
        else if (nadd.gt.0) then
!       setup first
          bstart=nadd
          if (nab.lt.(limb-bstart)) then
            limb=bstart+nab
          end if
!       part of the cycle over b need to be skipped
          do b=1,nadd
!           blank read and nothing more
            call reajalovy (lunab, nlength, wrk(possn0))
          end do
          nadd=0
!
        else
!       nabadd musi byt 0, inak nieje dobre
!       no part of the cycle over b need to be skipped
!       only setup
          bstart=0
          if (nab.lt.limb) then
            limb=nab
          end if
        end if
!
!
!*      define conditions for reading to stack
!       bstart - lower limit of b-sumation
!       bstop  - upper limit of b-sumation
!       nabnow - actual size af stack (if=1, no stacking)
!
!*.1    initialization
!Stare  bstart=0
!
!*.1    set unstacking conditions
4       bstart=bstart+1
        bstop=bstart
        nabnow=1
!
!*.2    without stacking if: 1)a=b; 2)a or b \in S; 3) nabnow.le.3
!       (N.B. it was found inconvenient to execute stacking if
!             size of stack is less then 4,or even more)
!
!       ad 1)
        if ((bstart.eq.a).and.(syma.eq.symb)) goto 6
!       ad 2)
        if ((bstart.le.nsb).or.(a.le.nsa)) goto 6
!       ad 3)
        left=limb-bstart+1
        if (syma.eq.symb) left=left-1
        if (left.le.3) goto 6
!
!*.3    set stacking conditions for real stacking
        if (left.le.nabstack) then
          nabnow=left
        else
          nabnow=nabstack
        end if
        bstop=bstart+nabnow-1
        todo=1
        possab=possabstack
!
!*.     reading to N or to stack
6       if (nabnow.eq.1) then
          CALL rea (lunab, nlength, wrk(possn0))
        else
          do iab=1,nabnow
          CALL rea (lunab, nlength, wrk(possabstack+(iab-1)*nlength))
          if ((possabstack+iab*nlength).ge.possm10) then
          end if
          end do
!         qq=wrk(possabstack)
        end if
!
!
          do 8 b=bstart,bstop
!
!      III.3       read N
!Stare    CALL rea (lunab, nlength, wrk(possn0))
!
!      III.4       def aeqb
          IF ((syma.eq.symb).and.(a.eq.b)) THEN
            aeqb=1
          ELSE
            aeqb=0
          END IF
!
!      III.5       def key
          IF ((a.gt.nsa).and.(b.gt.nsb)) THEN
            key=4
          ELSE IF ((a.gt.nsa).and.(b.le.nsb)) THEN
            key=2
          ELSE IF ((a.le.nsa).and.(b.gt.nsb)) THEN
            key=3
          ELSE
            key=1
          END IF
!
!      def yes for T25.1
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
!
!
!           T25
!           T2n(ab,ij)aaaa   <- sum(e>f-aa) [ <ab||ef>aaaa . Tau(ef,ij)aaaa ]
!           T2n(ab,ij)bbbb   <- sum(e>f-bb) [ <ab||ef>bbbb . Tau(ef,ij)bbbb ]
!           T2n(a,b,i,j)abab <- sum(e,f-ab) [ <ab||ef>abab . Tau(ef,ij)abab ]
!
        if (nabnow.eq.1) then
!       AB stacking is not active
!
!      III.6       get
!           -> M1 _a,_b(ef)aaaa  = <ab||ef>
!           -> M2 _a,_b(ef)bbbb  = <ab||ef>
!           -> M3 _a,_b(e,f)abab = <ab||ef>
!           -> M4 _b,_a(e,f)abab = <ba||fe>
!           free: H1-H4 (V4 reserved for stacking)
          CALL unpackab1 (wrk,wrksize,                                  &
     &     mapdn, mapin, mapdm1, mapim1, mapdm2, mapim2,                &
     &     mapdm3, mapim3, mapdm4, mapim4, ssn, key, aeqb)
!
!
          IF (yes.eq.1) THEN
!      T25.1.1      H1 (ij) = V1(ij,ef) . M1 (ef)
            CALL mult (wrk,wrksize,                                     &
     &       4, 2, 2, 2, mapdv1, mapiv1, 1, mapdm1, mapim1,             &
     &       ssn, mapdh1, mapih1,ssh1,possh10, rc)
!      T25.1.2      T2n(ab,ij)aaaa <- 1.0d0 . H1(ij)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a-nsa, b-nsb, syma, symb, 1.0d0,               &
     &       mapdh1, ssn, mapdt21, mapit21, 1, rc)
          END IF
!
          IF (aeqb.eq.0) THEN
!      T25.2.1      H1 (ij) = V2(ij,ef) . M2 (ef)
            CALL mult (wrk,wrksize,                                     &
     &       4, 2, 2, 2, mapdv2, mapiv2, 1, mapdm2, mapim2,             &
     &       ssn, mapdh1, mapih1,ssh1,possh10, rc)
!      T25.2.2      T2n(ab,ij)bbbb <- 1.0d0 . H1(ij)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a, b, syma, symb, 1.0d0, mapdh1, ssn,          &
     &       mapdt22, mapit22, 1, rc)
          END IF
!
          IF ((key.eq.2).or.(key.eq.4)) THEN
!      T25.3.1      H1 (i,j) = V3(i,j,e,f) . M3 (e,f)
            CALL mult (wrk,wrksize,                                     &
     &       4, 2, 2, 2, mapdv3, mapiv3, 1, mapdm3, mapim3,             &
     &       ssn, mapdh1, mapih1,ssh1,possh10, rc)
!      T25.3.2      T2n(a,b,i,j)abab <- 1.0d0 . H1(i,j)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a-nsa, b, syma, symb, 1.0d0, mapdh1,           &
     &       ssn, mapdt23, mapit23, 1, rc)
          END IF
!
          IF (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) THEN
!      T25.3.3      H1 (i,j) = V3(i,j,e,f) . M4 (e,f)
            CALL mult (wrk,wrksize,                                     &
     &       4, 2, 2, 2, mapdv3, mapiv3, 1, mapdm4, mapim4,             &
     &       ssn, mapdh1, mapih1,ssh1,possh10, rc)
!
!      T25.3.4      T2n(a,b,i,j)abab <- 1.0d0 . H1(i,j)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, b-nsb, a, symb, syma, 1.0d0, mapdh1,           &
     &       ssn, mapdt23, mapit23, 1, rc)
          END IF
!
!
        else
!       AB stacking is active
!
          if (todo.eq.1) then
!         once it need to be done
          todo=0
!
!      III.6       get
!           -> M1 _a,_b(ef)aaaa  = <ab||ef>
!           -> M2 _a,_b(ef)bbbb  = <ab||ef>
!           -> M3 _a,_b(e,f)abab = <ab||ef>
!           -> M4 _b,_a(e,f)abab = <ba||fe>
!           free: H1-H4 (V4 reserved for stacking)
!
!
!       def mapd and mapi for V4 _a,(ef,Bp)
       call grc0stack (nabnow,1,3,3,3,0,ssn,possv40,posst,mapdv4,mapiv4)
!           -> V4 _a,_b(ef)aaaa  = <ab||ef>
       call unpackab3 (wrk,wrksize,                                     &
     &      mapdn,mapin,mapdv4,mapiv4,ssn,nabnow,possabstack,nlength,1)
!
!      T25.1.1      H2 (ij,Bp) = V1(ij,ef) . V4 (ef,Bp)
        call multstack (wrk,wrksize,                                    &
     &                  mapdv1,mapdv4,mapdh2,mapiv1,mapiv4,mapih2,      &
     &                  1,ssn,possh20,nabnow)

!       def mapd and mapi for H1 _a,_b(ij)aa
        call grc0 (2,1,1,1,0,0,ssn,possh10,posst,mapdh1,mapih1)
        do bb=bstart,bstop
!      T25.1.2      ext H1(ij) <- H2(ij,_Bb) (only data transfer)
        call extstack (wrk,wrksize,                                     &
     &                 mapdh1,mapdh2,(bb-bstart+1),nabnow)
!      T25.1.3      T2n(ab,ij)aaaa <- 1.0d0 . H1(ij)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a-nsa, bb-nsb, syma, symb, 1.0d0,              &
     &       mapdh1, ssn, mapdt21, mapit21, 1, rc)
        end do
!
!
!
!       def mapd and mapi for V4 _a,(ef,Bp)
       call grc0stack (nabnow,1,4,4,4,0,ssn,possv40,posst,mapdv4,mapiv4)
!           -> V4 _a,_b(ef)bbbb  = <ab||ef>
       call unpackab3 (wrk,wrksize,                                     &
     &      mapdn,mapin,mapdv4,mapiv4,ssn,nabnow,possabstack,nlength,2)
!
!      T25.2.1      H2 (ij,Bp) = V2(ij,ef) . V4 (ef,Bp)
        call multstack (wrk,wrksize,                                    &
     &                  mapdv2,mapdv4,mapdh2,mapiv2,mapiv4,mapih2,      &
     &                  1,ssn,possh20,nabnow)
!
!       def mapd and mapi for H1 _a,_b(ij)bb
        call grc0 (2,1,2,2,0,0,ssn,possh10,posst,mapdh1,mapih1)
        do bb=bstart,bstop
!      T25.2.2      ext H1(ij) <- H2(ij,Bb)
        call extstack (wrk,wrksize,                                     &
     &                 mapdh1,mapdh2,(bb-bstart+1),nabnow)
!      T25.2.3      T2n(ab,ij)bbbb <- 1.0d0 . H1(ij)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a, bb, syma, symb, 1.0d0, mapdh1, ssn,         &
     &       mapdt22, mapit22, 1, rc)
        end do
!
!
!
!       def mapd and mapi for V4 _a,(ef,Bp)
       call grc0stack (nabnow,0,3,4,4,0,ssn,possv40,posst,mapdv4,mapiv4)
!           -> V4 _a,_b(e,f)abab = <ab||ef>
       call unpackab3 (wrk,wrksize,                                     &
     &      mapdn,mapin,mapdv4,mapiv4,ssn,nabnow,possabstack,nlength,3)
!
!      T25.3.1      H2 (i,j,Bp) = V3(i,j,e,f) . V4 (e,f,Bp)
        call multstack (wrk,wrksize,                                    &
     &                  mapdv3,mapdv4,mapdh2,mapiv3,mapiv4,mapih2,      &
     &                  1,ssn,possh20,nabnow)
!
!       def mapd and mapi for H1 _a,_b(ij)ab
        call grc0 (2,0,1,2,0,0,ssn,possh10,posst,mapdh1,mapih1)
        do bb=bstart,bstop
!      T25.3.2      ext H1(i,j) <- H2(ij,Bb)
        call extstack (wrk,wrksize,                                     &
     &                 mapdh1,mapdh2,(bb-bstart+1),nabnow)
!      T25.3.3      T2n(a,b,i,j)abab <- 1.0d0 . H1(i,j)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a-nsa, bb, syma, symb, 1.0d0, mapdh1,          &
     &       ssn, mapdt23, mapit23, 1, rc)
        end do
!
!
!
!       def mapd and mapi for V4 _a,(ef,Bp)
       call grc0stack (nabnow,0,3,4,3,0,ssn,possv40,posst,mapdv4,mapiv4)
!           -> V4 _b,_a(e,f)abab = <ba||fe>
       call unpackab3 (wrk,wrksize,                                     &
     &      mapdn,mapin,mapdv4,mapiv4,ssn,nabnow,possabstack,nlength,4)
!
!      T25.3.4      H2 (i,j,Bp) = V3(i,j,e,f) . V4 (e,f,Bp)
        call multstack (wrk,wrksize,                                    &
     &                  mapdv3,mapdv4,mapdh2,mapiv3,mapiv4,mapih2,      &
     &                  1,ssn,possh20,nabnow)
!
!       def mapd and mapi for H1 _a,_b(ij)ab
        call grc0 (2,0,1,2,0,0,ssn,possh10,posst,mapdh1,mapih1)
        do bb=bstart,bstop
!      T25.3.5      ext H1(i,j) <- H2(i,j,Bb)
        call extstack (wrk,wrksize,                                     &
     &                 mapdh1,mapdh2,(bb-bstart+1),nabnow)
!      T25.3.6      T2n(a,b,i,j)abab <- 1.0d0 . H1(i,j)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, bb-nsb, a, symb, syma, 1.0d0, mapdh1,          &
     &       ssn, mapdt23, mapit23, 1, rc)
        end do
!
          end if
        end if
!
!     T28
!      Q(ab,k,l)aaaa   <= sum(e-a)    [ <ab||ke>aaaa . T1o(e,l)aa ]
!     T2n(ab,ij)aaaa   <- Q(ab,i,j)aaaa - Q(ab,j,i)aaaa
!      Q(ab,k,l)bbbb   <= sum(e-b)    [ <ab||ke>bbbb . T1o(e,l)bb ]
!     T2n(ab,ij)bbbb   <- Q(ab,i,j)bbbb - Q(ab,j,i)bbbb
!     T2n(a,b,i,j)abab <- - sum(e-a)  [ <ab||je>abba . T1o(e,i)aa ]
!        <- + sum(e-b)  [ <ab||ie>abab . T1o(e,j)bb ]
!
!
!       put appropriate data from AB_stack to N mediate, if needed
        if (nabnow.gt.1) then
!       transfer data
          do nhelp=0,nlength-1
          wrk(possn0+nhelp)=wrk(possab+nhelp)
          end do
!       upgrade address
          possab=possab+nlength
        end if
!
!
!      III.7       get
!           -> M1 _a,_b(j,e)aaaa = <ab||je> (map is R1)
!           -> M2 _a,_b(j,e)bbbb = <ab||je> (map is R2)
!           -> M3 _a,_b(j,e)abba = <ab||je> (map is R3)
!           -> M4 _b,_a(j,e)abba = <ba||je> (map is R4)
!           -> H1 _a,_b(j,e)abab = <ab||je> (map is R5)
!           -> H2 _b,_a(j,e)abab = <ba||je> (map is R6)
!           free: H3,H4 (V4 reserved for stacking)
!
       CALL unpackab2 (wrk,wrksize,                                     &
     & mapdn,mapin,mapdr1,mapir1,mapdr2,mapir2,mapdr3,mapir3,           &
     & mapdr4,mapir4,mapdr5,mapir5,mapdr6,mapir6,ssn,key,aeqb)
!
!
!     T28
!      Q(ab,k,l)aaaa   <= sum(e-a)    [ <ab||ke>aaaa . T1o(e,l)aa ]
!     T2n(ab,ij)aaaa   <- Q(ab,i,j)aaaa - Q(ab,j,i)aaaa
!      Q(ab,k,l)bbbb   <= sum(e-b)    [ <ab||ke>bbbb . T1o(e,l)bb ]
!     T2n(ab,ij)bbbb   <- Q(ab,i,j)bbbb - Q(ab,j,i)bbbb
!     T2n(a,b,i,j)abab <- - sum(e-a)  [ <ab||je>abba . T1o(e,i)aa ]
!        <- + sum(e-b)  [ <ab||ie>abab . T1o(e,j)bb ]
!
!
!      def yes for T28.1
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
!
!
          IF (yes.eq.1) THEN
!      T28.1.1      H3(k,l) = M1(k,e) . T1o(e,l)aa
            CALL mult (wrk,wrksize,                                     &
     &       2, 2, 2, 1, mapdr1, mapir1, ssn, mapdt11,                  &
     &       mapit11, 1, mapdh3, mapih3,ssh3,possh30, rc)
!      T28.1.2      H4(kl)  = H3(k,l)-H3(l,k)
            CALL fack (wrk,wrksize,                                     &
     &       2, 1, mapdh3, ssh3, mapih3, mapdh4, mapih4,                &
     &       possh40, rc)
!      T28.1.3      T2n(ab,ij)aaaa <- 1.0d0 H4(ij)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a-nsa, b-nsb, syma, symb, 1.0d0,               &
     &       mapdh4, ssn, mapdt21, mapit21, 1, rc)
          END IF
!
          IF (aeqb.eq.0) THEN
!      T28.2.1      H3(k,l) = M2(k,e) . T1o(e,l)bb
            CALL mult (wrk,wrksize,                                     &
     &       2, 2, 2, 1, mapdr2, mapir2, ssn, mapdt12,                  &
     &       mapit12, 1, mapdh3, mapih3,ssh3,possh30, rc)
!      T28.2.2      H4(kl)  = H3(k,l)-H3(l,k)
            CALL fack (wrk,wrksize,                                     &
     &       2, 1, mapdh3, ssh3, mapih3, mapdh4, mapih4,                &
     &       possh40, rc)
!      T28.2.3      T2n(ab,ij)bbbb <- 1.0d0 H4(ij)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a, b, syma, symb, 1.0d0, mapdh4, ssn,          &
     &       mapdt22, mapit22, 1, rc)
          END IF
!
          IF ((key.eq.2).or.(key.eq.4)) THEN
!      T28.3.1      H3(j,i)  = M3(j,e) . T1o(e,i)aa
            CALL mult (wrk,wrksize,                                     &
     &       2, 2, 2, 1, mapdr3, mapir3, ssn, mapdt11,                  &
     &       mapit11, 1, mapdh3, mapih3,ssh3,possh30, rc)
!      T28.3.2      H4(i,j)  = H3(j,i)
            CALL map (wrk,wrksize,                                      &
     &       2, 2, 1, 0, 0, mapdh3, mapih3, ssn, mapdh4,                &
     &       mapih4, possh40, posst, rc)
!      T28.3.3      T2n(a,b,i,j)abab <- -1.0d0 H4(i,j)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a-nsa, b, syma, symb, -1.0d0, mapdh4,          &
     &       ssn, mapdt23, mapit23, 1, rc)
          END IF
!
          IF (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) THEN
!      T28.3.4      H3(j,i)  = M4(j,e) . T1o(e,i)aa
            CALL mult (wrk,wrksize,                                     &
     &       2, 2, 2, 1, mapdr4, mapir4, ssn, mapdt11,                  &
     &       mapit11, 1, mapdh3, mapih3,ssh3,possh30, rc)
!      T28.3.5      H4(i,j)  = H3(j,i)
            CALL map (wrk,wrksize,                                      &
     &       2, 2, 1, 0, 0, mapdh3, mapih3, ssn, mapdh4,                &
     &       mapih4, possh40, posst, rc)
!      T28.3.6      T2n(a,b,i,j)abab <- -1.0d0 H4(i,j)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, b-nsb, a, symb, syma, -1.0d0, mapdh4,          &
     &       ssn, mapdt23, mapit23, 1, rc)
          END IF
!
          IF ((key.eq.2).or.(key.eq.4)) THEN
!      T28.3.7      H3(i,j)  = H1(i,e) . T1o(e,j)bb
            CALL mult (wrk,wrksize,                                     &
     &       2, 2, 2, 1, mapdr5, mapir5, ssn, mapdt12,                  &
     &       mapit12, 1, mapdh3, mapih3,ssh3,possh30, rc)
!      T28.3.8      T2n(a,b,i,j)abab <- 1.0d0 H3(i,j)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, a-nsa, b, syma, symb, 1.0d0, mapdh3,           &
     &       ssn, mapdt23, mapit23, 1, rc)
          END IF
!
          IF (((key.eq.3).or.(key.eq.4)).and.(aeqb.eq.0)) THEN
!      T28.3.9      H3(j,i)  = H2(i,e) . T1o(e,j)bb
            CALL mult (wrk,wrksize,                                     &
     &       2, 2, 2, 1, mapdr6, mapir6, ssn, mapdt12,                  &
     &       mapit12, 1, mapdh3, mapih3,ssh3,possh30, rc)
!      T28.3.10     T2n(a,b,i,j)abab <- 1.0d0 H3(i,j)
            CALL add (wrk,wrksize,                                      &
     &       2, 4, 2, 5, b-nsb, a, symb, syma, 1.0d0, mapdh3,           &
     &       ssn, mapdt23, mapit23, 1, rc)
          END IF
!
!
  8     CONTINUE
!
!par
        nab=nab-(bstop-bstart+1)
        if (nab.eq.0) then
        goto 30
        end if
!
        if (bstop.lt.limb) then
          bstart=bstop
          goto 4
        end if
 10     CONTINUE
!
 21   CONTINUE
 20   CONTINUE
!
!      IV   close lunab
30    CALL filemanager (3, lunab, rc)
!
      RETURN
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(niter)
      END
