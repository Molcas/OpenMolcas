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
        kam=-1 ! dummy initialize
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
