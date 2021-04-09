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
        subroutine IniReord_t3(NaGrp,wrksize)
c
c       nacitanie vsupu a inicializacia premnennych
c       a tlac primitivnej hlavicky pre Reord procesz
c
#ifdef _MOLCAS_MPP_
        use Para_Info, only: MyRank, nProcs
#endif
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
cmp
#include "cholesky.fh"
#include "ccsd_t3compat.fh"
cmp
c
        integer NaGrp
        integer wrksize
cmp!
        integer nOrb(8),nOcc(8)
        integer ndelvirt

        integer LuSpool
        character*80 LINE
cmp
        integer rc
        real*8 FracMem
        character*3 msg

#ifdef _MOLCAS_MPP_
        integer jal1, jal2
#endif
cmp

c setup defaults

        Call Get_iArray('nOrb',nOrb,1)
        Call Get_iArray('nIsh',nOcc,1)

c
        no = nOcc(1)
        nv = nOrb(1) - nOcc(1)
c
        FracMem=0.0d0
        Call Cho_X_init(rc,FracMem) ! initialize cholesky info
c
c       take local # of Cholesky Vectors on this node
#ifdef _MOLCAS_MPP_
c
        do jal1=0,Nprocs-1
          NChLoc(jal1)=0
        end do
c
        NChLoc(MyRank)=NumCho(1)

        call gaigop (NChLoc(0),NProcs,'+')
c
        jal2=0
        do jal1=0,NProcs-1
          jal2=jal2+NChLoc(jal1)
        end do

        nc = jal2
#else
        nc = NumCho(1)
#endif

        Call Cho_X_final(rc)

        ndelvirt = 0
        LunAux = 13
        mhkey = 1
        generkey = 1
cmp!        NaGrp = 1
        Call get_iScalar('CHCCLarge',NaGrp)
        restkey = 0
        printkey = 1

c t3 specific keywords

        gen_files   = .True.
        run_triples = .True.
        t3_starta   = -1
        t3_stopa    = -1
        t3_startb   = -1
        t3_stopb    = -1
c
cmp!    read input file
c
      LuSpool = 17
      Call SpoolInp(LuSpool)
      Rewind(LuSpool)
 5      Read(LuSpool,'(A80)') LINE
       CALL UPCASE(LINE)
       IF( INDEX(LINE,'&CHT3') .EQ. 0 ) GOTO 5
 6     Read(LuSpool,'(A80)') LINE
       IF(LINE(1:1).EQ.'*') GOTO 6
       CALL UPCASE(LINE)
c
       IF (LINE(1:4).EQ.'TITL') THEN
       Read(LuSpool,*)

       ELSE IF (LINE(1:4).EQ.'FROZ') THEN ! FROZen
       Read(LuSpool,*) nfr
             if ((nfr.lt.0).or.(nfr.ge.no)) then
               write (6,*)
               write (6,*) 'Ilegal value for FROZen keyword : ',
     &                      nfr
               call abend()
             end if
             no = no - nfr

       ELSE IF (LINE(1:4).EQ.'DELE') THEN ! DELEted
       Read(LuSpool,*) ndelvirt
             if ((ndelvirt.lt.0).or.(ndelvirt.ge.nv)) then
               write (6,*)
               write (6,*) 'Ilegal value for DELEted keyword : ',
     &                      ndelvirt
               call abend()
             end if
             nv = nv - ndelvirt

cmp!       ELSE IF (LINE(1:4).EQ.'LARG') THEN ! LARGegroup
cmp!       Read(LuSpool,*) NaGrp
cmp!        if ((NaGrp.lt.1).or.(NaGrp.gt.32)) then
cmp!         write (6,*)
cmp!         write (6,*) 'Ilegal value for LARGegroup keyword : ',
cmp!     &                NaGrp
cmp!         write (6,*) 'Large segmentation must be -le 32'
cmp!         call abend()
cmp!        end if

cmp!       ELSE IF (LINE(1:4).EQ.'LUNA') THEN  ... toto sa nikdy nevyuzivalo
cmp!       Read(LuSpool,*) LunAux

       ELSE IF (LINE(1:4).EQ.'MHKE') THEN ! MHKEy
       Read(LuSpool,*) mhkey
           if ((mhkey.lt.0).or.(mhkey.gt.2)) then
              mhkey=1
              write(6,*)
              write(6,*) ' Warning!!! ',
     &                   ' MHKEy out of range, changed to 1'
           end if

       ELSE IF (LINE(1:4).EQ.'REST') THEN ! RESTart
            restkey = 1
            write (6,*)
            write (6,*) 'RESTart option is temporary disabled'
            write (6,*) 'No Restart possible (... yet).'
            call abend()

       ELSE IF (LINE(1:4).EQ.'PRIN') THEN ! PRINtkey
       Read(LuSpool,*) printkey
          if (((printkey.lt.0).or.(printkey.gt.10)).or.
     & ((printkey.gt.2).and.(printkey.lt.10))) then

            write (6,*)
            write (6,*) 'Ilegal value of the PRINtkey keyword: ',
     &                   printkey
            write (6,*) ' Use: 1  (Minimal) '
            write (6,*) '      2  (Minimal + Timings)'
            write (6,*) '      10 (Debug) '
            call abend()
          end if

       ELSE IF (LINE(1:4).EQ.'NOGE') THEN ! NOGEnerate
            gen_files = .False.

       ELSE IF (LINE(1:4).EQ.'NOTR') THEN ! NOTRiples
            run_triples = .False.

       ELSE IF (LINE(1:4).EQ.'ALOO') THEN ! ALOOp
       Read(LuSpool,*) t3_starta, t3_stopa
            if ((t3_starta.lt.-1).or.(t3_stopa.lt.-1)) then
               write (6,*) 'ALOOp values can be either: '
               write (6,*) '"-1" : indicating normal run, or'
               write (6,*) 'positive numbers!'
               call abend()
            end if

       ELSE IF (LINE(1:4).EQ.'BLOO') THEN ! BLOOp
       Read(LuSpool,*) t3_startb, t3_stopb
            if ((t3_startb.lt.-1).or.(t3_stopb.lt.-1)) then
               write (6,*) 'BLOOp values can be either: '
               write (6,*) '"-1" : indicating normal run, or'
               write (6,*) 'positive numbers!'
               call abend()
            end if

       ELSE IF (LINE(1:4).EQ.'END ') THEN
       GOTO 7
       END IF
       GOTO 6
7      CONTINUE

       Call Close_LuSpool(LuSpool)

c! take care of the cholesky vectors segmentation
c! to lead to < 100 blocks

cmp checks
        if (t3_starta.gt.t3_stopa) then
          write (6,*) 'Mismatch in input : '
          write (6,*) 'T3_STARTA = ',t3_starta
          write (6,*) 'T3_STOPA = ',t3_stopa
          call abend()
        end if
c
        if (t3_startb.gt.t3_stopb) then
          write (6,*) 'Mismatch in input : '
          write (6,*) 'T3_STARTB = ',t3_startb
          write (6,*) 'T3_STOPB = ',t3_stopb
          call abend()
        end if
c
        if ((t3_starta.lt.0).and.(t3_stopa.gt.0)) then
          write (6,*) 'Mismatch in input : '
          write (6,*) 'T3_STARTA = ',t3_starta
          write (6,*) 'T3_STOPA = ',t3_stopa
          call abend()
        end if
c
        if ((t3_startb.lt.0).and.(t3_stopb.gt.0)) then
          write (6,*) 'Mismatch in input : '
          write (6,*) 'T3_STARTB = ',t3_startb
          write (6,*) 'T3_STOPB = ',t3_stopb
          call abend()
        end if
c
cmp!        if ((t3_starta.gt.0).and.(t3_startb.lt.0)) then
cmp!          write (6,*) 'This restart combination not implemented'
cmp!          write (6,*) 'T3_STARTA = ',t3_starta
cmp!          write (6,*) 'T3_STARTB = ',t3_startb
cmp!          call abend
cmp!        end if
c
c2      tlac hlavicky
        write (6,*)
        write (6,*) '    Cholesky Based Closed-Shell (T) code'
        write (6,*)
      write (6,*) '--------------------------------------------------'

        write (6,'(A,i9)') ' Frozen Orbitals                   : ',
     & nfr
        write (6,'(A,i9)') ' Occupied Orbitals                 : ',
     & no
        write (6,'(A,i9)') ' Virtual Orbitals                  : ',
     & nv
        write (6,'(A,i9)') ' Total number of Cholesky Vectors  : ',
     & nc

      write (6,*) '--------------------------------------------------'

        write (6,'(A,i9)') ' Large Virtual Segmentation        : ',
     & NaGrp

      write (6,*) '--------------------------------------------------'

        msg = 'No'
        if (gen_files) msg = 'Yes'

        write (6,'(A,A5)') ' Generate Triples Scratch Files?        : ',
     & msg

        msg = 'No'
        if (.not.run_triples) msg = 'Yes'

        write (6,'(A,A5)') ' Stop after Scratch Files generation?   : ',
     & msg

      write (6,*) '--------------------------------------------------'

        if (t3_starta.eq.-1) then
          write (6,'(A,i4)') ' Calculating full loop A                '
        else
           write (6,'(A,i4)')
     & ' VO index triplet to start with in loop A : ',t3_starta
           write (6,'(A,i4)')
     & ' VO index triplet to stop  at   in loop A : ',t3_stopa
        end if

        if (t3_starta.eq.-1) then
          write (6,'(A,i4)') ' Calculating full loop B                '
        else
           write (6,'(A,i4)')
     & ' VO index triplet to start with in loop B : ',t3_startb
           write (6,'(A,i4)')
     & ' VO index triplet to stop  at   in loop B : ',t3_stopb
        end if

      write (6,*) '--------------------------------------------------'

        write (6,'(A,i9)') ' Lun Number for Aux. Matrixes      : ',
     & LunAux
        write (6,'(A,i9)') ' BLAS/FTN Matrix Handling          : ',
     & mhkey

        msg = 'No'
        if (restkey.eq.1) msg = 'Yes'

        write (6,'(A,A10)') ' Start from RstFil ?               : ',
     & msg
        write (6,'(A,i9)') ' Print level                       : ',
     & printkey

      write (6,*) '--------------------------------------------------'
        write (6,*)
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_integer(wrksize)
        end
