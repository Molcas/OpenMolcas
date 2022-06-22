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
        subroutine IniReord(NaGrp,NaSGrp,NchBlk,LunAux,wrksize)
c
c       nacitanie vsupu a inicializacia premnennych
c       a tlac primitivnej hlavicky pre Reord procesz
c
#ifdef _MOLCAS_MPP_
        use Para_Info, only: nProcs
#endif
        implicit none
#include "chcc1.fh"
#include "chcc_reord.fh"
cmp!
#include "parcc.fh"
cmp!

c
        integer NaGrp,NaSGrp,NchBlk
        integer LunAux,wrksize
cmp!
        integer nOrb(8),nOcc(8),nFro(8),nDel(8)
        integer intkey1,intkey2
        integer ndelvirt

        integer LuSpool
        character*80 LINE

#ifdef _MOLCAS_MPP_
        integer jal1
#endif
        integer NChLoc_min,NChLoc_max,NchBlk_tmp

        character*3 msg
cmp!

c setup defaults

        Call Get_iArray('nBas',nOrb,1) ! must read always nBas!!
        Call Get_iArray('nIsh',nOcc,1)
        Call Get_iArray('nFroPT',nFro,1) ! = 'nFro' in previous step
        Call Get_iArray('nDelPT',nDel,1) ! = 'nDel' in previous step
c
#ifdef _MOLCAS_MPP_
cmp     get min/max
        NChLoc_min=nc
        NChLoc_max=0
        do jal1=0,Nprocs-1
          if (NChLoc(jal1).le.NChLoc_min) NChLoc_min = NChLoc(jal1)
          if (NChLoc(jal1).ge.NChLoc_max) NChLoc_max = NChLoc(jal1)
        end do

cmp     calc reasonable starting value (200-300)
        if (NChLoc_min.ne.NChLoc_max) then
          NChBlk = int(NChLoc_min/2)
        else
          NChBlk = NChLoc_min
        end if

        if (NChBlk.ge.300) then
          NChBlk=min(200,int(NChBlk/2))
        end if

cmp     fix num of ChV blocks to be less then 100
        if (int(NChLoc_max/NChBlk).ge.100)
     & NChBlk = int(NChLoc_max/100) - 1
#else
        NChLoc_max=nc
cmp     calc reasonable starting value (200-300)
        if (nc.ge.300) then
          NChBlk=min(200,nc/2)
        else
          NChBlk=nc
        endif

        NChLoc_min=NChBlk
cmp     fix num of ChV blocks to be less then 100
        if (int(nc/NChBlk).ge.100)
     & NChBlk = int(nc/100) - 1
#endif

        nfr= nFro(1)
        no = nOcc(1)-nFro(1)
        ndelvirt = nDel(1)

        nv = nOrb(1)-nDel(1)-nOcc(1) ! nOrb defined as = to nBas, right!

        LunAux = 13
        mhkey = 1
        generkey = 1
        intkey1 = 0
        intkey2 = 0

        NaGrp = 0
        NaSGrp = 0
        W34DistKey = 1
        JoinLkey = 2 ! toto este nemam domyslene
        restkey = 0
        conv = 1.0d-6
        printkey = 1
        maxiter = 40

c
cmp!    read input file
c
      LuSpool = 17
      Call SpoolInp(LuSpool)
      Rewind(LuSpool)
 5      Read(LuSpool,'(A80)') LINE
       CALL UPCASE(LINE)
       IF( INDEX(LINE,'&CHCC') .EQ. 0 ) GOTO 5
 6     Read(LuSpool,'(A80)') LINE
       IF(LINE(1:1).EQ.'*') GOTO 6
       CALL UPCASE(LINE)
c
       IF (LINE(1:4).EQ.'TITL') THEN
       Read(LuSpool,*)

       ELSE IF (LINE(1:4).EQ.'FROZ') THEN
       Read(LuSpool,*) nfr
             if ((nfr.lt.0).or.(nfr.ge.no)) then
               write (6,*)
               write (6,*) 'Ilegal value for FROZen keyword : ',
     &                      nfr
               call abend()
             end if
             no = no + nFro(1) - nfr

       ELSE IF (LINE(1:4).EQ.'DELE') THEN
       Read(LuSpool,*) ndelvirt
             if ((ndelvirt.lt.0).or.(ndelvirt.gt.nv)) then
               write (6,*)
               write (6,*) 'Ilegal value for DELETED keyword : ',
     &                      ndelvirt
               call abend()
             end if
             nv = nv + nDel(1) - ndelvirt

       ELSE IF (LINE(1:4).EQ.'LARG') THEN
       Read(LuSpool,*) NaGrp
           if ((NaGrp.lt.0).or.(NaGrp.gt.maxGrp)) then
               write (6,*)
               write (6,*) 'Ilegal value for LARGE keyword : ',
     &                      NaGrp
               write (6,*) 'Large segmentation must be -le 32'
               call abend()
           end if

       ELSE IF (LINE(1:4).EQ.'SMAL') THEN
       Read(LuSpool,*) NaSGrp
           if ((NaSGrp.lt.0).or.(NaSGrp.gt.8)) then
               write (6,*)
               write (6,*) 'Ilegal value for SMALL keyword : ',
     &                      NaSGrp
               write (6,*) 'Small segmentation must be -le 8'
               call abend()
           end if

c          large == 0, small != 0 => quit
           if ((NaGrp.eq.0).and.(NaSGrp.ne.0)) then
               write (6,*)
               write (6,*) 'Small segmentation must be specified'
               write (6,*) 'with large segmentation, or both can'
               write (6,*) 'be left unspecified'
               call abend()
           end if

           if (NaGrp.ne.0) then
c             large != 0, small == 0 => small = 1
              if (NaSGrp.eq.0) then
                 NaSGrp = 1
              endif

c             large * small <= 64
              if ((NaGrp*NaSGrp).gt.maxSGrp) then
                write (6,*)
                write (6,*) 'Product of Large and Small segmen-'
                write (6,*) 'tation must be less or equal to 64'
                call abend()
              endif
           end if

       ELSE IF (LINE(1:4).EQ.'CHSE') THEN
       Read(LuSpool,*) NchBlk_tmp
           if ((NchBlk_tmp.lt.1).or.(NchBlk_tmp.gt.NChLoc_min)) then
               write (6,*)
               write (6,*) 'Ilegal value for CHSegment keyword  : ',
     &                      NchBlk_tmp
               write (6,*) 'Reseting to a reasonable value for    '
               write (6,*) 'this system :                         ',
     & NchBlk
           else if (int(NChLoc_max/NchBlk_tmp).ge.100) then
               write (6,*) 'Number of block of the MO Cholesky vector'
               write (6,*) 'exceeded the limit. Increasing value of  '
               write (6,*) 'the CHSEgmentation keyword to : ',
     & NchBlk
           else
              NchBlk = NchBlk_tmp
           end if

cmp!       ELSE IF (LINE(1:4).EQ.'LUNA') THEN  ... toto sa nikdy nevyuzivalo
cmp!       Read(LuSpool,*) LunAux

       ELSE IF (LINE(1:4).EQ.'MHKE') THEN
       Read(LuSpool,*) mhkey
           if ((mhkey.lt.0).or.(mhkey.gt.2)) then
              mhkey=1
              write(6,*)
              write(6,*) ' Warning!!!',
     &                   '  Matrix handling key out of range'
              write(6,*) ' parameter mhkey changed to 1'
           end if

       ELSE IF (LINE(1:4).EQ.'NOGE') THEN
            generkey = 0

       ELSE IF (LINE(1:4).EQ.'ONTH') THEN
            intkey1 = 1

       ELSE IF (LINE(1:4).EQ.'PREC') THEN
            intkey2 = 1

       ELSE IF (LINE(1:4).EQ.'NODI') THEN
            W34DistKey = 0

       ELSE IF (LINE(1:4).EQ.'JOIN') THEN
       Read(LuSpool,*) JoinLkey
            if ((JoinLkey.lt.0).or.(JoinLkey.gt.3)) then
               write (6,*)
               write (6,*) 'Ilegal value for Join keyword : ',
     &                      JoinLkey
               write (6,*) 'Use one of 0, 1, 2, 3'
               write (6,*) 'For details, see the manual ...'
               call abend()
            end if


       ELSE IF (LINE(1:4).EQ.'MAXI') THEN
       Read(LuSpool,*) maxiter
            if (maxiter.le.0) then
               write (6,*)
               write (6,*) 'Ilegal value of the MAXITER keyword: ',
     &                      maxiter
               write (6,*) 'Use integer > 0'
               call abend()
            end if

       ELSE IF (LINE(1:4).EQ.'REST') THEN
            restkey = 1
            write (6,*)
            write (6,*) 'This option is temporary disabled'
            write (6,*) 'No Restart possible (... yet).'
            call abend()

       ELSE IF (LINE(1:4).EQ.'THRE') THEN
       Read(LuSpool,*) conv

       ELSE IF (LINE(1:4).EQ.'PRIN') THEN
       Read(LuSpool,*) printkey
            if (((printkey.lt.0).or.(printkey.gt.10)).or.
     & ((printkey.gt.2).and.(printkey.lt.10))) then

               write (6,*)
               write (6,*) 'Ilegal value of the PRINT keyword: ',
     &                      printkey
               write (6,*) ' Use: 1  (Minimal) '
               write (6,*) '      2  (Minimal + Timings)'
               write (6,*) '      10 (Debug) '
               call abend()
            end if

       ELSE IF (LINE(1:4).EQ.'END ') THEN
       GOTO 7
       END IF
       GOTO 6
7      CONTINUE

       Call Close_LuSpool(LuSpool)

c! take care of the algorithm keyword
        if (intkey1.eq.intkey2) then
           if (intkey1.eq.0) then
              write (6,*)
              write (6,*) 'None of OnTheFly/PreCalculate'
              write (6,*) 'algorithm was selected. Using'
              write (6,*) 'default: PreCalculate (1)'
              intkey = 1
           else
              write (6,*)
              write (6,*) 'OnTheFly and PreCalculate keywords'
              write (6,*) 'are mutually exclusive'
              call abend()
           end if
        else
           if (intkey1.eq.1) then
               intkey = 0
           else
               intkey = 1
           end if
        end if

c2      tlac hlavicky
        write (6,*)
        write (6,*) '    Cholesky Based Closed-Shell CCSD code'
cmp!        write (6,*) ' Dedicated to the memory of Boris Jeltzin'
        write (6,*)
      write (6,*) '---------------------------------------------------'

        write (6,'(A,i9)') ' Frozen Orbitals                   : ',
     & nfr
        write (6,'(A,i9)') ' Occupied Orbitals                 : ',
     & no
        write (6,'(A,i9)') ' Virtual Orbitals                  : ',
     & nv
        write (6,'(A,i9)') ' Total number of Cholesky Vectors  : ',
     & nc

      write (6,*) '---------------------------------------------------'

        if (NaGrp.ne.0) then
          write (6,'(A,i9)') ' Large Virtual Segmentation        : ',
     & NaGrp
        else
          write (6,'(A,A9)') ' Large Virtual Segmentation        : ',
     & ' auto'
        end if

        if (NaSGrp.ne.0) then
          write (6,'(A,i9)') ' Small Virtual Segmentation        : ',
     & NaSGrp
        else
          write (6,'(A,A9)') ' Small Vectors Segmentation        : ',
     & ' auto'
        end if

        write (6,'(A,i9)') ' Cholesky Vectors Segmentation     : ',
     & NchBlk

      write (6,*) '---------------------------------------------------'

        msg = 'No'
        if (generkey.eq.1) msg = 'Yes'

        write (6,'(A,A4)') ' Generate Scratch Files?                : ',
     & msg
        write (6,'(A,i4)') ' Precalculate (1) / On-the-Fly (0) Alg. : ',
     & intkey
        write (6,'(A,i4)') ' 3 and 4-ext. MO integrals distribute?  : ',
     & W34DistKey
        write (6,'(A,i4)') ' Parallel Join of varios MO integrals   : ',
     & JoinLkey

      write (6,*) '---------------------------------------------------'

        write (6,'(A,E9.2)') ' Convergence Threshold             : ',
     & conv
        write (6,'(A,i9)') ' Maximum number of Iterations      : ',
     & maxiter

      write (6,*) '---------------------------------------------------'

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

      write (6,*) '---------------------------------------------------'
        write (6,*)

c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(wrksize)
        end
