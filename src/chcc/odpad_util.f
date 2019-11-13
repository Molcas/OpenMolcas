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
c       docasne utilitky
c       rea_chcc
c       wri_chcc
c       IniReord
c       RNFill
c       UrobTau
c       UrobChV
c       UrobInt
c       ReaW
c         rea1
c       GetX
c       SaveX
c       Exp1
c       Exp2
c       Exp2i
c       Exp4
c       UrobL1
c       UrobL2
c       UrobI1
c       UrobI2
c       UrobI3
c       UrobT2
c        TimeDelay
c
c       ----------------
c
        subroutine rea_chcc (lun,length,A)
c
c       nacitane bloku dat
c
        implicit none
        integer lun,length
        real*8 A(1:length)
c
        read (lun) A
c
        return
        end
c
c       ----------------
c
        subroutine wri_chcc (lun,length,A)
c
c       zapisanie bloku dat
c
        implicit none
        integer lun,length
        real*8 A(1:length)
c
        write (lun) A
c
        return
        end
c
c       ----------------
c
        subroutine IniReord(NaGrp,NaSGrp,NchBlk,LunAux,wrksize)
c
c       nacitanie vsupu a inicializacia premnennych
c       a tlac primitivnej hlavicky pre Reord procesz
c
        implicit none
#include "chcc1.fh"
#include "chcc_reord.fh"
cmp!
#include "chcc_parcc.fh"
#include "para_info.fh"
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
        character*80 TITLE

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
       TITLE=' '
 6     Read(LuSpool,'(A80)') LINE
       IF(LINE(1:1).EQ.'*') GOTO 6
       CALL UPCASE(LINE)
c
       IF (LINE(1:4).EQ.'TITL') THEN
       Read(LuSpool,'(A72)') TITLE

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
c
c --------------------------------
c
        subroutine RNFill (length,A,c)
c
c       fill an array with random numbers in interval (-c,c)
c
        implicit none
        integer length
        real*8 c
        real*8 A(1:length)
c
c       help variables
        integer i
c
c
        do i=1,length
c          A(i)=c*(srand()-0.5d0)
           A(i)=(1.0d-7)*i
        end do
c
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_real(c)
        end
c
c       ----------------
c
        subroutine UrobTau (Tau,NaGrp,LunTau)
c
c       vyraba subor LunTau so simulovanymi Tau amplitudami
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
c@@     include 'o2v4.fh' stare
#include "o3v3.fh"
c
        integer NaGrp,LunTau
        real*8 Tau(1)
c
c       help variables
        integer length
        integer aGrp,bGrp
c
c
c1      cycle over Groups a,b
c
        do aGrp=1,NaGrp
        do bGrp=1,aGrp
c
c1.1    def legth
c
        if (aGrp.eq.bGrp) then
c       groups of a and b are equal, reading for a'>=b'
          length=no*no*DimGrpv(aGrp)*(DimGrpv(bGrp)+1)/2
        else
c       aGrp>bGrp, reading for a',b' in given groups
          length=no*no*DimGrpv(aGrp)*DimGrpv(bGrp)
        end if
c
c
c1.2    full Tau with random numbers
c
        call RNFill (length,Tau(1),1.0d-2)
c
c1.3    read block
c
        write (6,*) aGrp,bGrp,length
        call wri_chcc (lunTau,length,Tau(1))
c
        end do
        end do
c
c
c2      rewind file
        rewind (LunTau)
c
        return
        end
c
c       ----------------
c
        subroutine UrobChV (L2,NaGrp,NbeGrp,LunAux)
c
c       vyraba fily so simulovanymi Cholesky vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
c@@     include 'o2v4.fh' stare
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,NbeGrp,LunAux
        real*8 L2(1)
c
c       help variables
        integer aGrp,beGrp,len
        real*8 schem
c
c1      cycle over a,be Groups
c
        do aGrp=1,NaGrp
        do beGrp=1,NbeGrp
c
c1.1      def length
          len=nc*DimGrpv(aGrp)*DimGrpv(beGrp)
c
c1.2      full L2 with random numbers
          schem=1.0d-2
          call RNFill (len,L2(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=L2Name(aGrp,beGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,L2Name(aGrp,beGrp))
c
c1.4      write L2 into proper file
          write (6,*) aGrp,beGrp,len
          call wri_chcc (LunAux,len,L2(1))
c
        end do
        end do
c
c
        return
        end
c
c       ----------------
c
        subroutine UrobInt (W,maxa,maxbe,n,lun)
c
c       vyraba subor LunInt so simulovanymi (a"be"|b"ga") integralmi
c       kde je za sebou N blokov
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer maxa,maxbe,n,lun
        real*8 W(1)
c
c       help variables
        integer length,i
c
c       open (unit=Lun,file='IntFil',form='unformatted')
c
c1      def legth
c
        length=maxa*maxa*maxbe*maxbe
c
c2      cycle over N
c
        do i=1,n
c
c2.1    full W with random numbers
c
        call RNFill (length,W(1),1.0d-2)
c
c2.2    write block
c
        write (6,*) 'Vint',i,length
        call wri_chcc (lun,length,W(1))
c
        end do
c
c
        rewind(lun)
c
        return
        end
c
c       -----------
c
        subroutine ReaW (W,aSGrp,beSGrp,bSGrp,gaSGrp,LunInt)
c
c       simuluje citanie VVVV integralov
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        real*8 W(1)
        integer aSGrp,beSGrp,bSGrp,gaSGrp,LunInt
c
c       help variables
        integer dim
c
        dim=DimSGrpa(aSGrp)*DimSGrpbe(beSGrp)*
     c      DimSGrpa(bSGrp)*DimSGrpbe(gaSGrp)
c
        call rea1 (LunInt,dim,W(1))
c
        return
        end
c
c       -----------
c
        subroutine rea1 (lun,length,A)
c
c       nacitane bloku dat z Sekvencneho suboru, ak koniec tak
c       rewind a citanie znova od zaciatku
c       Toto je ozaj ditry
c
        implicit none
        integer lun,length
        real*8 A(1:length)
c
        read (lun,end=99) A
        return
c
99      rewind(lun)
        read (lun) A
c
        end
c
c       -----------
c
        subroutine GetX (X,length,Lun,LunName,keyopen,keyclose)
c
c       this routine do
c       1) keyopen = 0 - nothing (i.e) file is opened
c                    1 - open LunName file with Lun
c                    2 - rewind Lun file
c                    3 - open LunName file with Lun with ACCESS='append'
c       2) read X  of dimension length
c       3) keyclose= 0 - nothing
c                    1 - close Lun file
c
c
        implicit none
        integer length,Lun,keyopen,keyclose
        real*8 X(1)
        character*6 LunName
c
c1
        if (keyopen.eq.1) then
*         open (unit=Lun,file=LunName,form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
        else if (keyopen.eq.2) then
          rewind(Lun)
        else if (keyopen.eq.3) then

cmp!          open (unit=Lun,file=LunName,form='unformatted',
cmp!     c          ACCESS='append')

          Call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
          call append_file_u(Lun)

        end if
c
c2
        call rea_chcc (Lun,length,X(1))
c
c3
        if (keyclose.eq.1) then
          close (Lun)
        end if
c
        return
        end
c
c       -----------
c
        subroutine SaveX (X,length,Lun,LunName,keyopen,keyclose)
c
c       this routine do
c       1) keyopen = 1 - open LunName file with Lun
c                    2 - rewind Lun file
c                    3 - open LunName file with Lun with ACCESS='append'
c                 else - nothing (i.e) file is opened
c       2) write X  of dimension length
c       3) keyclose= 1 - close Lun file
c                 else - nothing
c
c
        implicit none
        integer length,Lun,keyopen,keyclose
        real*8 X(1)
        character*6 LunName
c
c1
        if (keyopen.eq.1) then
*         open (unit=Lun,file=LunName,form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
        else if (keyopen.eq.2) then
          rewind(Lun)
        else if (keyopen.eq.3) then
cmp!          open (unit=Lun,file=LunName,form='unformatted',
cmp!     c          ACCESS='append')

          Call MOLCAS_BinaryOpen_Vanilla(Lun,LunName)
          call append_file_u(Lun)

        end if
c
c2
        call wri_chcc (Lun,length,X(1))
c
c3
        if (keyclose.eq.1) then
          close (Lun)
        end if
c
        return
        end
c
c       -----------
c
        subroutine Exp1 (A,B,dimi,dim2,dime)
c
c       this routine do:
c       expand A(i,pq) -> B(i,p,q)
c       for matrices, where A(i,pq)=A(i,qp)
c
c       parameter description:
c       A       - input matrix (I)
c       B       - outpun matrix (O)
c       dimi    - first dimension in A,B (I)
c       dim2    - # of pq (I)
c       dimf    - # of p (also q) (I)
c
        implicit none
        integer dimi,dim2,dime
        real*8 A(1:dimi,1:dim2)
        real*8 B(1:dimi,1:dime,1:dime)
c
c       help variables
        integer i,p,q,pq
c
        pq=0
        do p=1,dime
        do q=1,p
        pq=pq+1
c
          do i=1,dimi
            B(i,p,q)=A(i,pq)
          end do
c
          do i=1,dimi
            B(i,q,p)=A(i,pq)
          end do
c
        end do
        end do
c
        return
        end
c
c       -----------
c
        subroutine Exp2 (A,B,dima,dimb,dim2,dime)
c
c       this routine do:
c       expand A(a,b,pq) -> B(a,b,p,q)
c       for matrices, where A(a,b,pq)=A(a,b,qp)
c
c       parameter description:
c       A       - input matrix (I)
c       B       - outpun matrix (O)
c       dima,b  - first dimension in A,B (I)
c       dim2    - # of pq (I)
c       dimf    - # of p (also q) (I)
c
        implicit none
        integer dima,dimb,dim2,dime
        real*8 A(1:dima,1:dimb,1:dim2)
        real*8 B(1:dima,1:dimb,1:dime,1:dime)
c
c       help variables
        integer a1,b1,p,q,pq
c
        pq=0
        do p=1,dime
        do q=1,p
        pq=pq+1
c
          do b1=1,dimb
          do a1=1,dima
            B(a1,b1,p,q)=A(a1,b1,pq)
          end do
          end do
c
          do b1=1,dimb
          do a1=1,dima
            B(a1,b1,q,p)=A(a1,b1,pq)
          end do
          end do
c
        end do
        end do
c
        return
        end
c
c       -----------
c
        subroutine Exp2i (A,B,dima,dimb,dim2,dime)
c
c       this routine do:
c       expand A(a,b,pq) -> B(b,a,p,q)
c       for matrices, where A(a,b,pq)=A(a,b,qp)
c
c       parameter description:
c       A       - input matrix (I)
c       B       - outpun matrix (O)
c       dima,b  - first dimension in A,B (I)
c       dim2    - # of pq (I)
c       dimf    - # of p (also q) (I)
c
        implicit none
        integer dima,dimb,dim2,dime
        real*8 A(1:dima,1:dimb,1:dim2)
        real*8 B(1:dimb,1:dima,1:dime,1:dime)
c
c       help variables
        integer a1,b1,p,q,pq
c
        pq=0
        do p=1,dime
        do q=1,p
        pq=pq+1
c
          do a1=1,dima
          do b1=1,dimb
            B(b1,a1,p,q)=A(a1,b1,pq)
          end do
          end do
c
          do a1=1,dima
          do b1=1,dimb
            B(b1,a1,q,p)=A(a1,b1,pq)
          end do
          end do
c
        end do
        end do
c
        return
        end
c
c       -----------
c
        subroutine Exp4 (A,B,dima2,dima,dimp2,dimp)
c
c       this routine do:
c       expand A(ab,pq) -> B(b,a,p,q)
c       for matrices, where A(ab,pq)=A(ab,qp)=A(ba,pq)=A(ba,qp)
c        N.B. kvajto odflaknute
c
c       parameter description:
c       A       - input matrix (I)
c       B       - outpun matrix (O)
c
        implicit none
        integer dima2,dima,dimp2,dimp
        real*8 A(1:dima2,1:dimp2)
        real*8 B(1:dima,1:dima,1:dimp,1:dimp)
c
c       help variables
        integer a1,b1,ab,p,q,pq
c
        pq=0
        do p=1,dimp
        do q=1,p
        pq=pq+1
c
          ab=0
          do a1=1,dima
          do b1=1,a1
          ab=ab+1
            B(a1,b1,p,q)=A(ab,pq)
            B(a1,b1,q,p)=A(ab,pq)
            B(b1,a1,p,q)=A(ab,pq)
            B(b1,a1,q,p)=A(ab,pq)
          end do
          end do
c
        end do
        end do
c
        return
        end
c
c       ----------------
c
        subroutine UrobL1 (L1,NaGrp,LunAux)
c
c       vyraba fily so simulovanymi L1       vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,LunAux
        real*8 L1(1)
c
c       help variables
        integer aGrp,len
        real*8 schem
c
c1      cycle over a,be Groups
c
        do aGrp=1,NaGrp
c
c1.1      def length
          len=nc*DimGrpv(aGrp)*no
c
c1.2      full L1 with random numbers
          schem=1.0d-2
          call RNFill (len,L1(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=L1Name(aGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,L1Name(aGrp))
c
c1.4      write L1 into proper file
          write (6,*) aGrp,len
          call wri_chcc (LunAux,len,L1(1))
c
        close (LunAux)
c
        end do
c
c
        return
        end
c
c       ----------------
c
        subroutine UrobL2 (L2,NaGrp,NbeGrp,LunAux)
c
c       vyraba fily so simulovanymi Cholesky vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,NbeGrp,LunAux
        real*8 L2(1)
c
c       help variables
        integer aGrp,beGrp,len
        real*8 schem
c
c1      cycle over a,be Groups
c
        do aGrp=1,NaGrp
        do beGrp=1,NbeGrp
c
c1.1      def length
          len=nc*DimGrpv(aGrp)*DimGrpv(beGrp)
c
c1.2      full L2 with random numbers
          schem=1.0d-2
          call RNFill (len,L2(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=L2Name(aGrp,beGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,L2Name(aGrp,beGrp))
c
c1.4      write L2 into proper file
          write (6,*) aGrp,beGrp,len
          call wri_chcc (LunAux,len,L2(1))
c
          close (LunAux)
c
        end do
        end do
c
c
        return
        end
c
c       ----------------
c
        subroutine UrobI1 (I1,NaGrp,LunAux)
c
c       vyraba fily so simulovanymi I1       vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,LunAux
        real*8 I1(1)
c
c       help variables
        integer aGrp,len
        real*8 schem
c
c1      cycle over a,be Groups
c
        do aGrp=1,NaGrp
c
c1.1      def length
          len=no*DimGrpv(aGrp)*no*(no+1)/2
c
c1.2      full I1 with random numbers
          schem=1.0d-2
          call RNFill (len,I1(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=I1Name(aGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,I1Name(aGrp))
c
c1.4      write I1 into proper file
          write (6,*) aGrp,len
          call wri_chcc (LunAux,len,I1(1))
c
        close (LunAux)
c
        end do
c
c
        return
        end
c
c       ----------------
c
        subroutine UrobI2 (I2,NaGrp,NbeGrp,LunAux)
c
c       vyraba fily so simulovanymi I2       vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,NbeGrp,LunAux
        real*8 I2(1)
c
c       help variables
        integer aGrp,beGrp,len
        real*8 schem
c
c1      cycle over a,be Groups
c
        do aGrp=1,NaGrp
        do beGrp=1,NbeGrp
c
c1.1      def length
          len=no*no*DimGrpv(aGrp)*DimGrpv(beGrp)
c
c1.2      full I2 with random numbers
          schem=1.0d-2
          call RNFill (len,I2(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=I2Name(aGrp,beGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,I2Name(aGrp,beGrp))
c
c1.4      write I2 into proper file
          write (6,*) aGrp,beGrp,len
          call wri_chcc (LunAux,len,I2(1))
c
          close (LunAux)
c
        end do
        end do
c
c
        return
        end
c
c       ----------------
c
        subroutine UrobI3 (I3,NaGrp,NbeGrp,LunAux)
c
c       vyraba fily so simulovanymi I3       vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,NbeGrp,LunAux
        real*8 I3(1)
c
c       help variables
        integer aGrp,beGrp,len
        real*8 schem
c
c1      cycle over a,be Groups
c
        do aGrp=1,NaGrp
        do beGrp=1,NbeGrp
c
c1.1      def length
          if (aGrp.eq.beGrp) then
          len=no*(no+1)*DimGrpv(aGrp)*(DimGrpv(beGrp)+1)/4
          else
          len=no*(no+1)*DimGrpv(aGrp)*DimGrpv(beGrp)/2
          end if
c
c1.2      full I3 with random numbers
          schem=1.0d-2
          call RNFill (len,I3(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=I3Name(aGrp,beGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,I3Name(aGrp,beGrp))
c
c1.4      write I3 into proper file
          write (6,*) aGrp,beGrp,len
          call wri_chcc (LunAux,len,I3(1))
c
          close (LunAux)
c
        end do
        end do
c
c
        return
        end
c
c       ----------------
c
        subroutine UrobT2 (T2,NaGrp,NbeGrp,LunAux)
c
c       vyraba fily so simulovanymi T2       vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,NbeGrp,LunAux
        real*8 T2(1)
c
c       help variables
        integer aGrp,beGrp,len
        real*8 schem
c
c1      cycle over a,be Groups
c
        do aGrp=1,NaGrp
        do beGrp=1,NbeGrp
c
c1.1      def length
          len=no*(no+1)*DimGrpv(aGrp)*DimGrpv(beGrp)/2
c
c1.2      full T2 with random numbers
          schem=1.0d-2
          call RNFill (len,T2(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=T2Name(aGrp,beGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,T2Name(aGrp,beGrp))
c
c1.4      write T2 into proper file
          write (6,*) aGrp,beGrp,len
          call wri_chcc (LunAux,len,T2(1))
c
          close (LunAux)
c
        end do
        end do
c
c
        return
        end
c
c       ----------------
c
cmp!        subroutine TimeDelay (LunAux)
cmp!c
cmp!c        this routine do:
cmp!c        only on myRank=0
cmp!c        1) Read nTime, keyRep from ClcFil
cmp!c        2) wait nTime (sec)
cmp!c        3) for keyRep - 1 - goto 1
cmp!c                       0 - continue
cmp!c        on all nodes
cmp!c        4) allreduce Blnk (useless) in order to catch
cmp!c           time delay on all nodes
cmp!c
cmp!c
cmp!        implicit none
cmp!        integer LunAux
cmp!        include 'para_info.fh'
cmp!        real*8 Blnk(1)
cmp!        integer nTime,keyRep
cmp!        logical ClcFilExist
cmp!c
cmp!c1        Delaying on master node
cmp!c
cmp!        if (myRank.eq.0) then
cmp!          inquire (file='ClcFil',exist=ClcFilExist)
cmp!          if (ClcFilExist) then
cmp!1            open (unit=LunAux,file='ClcFil')
cmp!            read (LunAux,*) nTime,keyRep
cmp!            close (LunAux)
cmp!c
cmp!            write (6,*) ' Sleep for :',nTime, keyRep
cmp!            call sleep (nTime)
cmp!            if (keyRep.eq.1) then
cmp!                  goto 1
cmp!            end if
cmp!          end if
cmp!        end if
cmp!c
cmp!c2        Defining Blnk (useless)
cmp!c
cmp!        Blnk(1)=1.0d0*myRank
cmp!c
cmp!c3        Allreduce Blnk (i.e. waiting on all nodes)
cmp!c
cmp!#ifdef _MOLCAS_MPP_
cmp!c##     Synchronizacny bod:
cmp!c       Allreduce Blank
cmp!        call gadgop (Blnk(1),1,'+')
cmp!#endif
cmp!c
cmp!        return
cmp!        end
