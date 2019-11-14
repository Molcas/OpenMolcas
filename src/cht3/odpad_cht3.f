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
c DefParReord_t3
c DefParReordHlp1
c DefParReordHlp2
c IniReord_t3
c GetX_t3
c cht3_rea
c gather_t2
c grow_t2eq
c grow_t2neq
c grow_l1
c calc_MP2
c generate_juzekOE
c defcommon
c
c
c
c
c
c
c
c
c       ------------------------------------
c
        subroutine DefParReord_t3 (NaGrpR,maxdim)

c
c       This routine do:
c       define parameters in cht3_reord.fh using NaGrpR,maxdim
c
c       I/O parameter description:
c       NxGrpR   - # of groups in a (=b) set (I)
c       maxdim   - # maximal dimension of a (=b) Groups(O)
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
c
        integer NaGrpR,maxdim
c
c       help variables
c
        real*8 rdim
        integer i,j
        integer Up(1:MaxGrp),Low(1:MaxGrp)
c
c
c1      define parameters of Groups of a set
c
        rdim=1.0d0*nv/(1.0d0*NaGrpR)
c
        do i=1,NaGrpR
c
           if (i.eq.1) then
             Up(i)=int(rdim*i)
             Low(i)=1
           else if (i.eq.NaGrpR) then
             Up(i)=nv
             Low(i)=Up(i-1)+1
           else
             Up(i)=int(rdim*i)
             Low(i)=Up(i-1)+1
           end if
c
           DimGrpaR(i)=(Up(i)-Low(i))+1
c
        end do
c
c
c2      find maximal dimensions of a'
c
        maxdim=DimGrpaR(1)
        do i=1,NaGrpR
          if (DimGrpaR(i).gt.maxdim) then
          maxdim=DimGrpaR(i)
          end if
        end do
c
c
c3.1    def L2Name, T2Name, I2Name,I3Name
c
        do i=1,MaxGrp
        do j=1,MaxGrp
          call DefParReordHlp1(i,j,'L2',L2Name(i,j))
          call DefParReordHlp1(i,j,'T2',T2Name(i,j))
          call DefParReordHlp1(i,j,'I2',I2Name(i,j))
          call DefParReordHlp1(i,j,'I3',I3Name(i,j))
        end do
        end do
c
c3.2    def L1Name,I1Name
c
        do i=1,MaxGrp
          call DefParReordHlp2 (i,'L1vc',L1Name(i))
          call DefParReordHlp2 (i,'I1in',I1Name(i))
        end do
c
c3.3    def L0Name,I0Name
c
        L0Name='L0vctr'
        I0Name='I0intg'
c
        return
        end
c
c       ---------------------
c
        subroutine DefParReordHlp1 (i,j,Schem,Nomen)
c
c       help routine to DefPar..., producing names of Disc files
c       ex: Schem='XY', i=1, j=3  ->  Nomen='XY0103'
c        N.B. suspendovana rutina
c
        implicit none
        integer i,j
        character*2 Schem
        character*6 Nomen
c
c       help variables
        character*1 Chr(1:6)
        character*2 digit(1:64)
        character*2 ichr,jchr
        character*2 baza
        character*6 meno
c
        equivalence (Chr(1),meno)
        equivalence (Chr(1),baza)
        equivalence (Chr(3),ichr)
        equivalence (Chr(5),jchr)
c
c
c        quite a porno this piece
        digit(1)='01'
        digit(2)='02'
        digit(3)='03'
        digit(4)='04'
        digit(5)='05'
        digit(6)='06'
        digit(7)='07'
        digit(8)='08'
        digit(9)='09'
        digit(10)='10'
        digit(11)='11'
        digit(12)='12'
        digit(13)='13'
        digit(14)='14'
        digit(15)='15'
        digit(16)='16'
        digit(17)='17'
        digit(18)='18'
        digit(19)='19'
        digit(20)='20'
        digit(21)='21'
        digit(22)='22'
        digit(23)='23'
        digit(24)='24'
        digit(25)='25'
        digit(26)='26'
        digit(27)='27'
        digit(28)='28'
        digit(29)='29'
        digit(30)='30'
        digit(31)='31'
        digit(32)='32'
        digit(33)='33'
        digit(34)='34'
        digit(35)='35'
        digit(36)='36'
        digit(37)='37'
        digit(38)='38'
        digit(39)='39'
        digit(40)='40'
        digit(41)='41'
        digit(42)='42'
        digit(43)='43'
        digit(44)='44'
        digit(45)='45'
        digit(46)='46'
        digit(47)='47'
        digit(48)='48'
        digit(49)='49'
        digit(50)='50'
        digit(51)='51'
        digit(52)='52'
        digit(53)='53'
        digit(54)='54'
        digit(55)='55'
        digit(56)='56'
        digit(57)='57'
        digit(58)='58'
        digit(59)='59'
        digit(60)='60'
        digit(61)='61'
        digit(62)='62'
        digit(63)='63'
        digit(64)='64'
c
c
        baza=Schem
        ichr=digit(i)
        jchr=digit(j)
        Nomen=meno
c
        return
        end
c
c       ---------------------
c
        subroutine DefParReordHlp2 (i,Schem,Nomen)
c
c       help routine to DefParo2v4, producing names of Disc files
c       ex: Schem='XYZQ', i=1  ->  Nomen='XYZQ01'
c
        implicit none
        integer i
        character*4 Schem
        character*6 Nomen
c
c       help variables
        character*1 Chr(1:6)
        character*2 digit(1:64)
        character*2 ichr
        character*4 baza
        character*6 meno
c
        equivalence (Chr(1),meno)
        equivalence (Chr(1),baza)
        equivalence (Chr(5),ichr)
c
c
        digit(1)='01'
        digit(2)='02'
        digit(3)='03'
        digit(4)='04'
        digit(5)='05'
        digit(6)='06'
        digit(7)='07'
        digit(8)='08'
        digit(9)='09'
        digit(10)='10'
        digit(11)='11'
        digit(12)='12'
        digit(13)='13'
        digit(14)='14'
        digit(15)='15'
        digit(16)='16'
        digit(17)='17'
        digit(18)='18'
        digit(19)='19'
        digit(20)='20'
        digit(21)='21'
        digit(22)='22'
        digit(23)='23'
        digit(24)='24'
        digit(25)='25'
        digit(26)='26'
        digit(27)='27'
        digit(28)='28'
        digit(29)='29'
        digit(30)='30'
        digit(31)='31'
        digit(32)='32'
        digit(33)='33'
        digit(34)='34'
        digit(35)='35'
        digit(36)='36'
        digit(37)='37'
        digit(38)='38'
        digit(39)='39'
        digit(40)='40'
        digit(41)='41'
        digit(42)='42'
        digit(43)='43'
        digit(44)='44'
        digit(45)='45'
        digit(46)='46'
        digit(47)='47'
        digit(48)='48'
        digit(49)='49'
        digit(50)='50'
        digit(51)='51'
        digit(52)='52'
        digit(53)='53'
        digit(54)='54'
        digit(55)='55'
        digit(56)='56'
        digit(57)='57'
        digit(58)='58'
        digit(59)='59'
        digit(60)='60'
        digit(61)='61'
        digit(62)='62'
        digit(63)='63'
        digit(64)='64'
c
        baza=Schem
        ichr=digit(i)
        Nomen=meno
c
        return
        end
c
c       ----------------
c
        subroutine IniReord_t3(NaGrp,wrksize)
c
c       nacitanie vsupu a inicializacia premnennych
c       a tlac primitivnej hlavicky pre Reord procesz
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
cmp
#include "cholesky.fh"
#include "para_info.fh"
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
        character*80 TITLE
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
       TITLE=' '
 6     Read(LuSpool,'(A80)') LINE
       IF(LINE(1:1).EQ.'*') GOTO 6
       CALL UPCASE(LINE)
c
       IF (LINE(1:4).EQ.'TITL') THEN
       Read(LuSpool,'(A72)') TITLE

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
c
c       -----------
c
        subroutine GetX_t3 (X,length,Lun,LunName,keyopen,keyclose)
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
        call cht3_rea (Lun,length,X(1))
c
c3
        if (keyclose.eq.1) then
          close (Lun)
        end if
c
        return
        end
c
c       ----------------
c
        subroutine cht3_rea (lun,length,A)
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
c       -----------
c

        subroutine gather_t2(t2,t2_tmp,tmp)
c
c temporary routine. In future T2 block structure will be merged
c                    with the block structure of the (T) code
c
c this routine do :
c
c cycle through T2XY files and gather them into on T2 array
c T2(nv_beta,nv_alpha,no_beta,no_alpha)
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        integer a,b,dima,dimb
cmp        integer lenght
        integer lenght
        integer lasta,lastb
c
        real*8 t2(*),tmp(*),t2_tmp(*)
        integer a_tmp,b_tmp
c
cmp        write (6,*)
cmp        write (6,*) '------ DimGrpaR ------'
cmp        write (6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
cmp        write (6,*)
c
        do a=1,NvGrp
        do b=1,a
c
cmp@@        dima=nv/NvGrp
        dima=DimGrpaR(a)
        dimb=DimGrpaR(b)
c
cmp        write (6,'(A,i3,i3,2x,A6)') 'a,b,T2Name(a,b) ',a,b,T2Name(a,b)
c
         if (a.eq.b) then  ! a=b
c open the pertinent file
c
cmp@@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
        dimb=dima
c
cmp         write (6,*) 'dima = ',dima
         lenght=(dima*(dima+1)*no*no)/2
cmp         write (6,*) 'lenght = ',lenght
cmp!         write (6,*) 'file size (g77) = ',16+lenght*8
cmp         write (6,*) 'file size (ifort) = ',8+lenght*8
c
        call GetX_t3 (tmp,lenght,LunAux,T2Name(a,b),1,1)
c
         else ! a>b
c open the pertinent file
c
cmp@@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
cmp@@         if (b.eq.NvGrp) dimb=nv-((NvGrp-1)*dimb)
c
cmp         write (6,*) 'dima = ',dima
cmp         write (6,*) 'dimb = ',dimb
         lenght=dima*dimb*no*no
cmp         write (6,*) 'lenght = ',lenght
cmp!         write (6,*) 'file size (g77) = ',16+lenght*8
cmp         write (6,*) 'file size (ifort) = ',8+lenght*8
c
        call GetX_t3 (tmp,lenght,LunAux,T2Name(a,b),1,1)
c
         end if
c
c add its contents to the t2 array
c
cmp@@        lasta=(a-1)*(nv/NvGrp)
cmp@@        lastb=(b-1)*(nv/NvGrp)
cmp@@
          lasta=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.1) then
          do b_tmp=1,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
cmp@@
c
cmp        write (6,'(A,2(i5,2x),2(i3,x))') 'lasta, lastb, a, b = ',
cmp     & lasta,lastb,a,b
c
        if (a.eq.b) then ! expand and map
c expand and map l2_1 (a',b',m) <- tmp (m,ab')
        call expand4_12 (tmp,t2_tmp,dima,no,no)
        call grow_t2neq(t2,t2_tmp,dima,dimb,nv,no,
     & lasta,lastb)
        else
        call grow_t2neq(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb)
        end if
c
        end do
        end do
c
        return
        end
c
c ---------------
c
        subroutine grow_t2neq(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb)
c
c this routine do :
c
c grow amplitude file t2(a,b,i,j) by the segment in tmp
c for case sa != sb
c
        implicit none
        integer a,b,dima,dimb,nv,no,i,j
        integer lasta,lastb
        real*8 t2(1:nv,1:nv,1:no,1:no)
        real*8 tmp(1:dima,1:dimb,1:no,1:no)
c
cmp        write (6,*) 'grow_t2neq dima , dimb  ',dima,dimb
cmp        write (6,*) 'grow_t2neq lasta, lastb ',lasta,lastb
cmp        write (6,*) 'grow_t2neq no           ',no
c
c?        if (lasta.eq.lastb) then
c?        do j=1,no
c?        do i=1,no
c?        do a=1,dima
c?        do b=1,a
c?        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
c?cmp        if (a.ne.b) t2(lastb+b,lasta+a,j,i)=-1.0d0*tmp(a,b,j,i)
c?        if (a.ne.b) t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
c?        end do
c?        end do
c?        end do
c?        end do
c
c?        else
c
        do j=1,no
        do i=1,no
        do b=1,dimb
        do a=1,dima
        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
cmp        t2(lastb+b,lasta+a,j,i)=-1.0d0*tmp(b,a,j,i)
        t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
        end do
        end do
        end do
        end do
c
c?        end if
c
        return
        end
c
c ---------------
c
        subroutine grow_l1(l1,tmp,dima,nc,no,nv,last)
c
c this routine do :
c
c grow Cholesky vectors L1(m,i,a) by the segment in tmp
c
        implicit none
        integer a,dima,nc,nv,no,i,m,last
        real*8 l1(1:nc,1:no,1:nv)
        real*8 tmp(1:nc,1:no,1:dima)
c
cmp        write (6,*) 'grow_l1i ',dima
c
        do a=1,dima
        do i=1,no
        do m=1,nc
        l1(m,i,last+a)=tmp(m,i,a)
        end do
        end do
        end do
c
c
        return
        end
c
c --------------------
c
        subroutine calc_MP2 (w,e,no,nv)
c
c this is primitive checking routine to calculate 2nd order energy
c
        implicit none
        integer i,j,a,b,no,nv
        real*8 e(1:(no+nv)),w(1:nv,1:no,1:nv,1:no)
        real*8 e2,integral,denom
c
        e2=0.0d0
c
        do j=1,no
        do i=1,no
        do b=1,nv
        do a=1,nv
c
        denom=e(no+a)+e(no+b)-e(i)-e(j)
cmp!        write (6,'(4(i3,2x),A,3(f17.10,2x))') a,i,b,j,'w1, w2, denom ',
cmp!     & w(a,i,b,j),w(a,j,b,i),denom

        integral=(-1.0d0)*w(a,i,b,j)*(2.0d0*w(a,i,b,j)+
     & (-1.0d0)*w(a,j,b,i))

c!      write (6,*) integral

        e2=e2+(integral/denom)
c
        end do
        end do
        end do
        end do
c
        write (6,*) 'Druhy rad je asi = ',e2
c
        return
        end
c
c ---------
c
        subroutine generate_juzekOE (oe,oeh,oep,no,nv)
c
c this routine do :
c
c modifies standar oe record to one used by DIRCC
c
c oeh = ( oe_occ(alpha) ... oe_occ(beta) )   (alpha=beta for this
c                                             implementation)
c
c oep = ( oe_virt(alpha) ... oe_virt(beta) ) (alpha=beta for this
c                                             implementation)
        implicit none
        integer i,no,nv
c
        real*8 oe(*),oeh(*),oep(*)
c
c
        do i=1,no
        oeh(i)=oe(i)
        oeh(i+no)=oe(i)
        end do
c
        do i=1,nv
        oep(i)=oe(no+i)
        oep(i+nv)=oe(no+i)
        end do
c
        return
        end
c
c ---------------
c
        subroutine defcommon (nfr,no,nv)
c
c this routine do :
c
c define commons needed in DIRCC routines
c
        implicit none
        integer nfr,nv,no
        integer NOAB(2),NNOAB(3),NUAB(2),NNUAB(3)
        character*1 ich(3)
        integer IT,ITLAST,NBF,NOMX,NU,NOO,MX2,NNO,NNU,NUO,NSO
        integer me,nprocs
c
        COMMON/UHF/NOAB,NNOAB,NUAB,NNUAB,ICH
        COMMON/PARAM/IT,ITLAST,NBF,NOMX,NU,MX2,NNO,NNU,NUO,NSO
        common /my_mpi_world_com/ me, nprocs
c
cmp!        include 'task_info_inc'
cmp!        include 'ws_conn_inc'
c
c       logical llmpi
c
        integer IOPT
        COMMON/IOIND/IOPT(96)
c
c ----  UHF -----
c
        noab(1)=no
        noab(2)=no
c
        nnoab(1)=noab(1)*(noab(1)-1)/2
        nnoab(2)=noab(2)*(noab(2)-1)/2
        nnoab(3)=noab(1)*noab(2)
c
        nuab(1)=nv
        nuab(2)=nv
c
        nnuab(1)=(nuab(1)*(nuab(1)-1))/2
        nnuab(2)=(nuab(2)*(nuab(2)-1))/2
        nnuab(3)=nuab(1)*nuab(2)
c
        ich(1)="A"
        ich(2)="B"
        ich(3)="C"
c
c ----  PARAM ----
c
        nso=0
c??????????????????
        it=1
        itlast=1
        nbf=nfr+no+nv
        nomx=nbf
        noo=no
        nu=nv
        mx2=(nbf*(nbf+1))/2
        nno=(no*(no+1))/2
        nnu=(nu*(nu+1))/2
        nuo=no*nu
c
c ------ my_mpi_world_com --------
c
c zatial pre sekvencny chod
c
cmp!      me=0
cmp!      nprocs=1
cmp!      llmpi=.false.
cmp!      nws=1
cmp!      iws(1)=1
cmp!      lws(1)=.true.
c
c ------ pre Get3DM -----------
c
c ------    IOPT    -----------
c
        IOPT(14)=6
        IOPT(30)=0
        IOPT(76)=0
        IOPT(93)=0
        IOPT(93)=64
        IOPT(95)=0
c
C  sets IOPT(27) to an extreme number to force one file - temporary !!!  preskumat !!!
C  in prder to be compatible with RHF i=2^31-1
      iopt(27)=2147483647
c
c
c
c zatial tolko
c
        return
        end
c
c ---------------
c
        subroutine gen_vvvo(occ_ind,w3,l1_1,l2_1,tmp)
c
c this routine do
c
c regenerate VVVo integrals from cholesky vectors
c
c -------------------
c
c structure of the cholesky vector files :
c
c       L1(m,I ,A') L1vcxx xx - Group of A'
c
c       L2(m,A'B')  L2xxyy xx - Group of A', A'>=B'
c                          yy - Group of B'
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        integer a,b,c,dima,dimb,dimc,occ_ind
        integer lenght
        integer lasta,lastb,lastc
c
        real*8 w3(1:(nv*(nv+1))/2,1:nv)
        real*8 tmp(*),l1_1(*),l2_1(*)
c
        integer a_tmp,b_tmp,c_tmp
c
c algoritmus je dobry ak maxdim > no
c inak treba vymenit citanie L1 za L2
c
c dalo by sa to urobit podstatne lepsie, kedby
c dircc nevyzadoval VVV ako (ab,c) ale ako (a,b,c)
c
c mozno urob sort L1i (m,c') <- L1(m,i,c')
c ---
c
c1         loop over a'
c
        do a=1,NvGrp
c
c2        loop over b'
c
        do b=1,a
c
c2.1        read L2(m,a',b')
c
         if (a.eq.b) then  ! a=b
c open the pertinent file
c
cmp@        dima=nv/NvGrp
cmp@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
        dima=DimGrpaR(a)
c
cmp!         write (6,*) 'dima = ',dima
         dimb=dima
         lenght=(dima*(dima+1)*nc)/2
cmp!         write (6,*) 'lenght L2Name(a,b) = ',L2Name(a,b),lenght
cmp!     write (6,*) 'file size (g77) = ',16+lenght*8
cmp!         write (6,*) 'file size L2Name(a,b) (ifort) = ',
cmp!     & L2Name(a,b),8+lenght*8
c
        call GetX_t3 (tmp,lenght,LunAux,L2Name(a,b),1,1)
c
         else ! a>b
c open the pertinent file
c
cmp@@        dima=nv/NvGrp
cmp@@        dimb=nv/NvGrp
cmp@@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
cmp@@         if (b.eq.NvGrp) dimb=nv-((NvGrp-1)*dimb)
        dima=DimGrpaR(a)
        dimb=DimGrpaR(b)
c
cmp!         write (6,*) 'dima, dimb = ',dima,dimb
         lenght=dima*dimb*nc
cmp!         write (6,*) 'lenght L2Name(a,b) = ',L2Name(a,b),lenght
cmp!     write (6,*) 'file size (g77) = ',16+lenght*8
cmp!         write (6,*) 'file size L2Name(a,b) (ifort) = ',
cmp!     & L2Name(a,b),8+lenght*8
c
        call GetX_t3 (tmp,lenght,LunAux,L2Name(a,b),1,1)
        end if
c
c2.2        map  L2_1(a',b',m) <- tmp(m,a',b')
c
        if (a.eq.b) then ! expand and map
c expand and map l2_1 (a',b',m) <- tmp (m,ab')
        call exMap3_231 (tmp,l2_1,nc,dima)
        else
        call Map3_231_t3 (tmp,l2_1,nc,dima,dimb)
        end if
c
c3         loop over c'
c
        do c=1,NvGrp
c
c3.1        read L1(m,i,c')
c
cmp@@        dimc=nv/NvGrp
cmp@@         if (c.eq.NvGrp) dimc=nv-((NvGrp-1)*dimc)
        dimc=DimGrpaR(c)
c
cmp!         write (6,*) 'dimc = ',dimc
         lenght=nc*no*dimc
cmp!         write (6,*) 'lenght L1Name(c) = ',L1Name(c),lenght
cmp!         write (6,*) 'file size L1Name(c) (ifort) = ',
cmp!     & L1Name(c),8+8*lenght
c
        call GetX_t3 (tmp,lenght,LunAux,L1Name(c),1,1)
c
c3.2        extract l1_1 (m,c')_i <- tmp (m,i,c')
c toto by sa dalo nahradit mapovanim
c
        call ext_o_32 (tmp,l1_1,nc,no,dimc,occ_ind)
c
c3.2.1        zero tmp
c
        call zeroma(tmp,1,dima*dimb*dimc)
c
c3.3         mult tmp (a',b',c') <- L2_1 (a',b',m) l1_1 (m,c')
c
               call mc0c1a3b
     & (dima*dimb,nc,nc,dimc,dima*dimb,dimc,
     & dima*dimb,nc,dimc,l2_1,l1_1,tmp)
C
c3.4        add W(ab,c) <- tmp (a'b'c')
c
cmp@@        lasta=(a-1)*(nv/NvGrp)
cmp@@        lastb=(b-1)*(nv/NvGrp)
cmp@@        lastc=(c-1)*(nv/NvGrp)
c
          lasta=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.1) then
          do b_tmp=1,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
c
          lastc=0
        if (c.gt.1) then
          do c_tmp=1,c-1
          lastc=lastc+DimGrpaR(c_tmp)
          end do
        end if
c
cmp!        write (6,'(A,3(i4),2x,3(i4))') 'lasta, lastb, lastc = ',
cmp!     & lasta,lastb,lastc,a,b,c

c sme v gen_vvvo
        call grow_w3 (w3,tmp,
     & nv,nv,dima,dimb,dimc,lasta,lastb,lastc)
c
c3.5        end loop over c'
        end do
c4        end loop over b'
        end do
c5        end loop over a'
        end do
c
        return
        end
c
c --------------------
c
        subroutine pack32_23 (AA,BB,d1,d2)
c
c this routine do :
c
c A(a,b,c) -> B(a,bc)  bc : b>=c,  dimb must eq dimc
c
        implicit none
        integer d1,d2,a,b,c,bc
        real*8 AA(d1,d2,d2),BB(d1,(d2*(d2+1)/2))
c
        bc=0
        do b=1,d2
        do c=1,b
        do a=1,d1
c
        bc=bc+1
        BB(a,bc)=AA(a,b,c)
        end do
        end do
        end do
c
        return
        end
c
c
c --------------------
c
        subroutine pack23_23 (AA,BB,d1,d2)
c
c this routine do :
c
c A(a,bc) -> B(a,b,c)  bc : b>=c,  dimb must eq dimc
c
        implicit none
        integer d1,d2,a,b,c,bc
        real*8 BB(d1,d2,d2),AA(d1,(d2*(d2+1)/2))
c
        bc=0
        do b=1,d2
        do c=1,b
        do a=1,d1
c
        bc=bc+1
        BB(a,b,c)=AA(a,bc)
        BB(a,c,b)=AA(a,bc)
        end do
        end do
        end do
c
        return
        end
c
c -----------------------
c
        subroutine exMap3_231 (A,B,d1,d2)
c
c this routine do :
c
c A (a,bc) -> B(b,c,a)
c
        implicit none
        integer d1,d2,i1,i2,i3,i23
        real*8 A(1:d1,1:(d2*(d2+1))/2)
        real*8 B(1:d2,1:d2,1:d1)
c
        i23=0
        do i2=1,d2
        do i3=1,i2
        i23=i23+1
        do i1=1,d1
c
        B(i2,i3,i1)=A(i1,i23)
        B(i3,i2,i1)=A(i1,i23)
c
        end do
        end do
        end do
c
        return
        end
c
c -------------------------
c
c3.2        extract l1_1 (m,c')_i <- tmp (m,i,c')
c
        subroutine ext_o_32 (A,B,nc,no,dima,occ_ind)
c
c this routine do :
c
c extract B (m,a')_i <- A (m,i,a')
c
        implicit none
        integer i1,i2,occ_ind,dima,nc,no
        real*8 A(1:nc,1:no,1:dima),B(1:nc,1:dima)
c
        do i2=1,dima
        do i1=1,nc
c
        B(i1,i2)=A(i1,occ_ind,i2)
c
        end do
        end do
c
        return
        end
c
c -----------------------------
c
        subroutine grow_w3_old (w3,AA,nv,d2,dima,dimb,dimc,
     & lasta,lastb,lastc)
c
c this routine do :
c
c add the block contribution AA(a',b',c') to w3(a>=b,c)
c
        implicit none
c
        integer a,b,c,dima,dimb,dimc,lasta,lastb,lastc,ab,nv
        integer a_point,b_point
        integer d2
        real*8 w3(1:(nv*(nv+1))/2,1:d2)
        real*8 AA(1:dima,1:dimb,1:dimc)
        integer a_old,b_old
c
        if ((dima.eq.0).or.(dimb.eq.0)) then
        write (6,*) 'dima, dimb = ',dima,dimb
        write (6,*) 'zle je'
        call abend()
        end if
c
        a_point=0
        b_point=0
        ab=0
        write (6,'(A,3(i5))') 'lasta, lastb, lastc = ',lasta,lastb,lastc
        write (6,'(A,2(i5))') 'dima, dimb          = ',dima,dimb
c
        a_old=0
        b_old=0
c
        do a=1,nv
        b_point=0
        do b=1,a
        ab=ab+1
        if ((a.ge.(lasta+1)).and.(a.le.(lasta+dima))) then
c
        if (a.ne.a_old) then
        a_point=a_point+1
        a_old=a
        end if
c
        if ((b.ge.max(1,lastb+1)).and.
     & (b.le.min(a,lastb+dimb))) then
c
!        write (6,*) 'b, b_old = ',b,b_old
        if ((b.ne.b_old).or.(b.eq.max(1,lastb+1))) then
!        write (6,*) 'wft'
        b_point=b_point+1
        b_old=b
        end if
c
c!        if (lastc.eq.0) write (6,'(A,5(i5))') 'ab, a, b, a_point, b_point = ',
c!     & ab,a,b,a_point,b_point
        do c=1,dimc
        w3(ab,lastc+c)=AA(a_point,b_point,c)
        end do
c
        end if
        end if
c
        end do
        end do
c
        return
        end
c
c -----------------------------
c
        subroutine grow_w3 (w3,AA,nv,d2,dima,dimb,dimc,
     & lasta,lastb,lastc)
c
c this routine do :
c
c add the block contribution AA(a',b',c') to w3(a>=b,c)
c
        implicit none
c
        integer a,b,c,dima,dimb,dimc,lasta,lastb,lastc,ab,nv
        integer a_point,b_point
        integer d2
        real*8 w3(1:(nv*(nv+1))/2,1:d2)
        real*8 AA(1:dima,1:dimb,1:dimc)
        integer a_old,b_old
c
        if ((dima.eq.0).or.(dimb.eq.0)) then
        write (6,*) 'dima, dimb = ',dima,dimb
        write (6,*) 'zle je'
        call abend()
        end if
c
        a_point=0
        b_point=0
        ab=0
cmp        write (6,'(A,3(i5))') 'lasta, lastb, lastc = ',lasta,lastb,lastc
cmp        write (6,'(A,2(i5))') 'dima, dimb          = ',dima,dimb
c
        a_old=0
        b_old=0
c
        do a=1,nv
        b_point=0
        do b=1,a
        ab=ab+1
        if ((a.ge.(lasta+1)).and.(a.le.(lasta+dima))) then
c
        if (a.ne.a_old) then
        a_point=a_point+1
        a_old=a
        end if
c
        if ((b.ge.max(1,lastb+1)).and.
     & (b.le.min(a,lastb+dimb))) then
c
!        write (6,*) 'b, b_old = ',b,b_old
        if ((b.ne.b_old).or.(b.eq.max(1,lastb+1))) then
!        write (6,*) 'wft'
        b_point=b_point+1
        b_old=b
        end if
c
cmp        if (lastc.eq.0) write (6,'(A,5(i5))') 'ab, a, b, a_point, b_point = ',
cmp     & ab,a,b,a_point,b_point
        do c=1,dimc
        w3(ab,lastc+c)=AA(a_point,b_point,c)
        end do
c
        end if
        end if
c
        end do
        end do
c
        return
        end
c
c --------------------
c
        subroutine pack32_12 (AA,BB,d1,d2)
c
c this routine do :
c
c A(a,b,c) -> B(ab,c)  bc : b>=c,  dimb must eq dimc
c
        implicit none
        integer d1,d2,a,b,c,ab
        real*8 AA(d1,d1,d2),BB(1:(d1*(d1+1)/2),1:d2)
c
        do c=1,d2
        ab=0
        do a=1,d1
        do b=1,a
        ab=ab+0
c
        BB(ab,c)=AA(a,b,c)
        end do
        end do
        end do
c
        return
        end
c
c --------------------
c
        subroutine pack43_34 (AA,BB,d1,d2,d3)
c
c this routine do :
c
c A(a,b,c,d) -> B(a,b,cd)  c>=d
c
        implicit none
        integer d1,d2,d3,a,b,c,d,cd
        real*8 AA(d1,d2,d3,d3),BB(d1,d2,(d3*(d3+1)/2))
c
        cd=0
        do c=1,d3
        do d=1,c
        cd=cd+1
        do b=1,d2
        do a=1,d1
c
        BB(a,b,cd)=AA(a,b,c,d)
        end do
        end do
        end do
        end do
c
        return
        end
c
c -------------------
c
        subroutine gen_oovo (w,l0,l1,tmp)
c
c this routine genetates (ij,a,k) integrals from
c blocked MO cholesky vectors
c
c --------
c
c       L0(m,IJ)    L0vctr  I>=J
c       L1(m,I ,A') L1vcxx xx - Group of A'
c

        implicit none
c
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        real*8 tmp(*),l0(*),l1(*),w(*)
        integer a,dima,lenght,last
c
        integer a_tmp
c
c1        read tmp(m,IJ)
c
        lenght=nc*(no*(no+1))/2
c
cmp!        write (6,'(A,A6)') 'L0vcrt ','L0vcrt'
cmp!        write (6,*) 'lenght = ',lenght
cmp!        write (6,*) 'file size (ifort) = ',8+8*lenght
c
        call GetX_t3 (tmp,lenght,LunAux,'L0vctr',1,1)
c2        map l0(IJ,m)   <- tmp(m,IJ)
c
        call Map2_21_t3 (tmp,l0,nc,(no*(no+1)/2))
c
c3        loop over A'
c
        do a=1,NvGrp
cmp@@        dima=nv/NvGrp
        dima=DimGrpaR(a)
c
c4        read tmp(m,I,A')
c
cmp!        write (6,'(A,i3,2x,A6)') 'a,L1Name(a) ',a,L1Name(a)
c
cmp@@         if (a.eq.NvGrp) dima=nv-((NvGrp-1)*dima)
c
cmp!         write (6,*) 'dima = ',dima
         lenght=nc*no*dima
cmp!         write (6,*) 'lenght = ',lenght
cmp!         write (6,*) 'file size (ifort) = ',8+8*lenght
c
        call GetX_t3 (tmp,lenght,LunAux,L1Name(a),1,1)
c
c5        grow l1(m,I,A)
c
cmp        last=(a-1)*(nv/NvGrp)
c
          last=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          last=last+DimGrpaR(a_tmp)
          end do
        end if
c
        call grow_l1(l1,tmp,dima,nc,no,nv,last)
c
c6        end loop over A'
c
        end do
c
c7        map tmp(m,A,I) <- l1(m,I,A)
c
        call Map3_132_t3 (l1,tmp,nc,no,nv)
c
c7.1        zero w
c
        call zeroma (w,1,((no*(no+1))/2)*nv*no)
c
c8        mult w(IJ,A,I)  <- l0(IJ,m) . tmp(m,A,I)
c
        call mc0c1a3b (
     & (no*(no+1))/2,nc,
     &  nc,nv*no,
     & (no*(no+1))/2,nv*no,
     & (no*(no+1))/2,nc,nv*no,l0,tmp,w)
c
        return
        end
c
c ----------------------
c
        subroutine grow_l2(A,B,nc,nv,dima,dimb,lasta,lastb)
c
c this routine do :
c
c grow A(A,B,m) from the blocked cholesky vectors
c B(a',b',m)
c
        implicit none
        integer i1,i2,i3,dima,dimb,nc,nv
        integer lasta,lastb
        real*8 A(nv,nv,nc),B(dima,dimb,nc)
c
        do i3=1,nc
        do i1=1,dima
        do i2=1,dimb
        A(lasta+i1,lastb+i2,i3)=B(i1,i2,i3)
        A(lastb+i2,lasta+i1,i3)=B(i1,i2,i3)
        end do
        end do
        end do
c
        return
        end

c
c ----------------------
c
        subroutine expand4_12 (AA,BB,d1,d2,d3)
c
c this routine do :
c
c A(ab,i,j) -> A(a,b,i,j)
c
        implicit none
        integer d1,d2,d3,a,b,i,j,ab
        real*8 AA(1:(d1*(d1+1))/2,d2,d3),BB(1:d1,1:d1,1:d2,1:d3)
c
        ab=0
        do a=1,d1
        do b=1,a
        ab=ab+1
        do i=1,d2
        do j=1,d3
        BB(a,b,i,j)=AA(ab,i,j)
        if (a.ne.b) BB(b,a,j,i)=AA(ab,i,j)
        end do
        end do
        end do
        end do
c
        return
        end
c
c       -----------
c
        subroutine GetRest_t3 (t1,t1_tmp,E2old)
c
c        this file read 1) T1o
c                      2) E1old,E2old,niter
c        from RstFil file
c
        implicit none
#include "cht3_ccsd1.fh"
#include "ccsd_t3compat.fh"
        integer niter
        real*8 E1old,E2old
        real*8 t1(*),t1_tmp(*)
c
c        help variables
        integer length,i
c
*       open (unit=LunAux,File='RstFil',form='unformatted')
        Call MOLCAS_BinaryOpen_Vanilla(LunAux,'RstFil')
        length=nv*no
cmp        write (*,*) 'no, nv, length = ',no,nv,length
        call cht3_rea (LunAux,length,t1)
c
        call transp (t1,t1_tmp,nv,no)
c
        do i=1,length
        t1(i+length)=t1_tmp(i)
        t1(i)=t1_tmp(i)
        end do
c
c
        read (LunAux) E1old,E2old,niter

        if (printkey.gt.1) then
        write (6,'(A,2(f15.12,1x))') 'Results from CCSD : E1, E2 ',
     & E1old,E2old
        end if

        close (LunAux)
c
c
        return
        end
c
c ----------------------
c
        subroutine expand3_23 (A,B,dim1,dim2)
c
c this routine do :
c
c AA(a,bc), b>=c -> BB(a,b,c)
c

        implicit none
        integer a,b,c,bc,dim1,dim2
        real*8 AA(1:dim1,1:(dim2*(dim2+1))/2)
        real*8 BB(1:dim1,1:dim2,1:dim2)
c
        bc=0
        do b=1,dim2
        do c=1,b
        bc=bc+1
        do a=1,dim1
c
        BB(a,b,c)=AA(a,bc)
        BB(a,c,b)=AA(a,bc)
c
        end do
        end do
        end do
c
        return
        end
c
c --------------
c
        subroutine transp (AA,BB, dim1,dim2)
c
c this routine do :
c
c AA(a,b) => BB(b,a)
c
        implicit none
        integer dim1,dim2,i,j
        real*8 AA(1:dim1,1:dim2),BB(1:dim2,1:dim1)
c
        do i=1,dim1
        do j=1,dim2
        BB(j,i)=AA(i,j)
        end do
        end do
c
        return
        end
c
c ----------------------
c
        subroutine gen_vvoo (w,l1,tmp,l2)
c
c this routine do
c
c regenerate (ab,ij) integrals from blocked
c MO cholesky vectors
c
c <vv|oo> = (vo|vo)
c
c --------
c
c       L1(m,I,A')
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        real*8 tmp(*),l1(*),w(*),l2(*)
        integer a,b,dima,dimb,lenght,lasta,lastb
        integer a_tmp,b_tmp
c
        do a=1,NvGrp
c
c1        read tmp(m,I,A')
c
        dima=DimGrpaR(a)
         lenght=nc*no*dima
         call GetX_t3 (tmp,lenght,LunAux,L1Name(a),1,1)
c
c5        map l1 (A',I,m) <- tmp (m,I,A')
c
         call Map3_321_t3 (tmp,l1,nc,no,dima)
c
c  ----- read tmp(m,I,B')
c
        do b=1,a
c
        dimb=DimGrpaR(b)
         lenght=nc*no*dimb
         call GetX_t3 (tmp,lenght,LunAux,L1Name(b),1,1)
c
c4        map l2 (m,B',I) <- tmp (m,I,B')
c
        call Map3_132_t3 (tmp,l2,nc,no,dimb)
c
c        zero tmp
c
        call zeroma (tmp,1,dima*no*dimb*no)
c
c7      mult tmp(A',I,B',J) <- l1 (A',I,m) . l2(m,B',J)
c
        call mc0c1a3b (
     & dima*no,nc,nc,dimb*no,
     & dima*no,dimb*no,
     & dima*no,nc,dimb*no,l1,l2,tmp)
c
          lasta=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.1) then
          do b_tmp=1,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
c
c        write (6,'(A,4(i4,x))') 'BB1 dima, dimb, lasta, lastb ',
c     & dima,dimb,lasta,lastb
        call grow_vvoo(w,tmp,no,nv,dima,dimb,lasta,lastb)
c
        if (a.ne.b) then
c
        call Map4_3412_t3 (tmp,l2,dima,no,dimb,no)
c
c        write (6,'(A,4(i4,x))') 'BB2 dima, dimb, lasta, lastb ',
c     & dimb,dima,lastb,lasta
        call grow_vvoo(w,l2,no,nv,dimb,dima,lastb,lasta)
c
        end if
c
c3        end loop over B'
c
        end do
c
c3        end loop over A'
c
        end do
c
        return
        end
c
c -----------------------
c
        subroutine grow_vvoo(A,B,no,nv,dima,dimb,lasta,lastb)
c
c this routine do :
c
c grow A(1324)/(vvoo) from the blocked cholesky vectors
c B(1234)/(vovo)
c
        implicit none
        integer i1,i2,i3,i4,dima,dimb,no,nv
        integer lasta,lastb
        real*8 A(1:nv,1:nv,1:no,1:no),B(1:dima,1:no,1:dimb,1:no)
c
c!        write (6,'(A,4(i3,x))') 'AA lasta, lastb, dima, dimb ',
c!     & lasta,lastb,dima,dimb
c
        do i4=1,no
        do i3=1,no
        do i1=1,dima
        do i2=1,dimb
        A(lasta+i1,lastb+i2,i3,i4)=B(i1,i3,i2,i4)
        end do
        end do
        end do
        end do
c
        return
        end

c
c ----------------------
c

        subroutine block_interf(
     &   ind1f,ind1l,
     &   ind2f,ind2l,
     &   b1f,b1l,nind_b1f,nind_b1l,
     &   b2f,b2l,nind_b2f,nind_b2l)
c
c this routine do :
c
c Interface between Palo and Juzek virtual orbitals blocking structure
c
c ind1f,ind1l,ind2f,ind2l - first and last absolute of the VO indexes of interess
c
c b1f,b1l   - first and last Palo's block which contain ind1
c b2f,b2l   - first and last Palo's block which contain ind2
c
c nind_b1f  - sum of VOs in blocks 1,2, ..., b1f-1
c nind_b2f  - sum of VOs in blocks 1,2, ..., b2f-1
c
c nind_b1l  - # of VOs in b1f before ind1
c nind_b2l  - # of VOs in b2f before ind2
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "ccsd_t3compat.fh"
c
        integer i,sum
        integer ind1f,ind1l,ind2f,ind2l
        integer b1f,b2f,b1l,b2l
        integer nind_b1f,nind_b2f,nind_b1l,nind_b2l
        logical found1,found2,found3,found4
c
c set b1f, b2f, b1l, b2l
c
          sum=0
          found1=.false.
          found2=.false.
          found3=.false.
          found4=.false.
c
        do i=1,NvGrp
          sum=sum+DimGrpaR(i)
c
          if ((ind1f.le.sum).and.(.not.found1)) then
            b1f=i
            found1=.true.
          end if
c
          if ((ind1l.le.sum).and.(.not.found2)) then
            b1l=i
            found2=.true.
          end if
c
          if ((ind2f.le.sum).and.(.not.found3)) then
            b2f=i
            found3=.true.
          end if
c
          if ((ind2l.le.sum).and.(.not.found4)) then
            b2l=i
            found4=.true.
          end if
c
        end do
c
cmp        write (*,*) 'b1f, b1l, b2f, b2l ',b1f,b1l,b2f,b2l
c
c set nind_b1f, nind_b1l
c
        if (b1f.gt.1) then
            sum=0
          do i=1,b1f-1
            sum=sum+DimGrpaR(i)
          end do
            nind_b1f=sum
        else
            nind_b1f=0
        end if
            nind_b1l=ind1f-nind_b1f-1
c
c set nind_b1f, nind_b1l
c
        if (b2f.gt.1) then
            sum=0
          do i=1,b2f-1
            sum=sum+DimGrpaR(i)
          end do
            nind_b2f=sum
        else
            nind_b2f=0
        end if
            nind_b2l=ind2f-nind_b2f-1
c
        return
        end
c
c       -----------
c
        subroutine gather_t2_blocked(
     & length1,length2,
     & ngaf,ngal,ngbf,ngbl,
     & t2,t2_tmp,tmp,switch)
c
c length1 = length of the 1st VO index (=< vblock)
c length2 = length of the 2nd VO index (=< vblock)
c
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        integer ngaf,ngal,ngbf,ngbl
        integer a,b,dima,dimb
cmp        integer lenght
        integer lenght
        integer lasta,lastb
        integer length1,length2
c
        real*8 t2(1:(length1*length2*no*no))
        real*8 tmp(1:(maxdim*maxdim*no*no))
        real*8 t2_tmp(1:(maxdim*maxdim*no*no))
        integer a_tmp,b_tmp
c
        logical sym
        logical switch
c
        sym=.false.
c
cmp        if (ngaf.eq.ngbf) sym=.true.
        if ((ngaf.eq.ngbf).and.(ngal.eq.ngbl)) sym=.true.
c
cmp        write (6,*)
cmp        write (6,*) '------ DimGrpaR ------'
cmp        write (6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
cmp        write (6,*)
c
        do a=ngaf,ngal
        do b=ngbf,min0(a,ngbl)
c
        dima=DimGrpaR(a)
        dimb=DimGrpaR(b)
c
         if (a.eq.b) then  ! a=b
          lenght=(dima*(dima+1)*no*no)/2
          call GetX_t3 (tmp,lenght,LunAux,T2Name(a,b),1,1)
         else ! a>b
          lenght=dima*dimb*no*no
          call GetX_t3 (tmp,lenght,LunAux,T2Name(a,b),1,1)
         end if
c
          lasta=0
        if (a.gt.ngaf) then
          do a_tmp=ngaf,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.ngbf) then
          do b_tmp=ngbf,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
c
cmp        write (6,'(A,2(i5,2x),2(i3,x))') 'lasta, lastb, a, b = ',
cmp     & lasta,lastb,a,b
c
        if (a.eq.b) then ! expand and map
        call expand4_12 (tmp,t2_tmp,dima,no,no)
        call grow_t2_blocked(t2,t2_tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b,sym,switch)
        else
        call grow_t2_blocked(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b,sym,switch)
        end if
c
        end do
        end do
c
        return
        end
c
c --------------------
c
        subroutine grow_t2_blocked(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,grpa,grpb,sym,switch)
c
c this routine do :
c
c
        implicit none
        integer a,b,dima,dimb,nv,no,i,j
        integer lasta,lastb
        integer grpa,grpb
        integer length1,length2
cmp        real*8 t2(1:nv,1:nv,1:no,1:no)
        real*8 t2(1:length1,1:length2,1:no,1:no)
        real*8 tmp(1:dima,1:dimb,1:no,1:no)
        logical sym
        logical switch
c
cmp        write (6,*) 'grow_t2neq dima , dimb  ',dima,dimb
cmp        write (6,*) 'grow_t2neq lasta, lastb ',lasta,lastb
cmp        write (6,*) 'grow_t2neq no           ',no
c
cmp        if (lasta.eq.lastb) then
c?        if (grpa.eq.grpb) then
c?        do j=1,no
c?        do i=1,no
c?        do a=1,dima
c?        do b=1,a
c?        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
c?        if (a.ne.b) t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
c?        end do
c?        end do
c?        end do
c?        end do
c
c?        else
c
        do j=1,no
        do i=1,no
        do b=1,dimb
        do a=1,dima
          if (.not.switch) then
        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
          else
cmp!        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,j,i)
        t2(lasta+a,lastb+b,i,j)=1.0d0*tmp(a,b,i,j)
          end if
cmpn        t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
c
        if (sym) then
        t2(lastb+b,lasta+a,j,i)=1.0d0*tmp(a,b,i,j)
        end if
c
        end do
        end do
        end do
        end do
c
c?        end if
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nv)
        call Unused_integer(grpa)
        call Unused_integer(grpb)
      end if
        end
c
c ---------------
c
        subroutine gather_t2anti_blocked(
     & length1,length2,
     & ngaf,ngal,ngbf,ngbl,
     & t2,t2_tmp,tmp)
c
c length1 = length of the 1st VO index (nv)
c length2 = length of the 2nd VO index (=< vblock)
c
c This routine generates T2 amplitudes in this form :
c
c T2 = t2(a,b,j<i) - t2(b,a,j<i) finaly : <a,b,(i<j) >  ;  a in nv, b in vblock
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        integer ngaf,ngal,ngbf,ngbl
        integer a,b,dima,dimb
        integer lenght
        integer lasta,lastb
        integer length1,length2
c
        real*8 t2(*),tmp(*),t2_tmp(*)
        integer a_tmp,b_tmp
c
        logical switch
        integer aa,bb
c
cmp        write (6,*)
cmp        write (6,*) '------ DimGrpaR ------'
cmp        write (6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
cmp        write (6,*)
c
        do a=1,NvGrp
        do b=ngbf,ngbl
c
          switch=.false.
        if (a.ge.b) then
          aa=a
          bb=b
        else
          aa=b
          bb=a
          switch=.true.
        end if
c
        dima=DimGrpaR(aa)
        dimb=DimGrpaR(bb)
c
         if (aa.eq.bb) then  ! aa=bb
          lenght=(dima*(dima+1)*no*no)/2
          call GetX_t3 (tmp,lenght,LunAux,T2Name(aa,bb),1,1)
         else ! aa>bb
          lenght=dima*dimb*no*no
          call GetX_t3 (tmp,lenght,LunAux,T2Name(aa,bb),1,1)
         end if
c
          lasta=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.ngbf) then
          do b_tmp=ngbf,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
c
cmp        write (6,'(A,2(i5,2x),L,2(i3,x))') 'lasta, lastb, switch, a, b = ',
cmp     & lasta,lastb,switch,a,b
c
        if (aa.eq.bb) then ! expand tmp
        call expand4_12 (tmp,t2_tmp,dima,no,no)
        call grow_t2anti_blocked1(t2,t2_tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b)
        else
         if (.not.switch) then
        call grow_t2anti_blocked1(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b)
         else
        call grow_t2anti_blocked2(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b)
        end if
        end if
c
        end do
        end do
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(ngaf)
        call Unused_integer(ngal)
      end if
        end
c
c --------------------
c
        subroutine grow_t2anti_blocked1(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,grpa,grpb)
c
c this routine do :
c
c
        implicit none
        integer a,b,dima,dimb,nv,no,i,j,ij
        integer lasta,lastb
        integer grpa,grpb
        integer length1,length2
cmp        real*8 t2(1:nv,1:nv,1:no,1:no)
        real*8 t2(1:length1,1:length2,1:(((no-1)*no)/2))
        real*8 tmp(1:dima,1:dimb,1:no,1:no)
c
cmp        write (6,*) 'lasta+dima, length1 ',lasta+dima, length1
cmp        write (6,*) 'lastb+dimb, length2 ',lastb+dimb, length2
c
cmp        write (6,'(A,2(i10,x),i3)') 'length1, length2, no ',length1, length2, no
cmp        write (6,'(A,4(i3,x))') 'dima, dimb, lasta, lastb',
cmp     & dima,dimb,lasta,lastb
c
        ij=0
        do i=2,no
        do j=1,i-1
        ij=ij+1
        do b=1,dimb
        do a=1,dima
c
        t2(lasta+a,lastb+b,ij)=tmp(a,b,i,j)+(-1.0d0*tmp(a,b,j,i))
c
        end do
        end do
        end do
        end do
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nv)
        call Unused_integer(grpa)
        call Unused_integer(grpb)
      end if
        end
c
c ---------------
c
        subroutine grow_t2anti_blocked2(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,grpa,grpb)
c
c this routine do :
c
c
        implicit none
        integer a,b,dima,dimb,nv,no,i,j,ij
        integer lasta,lastb
        integer grpa,grpb
        integer length1,length2
cmp        real*8 t2(1:nv,1:nv,1:no,1:no)
        real*8 t2(1:length1,1:length2,1:(((no-1)*no)/2))
        real*8 tmp(1:dima,1:dimb,1:no,1:no)
c
cmp        write (6,*) 'lasta+dima, length1 ',lasta+dima, length1
cmp        write (6,*) 'lastb+dimb, length2 ',lastb+dimb, length2
cmp        write (6,*) 'lasta, lastb ',lasta,lastb
cmp        write (6,*) 'dima, dimb ',dima,dimb
c
        ij=0
        do i=2,no
        do j=1,i-1
        ij=ij+1
        do b=1,dima
        do a=1,dimb
c
        t2(lasta+a,lastb+b,ij)=tmp(b,a,j,i)+(-1.0d0*tmp(b,a,i,j))
c
        end do
        end do
        end do
        end do
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nv)
        call Unused_integer(grpa)
        call Unused_integer(grpb)
      end if
        end
c
c ----------------------
c
        subroutine gen_vvoo_blocked (w,l1,tmp,l2,
     & length1,length2,ngaf,ngal,ngbf,ngbl)
c
c this routine do
c
c regenerate (ab,ij) integrals from blocked
c MO cholesky vectors
c
c <vv|oo>
c
c --------
c
c       L1(m,I,A')
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        real*8 tmp(*),l1(*),w(*),l2(*)
        integer a,b,dima,dimb,lenght,lasta,lastb
        integer a_tmp,b_tmp
c
        integer ngaf,ngal,ngbf,ngbl
        integer length1,length2
        logical sym
c
        sym=.false.
c
        if ((ngaf.eq.ngbf).and.(ngal.eq.ngbl)) sym=.true.
c
        do a=ngaf,ngal
c
c1        read tmp(m,I,A')
c
        dima=DimGrpaR(a)
        lenght=nc*no*dima
        call GetX_t3 (tmp,lenght,LunAux,L1Name(a),1,1)
c
c5        map l1 (A',I,m) <- tmp (m,I,A')
c
        call Map3_321_t3 (tmp,l1,nc,no,dima)
c
c  ----- read tmp(m,I,B')
c
        do b=ngbf,min0(a,ngbl)
c
        dimb=DimGrpaR(b)
        lenght=nc*no*dimb
        call GetX_t3 (tmp,lenght,LunAux,L1Name(b),1,1)
c
c4        map l2 (m,B',I) <- tmp (m,I,B')
c
        call Map3_132_t3 (tmp,l2,nc,no,dimb)
c
c        zero tmp
c
        call zeroma (tmp,1,dima*no*dimb*no)
c
c7      mult tmp(A',I,B',J) <- l1 (A',I,m) . l2(m,B',J)
c
        call mc0c1a3b (
     & dima*no,nc,nc,dimb*no,
     & dima*no,dimb*no,
     & dima*no,nc,dimb*no,l1,l2,tmp)
c
          lasta=0
        if (a.gt.ngaf) then
          do a_tmp=ngaf,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.ngbf) then
          do b_tmp=ngbf,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
c
c8        grow w(A',B',I,J) <- tmp(A',I,B',J)
c
cmp        write (6,'(A,6(i4,x))') 'BB1 dima, dimb, lasta, lastb ',
cmp     & dima,dimb,lasta,lastb,a,b
        call grow_vvoo_blocked(w,tmp,no,nv,dima,dimb,lasta,lastb,
     & length1,length2,a,b,sym)
c
c?        if (a.ne.b) then
c?c
c?        call Map4_3412_t3 (tmp,l2,dima,no,dimb,no)
c?c
c?c        write (6,'(A,4(i4,x))') 'BB2 dima, dimb, lasta, lastb ',
c?c     & dimb,dima,lastb,lasta
c?        call grow_vvoo_blocked(w,l2,no,nv,dimb,dima,lastb,lasta,
c?     & length1,length2,a,b)
c?c
c?        end if
c
c3        end loop over B'
c
        end do
c
c3        end loop over A'
c
        end do
c
        return
        end
c
c -----------------------
c
        subroutine grow_vvoo_blocked(AA,BB,no,nv,dima,dimb,lasta,lastb,
     & length1,length2,a,b,sym)
c
c this routine do :
c
c grow A(1324)/(vvoo) from the blocked cholesky vectors
c B(1234)/(vovo)
c
        implicit none
        integer i1,i2,i3,i4,dima,dimb,no,nv
        integer lasta,lastb
        integer length1,length2,a,b
        real*8 AA(1:length1,1:length2,1:no,1:no)
        real*8 BB(1:dima,1:no,1:dimb,1:no)
        logical sym
c
cmp        write (6,'(A,4(i3,x))') 'AA lasta, lastb, dima, dimb ',
cmp     & lasta,lastb,dima,dimb
cmp        write (6,'(A,2(i9,x))') 'chk_a ',lasta+dima,length1
cmp        write (6,'(A,2(i9,x))') 'chk_b ',lastb+dimb,length2
c
        do i4=1,no
        do i3=1,no
        do i1=1,dima
        do i2=1,dimb
        AA(lasta+i1,lastb+i2,i3,i4)=BB(i1,i3,i2,i4)
c
        if (sym) AA(lastb+i2,lasta+i1,i4,i3)=BB(i1,i3,i2,i4)
c
        end do
        end do
        end do
        end do
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nv)
        call Unused_integer(a)
        call Unused_integer(b)
      end if
        end

c
c ----------------------
c
        subroutine gather_t2_fblocked(
     & length1,length2,
     & ngaf,ngal,ngbf,ngbl,
     & t2,t2_tmp,tmp)
c
c length1 = length of the 1st VO index (nv)
c length2 = length of the 2nd VO index (=< vblock)
c
c This routine generates T2 amplitudes in this form :
c
c T2 = t2(a,b,j<i) - t2(b,a,j<i) finaly : <a,b,(i<j) >  ;  a in nv, b in vblock
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
#include "ccsd_t3compat.fh"
c
        integer ngaf,ngal,ngbf,ngbl
        integer a,b,dima,dimb
        integer lenght
        integer lasta,lastb
        integer length1,length2
c
        real*8 t2(*),tmp(*),t2_tmp(*)
        integer a_tmp,b_tmp
c
        logical switch
        integer aa,bb
c
cmp        write (6,*)
cmp        write (6,*) '------ DimGrpaR ------'
cmp        write (6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
cmp        write (6,*)
c
        do a=1,NvGrp
        do b=ngbf,ngbl
c
          switch=.false.
        if (a.ge.b) then
          aa=a
          bb=b
        else
          aa=b
          bb=a
          switch=.true.
        end if
c
        dima=DimGrpaR(aa)
        dimb=DimGrpaR(bb)
c
         if (aa.eq.bb) then  ! aa=bb
          lenght=(dima*(dima+1)*no*no)/2
          call GetX_t3 (tmp,lenght,LunAux,T2Name(aa,bb),1,1)
         else ! aa>bb
          lenght=dima*dimb*no*no
          call GetX_t3 (tmp,lenght,LunAux,T2Name(aa,bb),1,1)
         end if
c
          lasta=0
        if (a.gt.1) then
          do a_tmp=1,a-1
          lasta=lasta+DimGrpaR(a_tmp)
          end do
        end if
c
          lastb=0
        if (b.gt.ngbf) then
          do b_tmp=ngbf,b-1
          lastb=lastb+DimGrpaR(b_tmp)
          end do
        end if
c
cmp        write (6,'(A,2(i5,2x),L,2(i3,x))') 'lasta, lastb, switch, a, b = ',
cmp     & lasta,lastb,switch,a,b
c
        if (aa.eq.bb) then ! expand tmp
        call expand4_12 (tmp,t2_tmp,dima,no,no)
        call grow_t2_fblocked1(t2,t2_tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b)
        else
         if (.not.switch) then
        call grow_t2_fblocked1(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b)
         else
        call grow_t2_fblocked2(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,a,b)
        end if
        end if
c
        end do
        end do
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(ngaf)
        call Unused_integer(ngal)
      end if
        end
c
c --------------------
c
        subroutine grow_t2_fblocked1(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,grpa,grpb)
c
c this routine do :
c
c
        implicit none
        integer a,b,dima,dimb,nv,no,i,j
c       integer ij
        integer lasta,lastb
        integer grpa,grpb
        integer length1,length2
cmp        real*8 t2(1:nv,1:nv,1:no,1:no)
        real*8 t2(1:length1,1:length2,1:no,1:no)
        real*8 tmp(1:dima,1:dimb,1:no,1:no)
c
cmp        write (6,*) 'lasta+dima, length1 ',lasta+dima, length1
cmp        write (6,*) 'lastb+dimb, length2 ',lastb+dimb, length2
c
cmp        write (6,'(A,2(i10,x),i3)') 'length1, length2, no ',length1, length2, no
cmp        write (6,'(A,4(i3,x))') 'dima, dimb, lasta, lastb',
cmp     & dima,dimb,lasta,lastb
c
cmp        ij=0
cmp        do i=2,no
cmp        do j=1,i-1
cmp        ij=ij+1
        do i=1,no
        do j=1,no
        do b=1,dimb
        do a=1,dima
c
cmp        t2(lasta+a,lastb+b,ij)=tmp(a,b,i,j)+(-1.0d0*tmp(a,b,j,i))
        t2(lasta+a,lastb+b,i,j)=tmp(a,b,i,j)
c
        end do
        end do
        end do
        end do
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nv)
        call Unused_integer(grpa)
        call Unused_integer(grpb)
      end if
        end
c
c ---------------
c
c
        subroutine grow_t2_fblocked2(t2,tmp,dima,dimb,nv,no,
     & lasta,lastb,length1,length2,grpa,grpb)
c
c this routine do :
c
c
        implicit none
        integer a,b,dima,dimb,nv,no,i,j
c       integer ij
        integer lasta,lastb
        integer grpa,grpb
        integer length1,length2
cmp     real*8 t2(1:nv,1:nv,1:no,1:no)
        real*8 t2(1:length1,1:length2,1:no,1:no)
        real*8 tmp(1:dima,1:dimb,1:no,1:no)
c
cmp     write (6,*) 'lasta+dima, length1 ',lasta+dima, length1
cmp     write (6,*) 'lastb+dimb, length2 ',lastb+dimb, length2
cmp     write (6,*) 'lasta, lastb ',lasta,lastb
cmp     write (6,*) 'dima, dimb ',dima,dimb
c
cmp        ij=0
cmp        do i=2,no
cmp        do j=1,i-1
cmp        ij=ij+1
        do i=1,no
        do j=1,no
        do b=1,dima
        do a=1,dimb
c
cmp        t2(lasta+a,lastb+b,ij)=tmp(b,a,j,i)+(-1.0d0*tmp(b,a,i,j))
        t2(lasta+a,lastb+b,i,j)=tmp(b,a,j,i)
c
        end do
        end do
        end do
        end do
c
        return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(nv)
        call Unused_integer(grpa)
        call Unused_integer(grpb)
      end if
        end
c
c ---------------
c
        subroutine check_create_klvab_t3_mem (vblock)
c
c this routine finds the upper estimate of the memory
c requirements of the most demanding step in create_klvab_t3
c
        implicit none
c
#include "cht3_ccsd1.fh"
#include "ccsd_t3compat.fh"
        integer vblock,vblock_my
        integer mem,mem_trial,mem_avail
c
c.0 - calculate vblock_my
c
        call my_block (vblock,vblock_my)
c
        if (printkey.ge.10) then
        write (6,*)
        write (6,*) 'check_create_klvab_t3_mem '
        write (6,*)
        write (6,'(A,3(i5,1x))') 'nc,no,nv',nc,no,nv
        write (6,'(A,3(i5,1x))') 'maxdim,vblock,vblock_my',
     & maxdim,vblock,vblock_my
        end if
c
c.1 !create
        mem=vblock*vblock*(no+nv)+
     & nv*((nv*(nv+1))/2)+nv*nv+nc*maxdim+nc*maxdim*maxdim+
     & max(nc*maxdim*maxdim,nc*no*maxdim,maxdim*maxdim*maxdim)
c.2 !klvaa_vvvo
        mem_trial=vblock*vblock*(no+nv)+
     & (nv*(nv*(nv+1))/2)+nv*nv+vblock_my*vblock_my*no*no+
     & 2*maxdim*maxdim*no*no
c
        if (mem_trial.gt.mem) mem=mem_trial
c.3 !create
        mem_trial=vblock*vblock*(no+nv)+
     & (nv*(nv*(nv+1))/2)+nv*nv+vblock_my*vblock_my*no*no+
     & 2*maxdim*maxdim*no*no
c
        if (mem_trial.gt.mem) mem=mem_trial
c.4 !create
        mem_trial=no*no*vblock*(no+nv)+
     & no*nv*(no*(no+1)/2)+vblock*no*no+nc*(no*(no+1)/2)+
     & nc*no*nv+max(nc*((no*(no+1))/2),nc*no*maxdim,nc*no*nv)
c
        if (mem_trial.gt.mem) mem=mem_trial
c.5 !create
        mem_trial=no*no*vblock*(no+nv)+
     & no*nv*(no*(no+1)/2)+vblock*no*no+nv*vblock_my*no*no+
     & 2*maxdim*maxdim*no*no
c
        if (mem_trial.gt.mem) mem=mem_trial
c.6 !klvaa_oovo
        mem_trial=no*no*vblock*(no+nv)+
     & no*nv*(no*(no+1)/2)+vblock*no*no+
     & nv*vblock_my*(((no-1)*no)/2)+2*maxdim*maxdim*no*no
c
        if (mem_trial.gt.mem) mem=mem_trial
c.7 !klvaa_oovo
        mem_trial=(((no-1)*no)/2)*vblock*vblock+
     & vblock_my*vblock_my*no*no+nc*no*maxdim+
     & 2*max(nc*no*maxdim,maxdim*maxdim*no*no)
c
        if (mem_trial.gt.mem) mem=mem_trial
c.8 !klvaa_oovo
        mem_trial=no*no*vblock*vblock+
     & vblock_my*vblock_my*no*no+nc*no*maxdim+
     & 2*max(nc*no*maxdim,maxdim*maxdim*no*no)
c
        if (mem_trial.gt.mem) mem=mem_trial
c
           if (printkey.ge.10) then
              write (6,*)
              write (6,'(A,f10.1,A,f7.1,A,f3.1,A)')
     &       'Memory required for the reorg. step = ',
     &        (8.0d0*mem)/(1024),' kb ',
     &        (8.0d0*mem)/(1024*1024),' Mb ',
     &        (8.0d0*mem)/(1024*1024*1024),' Gb '
           end if
c
c - calculate available free memory
c
        Call GetMem('(T)','Max','Real',mem_avail,mem_avail)
c
           if (printkey.ge.10) then
              write (6,'(A,f10.1,A,f7.1,A,f3.1,A)')
     &       'Available memory                    = ',
     &        (8.0d0*mem_avail)/(1024),' kb ',
     &        (8.0d0*mem_avail)/(1024*1024),' Mb ',
     &        (8.0d0*mem_avail)/(1024*1024*1024),' Gb '
              write (6,*)
           end if
c
c - check, if mem fits
c
        if (mem_avail.lt.mem) then
          write (6,*) 'Not enough memory for the transformation step '
          call Abend()
        end if
c
        return
        end
c
c ----------------
c
        subroutine my_block (vblock,vblock_my)
c
c this subroutine calculates maximum overlap between juzek's
c vblock segmentation and palo's dimgrp
c
        implicit none
c
#include "cht3_ccsd1.fh"
#include "ccsd_t3compat.fh"
#include "cht3_reord.fh"
c
        integer i,j,i_tmp,i_f,i_l,poss
        integer vblock,vblock_my,vblock_my_tmp
        logical found
c
        vblock_my=0
        i_l=0
        i_f=0
c
        do i=1,nv,vblock
c
c - find initial possition of the i-th juzek's block
c
          poss=0
          found=.false.
        do j=1,NvGrp
          poss=poss+DimGrpaR(j)
          if ((i.le.poss).and.(.not.found)) then
            i_f=j
            found=.true.
cmp        write (6,'(A,3(i5,x))') 'i,i_f,poss     = ',
cmp     & i,i_f,poss
          end if
        end do
c
        if ((i+vblock-1).le.nv) then
          i_tmp=i+vblock-1
        else
          i_tmp=nv
        end if
cmp        write (6,'(A,2(i5,x))') 'i,i_tmp        = ',
cmp     & i,i_tmp
c
c - find terminal possition of the i-th juzek's block
c
          poss=0
          found=.false.
        do j=1,NvGrp
          poss=poss+DimGrpaR(j)
          if ((i_tmp.le.poss).and.(.not.found)) then
            i_l=j
            found=.true.
cmp        write (6,'(A,3(i5,x))') 'i_tmp,i_l,poss = ',
cmp     & i_tmp,i_l,poss
          end if
        end do
c
        vblock_my_tmp=0
        do j=i_f,i_l
        vblock_my_tmp=vblock_my_tmp+DimGrpaR(j)
        end do
c
        if (vblock_my_tmp.gt.vblock_my)
     & vblock_my=vblock_my_tmp
c
cmp        write (6,'(A,2(i5,x))') 'vblock_my_tmp, vblock_my',
cmp     & vblock_my_tmp, vblock_my
cmp        write (6,*)
        end do
c
        return
        end
c
c ----------------
c
        subroutine check_loops(nv,vblock,nla,nlb)
c
        integer nv,vblock,nla,nlb
        integer nuga,nga,ngb,ngc

        nuga=nv/vblock
cmp! pridavok
        if((nuga*vblock).lt.nv)nuga=nuga+1
c
        nla=0
        do nga=1,nuga
        do ngb=1,nga
        do ngc=1,ngb
        nla=nla+1
        end do
        end do
        end do
c
        nlb=0
        do nga=1,nuga
        do ngb=1,nga
        do ngc=1,nuga
        nlb=nlb+1
        end do
        end do
        end do
c
        return
        end
