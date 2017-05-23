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
* Copyright (C) 1993, Markus P. Fuelscher                              *
*               1993, Per Ake Malmqvist                                *
*               Pavel Neogrady                                         *
************************************************************************
       Subroutine RdInpPN(run_triples,run_sort)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     - read input                                                     *
*     - set defaults                                                   *
*                                                                      *
*     calling parameters: none                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and P.-AA. Malmqvist                              *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*     reduced by P.N.                                                  *
*                                                                      *
************************************************************************
       Implicit Real*8 (A-H,O-Z)


#include "SysDef.fh"
#include "ccsort.fh"
#include "reorg.fh"
#include "motra.fh"
*
       Real*8 Weights(mxRoot)
       Parameter ( nCmd =20 )
       Character*4 Command,Cmd(nCmd)
       Character*72  Line,Blank
       Data Cmd /'TITL','END ','CCSD','CCT ','CLOS',
     &           'OPEN','FROZ','DELE','PRIN','NOOP','IOKE','ZROF',
     &           'DENO','SHIF','ACCU','ADAP','EXTR','TRIP','NOSO',
     &           'ITER'/
       Logical run_triples,run_sort

*
*---  Initialize -------------------------------------------------------*
       Do i=1,72
       Blank(i:i)=' '
       End Do
       LROOT=0
       MAXIT=0
       CONV=1.0D-06
       THRSHN=1.0D-10
       THRSHS=1.0D-08
       THRSHF=0.05D0
       ORBIN='DEFAULT '
       THRENE=1.50D0
       ORBIT='DEFAULT '
       THROCC=0.0D0
       FOCKTYPE='STANDARD'
       HZERO='STANDARD'
       METHOD='CONJ'
       IFJAC=0
       RFpert=.false.
       CITHR=0.05D00
       NTIT=0
c
       lunsta=21
       luna1=22
       luna2=23
       luna3=24
       luna4=25
       lunab=50
       lunt3=26
       lunda1=9
       lunda2=10
c
*---  Open JOBIPH and LUONEM files ------------------------------------*
*
*.    Job interface
       JOBIPH=15
*.    Job interface
       CALL DANAME(JOBIPH,'JOBIPH')
*
*---  Read input from JOBIPH file -------------------------------------*
       IAD15=0
       CALL iDAFILE(JOBIPH,2,IADR15,15,IAD15)
cDIVNUO
       IAD15=IADR15(1)
       CALL WR_RASSCF_Info(JOBIPH,2,iAd15,
     &                     NACTEL,ISPIN,NSYM,LSYM,
     &                     NFRO,NISH,NASH,NDEL,NBAS,8,
     &                     NAME,4*2*MXORB,
     &                     NCONF,HEADER,2*72,
     &                     TITLE,4*18*MXTIT,POTNUC,
     &                     LROOTS,NROOTS,IROOT,MXROOT,NRAS1,
     &                     NRAS2,NRAS3,NHOLE1,NELE3,IPT2,Weights)
       SPIN=DBLE(ISPIN-1)/2
c
c     define defaults for REORG
c
       ntAsh = 0
       do isym = 1,nSym
         ntAsh = ntAsh+nAsh(iSym)
       end do
       cckey=1
       t3key=1
       clopkey=1
       if ( ntAsh.eq.0 ) clopkey=2
       do nhelp=1,nsym
         ndelr(nhelp)=NDEL(nhelp)
         nfror(nhelp)=NFRO(nhelp)
       end do
       maxsize=0
CGG       fullprint=0
       noop=0
       iokey=1
       zrkey=1
       run_triples=.true.
       run_sort=.true.
*
*---  Read input from TRAONE file -------------------------------------*
      Call RdTraOne
      Do iSym=1,nSym
        nFror(iSym) = nFroX(iSym)
        nDelr(iSym) = nDelX(iSym)
      End Do
*
*---  Read input from LuSpool -----------------------------------------*
      LuSpool = 17
      Call SpoolInp(LuSpool)
      Rewind(LuSpool)
       Command='&REO'
       Call RdNlst(LuSpool,'CCSDT')
 10    Read(LuSpool,'(A)',End=9910) Line
       If ( Line(1:1).eq.'*' ) Goto 10
       If ( Line.eq.Blank ) Goto 10
       Command=Line(1:4)
       Call UpCase(Command)
       jCmd=0
       Do iCmd=1,nCmd
       If ( Command.eq.Cmd(iCmd) ) jCmd=iCmd
       End Do
       If ( jCmd.eq.0 ) Goto 9930
 20     Goto (100,2100,200,300,400,500,600,700,900,1000,1100,1200,
     & 1300,1400,1500,1600,1700,1800,1900,2000) jCmd
*
*---  process TITLE    command ----------------------------------------*
 100   continue
       Read(LuSpool,'(A)',End=9910) Line
       Command=Line(1:4)
       Call UpCase(Command)
       If ( Command(1:1).eq.'*' ) Goto 100
       jCmd=0
       Do iCmd=1,nCmd
       If ( Command.eq.Cmd(iCmd) ) jCmd=iCmd
       End Do
       If ( jCmd.ne.0 ) Goto 20
       if (nTit.ge.mxTit) Goto 100
       nTit=nTit+1
       If ( nTit.le.mxTit ) Read (Line,'(18A4)') (Title(nTit,i),i=1,18)
       Goto 100
c---  CCSD command ------------
 200   continue
       cckey=1
       t3key=0
       run_triples=.false.
       goto 100
c---  CCT  command ------------
 300   continue
       cckey=1
       t3key=1
       run_triples=.true.
       goto 100
c---  CLOSed  command ------------
 400   continue
       clopkey=2
       goto 100
c---  OPEN  command ------------
 500   continue
       clopkey=1
       goto 100
c---  FROZen  command ------------
 600   continue
       read (LuSpool,*) (nfror(nhelp),nhelp=1,nsym)
       goto 100
c---  DELEte  command ------------
 700   continue
       read (LuSpool,*) (ndelr(nhelp),nhelp=1,nsym)
       goto 100
*
c---  PRINt  command ------------
 900   continue
       read (LuSpool,*) fullprint
       goto 100
*
c---  NOOPeration  command ------------
 1000  continue
       noop=1
       goto 100
*
c---  IOKEy  command ------------
 1100  continue
       read (LuSpool,*) iokey
       if ((iokey.lt.1).or.(iokey.gt.2)) then
       iokey=2
       end if
       goto 100
*
c---  ZROFf  command ------------
 1200  continue
       zrkey=0
       goto 100

c---  DENO command ---------
1300  continue
      goto 100
*
c---  SHIF command ---------
1400  continue
      goto 100
c---  ACCU command ---------
1500  continue
      goto 100
c---  ADAP command ---------
1600  continue
      goto 100
c---  EXTR command ---------
1700  continue
      goto 100
c---  TRIP command ---------
1800  continue
      goto 100
c---  NOSOrt  command ------------
1900   continue
       run_sort=.false.
       goto 100
c---  ITER command ---------
2000  continue
      goto 100
*
*
*---  The end of the input section, complete input processing ---------*
 2100  continue
       NFROT=0
       NISHT=0
       NASHT=0
       NRAS1T=0
       NRAS2T=0
       NRAS3T=0
       NOSHT=0
       NSSHT=0
       NDELT=0
       NORBT=0
       NBAST=0
       NBAS2=0
       NORB1=0
       DO ISYM=1,NSYM
         NIES(ISYM)=NISHT
         NAES(ISYM)=NASHT
         NSES(ISYM)=NSSHT
         NOSH(ISYM)=NISH(ISYM)+NASH(ISYM)
         NSSH(ISYM)=NBAS(ISYM)-NFRO(ISYM)-NOSH(ISYM)-NDEL(ISYM)
         NORB(ISYM)=NOSH(ISYM)+NSSH(ISYM)
         NORBT=NORBT+NORB(ISYM)
         NBAS2=NBAS2+NBAS(ISYM)**2
         NORB1=NORB1+(NORB(ISYM)**2+NORB(ISYM))/2
         NFROT=NFROT+NFRO(ISYM)
         NISHT=NISHT+NISH(ISYM)
         NASHT=NASHT+NASH(ISYM)
         NOSHT=NOSHT+NOSH(ISYM)
         NRAS1T=NRAS1T+NRAS1(ISYM)
         NRAS2T=NRAS2T+NRAS2(ISYM)
         NRAS3T=NRAS3T+NRAS3(ISYM)
         NSSHT=NSSHT+NSSH(ISYM)
         NDELT=NDELT+NDEL(ISYM)
         NBAST=NBAST+NBAS(ISYM)
       END DO
*
*---  Identify the wave function type ---------------------------------*
       ISCF=0
       IF(NASHT.EQ.0) ISCF=1
       IF(NACTEL.EQ.2*NASHT) ISCF=1
       if ( iSpin.gt.1 )  then
         IF((ISPIN.EQ.NACTEL+1).AND.(NACTEL.EQ.NASHT)) ISCF=2
       end if
c
c     test agreement between REORG input and JOBIPH
c       if (clopkey.ne.(3-ISCF)) then
c         write(6,*) ' Diference in closed/open specification'
c         write(6,*) ' Plaese, correct REORG input file'
c         write(6,*)   clopkey,ISCF
c         Call Quit(16)
c       end if
*
c      Should not be necessary anymore, let's just hardwire it in.
c      This will save a keyword
       clopkey = 3-ISCF
*---  Identify the reference function ---------------------------------*
       IF(ISCF.GT.0) THEN
       LROOT=1
       NROOTS=1
       IROOT(1)=1
       END IF
       IF(LROOT.EQ.0 .AND. NROOTS.EQ.1) LROOT=IROOT(1)
*
*---  Create the symmetry multiplication table ------------------------*
       MUL(1,1)=1
       M=1
       DO N=1,3
         DO I=1,M
           DO J=1,M
             MUL(I+M,J)=M+MUL(I,J)
             MUL(I,J+M)=MUL(I+M,J)
             MUL(I+M,J+M)=MUL(I,J)
           END DO
         END DO
         M=2*M
       END DO
       Call Close_LuSpool(LuSpool)

*
*
*---  Exit ------------------------------------------------------------*
       Return
*
*---  Error exits -----------------------------------------------------*
 9910  Write(6,*)
       Write(6,*) ' *** input error ***'
       Write(6,*) ' hitting end of file mark'
       Write(6,*)
       Call Quit_OnUserError()
 9930  Write(6,*)
       Write(6,*) ' *** input error ***'
       Write(6,*) ' unknown input'
       Write(6,*) ' line: ',Line
       Write(6,*)
       Call Quit_OnUserError()
       End
