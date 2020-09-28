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
* Copyright (C) 1994, Per Ake Malmqvist                                *
*               1995,1996, Pavel Neogrady                              *
************************************************************************
       SUBROUTINE REORG(run_triples,IRETURN)
       IMPLICIT REAL*8 (A-H,O-Z)
       Logical run_triples
#include "reorg.fh"
#include "WrkSpc.fh"
       fullprint=0
       If (iPrintLevel(-1).LE.0) fullprint=-1
       Call GetMem('FIRAS','ALLO','REAL',ipFIRAS,mbas*mbas)
       Call GetMem('FI',   'ALLO','REAL',ipFI,mbas*mbas)
       call REORG_(Work(ipFIRAS),Work(ipFI),run_triples,IRETURN)
       Call GetMem('FIRAS','FREE','REAL',ipFIRAS,mbas*mbas)
       Call GetMem('FI',   'FREE','REAL',ipFI,mbas*mbas)
       return
       end
       SUBROUTINE REORG_(FIRAS,FI,run_triples,IRETURN)

*----------------------------------------------------------------------*
*     1994  PER-AAKE MALMQUIST                                         *
*     DEPARTMENT OF THEORETICAL CHEMISTRY                              *
*     UNIVERSITY OF LUND, SWEDEN                                       *
*                                                                      *
* modified by P.N. (Dec. 1995)                                         *
* all allocation memory routines removed by P.N. (8.03.1996)           *
*----------------------------------------------------------------------*

C
C FILES USED:
C     TRAINT    2 electron MO INTEGRALS
C     JOBIPH    THE JOB-INTERFACE FILE AS PRODUCED BY THE RASSCF
C               PROGRAM
C INPUT
C     AT PRESENT:  THE INPUT FILE IS SEARCHED FOR THE STRING
C     '&REORG '. Input is only TITLE
C
C LIMITATIONS
c     Like in CASPT2
C
************************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"
#include "ccsort.fh"
#include "reorg.fh"
#include "intgrl.fh"
#include "WrkSpc.fh"

       real*8 FIRAS(1:mbas*mbas)
       real*8 EPSRAS(1:mbas)
       real*8 FI(1:mbas*mbas)
       real*8 EPS(1:mbas)
       real*8 ene(mxRoot,mxIter)
c       real*8 Ene(25,50)
c
       INTEGER NOIPSB(106)
       Logical run_triples,run_sort
c
c
C PRINT THE PROGRAM HEADER
       call ccsort_helloPN

c*    READ AND ECHO INPUT DATA, READ JOBIPH, PRINT INPUT DATA,
       CALL RDINPPN(run_triples,run_sort)
       If (fullprint.GE.0) CALL PRINPPN
       CALL CHKINP_CCSORT
c
       if (run_sort) then
c
c*    read FI from JOBIPH
       ntot3=0
       ntot2=0
       do i=1,nsym
       ntot3=ntot3+(norb(i)*(norb(i)+1))/2
       ntot2=ntot2+norb(i)
       end do
c
c*    pick the total energy from the JOBIPH file
c
       iad15=iadr15(6)
       lad15=mxroot*mxiter
       Call dDaFile(JOBIPH,2,Ene,lad15,iad15)
       EScf=0.0d0
       i=1
c
c*    take the last non-zero energy stored
c
       Do While ((Ene(LROOT,i).ne.0.0D0) .and. (i.le.mxIter))
         Escf = Ene(LROOT,i)
         i=i+1
       End Do
       If (fullprint.GE.0) then
          write(6,*)
          write(6,'(6X,A,F16.8)') 'SCF energy:',Escf
          write(6,'(6X,A)')       '-----------'
          write(6,*)
       EndIf
c
c*    get fi from previous RASSCF
c
       iad15=iadr15(10)
       call ddafile(JOBIPH,2,firas(1),ntot3,iad15)
c
c*    get eps from previous RASSCF
c
       iad15=iadr15(11)
       call ddafile(JOBIPH,2,epsras(1),ntot2,iad15)
c
c*    reduce fi,eps and update n's
       call mod1 (nsym,nfro,nish,nash,nssh,ndel,norb,nfror,ndelr,
     &            firas,fi,epsras,eps)
c
c*    def diagonal Fok for closed shell
c
       if (clopkey.eq.2) then
         call mod2 (nsym,nish,nash,nssh,norb,fi,eps)
       end if
c
c*    define noa,nob,nva,nvb
c
       do i=1,nsym
         noa(i)=nish(i)+nash(i)
         nob(i)=nish(i)
         nva(i)=nssh(i)
         nvb(i)=nssh(i)+nash(i)
       end do
c
       if (nsym.lt.8) then
         do i=1+nsym,8
           noa(i)=0
           nob(i)=0
           nva(i)=0
           nvb(i)=0
         end do
       end if
c
       if (fullprint.gt.1) then
         write(6,*)
         write(6,'(6X,A)') 'Diagonal Fock matrix elements and '//
     &                      'orbital energies:'
         write(6,'(6X,A)') '----------------------------------'//
     &                      '-----------------'
         write(6,*)
         write(6,'(6X,A)') '----------------------------------------'
         write(6,'(6X,A)') '   i      F(i,i)           eps(i)       '
         write(6,'(6X,A)') '----------------------------------------'
         ij=0
         do i=1,norb(1)
           do j=1,i
             ij=ij+1
             if (i.eq.j) then
               write(6,'(6X,I4,2F18.10)') i,fi(ij),eps(i)
             end if
           end do
         end do
         write(6,'(6X,A)') '----------------------------------------'
         write(6,*)
       end if
c
c*    prepair adress (stupid)
c
CFUE
CFUE   The unit number of the transformed two electron integrals
CFUE   must be 40, 50, 60, 70, 80 or 90. Any other number will
CFUE   not be compatible with the I/O driver in MOLCAS.
CFUE
       LUINTM=40
cJR       call DANAME (LUINTM,'TRAINT')
       call DANAME_MF (LUINTM,'TRAINT')
       call mkadress (NOIPSB)
c
c*    open TRAINT and call action
c
       Call GetMem('FOKA','ALLO','REAL',ipFOKA,(mbas**2+mbas)/2)
       Call GetMem('FOKB','ALLO','REAL',ipFOKB,(mbas**2+mbas)/2)

       call action (Work(ipFOKA),Work(ipFOKB),fi,eps)
       Call GetMem('FOKA','FREE','REAL',ipFOKA,(mbas**2+mbas)/2)
       Call GetMem('FOKB','FREE','REAL',ipFOKB,(mbas**2+mbas)/2)
c
c      close files
c
       call daclos(luintm)
       call daclos(jobiph)
c
       else
c      case, when SORT was skipped
       write (6,*) ' SORT part was skipped'
       write (6,*) ' Input parameters are from last actual run of SORT'
       end if
c
       ireturn=0
       return
       END
c
c     ---------------------------
c
       subroutine mod1 (nsym,nfro,nish,nash,nssh,ndel,norb,nfror,ndelr,
     & firas,fi,epsras,eps)
c
c     this routine do:
c     1) reduce firas, epsras if nfror>nfro, ndelr>ndel
c     2) redefine nfro,nish,nash,nssh,ndel,norb to proper ones
c
       integer nsym
       integer nfro(1:8)
       integer ndel(1:8)
       integer nish(1:8)
       integer nash(1:8)
       integer nssh(1:8)
       integer norb(1:8)
       integer nfror(1:8)
       integer ndelr(1:8)
       real*8 fi(*)
       real*8 firas(*)
       real*8 eps(*)
       real*8 epsras(*)
c
c     help variables
c
       integer p,q,pqras,pqnew,pras,pnew,isym
       integer ndf,ndd,nup,nlow

c
c1    reduce fi
c
       pqras=0
       pqnew=0
       do 100 isym=1,nsym
c
       ndf=nfror(isym)-nfro(isym)
       ndd=ndelr(isym)-ndel(isym)
       nlow=ndf+1
       nup=norb(isym)-ndd
c
       do 50 p=1,norb(isym)
       do 51 q=1,p
       pqras=pqras+1
c
       if ((p.ge.nlow).and.(p.le.nup)) then
       if ((q.ge.nlow).and.(q.le.nup)) then
       pqnew=pqnew+1
       fi(pqnew)=firas(pqras)
       end if
       end if
c
 51     continue
 50     continue
c
 100    continue
c
c2    reduce eps
c
       pras=0
       pnew=0
       do 200 isym=1,nsym
c
       ndf=nfror(isym)-nfro(isym)
       ndd=ndelr(isym)-ndel(isym)
       nlow=ndf+1
       nup=norb(isym)-ndd
c
       do 150 p=1,norb(isym)
       pras=pras+1
c
       if ((p.ge.nlow).and.(p.le.nup)) then
       pnew=pnew+1
       eps(pnew)=epsras(pras)
       end if
c
 150    continue
c
 200    continue
c
c3    define new nfro,nish,nash,nssh,ndel,norb
c
       do 300 isym=1,nsym
       nash(isym)=nash(isym)
       nish(isym)=nish(isym)-nfror(isym)+nfro(isym)
       nssh(isym)=nssh(isym)-ndelr(isym)+ndel(isym)
       norb(isym)=norb(isym)-nfror(isym)+nfro(isym)-ndelr(isym)
     &           +ndel(isym)
       nfro(isym)=nfror(isym)
 300    continue
c
c
       return
       end
c
c     ---------------------------
c
       subroutine mod2 (nsym,nish,nash,nssh,norb,fi,eps)
c
c     this routine define fi(p,q) = delta(p,q).eps(p)
c     and redefine nish=nish+nash, nash=0
c
c     this is suitable for closed shell case
c
       integer nsym
       integer nish(1:8)
       integer nash(1:8)
       integer nssh(1:8)
       integer norb(1:8)
       real*8 fi(*)
       real*8 eps(*)
c
c     help variables
c
       integer p,q,isym,pq,padd
c
c1    redefine foki
c
       pq=0
       padd=0
       do 200 isym=1,nsym
c
       do 100 p=1,norb(isym)
       do 101 q=1,p
       pq=pq+1
       if (p.eq.q) then
       fi(pq)=eps(padd+p)
       else
       fi(pq)=0.0d0
       end if
 101    continue
 100    continue
c
       padd=padd+norb(isym)
 200    continue
c
c2    redefine n's
c
       do 300 isym=1,nsym
       nish(isym)=nish(isym)+nash(isym)
       nash(isym)=0
 300    continue
c
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(nssh)
       end
