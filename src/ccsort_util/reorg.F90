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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               1995,1996, Pavel Neogrady                              *
!***********************************************************************
       SUBROUTINE REORG(run_triples,IRETURN)
       IMPLICIT REAL*8 (A-H,O-Z)
       Logical run_triples
#include "reorg.fh"
#include "stdalloc.fh"
       real*8, allocatable :: FIRAS(:),FI(:)
       fullprint=0
       If (iPrintLevel(-1).LE.0) fullprint=-1
       call mma_Allocate(FIRAS,mbas*mbas,Label='FIRAS')
       call mma_Allocate(FI,mbas*mbas,Label='FI')
       call REORG_(FIRAS,FI,run_triples,IRETURN)
       call mma_Deallocate(FIRAS)
       call mma_Deallocate(FI)
       return
       end
       SUBROUTINE REORG_(FIRAS,FI,run_triples,IRETURN)

!----------------------------------------------------------------------*
!     1994  PER-AAKE MALMQUIST                                         *
!     DEPARTMENT OF THEORETICAL CHEMISTRY                              *
!     UNIVERSITY OF LUND, SWEDEN                                       *
!                                                                      *
! modified by P.N. (Dec. 1995)                                         *
! all allocation memory routines removed by P.N. (8.03.1996)           *
!----------------------------------------------------------------------*

!
! FILES USED:
!     TRAINT    2 electron MO INTEGRALS
!     JOBIPH    THE JOB-INTERFACE FILE AS PRODUCED BY THE RASSCF
!               PROGRAM
! INPUT
!     AT PRESENT:  THE INPUT FILE IS SEARCHED FOR THE STRING
!     '&REORG '. Input is only TITLE
!
! LIMITATIONS
!     Like in CASPT2
!
!***********************************************************************
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
!       real*8 Ene(25,50)
!
       INTEGER NOIPSB(106)
       Logical run_triples,run_sort
!
!
! PRINT THE PROGRAM HEADER
       call ccsort_helloPN

!*    READ AND ECHO INPUT DATA, READ JOBIPH, PRINT INPUT DATA,
       CALL RDINPPN(run_triples,run_sort)
       If (fullprint.GE.0) CALL PRINPPN
       CALL CHKINP_CCSORT
!
       if (run_sort) then
!
!*    read FI from JOBIPH
       ntot3=0
       ntot2=0
       do i=1,nsym
       ntot3=ntot3+(norb(i)*(norb(i)+1))/2
       ntot2=ntot2+norb(i)
       end do
!
!*    pick the total energy from the JOBIPH file
!
       iad15=iadr15(6)
       lad15=mxroot*mxiter
       Call dDaFile(JOBIPH,2,Ene,lad15,iad15)
       EScf=0.0d0
       i=1
!
!*    take the last non-zero energy stored
!
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
!
!*    get fi from previous RASSCF
!
       iad15=iadr15(10)
       call ddafile(JOBIPH,2,firas(1),ntot3,iad15)
!
!*    get eps from previous RASSCF
!
       iad15=iadr15(11)
       call ddafile(JOBIPH,2,epsras(1),ntot2,iad15)
!
!*    reduce fi,eps and update n's
       call mod1 (nsym,nfro,nish,nash,nssh,ndel,norb,nfror,ndelr,       &
     &            firas,fi,epsras,eps)
!
!*    def diagonal Fok for closed shell
!
       if (clopkey.eq.2) then
         call mod2 (nsym,nish,nash,nssh,norb,fi,eps)
       end if
!
!*    define noa,nob,nva,nvb
!
       do i=1,nsym
         noa(i)=nish(i)+nash(i)
         nob(i)=nish(i)
         nva(i)=nssh(i)
         nvb(i)=nssh(i)+nash(i)
       end do
!
       if (nsym.lt.8) then
         do i=1+nsym,8
           noa(i)=0
           nob(i)=0
           nva(i)=0
           nvb(i)=0
         end do
       end if
!
       if (fullprint.gt.1) then
         write(6,*)
         write(6,'(6X,A)') 'Diagonal Fock matrix elements and '//       &
     &                      'orbital energies:'
         write(6,'(6X,A)') '----------------------------------'//       &
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
!
!*    prepair adress (stupid)
!
!FUE
!FUE   The unit number of the transformed two electron integrals
!FUE   must be 40, 50, 60, 70, 80 or 90. Any other number will
!FUE   not be compatible with the I/O driver in MOLCAS.
!FUE
       LUINTM=40
!JR       call DANAME (LUINTM,'TRAINT')
       call DANAME_MF (LUINTM,'TRAINT')
       call mkadress (NOIPSB)
!
!*    open TRAINT and call action
!
       Call GetMem('FOKA','ALLO','REAL',ipFOKA,(mbas**2+mbas)/2)
       Call GetMem('FOKB','ALLO','REAL',ipFOKB,(mbas**2+mbas)/2)

       call action (Work(ipFOKA),Work(ipFOKB),fi,eps)
       Call GetMem('FOKA','FREE','REAL',ipFOKA,(mbas**2+mbas)/2)
       Call GetMem('FOKB','FREE','REAL',ipFOKB,(mbas**2+mbas)/2)
!
!      close files
!
       call daclos(luintm)
       call daclos(jobiph)
!
       else
!      case, when SORT was skipped
       write (6,*) ' SORT part was skipped'
       write (6,*) ' Input parameters are from last actual run of SORT'
       end if
!
       ireturn=0
       return
       END
