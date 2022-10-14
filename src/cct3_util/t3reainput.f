!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine t3reainput
!
!     this routine do:
!     1) read INPDAT file, produced by REORG with mul,nsym,noa,nob,nva,nvb,norb,eps
!     2) read input file for NIT3 to read (parameters transported through cmm common)

!   ####################
!   Due to the merging of CC input files to one, to avoid conflicts
!   Denominators in CCT3 has become T3Denominators
!   and Shift has become T3Shift
!                           (JR) Lund 2003
!     ! title   - jobtitle
!     1-ntit rows with 72 characters
!     no default
!     ! ntit    - number of rowws in jobtitle
!     ntit is limited to 10 CGG From now 1
!     default (1)
!     ! typt3   - type of T3 cpntribution
!     0 - CCSD
!     1 - CCSD+T(CCSD)
!     2 - CCSD(T) Ragh
!     3 - CCSD(T) Bart
!     default (3)
!     ! typden  - type of denominator (division of fok)
!     0 - diagonal elements
!     1 - average of faa and fbb
!     2 - orbital energies
!     default (0)
!     ! keysa   - Spin adaptation key
!     0 - no adaptation
!     1 - T2 DDVV adaptation
!     2 - T2 DDVV + T1 DV adaptation
!     3 - full T1 and T2 adaptation (only for doublets)
!     4 - full T2 adaptation without SDVS (only for doublets)
!     (default=0)
!     ! filerst - name for CCSD results containing file
!     (default=RSTART)
!     ! mchntyp - type of machine in matrix multiplication
!     1 - C=A*B is faster or comparable with C=AT*B
!     2 - C=AT*B is faster
!     (default=1)
!     ! slim    - limitation for usieng C=AT*B
!     no default (suitable=2.0d0)
!     ! shifhto - shift for occupied
!     (default=0.0)
!     ! shifhtv - shift for virtuals
!     (default=0.0)
!     ! maxspace - maximal allowed work space
!     (default=0 - unlimited)
!     ! fullprint - level of printing control key
!     (default=0)
!     ! noop - No Operation key
!     (default=no)
!     & iokey - I/O control key
!       1 - Fortran I/O system
!       2 - MOLCAS DA IO system
!     (default=2)
!     & mhkey - Matrix handling control key
!       1 - ESSL routines
!       2 - Fortran I/O system
!     (default=1)
!     .....   - can be added
!
!     3) initialize nshf
!
!
#include "t31.fh"
!
!     help variables
!
       character*80 LINE
       integer nhelp
!
!1    read INPDAT
!
!       open (unit=1,file='INPDAT',form='unformatted')
       call molcas_binaryopen_vanilla(1,'INPDAT')
       read (1) nactel,ispin,nsym,lsym,mmul,noa,nob,nva,nvb,norb,eps
       close (1)
!
!2    def dimm
!
       do 10 nhelp=1,nsym
!
       dimm(1,nhelp)=noa(nhelp)
       dimm(2,nhelp)=nob(nhelp)
       dimm(3,nhelp)=nva(nhelp)
       dimm(4,nhelp)=nvb(nhelp)
       dimm(5,nhelp)=norb(nhelp)
!
 10     continue
!
!3    define nshf
!
       do 20 nhelp=1,maxorb
       nshf(nhelp)=(nhelp-1)*(nhelp-2)/2
 20     continue
!
!
!4    define defaults
!
       typt3=3
       ntit=1
       typden=0
       keysa=0
       filerst='RSTART'
       mchntyp=1
       slim=1.0d0
       shifto=0.0d0
       shiftv=0.0d0
       maxspace=0
!GG       fullprint=0
       noop=0
       iokey=1
       mhkey=1
       ijsegkey=0
         symimin=1
         symjmin=1
         symimax=nsym
         symjmax=nsym
         imin=0
         jmin=0
         imax=0
         jmax=0
!
!5    read input file
!
      LuSpool = 17
      Call SpoolInp(LuSpool)
      Rewind(LuSpool)
      TITLE=' '
 5      Read(LuSpool,'(A80)') LINE
       CALL UPCASE(LINE)
       IF( INDEX(LINE,'&CCSDT') .EQ. 0 ) GOTO 5
       NTIT=1
 6      Read(LuSpool,'(A80)') LINE
       IF(LINE(1:1).EQ.'*') GOTO 6
       CALL UPCASE(LINE)
!
       IF (LINE(1:4).EQ.'TITL') THEN
       Read(LuSpool,'(A72)') TITLE
       ELSE IF (LINE(1:4).EQ.'TRIP') THEN
       Read(LuSpool,*) typt3
       ELSE IF (LINE(1:4).EQ.'T3DE') THEN
       Read(LuSpool,*) typden
       ELSE IF (LINE(1:4).EQ.'ADAP') THEN
       Read(LuSpool,*) keysa
       if ((keysa.gt.4).or.(keysa.lt.0)) then
       keysa=0
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, keysa was changed to 0'
         end if
       end if
       if ((keysa.ne.0).and.(typden.eq.0)) then
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, typden is incompatible with SA'
         end if
       end if
       ELSE IF (LINE(1:4).EQ.'LOAD') THEN
       Read(LuSpool,*) filerst
       ELSE IF (LINE(1:4).EQ.'MACH') THEN
       Read(LuSpool,*) mchntyp,slim
       if ((mchntyp.lt.1).or.(mchntyp.gt.2)) then
       mchntyp=1
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, mchntyp was changed to 1'
         end if
       end if
       ELSE IF (LINE(1:4).EQ.'T3SH') THEN
       Read(LuSpool,*) shifto,shiftv
       ELSE IF (LINE(1:4).EQ.'PRIN') THEN
       Read(LuSpool,*) fullprint
       ELSE IF (LINE(1:4).EQ.'NOOP') THEN
       noop=1
       ELSE IF (LINE(1:4).EQ.'IOKE') THEN
       Read(LuSpool,*) iokey
       if ((iokey.lt.0).or.(iokey.gt.2)) then
       iokey=2
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, I/O key out of range'
         write(6,*) ' parameter iokey changed to 2'
         end if
       end if
       ELSE IF (LINE(1:4).EQ.'MHKE') THEN
       Read(LuSpool,*) mhkey
       if ((mhkey.lt.0).or.(mhkey.gt.2)) then
       mhkey=1
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, Matrix handling key is out of range'
         write(6,*) ' parameter iokey changed to 1'
         end if
       end if
       ELSE IF (LINE(1:4).EQ.'IJSE') THEN
       ijsegkey=1
       Read(LuSpool,*) symimin,imin,symjmin,jmin,                       &
     &                 symimax,imax,symjmax,jmax
       ELSE IF (LINE(1:4).EQ.'END ') THEN
       GOTO 7
       END IF
       GOTO 6
 7      CONTINUE
!
       Call Close_LuSpool(LuSpool)
       return
       end
