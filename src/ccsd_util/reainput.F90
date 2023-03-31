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
       subroutine reainput
!
!     this routine do:
!     1) read INPDAT file, produced by REORG with mul,nsym,noa,nob,nva,nvb,norb,eps
!     2) read input file for CCSD to read (parameters transported through cmm common)
!     ! title   - jobtitle
!     1-ntit rows with 72 characters
!     no default
!     ! ntit    - number of rows in jobtitle
!     ntit is limited to 10 !!! NOT ANY MORE !!! ntit is limited to 1
!     default (1)
!     ! maxiter - maximum number of iterations
!     default (10)
!     ! typt3   - type of T3 cpntribution
!     0 - CCSD
!     1 - CCSD+T(CCSD)
!     2 - CCSD(T)
!     default (0)
!     ! typden  - type of denominator (division of fok)
!     0 - diagonal elements
!     1 - average of faa and fbb
!     2 - orbital energies
!     default (0)
!     ! firstext- first iteration wnere extrapolation is used
!     (default-no) less than cycext is not possible
!     ! cycext  - cycle of extrapolation
!     (default-no) limited to 2-4
!     ! ccnonv  - energy convergence criterion
!     (default=1.0d-6)
!     ! keysa   - Spin adaptation key
!     0 - no adaptation
!     1 - T2 DDVV adaptation
!     2 - T2 DDVV + T1 DV adaptation
!     3 - full T1 and T2 adaptation (only for doublets)
!     4 - full T2 adaptation without SDVS (only for doublets)
!     (default=0)
!     ! keyrst  - restart key
!     0 - no saving restart informations (amplitudes)
!     1 - save amplitudes
!     2 - start from previous informations
!     (default=1)
!     ! filerst - name for restart informations file
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
!     ! fullprint - level of printing contrlo key
!     (default=0)
!     & noop  - no operation key
!     (default=no)
!     & iokey - I/O control key
!       1 - Fortran I/O system
!       2 - MOLCAS DA IO system
!     (default=2)
!     & mhkey - Matrix handling control key
!       1 - ESSL routines
!       2 - Fortran I/O system
!     & noccsd - key to supress CCSD run
!     (default=no)c     (default=1)
!     .....   - can be added
!
#include "ccsd1.fh"
!
!     help variables
!
       character*80 LINE
       integer nhelp
       integer f_iostat,f_recl
       logical is_error
!
!1    read INPDAT
!
       call molcas_open_ext2(1,'INPDAT','sequential','unformatted',     &
     &             f_iostat,.false.,f_recl,'unknown',is_error)
!       open (unit=1,file='INPDAT',form='unformatted')
       read (1) nactel,ispin,nsym,lsym,mmul,                            &
     &          noa,nob,nva,nvb,norb,eps,Escf
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
!4    define defaults
!
       maxiter=30
       typt3=0
       ntit=1
       typden=0
       yesext=0
       firstext=0
       cycext=0
       ccconv=1.0d-7
       keysa=0
       keyrst=1
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
       noccsd=0
!
!5    read input file
!
      LuSpool = 17
      Call SpoolInp(LuSpool)
      Rewind(LuSpool)
 5      Read(LuSpool,'(A80)') LINE
       CALL UPCASE(LINE)
       IF( INDEX(LINE,'&CCSDT') .EQ. 0 ) GOTO 5
       TITLE=' '
 6     Read(LuSpool,'(A80)') LINE
       IF(LINE(1:1).EQ.'*') GOTO 6
       CALL UPCASE(LINE)
!
       IF (LINE(1:4).EQ.'TITL') THEN
       Read(LuSpool,'(A72)') TITLE
       ELSE IF (LINE(1:4).EQ.'ITER') THEN
       Read(LuSpool,*) maxiter
       ELSE IF (LINE(1:4).EQ.'DENO') THEN
       Read(LuSpool,*) typden
       if ((typden.lt.0).or.(typden.gt.2)) then
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, Invalid type of denominators'
         write(6,*) ' parameter typden changed to 2'
         end if
       typden=2
       end if
       ELSE IF (LINE(1:4).EQ.'EXTR') THEN
       yesext=1
       Read(LuSpool,*) firstext,cycext
       if ((cycext.lt.2).or.(cycext.gt.4)) then
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, Size of DIIS procedure out of range'
         write(6,*) ' parameter cycext changed to 4'
         end if
       cycext=4
       end if
       if (firstext.lt.cycext) then
       firstext=cycext
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, First DIIS iteration is smaller '//   &
     &              'then DIIS size'
         write(6,*) ' parameter firstext was changed to:',firstext
         end if
       end if
       ELSE IF (LINE(1:4).EQ.'ACCU') THEN
       Read(LuSpool,*) ccconv
       ELSE IF (LINE(1:4).EQ.'ADAP') THEN
       Read(LuSpool,*) keysa
       if ((keysa.gt.4).or.(keysa.lt.0)) then
       keysa=0
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, Adaptation key out of range'
         write(6,*) ' parameter keysa changed to 0'
         end if
       end if
       if ((keysa.ne.0).and.(typden.eq.0)) then
       typden=2
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, typden is incompatible with SA'
         write(6,*) ' type of denominators changed to 2 - Orb. energies'
         end if
       end if
       ELSE IF (LINE(1:4).EQ.'REST') THEN
       Read(LuSpool,*) keyrst
       if ((keyrst.lt.0).or.(keyrst.gt.2)) then
       keyrst=1
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, Restart key out of range'
         write(6,*) ' parameter keyrst changed to 1'
         end if
       end if
       Read(LuSpool,*) filerst
       ELSE IF (LINE(1:4).EQ.'MACH') THEN
       Read(LuSpool,*) mchntyp,slim
       if ((mchntyp.lt.1).or.(mchntyp.gt.2)) then
       mchntyp=1
         if (fullprint.ge.0) then
         write(6,*) ' Warning!!!, Machinetype out of range'
         write(6,*) ' parameter mchtyp changed to 1'
         end if
       end if
       ELSE IF (LINE(1:4).EQ.'SHIF') THEN
       Read(LuSpool,*) shifto,shiftv
       ELSE IF (LINE(1:4).EQ.'PRIN') THEN
       Read(LuSpool,*) fullprint
       if ((fullprint.lt.0).or.(fullprint.gt.3)) then
       fullprint=0
       write(6,*) ' Warning!!!, Printing key out of range'
       write(6,*) ' parameter fullprint changed to 0'
       end if
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
         write(6,*) ' Warning!!!, Matrix handling key out of range'
         write(6,*) ' parameter mhkey changed to 1'
         end if
       end if
       ELSE IF (LINE(1:4).EQ.'NOSD') THEN
       noccsd=1
       ELSE IF (LINE(1:4).EQ.'END ') THEN
       GOTO 7
       END IF
       GOTO 6
 7      CONTINUE
!
       Call Close_LuSpool(LuSpool)
       return
       end
