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
c     this fil contains following routines:
c     initfiles
c     reainput
c     wrhead
c     divfok
c     fokunpck1
c     fokunpck2
c     fokunpck3
c     fokunpck4
c     fokunpck5
c     reaintsta
c
c     -----------------------
c
       subroutine initfiles (lenght,lenv,lenn)
c
c     this routine distribute work space WRK for required files
c     for fix mediates it defines also mapd and mapi, for help mediates
c     it estimates their lenght and distribute WRK (i.e. def poss0 parameters)
c     !N.B. This routine cannot run with +OP2 level
c
c       lenght  -  total length of all work space needed
c       lenv    -  length of V - type array
c       lenn    -  length of N - type array
c
#include "ccsd1.fh"
#include "ccsd2.fh"
c
       integer lenght,lenv,lenn
c
c     help variable
c
       integer posst,symp,symq,sympq,symr,syms
       integer lenghtv,lenghtm,lenghth,lenghtn
       integer maxnoa,maxnvb,maxnorb
       integer maxov(1:8)
c
c1    maps and possitions for fix mediated
c
c1.0  maps for DP - diagonal part
c     N.B. DP has one degree of freedom, while other 1 index has none
c     DP1 - dp(p)a
c     DP2 - dp(p)b
c
       do symp=1,nsym
       do symq=1,nsym
       do symr=1,nsym
       mapidp1(symp,symq,symr)=0
       mapidp2(symp,symq,symr)=0
       end do
       end do
       end do
c
       posst=1
c
       possdp10=posst
       mapddp1(0,1)=5
       mapddp1(0,2)=0
       mapddp1(0,3)=0
       mapddp1(0,4)=0
       mapddp1(0,5)=nsym
       mapddp1(0,6)=0
c
       do symp=1,nsym
       mapddp1(symp,1)=posst
       mapddp1(symp,2)=norb(symp)
       mapddp1(symp,3)=symp
       mapddp1(symp,4)=1
       mapddp1(symp,5)=1
       mapddp1(symp,6)=1
       mapidp1(symp,1,1)=symp
       posst=posst+norb(symp)
       end do
c
       possdp20=posst
       mapddp2(0,1)=5
       mapddp2(0,2)=0
       mapddp2(0,3)=0
       mapddp2(0,4)=0
       mapddp2(0,5)=nsym
       mapddp2(0,6)=0
c
       do symp=1,nsym
       mapddp2(symp,1)=posst
       mapddp2(symp,2)=norb(symp)
       mapddp2(symp,3)=symp
       mapddp2(symp,4)=1
       mapddp2(symp,5)=1
       mapddp2(symp,6)=1
       mapidp2(symp,1,1)=symp
       posst=posst+norb(symp)
       end do
c
c
c1.1  maps for T1
c     T11 - t1oaa(a,i)
c     T12 - t1obb(a,i)
c     T13 - t1naa(a,i)
c     T14 - t1nbb(a,i)
c
       posst110=posst
       call grc0 (2,0,3,1,0,0,1,
     & posst110,posst,mapdt11,mapit11)
       posst120=posst
       call grc0 (2,0,4,2,0,0,1,
     & posst120,posst,mapdt12,mapit12)
       posst130=posst
       call grc0 (2,0,3,1,0,0,1,
     & posst130,posst,mapdt13,mapit13)
       posst140=posst
       call grc0 (2,0,4,2,0,0,1,
     & posst140,posst,mapdt14,mapit14)
c
c
c1.2  maps for F1
c     F11 - FI(a,e)aa
c     F12 - FI(a,e)bb
c
       possf110=posst
       call grc0 (2,0,3,3,0,0,1,
     & possf110,posst,mapdf11,mapif11)
       possf120=posst
       call grc0 (2,0,4,4,0,0,1,
     & possf120,posst,mapdf12,mapif12)
c
c
c1.3  maps for F2
c     F21 - FII(m,i)aa
c     F22 - FII(m,i)bb
c
       possf210=posst
       call grc0 (2,0,1,1,0,0,1,
     & possf210,posst,mapdf21,mapif21)
       possf220=posst
       call grc0 (2,0,2,2,0,0,1,
     & possf220,posst,mapdf22,mapif22)
c
c
c1.4  maps for F3
c     F31 - FIII(e,m)aa
c     F32 - FIII(e,m)bb
c
       possf310=posst
       call grc0 (2,0,3,1,0,0,1,
     & possf310,posst,mapdf31,mapif31)
       possf320=posst
       call grc0 (2,0,4,2,0,0,1,
     & possf320,posst,mapdf32,mapif32)
c
c
c1.5  maps for FK
c     FK1 - f(a,b)aa
c     FK2 - f(a,b)bb
c     FK3 - f(a,i)aa
c     FK4 - f(a,i)bb
c     FK5 - f(i,j)aa
c     FK6 - f(i,j)bb
c
       possfk10=posst
       call grc0 (2,0,3,3,0,0,1,
     & possfk10,posst,mapdfk1,mapifk1)
       possfk20=posst
       call grc0 (2,0,4,4,0,0,1,
     & possfk20,posst,mapdfk2,mapifk2)
       possfk30=posst
       call grc0 (2,0,3,1,0,0,1,
     & possfk30,posst,mapdfk3,mapifk3)
       possfk40=posst
       call grc0 (2,0,4,2,0,0,1,
     & possfk40,posst,mapdfk4,mapifk4)
       possfk50=posst
       call grc0 (2,0,1,1,0,0,1,
     & possfk50,posst,mapdfk5,mapifk5)
       possfk60=posst
       call grc0 (2,0,2,2,0,0,1,
     & possfk60,posst,mapdfk6,mapifk6)
c
c
c1.6  maps for T2
c     T21 - t2n(ab,ij)aaaa
c     T22 - t2n(ab,ij)bbbb
c     T33 - t2n(a,b,i,j)abab
c
       posst210=posst
       call grc0 (4,4,3,3,1,1,1,
     & posst210,posst,mapdt21,mapit21)
       posst220=posst
       call grc0 (4,4,4,4,2,2,1,
     & posst220,posst,mapdt22,mapit22)
       posst230=posst
       call grc0 (4,0,3,4,1,2,1,
     & posst230,posst,mapdt23,mapit23)
c
c
c1.7  maps for W0
c     W01 - <mn||ij>aaaa
c     W02 - <mn||ij>bbbb
c     W03 - <mn||ij>abab
c
       possw010=posst
       call grc0 (4,4,1,1,1,1,1,
     & possw010,posst,mapdw01,mapiw01)
       possw020=posst
       call grc0 (4,4,2,2,2,2,1,
     & possw020,posst,mapdw02,mapiw02)
       possw030=posst
       call grc0 (4,0,1,2,1,2,1,
     & possw030,posst,mapdw03,mapiw03)
c
c
c1.8  maps for W1
c     W11 - <ie||mn>aaaa
c     W12 - <ie||mn>bbbb
c     W13 - <ie||mn>abab
c     W14 - <ie||mn>baab
c
       possw110=posst
       call grc0 (4,3,1,3,1,1,1,
     & possw110,posst,mapdw11,mapiw11)
       possw120=posst
       call grc0 (4,3,2,4,2,2,1,
     & possw120,posst,mapdw12,mapiw12)
       possw130=posst
       call grc0 (4,0,1,4,1,2,1,
     & possw130,posst,mapdw13,mapiw13)
       possw140=posst
       call grc0 (4,0,2,3,1,2,1,
     & possw140,posst,mapdw14,mapiw14)
c
c
c2    for help files mapps are irrelevant,
c     here only estimation of maximal lenght is done to
c     define poss0 of help files
c     we have:
c     four V files - of vvoo type
c     four M files - of vvo  type
c     four H files - of voo  type
c     one  N file  - of nn   type
c
c2.*  def max{noa}, max{norb} ,max{nvb}, maxov(isym)=max{noa(isym),nvb(isym)}
c
       maxnoa=noa(1)
       maxnvb=nvb(1)
       maxnorb=norb(1)
       do 100 symp=1,nsym
       if (noa(symp).gt.maxnoa) then
       maxnoa=noa(symp)
       end if
       if (norb(symp).gt.maxnorb) then
       maxnorb=norb(symp)
       end if
       if (nvb(symp).gt.maxnvb) then
       maxnvb=nvb(symp)
       end if
       if (nvb(symp).gt.noa(symp)) then
       maxov(symp)=nvb(symp)
       else
       maxov(symp)=noa(symp)
       end if
 100    continue
c
c2.*  def lenghts of V,M,H and N fils
c
       lenghtv=0
       lenghtm=0
       lenghth=0
       lenghtn=0
c
       do 200 symp=1,nsym
c     symq is not known for N file
c     instead of norb(symr) maxnorb will be used so that reallenght<=lenght
       lenghtn=lenghtn+norb(symp)*maxnorb
       do 200 symq=1,nsym
       sympq=mmul(symp,symq)
c     symr is not known for M and H files
c     instead of noa(symr) maxnoa will be used so that reallenght<=lenght
       lenghtm=lenghtm+maxov(symp)*maxov(symq)*maxnoa
       lenghth=lenghth+maxov(symp)*noa(symq)*maxnoa
       do 200 symr=1,nsym
       syms=mmul(sympq,symr)
       lenghtv=lenghtv+maxov(symp)*maxov(symq)*noa(symr)*noa(syms)
 200    continue
c
c2.1  V - files
c
       possv10=posst
       posst=posst+lenghtv
       possv20=posst
       posst=posst+lenghtv
       possv30=posst
       posst=posst+lenghtv
       possv40=posst
       posst=posst+lenghtv
       lenv=lenghtv
c
c2.2  M - files
c
       possm10=posst
       posst=posst+lenghtm
       possm20=posst
       posst=posst+lenghtm
       possm30=posst
       posst=posst+lenghtm
       possm40=posst
       posst=posst+lenghtm
c
c2.3  H - files
c
       possh10=posst
       posst=posst+lenghth
       possh20=posst
       posst=posst+lenghth
       possh30=posst
       posst=posst+lenghth
       possh40=posst
       posst=posst+lenghth
c
c2.4  N,P - files
c
       possn0=posst
       posst=posst+lenghtn
       possp0=posst
       posst=posst+lenghtn
       lenn=lenghtn
c
c2.5  dedlare space for help matrix D in for matrix multiplication C=AT*B if
c     mchntyp=2
c
       if (mchntyp.eq.2) then
       possd0=posst
       if (maxnoa.le.maxnvb) then
       posst=posst+maxnoa*maxnoa*maxnvb*maxnvb
       else
       posst=posst+maxnoa*maxnoa*maxnoa*maxnoa
       end if
       end if
c
c2.6   def size of Work space
       lenght=posst-1
c
       return
       end
c
c     ----------------------
c
       subroutine reainput
c
c     this routine do:
c     1) read INPDAT file, produced by REORG with mul,nsym,noa,nob,nva,nvb,norb,eps
c     2) read input file for CCSD to read (parameters transported through cmm common)
c     ! title   - jobtitle
c     1-ntit rows with 72 characters
c     no default
c     ! ntit    - number of rows in jobtitle
c     ntit is limited to 10 !!! NOT ANY MORE !!! ntit is limited to 1
c     default (1)
c     ! maxiter - maximum number of iterations
c     default (10)
c     ! typt3   - type of T3 cpntribution
c     0 - CCSD
c     1 - CCSD+T(CCSD)
c     2 - CCSD(T)
c     default (0)
c     ! typden  - type of denominator (division of fok)
c     0 - diagonal elements
c     1 - average of faa and fbb
c     2 - orbital energies
c     default (0)
c     ! firstext- first iteration wnere extrapolation is used
c     (default-no) less than cycext is not possible
c     ! cycext  - cycle of extrapolation
c     (default-no) limited to 2-4
c     ! ccnonv  - energy convergence criterion
c     (default=1.0d-6)
c     ! keysa   - Spin adaptation key
c     0 - no adaptation
c     1 - T2 DDVV adaptation
c     2 - T2 DDVV + T1 DV adaptation
c     3 - full T1 and T2 adaptation (only for doublets)
c     4 - full T2 adaptation without SDVS (only for doublets)
c     (default=0)
c     ! keyrst  - restart key
c     0 - no saving restart informations (amplitudes)
c     1 - save amplitudes
c     2 - start from previous informations
c     (default=1)
c     ! filerst - name for restart informations file
c     (default=RSTART)
c     ! mchntyp - type of machine in matrix multiplication
c     1 - C=A*B is faster or comparable with C=AT*B
c     2 - C=AT*B is faster
c     (default=1)
c     ! slim    - limitation for usieng C=AT*B
c     no default (suitable=2.0d0)
c     ! shifhto - shift for occupied
c     (default=0.0)
c     ! shifhtv - shift for virtuals
c     (default=0.0)
c     ! maxspace - maximal allowed work space
c     (default=0 - unlimited)
c     ! fullprint - level of printing contrlo key
c     (default=0)
c     & noop  - no operation key
c     (default=no)
c     & iokey - I/O control key
c       1 - Fortran I/O system
c       2 - MOLCAS DA IO system
c     (default=2)
c     & mhkey - Matrix handling control key
c       1 - ESSL routines
c       2 - Fortran I/O system
c     & noccsd - key to supress CCSD run
c     (default=no)c     (default=1)
c     .....   - can be added
c
#include "ccsd1.fh"
c
c     help variables
c
       character*80 LINE
       integer nhelp
       integer f_iostat,f_recl
       logical is_error
c
c1    read INPDAT
c
       call molcas_open_ext2(1,'INPDAT','sequential','unformatted',
     &             f_iostat,.false.,f_recl,'unknown',is_error)
c       open (unit=1,file='INPDAT',form='unformatted')
       read (1) nactel,ispin,nsym,lsym,mmul,
     &          noa,nob,nva,nvb,norb,eps,Escf
       close (1)
c
c2    def dimm
c
       do 10 nhelp=1,nsym
c
       dimm(1,nhelp)=noa(nhelp)
       dimm(2,nhelp)=nob(nhelp)
       dimm(3,nhelp)=nva(nhelp)
       dimm(4,nhelp)=nvb(nhelp)
       dimm(5,nhelp)=norb(nhelp)
c
 10     continue
c
c3    define nshf
c
       do 20 nhelp=1,maxorb
       nshf(nhelp)=(nhelp-1)*(nhelp-2)/2
 20     continue

c
c4    define defaults
c
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
CGG       fullprint=0
       noop=0
       iokey=1
       mhkey=1
       noccsd=0
c
c5    read input file
c
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
c
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
         write(6,*) ' Warning!!!, First DIIS iteration is smaller '//
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
c
       Call Close_LuSpool(LuSpool)
       return
       end
c
c     ----------------------
c
c
c
c     ----------------------
c
       subroutine wrhead
c
c     this routine write head of ther output file
c     parameters transported through cmm common
c
#include "ccsd1.fh"
c
c     help variables
c
       integer nhelp,nhelp1,nhelp2
c
c1    write header and title
c
*      Call HelloPN
      if (noop.eq.1) then
      write(6,'(6X,A)') ' No Operation is required'
      write(6,'(6X,A)') ' Happy Landing'
      Call Finish(0)
      end if
c
      if ( title(1:1).ne.' ' ) then
          write(6,*)
          write(6,'(6X,112A1)') ('*',i=1,112)
          write(6,'(6X,A1,19X,A72,19X,A1)') '*',title,'*'
          write(6,'(6X,112A1)') ('*',i=1,112)
      end if
c
c2    write occupations and matrix multiplication tablaux
c
       Write(6,*)
       Write(6,*)
       Write(6,'(6X,A)')'Wave function specifications:'
       Write(6,'(6X,A)')'-----------------------------'
       Write(6,*)
       Write(6,'(6X,A,T45,I6)') 'Spin mutiplicity',ispin
       Write(6,'(6X,A,T45,I6)')   'State symmetry',lsym
       Write(6,*)
       Write(6,'(6X,A)')'Orbital specifications:'
       Write(6,'(6X,A)')'-----------------------'
       Write(6,*)
       Write(6,'(6X,A,T47,8I4)')
     & 'Symmetry species',
     & (nhelp,nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'Total no. of orbitals',
     & (norb(nhelp),nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'No. of occupied orbitals with alpha spin',
     & (noa(nhelp),nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'No. of occupied orbitals with beta spin',
     & (nob(nhelp),nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'No. of virtual orbitals with alpha spin',
     & (nva(nhelp),nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'No. of virtual orbitals with beta spin',
     & (nvb(nhelp),nhelp=1,nsym)
       Write(6,*)
           if (fullprint.gt.1) then
           write(6,*)
           write(6,16) nsym
 16        format (' NUMBER OF IRREPS             :',i3)
           write(6,17)
17         format (' MATRIX MULTIPLICATION TABLE')
           write(6,*)
           write(6,18) (nhelp,nhelp=1,nsym)
18         format (' IRREP #',9x,8(i3,2x))
           write(6,*)
           do 19 nhelp1=1,nsym
           write(6,20) nhelp1,(mmul(nhelp1,nhelp),nhelp=1,nsym)
19         continue
20         format (' IRREP #',i3,6x,8(i3,2x))
           write(6,*)
           end if
       Write(6,*)
       Write(6,'(6X,A)')'Methods and options:'
       Write(6,'(6X,A)')'--------------------'
       Write(6,*)
       Write(6,'(6X,A,T45,I3)') 'Max no. of iterations',maxiter
       if (typden.eq.0) then
         Write(6,'(6X,A,T45,A)') 'Type of denominators',
     &   'diagonal Fock matrix elements'
       else if (typden.eq.1) then
         Write(6,'(6X,A,T45,A)') 'Type of denominators',
     &   'spin averaged diagonal Fock matrix elements'
       else
         Write(6,'(6X,A,T45,A)') 'Type of denominators',
     &   'orbital energies'
       end if
       Write(6,'(6X,A,T45,F22.14)') 'energy convergence criterium',
     & ccconv
       if (yesext.ne.0) then
       end if
c
c6    write orbital energies per symmetry
c
       write(6,*)
       if (fullprint.gt.0) then
       write(6,51)
 51     format (' LIST OF ORBITAL ENERGIES')
       write(6,*)
c
       nhelp2=1
       do nhelp=1,nsym
       write(6,54) nhelp
 54     format (' IRREDUCIBLE REPRESENTATION NO:',i2)
       do nhelp1=1,norb(nhelp)
       write(6,55) nhelp1,eps(nhelp2)
 55     format (' ORBITAL NO:',i3,5x,f16.10)
       nhelp2=nhelp2+1
       end do
       end do
       write(6,*)
       end if
c
c6    extrapolation parameters
c
       if (yesext.eq.0) then
       write(6,61)
 61     format (' DIIS EXTRAPOLATION USED     : NO')
       else
        write(6,62)
 62      format (' DIIS EXTRAPOLATION USED     : YES')
        write(6,63) firstext
 63      format (' FIRST ITERATION OF EXT.     :',i3)
        write(6,64) cycext
 64      format (' EXTRAPOLATION CYCLE         :',i3)
       end if
       write(6,*)
c
c7    convergence criterion
c
CFUE       write(6,71) ccconv
CFUE 71     format (' ENERGY CONVERGENCE CRITER.  : ',d10.5)
c
c8    spin adaptation
c
       if (keysa.eq.0) then
       write(6,81)
       else if (keysa.eq.1) then
       write(6,82)
       else if (keysa.eq.2) then
       write(6,83)
       else if (keysa.eq.3) then
       write(6,84)
       else if (keysa.eq.4) then
       write(6,85)
       end if
 81     format (' SPIN ADAPTATION             : NONE ')
 82     format (' SPIN ADAPTATION             : T2 DDVV ')
 83     format (' SPIN ADAPTATION             : T2 DDVV + T1 DV ')
 84     format (' SPIN ADAPTATION             : T1 AND T2 FULL ')
 85     format (' SPIN ADAPTATION             : T2 FULL WITHOUT SDVS')
c
c9    restart status
       if (keyrst.eq.0) then
       write(6,91)
       else if (keyrst.eq.1) then
       write(6,92) filerst
       else if (keyrst.eq.2) then
       write(6,93) filerst
       end if
 91     format (' RESTART STATUS              : NONE ')
 92     format (' RST. INF. WILL BE SAVED IN  : ',a6)
 93     format (' RST. INF. WILL BE LOAD FROM : ',a6)
c
c10   write matrix multiplication performance
c
       if (mchntyp.eq.1) then
       write(6,101)
       else if (mchntyp.eq.2) then
       write(6,102) slim
       end if
 101    format (' PREFERENCE MATRIX MULT.     : NORMAL')
 102    format (' PREFERENCE MATRIX MULT.     : TRANSP ; LIMIT =',d12.5)
c
c11   write denominator shifts
       write(6,111) shifto
       write(6,112) shiftv
 111    format (' DENOMINATOR SHIFT FOR OCC.  : ',d12.5)
 112    format (' DENOMINATOR SHIFT FOR VIRT. : ',d12.5)
c
c12   write workspace
c
c      if (maxspace.eq.0) then
c      write(6,121)
c      else
c      write(6,122) maxspace
c      end if
c121    format (' MAXIMAL ALLOWED WORK SPACE  : UNLIMITED')
c122    format (' MAXIMAL ALLOWED WORK SPACE  : ',i10)
c
c13   level of printing
c
       if (fullprint.eq.0) then
       write(6,131)
 131    format (' LEVEL OF OUTPUT PRINTING    : MINIMAL')
       else if (fullprint.eq.1) then
       write(6,132)
 132    format (' LEVEL OF OUTPUT PRINTING    : MEDIUM')
       else if (fullprint.eq.2) then
       write(6,133)
 133    format (' LEVEL OF OUTPUT PRINTING    : MAXIMAL')
       else if (fullprint.eq.3) then
       write(6,134)
 134    format (' LEVEL OF OUTPUT PRINTING    : DEBUG')
       end if
c
c14   I/O handling
c
       if (iokey.eq.1) then
       write(6,141)
 141    format (' INPUT/OUTPUT HANDLING       : Standard SQ ')
       else
       write(6,142)
 142    format (' INPUT/OUTPUT HANDLING       : Molcas4  DA ')
       end if
c
c15    Matrix handling
c
       if (mhkey.eq.1) then
       write(6,151)
 151    format (' MATRIX OPERATIONS           : ESSL        ')
       else
       write(6,152)
 152    format (' MATRIX OPERATIONS           : Fortran code')
       end if
c
c*    can be added
c
       write(6,*)
       write(6,*)
c
       return
       end
c
c     -----------------------------------------
c
       subroutine divfok (wrk,wrksize,
     & mapdfa,mapifa,possfa0,mapdfb,mapifb,possfb0,
     & mapdfk1,mapifk1,possfk10,mapdfk2,mapifk2,possfk20,
     & mapdfk3,mapifk3,possfk30,mapdfk4,mapifk4,possfk40,
     & mapdfk5,mapifk5,possfk50,mapdfk6,mapifk6,possfk60,
     & mapddp1,mapidp1,possdp10,mapddp2,mapidp2,possdp20,rc)
c
c     this routine divide fok(p,q) -> fk(a,b) + fk(a,i) + f(i,j) + dp(p)
c     to diagonal part and rest
c
c     mapd and mapi for:
c     fa,fb - fok(p,q)aa,bb
c     fk1-6 - f(ab)aa,f(ab)bb,f(ai)aa,f(ai)bb,f(ij)aa,f(ij)bb
c     dp1,2 - diagonal part dp(p)a,b
c     rc    - return (error) code
c
c
#include "ccsd1.fh"
#include "wrk.fh"
       integer rc
c
c1    maps for FOKA,FOKB
c
       integer mapdfa(0:512,1:6)
       integer mapifa(1:8,1:8,1:8)
       integer possfa0
c
       integer mapdfb(0:512,1:6)
       integer mapifb(1:8,1:8,1:8)
       integer possfb0
c
c2    maps for FK
c     FK1 - f(a,b)aa
c     FK2 - f(a,b)bb
c     FK3 - f(a,i)aa
c     FK4 - f(a,i)bb
c     FK5 - f(i,j)aa
c     FK6 - f(i,j)bb
c
       integer mapdfk1(0:512,1:6)
       integer mapifk1(1:8,1:8,1:8)
       integer possfk10
c
       integer mapdfk2(0:512,1:6)
       integer mapifk2(1:8,1:8,1:8)
       integer possfk20
c
       integer mapdfk3(0:512,1:6)
       integer mapifk3(1:8,1:8,1:8)
       integer possfk30
c
       integer mapdfk4(0:512,1:6)
       integer mapifk4(1:8,1:8,1:8)
       integer possfk40
c
       integer mapdfk5(0:512,1:6)
       integer mapifk5(1:8,1:8,1:8)
       integer possfk50
c
       integer mapdfk6(0:512,1:6)
       integer mapifk6(1:8,1:8,1:8)
       integer possfk60
c
c
c3    maps for DP - diagonal part
c     DP1 - dp(p)a
c     DP2 - dp(p)b
c
       integer mapddp1(0:512,1:6)
       integer mapidp1(1:8,1:8,1:8)
       integer possdp10
c
       integer mapddp2(0:512,1:6)
       integer mapidp2(1:8,1:8,1:8)
       integer possdp20
c
c     help variables
c
       integer symp,rc1
       integer iifoka,iifokb,iifok,iifaa,iifai,iifii,iidpa,iidpb,iidp
       integer possfoka,possfokb,possfok,possfaa,possfai,possfii
       integer possdpa,possdpb,possdp
c
       rc=0
c
c1    define dp
c
       do 500 symp=1,nsym
c
       iidpa=mapidp1(symp,1,1)
       possdpa=mapddp1(iidpa,1)
       iidpb=mapidp2(symp,1,1)
       possdpb=mapddp2(iidpb,1)
       iifoka=mapifa(symp,1,1)
       possfoka=mapdfa(iifoka,1)
       iifokb=mapifb(symp,1,1)
       possfokb=mapdfb(iifokb,1)
c
       if (norb(symp).gt.0) then
       call fokunpck5 (symp,wrk(possfoka),wrk(possfokb),
     & wrk(possdpa),wrk(possdpb),norb(symp),rc1)
       end if
c
 500    continue
c
c2    define faa,fai,fii
c
       do 1000 symp=1,nsym
       if (norb(symp).eq.0) goto 1000
c
c2.1  alpha case
c
       iifok=mapifa(symp,1,1)
       iifaa=mapifk1(symp,1,1)
       iifai=mapifk3(symp,1,1)
       iifii=mapifk5(symp,1,1)
       iidp=mapidp1(symp,1,1)
c
       possfok=mapdfa(iifok,1)
       possfaa=mapdfk1(iifaa,1)
       possfai=mapdfk3(iifai,1)
       possfii=mapdfk5(iifii,1)
       possdp=mapddp1(iidp,1)
c
       call fokunpck1 (wrk(possfok),wrk(possdp),norb(symp))
       if (nva(symp).gt.0) then
       call fokunpck2 (wrk(possfok),wrk(possfaa),norb(symp),nva(symp),
     &                 noa(symp))
       end if
       if ((noa(symp)*nva(symp)).gt.0) then
       call fokunpck3 (wrk(possfok),wrk(possfai),norb(symp),nva(symp),
     &                 noa(symp))
       end if
       if (noa(symp).gt.0) then
       call fokunpck4 (wrk(possfok),wrk(possfii),norb(symp),noa(symp))
       end if
c
c2.2  alpha case
c
       iifok=mapifb(symp,1,1)
       iifaa=mapifk2(symp,1,1)
       iifai=mapifk4(symp,1,1)
       iifii=mapifk6(symp,1,1)
       iidp=mapidp2(symp,1,1)
c
       possfok=mapdfb(iifok,1)
       possfaa=mapdfk2(iifaa,1)
       possfai=mapdfk4(iifai,1)
       possfii=mapdfk6(iifii,1)
       possdp=mapddp2(iidp,1)
c
       call fokunpck1 (wrk(possfok),wrk(possdp),norb(symp))
       if (nvb(symp).gt.0) then
       call fokunpck2 (wrk(possfok),wrk(possfaa),norb(symp),nvb(symp),
     &                 nob(symp))
       end if
       if ((nob(symp)*nvb(symp)).gt.0) then
       call fokunpck3 (wrk(possfok),wrk(possfai),norb(symp),nvb(symp),
     &                 nob(symp))
       end if
       if (nob(symp).gt.0) then
       call fokunpck4 (wrk(possfok),wrk(possfii),norb(symp),nob(symp))
       end if
c
 1000   continue
c
       return
c Avoid unused argument warnings
       if (.false.) then
         call Unused_integer(possfa0)
         call Unused_integer(possfb0)
         call Unused_integer(possfk10)
         call Unused_integer(possfk20)
         call Unused_integer(possfk30)
         call Unused_integer(possfk40)
         call Unused_integer(possfk50)
         call Unused_integer(possfk60)
         call Unused_integer(possdp10)
         call Unused_integer(possdp20)
       end if
       end
c
c     -------------
c
       subroutine fokunpck1 (fok,dp,dimfok)
c
c     this routine do Fok = Fok - dp
c     fok    - Fok matrix (I/O)
c     dp     - Diagonal part vector (I)
c     dimfok - dimension for Fok matrix - norb (I)
c
       integer dimfok
c
       real*8 fok(1:dimfok,1:dimfok)
       real*8 dp(1:dimfok)
c
c     help variables
c
       integer p
c
c1    substract dp from Fok
       do 100 p=1,dimfok
       fok(p,p)=fok(p,p)-dp(p)
 100    continue
c
       return
       end
c
c     -------------
c
       subroutine fokunpck2 (fok,faa,dimfok,dimfa,shift)
c
c     this routine distribute (Fok - dp) to Faa
c     fok    - Fok matrix (I)
c     faa    - Faa matrix (O)
c     dimfok - dimension for Fok matrix - norb (I)
c     dimfa  - dimension of virtuals - nv (I)
c
       integer dimfok,dimfa,shift
c
       real*8 fok(1:dimfok,1:dimfok)
       real*8 faa(1:dimfa,1:dimfa)
c
c     help variables
c
       integer a,b
c
c1    distribute Fok to Faa
       do 200 b=1,dimfa
       do 200 a=1,dimfa
       faa(a,b)=fok(shift+a,shift+b)
 200    continue
c
       return
       end
c
c     -------------
c
       subroutine fokunpck3 (fok,fai,dimfok,dimfa,dimfi)
c
c     this routine distribute (Fok - dp) to Fai
c     fok    - Fok matrix (I)
c     fai    - Fai matrix (O)
c     dimfok - dimension for Fok matrix - norb (I)
c     dimfa  - dimension of virtuals - nv (I)
c     dimfi  - dimension of occupied - no (I)
c
       integer dimfok,dimfa,dimfi
c
       real*8 fok(1:dimfok,1:dimfok)
       real*8 fai(1:dimfa,1:dimfi)
c
c     help variables
c
       integer a,i
c
c1    distribute Fok to Fai
       do 300 i=1,dimfi
       do 300 a=1,dimfa
       fai(a,i)=fok(dimfi+a,i)
 300    continue
c
       return
       end
c
c     -------------
c
       subroutine fokunpck4 (fok,fii,dimfok,dimfi)
c
c     this routine distribute (Fok - dp) to Fii
c     fok    - Fok matrix (I)
c     fii    - Fii matrix (O)
c     dimfok - dimension for Fok matrix - norb (I)
c     dimfi  - dimension of occupied - no (I)
c
       integer dimfok,dimfi
c
       real*8 fok(1:dimfok,1:dimfok)
       real*8 fii(1:dimfi,1:dimfi)
c
c     help variables
c
       integer i,j
c
c1    distribute Fok to Fii
       do 400 j=1,dimfi
       do 400 i=1,dimfi
       fii(i,j)=fok(i,j)
 400    continue
c
       return
       end
c
c     -------------
c
       subroutine fokunpck5 (symp,foka,fokb,dpa,dpb,dimfok,rc)
c
c     this routine produce dpa,dpb from foka,fokb
c     for some cases
c     shifto,shiftv will be also added
c
c     symp   - symmtry of this block
c     foka   - Fok aa matrix (I)
c     fokb   - Fok bb matrix (I)
c     dpa    - Diagonal part alfa vector (O)
c     dpa    - Diagonal part beta vector (O)
c     dimfok - dimension for Fok matrix - norb (I)
c     rc     - return (error) code
c
#include "ccsd1.fh"
       integer symp,dimfok,rc
c
       real*8 foka(1:dimfok,1:dimfok)
       real*8 fokb(1:dimfok,1:dimfok)
       real*8 dpa(1:dimfok)
       real*8 dpb(1:dimfok)
c
c     help variables
c
       integer p,nhelp1,nhelp2
c
       rc=0
c
       if (typden.eq.0) then
c1    diagonal elements are required
c
       do 100 p=1,dimfok
       dpa(p)=foka(p,p)
       dpb(p)=fokb(p,p)
 100    continue
c
       else if (typden.eq.1) then
c2    (faa+fbb)/2 are required
c
       do 200 p=1,dimfok
       dpa(p)=(foka(p,p)+fokb(p,p))/2
       dpb(p)=dpa(p)
 200    continue
c
       else if (typden.eq.2) then
c3    orbital energies are required
c
c3.1  def shift
       if (symp.eq.1) then
       nhelp1=0
       else
       nhelp1=0
       do 300 nhelp2=1,symp-1
       nhelp1=nhelp1+norb(nhelp2)
 300    continue
       end if
c
c3.2  map oe to dp
       do 400 p=1,dimfok
       dpa(p)=eps(nhelp1+p)
       dpb(p)=eps(nhelp1+p)
 400    continue
c
       else
c     RC=1 : invalid key (NCI/Stup)
       rc=1
       end if
c
       if ((keysa.eq.3).or.(keysa.eq.4)) then
c     for full adaptation scheme only D and V orbitals are shifted
c
       do 501 p=1,nob(symp)
       dpa(p)=dpa(p)-shifto
       dpb(p)=dpb(p)-shifto
 501    continue
c
       do 502 p=1+noa(symp),norb(symp)
       dpa(p)=dpa(p)+shiftv
       dpb(p)=dpb(p)+shiftv
 502    continue
c
       else
c     for other schemes all orbitals are shifted
c
       do 511 p=1,noa(symp)
       dpa(p)=dpa(p)-shifto
 511    continue
c
       do 512 p=1,nob(symp)
       dpb(p)=dpb(p)-shifto
 512    continue
c
       do 513 p=1+noa(symp),norb(symp)
       dpa(p)=dpa(p)+shiftv
 513    continue
c
       do 514 p=1+nob(symp),norb(symp)
       dpb(p)=dpb(p)+shiftv
 514    continue
c
       end if
c
        if (fullprint.ge.2) then
        write (6,*) ' Diagonal part Dp aa, bb for irrep: ',symp
        do p=1,norb(symp)
        write (6,99) p,dpa(p),dpb(p)
99      format (2x,i4,2(f20.14,2x))
        end do
        end if
c
       return
       end
c
c     -------------------------------------
c
       subroutine reaintsta (wrk,wrksize)
c
c     this routine read integral file INTSTA (reorg), which contains
c     following integrals: foka,fokb,
c     <kl||ij>aaaa,<kl||ij>bbbb,<kl||ij>abab
c     <ka||ij>aaaa,<ka||ij>bbbb,<ka||ij>abab,<ka||ij>baab
c     <ab||ij>aaaa,<ab||ij>bbbb,<ab||ij>abab
c
c     two electron integrals are readed to their fix files,
c     foka,fokb are readed to N,P help files
c
c     use and destroy : V1-3, N,P
c
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "filemgr.fh"
#include "wrk.fh"
c
c     help variables
c
       integer lunsta,rc,f_recl,f_iostat
       logical is_error
c
c*    open INTSTA file
       lunsta=1
       if (iokey.eq.1) then
c      Fortran IO
       call molcas_open_ext2(lunsta,'INTSTA','sequential','unformatted',
     &                     f_iostat,.false.,f_recl,'unknown',is_error)
c       open (unit=lunsta,file='INTSTA',form='unformatted')
c
       else
c      MOLCAS IO
       call daname (lunsta,'INTSTA')
       daddr(lunsta)=0
       end if
c
c1    read foka to N
       call getmediate (wrk,wrksize,
     & lunsta,possn0,mapdn,mapin,rc)
c
c2    read fokb to P
       call getmediate (wrk,wrksize,
     & lunsta,possp0,mapdp,mapip,rc)
c
c
c3    read <kl||ij>aaaa to W01
       call getmediate (wrk,wrksize,
     & lunsta,possw010,mapdw01,mapiw01,rc)
c
c4    read <kl||ij>bbbb to W02
       call getmediate (wrk,wrksize,
     & lunsta,possw020,mapdw02,mapiw02,rc)
c
c5    read <kl||ij>abab to W03
       call getmediate (wrk,wrksize,
     & lunsta,possw030,mapdw03,mapiw03,rc)
c
c
c6    read <ie||mn>aaaa to W11
       call getmediate (wrk,wrksize,
     & lunsta,possw110,mapdw11,mapiw11,rc)
c
c7    read <ie||mn>bbbb to W12
       call getmediate (wrk,wrksize,
     & lunsta,possw120,mapdw12,mapiw12,rc)
c
c8    read <ie||mn>abab to W13
       call getmediate (wrk,wrksize,
     & lunsta,possw130,mapdw13,mapiw13,rc)
c
c9    read <ie||mn>baab to W14
       call getmediate (wrk,wrksize,
     & lunsta,possw140,mapdw14,mapiw14,rc)
c
c
c10   read <ab||ij>aaaa to V1
       call getmediate (wrk,wrksize,
     & lunsta,possv10,mapdv1,mapiv1,rc)
c
c11   read <ab||ij>bbbb to V2
       call getmediate (wrk,wrksize,
     & lunsta,possv20,mapdv2,mapiv2,rc)
c
c12   read <ab||ij>abab to V3
       call getmediate (wrk,wrksize,
     & lunsta,possv30,mapdv3,mapiv3,rc)
c
c*    close INTSTA file
c
       if (iokey.eq.1) then
c      Fortran IO
       close (lunsta)
c
       else
c      MOLCAS IO
       call daclos (lunsta)
       end if
c
       return
       end
c
c      -------------------------------------------------------------
c
