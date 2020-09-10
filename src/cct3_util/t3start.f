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
c     t3initfiles
c     t3reainput
c     t3wrhead
c     divfok
c     cct3_fokunpck1
c     cct3_fokunpck2
c     cct3_fokunpck3
c     cct3_fokunpck4
c     cct3_fokunpck5
c     t3reaintsta
c     t3reaccsd
c     t3wresult
c     t3gresult
c
c
c     -----------------------
c
       subroutine t3initfiles (lenght)
c
c     this routine distribute work space WRK for required files
c     for fix mediates it defines also mapd and mapi, for help mediates
c     it estimates their lenght and distribute WRK (i.e. def poss0 parameters)
c
c     lenght - overal requirements of work space (O)
c
c     !N.B. This routine cannot run with +OP2 level
c
       integer lenght
c
#include "t31.fh"
#include "t32.fh"
c
c     help variable
c
       integer posst,symp,symq,symr
       integer sizew,sizem,sizeh,sizen,sizel,sizer
       integer maxnoa,maxnvb,maxnorb
       integer nhelp1,nhelp2
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
c
       posst110=posst
       call cct3_grc0 (2,0,3,1,0,0,1,
     & posst110,posst,mapdt11,mapit11)
       posst120=posst
       call cct3_grc0 (2,0,4,2,0,0,1,
     & posst120,posst,mapdt12,mapit12)
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
       call cct3_grc0 (2,0,3,3,0,0,1,
     & possfk10,posst,mapdfk1,mapifk1)
       possfk20=posst
       call cct3_grc0 (2,0,4,4,0,0,1,
     & possfk20,posst,mapdfk2,mapifk2)
       possfk30=posst
       call cct3_grc0 (2,0,3,1,0,0,1,
     & possfk30,posst,mapdfk3,mapifk3)
       possfk40=posst
       call cct3_grc0 (2,0,4,2,0,0,1,
     & possfk40,posst,mapdfk4,mapifk4)
       possfk50=posst
       call cct3_grc0 (2,0,1,1,0,0,1,
     & possfk50,posst,mapdfk5,mapifk5)
       possfk60=posst
       call cct3_grc0 (2,0,2,2,0,0,1,
     & possfk60,posst,mapdfk6,mapifk6)
c
c
c1.6  maps for T2
c     T21 - t2o(ab,ij)aaaa
c     T22 - t2o(ab,ij)bbbb
c     T23 - t2o(a,b,i,j)abab
c
       posst210=posst
       call cct3_grc0 (4,4,3,3,1,1,1,
     & posst210,posst,mapdt21,mapit21)
       posst220=posst
       call cct3_grc0 (4,4,4,4,2,2,1,
     & posst220,posst,mapdt22,mapit22)
       posst230=posst
       call cct3_grc0 (4,0,3,4,1,2,1,
     & posst230,posst,mapdt23,mapit23)
c
c
c1.8  maps for W1
c     W11 - <ie||mn>aaaa
c     W12 - <ie||mn>bbbb
c     W13 - <ie||mn>abab
c     W14 - <ie||mn>baab
c
       possw110=posst
       call cct3_grc0 (4,3,1,3,1,1,1,
     & possw110,posst,mapdw11,mapiw11)
       possw120=posst
       call cct3_grc0 (4,3,2,4,2,2,1,
     & possw120,posst,mapdw12,mapiw12)
       possw130=posst
       call cct3_grc0 (4,0,1,4,1,2,1,
     & possw130,posst,mapdw13,mapiw13)
       possw140=posst
       call cct3_grc0 (4,0,2,3,1,2,1,
     & possw140,posst,mapdw14,mapiw14)
c
c
c1.9  maps for W2
c     W21 - <ab||ij>aaaa
c     W22 - <ab||ij>bbbb
c     W23 - <a,b|i,j>abab
c
       possw210=posst
       call cct3_grc0 (4,4,3,3,1,1,1,
     & possw210,posst,mapdw21,mapiw21)
       possw220=posst
       call cct3_grc0 (4,4,4,4,2,2,1,
     & possw220,posst,mapdw22,mapiw22)
       possw230=posst
       call cct3_grc0 (4,0,3,4,1,2,1,
     & possw230,posst,mapdw23,mapiw23)
c
c
c2    for help files mapps are irrelevant,
c     here only estimation of maximal lenght is done to
c     define poss0 of help files
c     we have:
c     2  W,V files - of vv2 type
c     2    L files - of vvv (vvo) type
c     3    R files - of vv2+ type
c     3    M files - of vv (vo)  type
c     3    H files - of v (o)  type
c     2  N,P files - of nn   type
c
c
       possw0=posst
c
c2.*  find maxsize of W,L,M,H
       sizew=0
       sizel=0
       sizer=0
       sizem=0
       sizeh=0
c
       do 50 nhelp1=1,nsym
c
c     W,V files
       call cct3_t3grc0 (3,2,4,4,4,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizew) then
       sizew=nhelp2
       end if
c
c     L files
       call cct3_t3grc0 (3,0,4,4,4,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizel) then
       sizel=nhelp2
       end if
       call cct3_t3grc0 (3,0,1,4,4,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizel) then
       sizel=nhelp2
       end if
c
c     R files
       call cct3_t3grc0 (3,8,4,4,4,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizer) then
       sizer=nhelp2
       end if
c
c     M files
       call cct3_t3grc0 (2,0,4,4,0,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizem) then
       sizem=nhelp2
       end if
       call cct3_t3grc0 (2,0,1,4,0,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizem) then
       sizem=nhelp2
       end if
c
c     H files
       call cct3_t3grc0 (1,0,4,0,0,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizeh) then
       sizeh=nhelp2
       end if
       call cct3_t3grc0 (1,0,1,0,0,0,nhelp1,possw0,posst,mapdw,mapiw)
       nhelp2=posst-possw0
       if (nhelp2.gt.sizeh) then
       sizeh=nhelp2
       end if
c
 50     continue
c
c
c2.*  def max{noa}, max{norb} ,max{nvb}
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
 100    continue
c
c
c2.*  def lenghts of N fils
c
       sizen=0
c
       do 200 symp=1,nsym
c     symq is not known for N file
c     instead of norb(symr) maxnorb will be used so that reallenght<=lenght
       sizen=sizen+norb(symp)*maxnorb
 200    continue
c
c2.1  W,V - files
c
c     possw0 is defined
       posst=possw0+sizew
       possv0=posst
       posst=posst+sizew
c
c2.2  L - files
c
       possl10=posst
       posst=posst+sizel
       possl20=posst
       posst=posst+sizel
c
c2.3  R - files
c
       possr10=posst
       posst=posst+sizer
       possr20=posst
       posst=posst+sizer
       possr30=posst
       posst=posst+sizer
c
c2.4  M - files
c
       possm10=posst
       posst=posst+sizem
       possm20=posst
       posst=posst+sizem
       possm30=posst
       posst=posst+sizem
c
c2.5  H - files
c
       possh10=posst
       posst=posst+sizeh
       possh20=posst
       posst=posst+sizeh
       possh30=posst
       posst=posst+sizeh
c
c2.6  N,P - files
c
       possn0=posst
       posst=posst+sizen
       possp0=posst
       posst=posst+sizen
c
c2.7  dedlare space for help matrix D in for matrix multiplication C=AT*B if
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
       lenght=posst-1
c
       return
       end
c
c     ----------------------
c
       subroutine t3reainput
c
c     this routine do:
c     1) read INPDAT file, produced by REORG with mul,nsym,noa,nob,nva,nvb,norb,eps
c     2) read input file for NIT3 to read (parameters transported through cmm common)

c   ####################
c   Due to the merging of CC input files to one, to avoid conflicts
c   Denominators in CCT3 has become T3Denominators
c   and Shift has become T3Shift
c                           (JR) Lund 2003
c     ! title   - jobtitle
c     1-ntit rows with 72 characters
c     no default
c     ! ntit    - number of rowws in jobtitle
c     ntit is limited to 10 CGG From now 1
c     default (1)
c     ! typt3   - type of T3 cpntribution
c     0 - CCSD
c     1 - CCSD+T(CCSD)
c     2 - CCSD(T) Ragh
c     3 - CCSD(T) Bart
c     default (3)
c     ! typden  - type of denominator (division of fok)
c     0 - diagonal elements
c     1 - average of faa and fbb
c     2 - orbital energies
c     default (0)
c     ! keysa   - Spin adaptation key
c     0 - no adaptation
c     1 - T2 DDVV adaptation
c     2 - T2 DDVV + T1 DV adaptation
c     3 - full T1 and T2 adaptation (only for doublets)
c     4 - full T2 adaptation without SDVS (only for doublets)
c     (default=0)
c     ! filerst - name for CCSD results containing file
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
c     ! fullprint - level of printing control key
c     (default=0)
c     ! noop - No Operation key
c     (default=no)
c     & iokey - I/O control key
c       1 - Fortran I/O system
c       2 - MOLCAS DA IO system
c     (default=2)
c     & mhkey - Matrix handling control key
c       1 - ESSL routines
c       2 - Fortran I/O system
c     (default=1)
c     .....   - can be added
c
c     3) initialize nshf
c
c
#include "t31.fh"
c
c     help variables
c
       character*80 LINE
       integer nhelp
c
c1    read INPDAT
c
c       open (unit=1,file='INPDAT',form='unformatted')
       call molcas_binaryopen_vanilla(1,'INPDAT')
       read (1) nactel,ispin,nsym,lsym,mmul,noa,nob,nva,nvb,norb,eps
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
c
c4    define defaults
c
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
CGG       fullprint=0
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
c
c5    read input file
c
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
c
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
       Read(LuSpool,*) symimin,imin,symjmin,jmin,
     c                 symimax,imax,symjmax,jmax
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
       subroutine t3wrhead
c
c     this routine write head of ther output file
c     parameters transported through cmm common
c
#include "t31.fh"
c
c     help variables
c
       integer nhelp,nhelp1,nhelp2
c
c1    write title
c
      if (noop.eq.1) then
      write(6,'(6X,A)') ' No operation is required'
      write(6,'(6X,A)') ' Happy Landing'
      Call Finish(0)
      end if
c
CGG       do 10 nhelp=1,ntit
CGG       write(6,9) title(nhelp)
CGG 10     continue
CGG 9      format (A72)
CGG       write(6,*)
c
c2    write occupations and matrix multiplication tablaux
c
       write(6,11) (norb(nhelp),nhelp=1,nsym)
 11     format (' NORB ',8(i3,2x))
       write(6,12) (noa(nhelp),nhelp=1,nsym)
 12     format (' NOA  ',8(i3,2x))
       write(6,13) (nob(nhelp),nhelp=1,nsym)
 13     format (' NOB  ',8(i3,2x))
       write(6,14) (nva(nhelp),nhelp=1,nsym)
 14     format (' NVA  ',8(i3,2x))
       write(6,15) (nvb(nhelp),nhelp=1,nsym)
 15     format (' NVB  ',8(i3,2x))
c
       if (fullprint.gt.1) then
       write(6,*)
       write(6,*)
       write(6,17)
 17     format (' MATRIX MULTIPLICATION TABLE')
       write(6,*)
       write(6,18) (nhelp,nhelp=1,nsym)
 18     format (' IRREP #',9x,8(i3,2x))
       write(6,*)
       do 19 nhelp1=1,nsym
       write(6,20) nhelp1,(mmul(nhelp1,nhelp),nhelp=1,nsym)
 19     continue
 20     format (' IRREP #',i3,6x,8(i3,2x))
       write(6,*)
       end if
c
       write(6,*)
       write(6,16) nsym
 16     format (' NUMBER OF IRREPS             :',i3)
       write(6,161) ispin
 161    format (' MULTIPLICITY                 :',i3)
       write(6,162) lsym
 162    format (' OVERALL SYMMETRY STATE       :',i3)
c
c4    write type of triples
c
       if (typt3.eq.0) then
       write(6,31)
 31     format (' METHOD                       : CCSD')
       else if (typt3.eq.1) then
       write(6,32)
 32     format (' METHOD                       : CCSD+T(CCSD)')
       else if (typt3.eq.2) then
       write(6,33)
 33     format (' METHOD                       : CCSD+T(CCSD)',
     &  '+<T3(1)WT1> = CCSD(T)' )
       else if (typt3.eq.3) then
       write(6,34)
 34     format (' METHOD                       : CCSD+T(CCSD)',
     &  '+<T3(1)WT1>+<T3(1)UT2> = CCSD(T)' )
       end if
c
c5    write type of fok division
c
       if (typden.eq.0) then
       write(6,41)
       else if (typden.eq.1) then
       write(6,42)
       else
       write(6,43)
       end if
 41     format (' TYPE OF DENOMINATOR          : DIAGONAL')
 43     format (' TYPE OF DENOMINATOR          : ORBITAL ENERGIES')
 42     format (' TYPE OF DENOMINATOR          : (FAA+FBB)/2')
c
c
c6    write orbital energies per symmetry
c
       if (fullprint.gt.0) then
       write(6,*)
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
c9    type file whre CCSD results are
       write(6,91) filerst
 91     format (' CCSD RESULTS LOAD FROM FILE : ',a6)
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
c121    format (' MAXIMAL ALLOWED WORK SPACE  : Unlimited')
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
       else
       write(6,133)
 133    format (' LEVEL OF OUTPUT PRINTING    : MAXIMAL')
       end if
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
       if (mhkey.eq.1) then
       write(6,151)
 151    format (' MATRIX HANDLING             : ESSL        ')
       else
       write(6,152)
 152    format (' MATRIX HANDLING             : Fortran code')
       end if
c
c20   Print parameters of IJ cycle segmentation
c
       if (ijsegkey.eq.0) then
       write(6,201)
 201    format (' IJ CYCLE SEGMENTED          : NO          ')
       else
       write(6,202)
       write(6,203) symimin
       write(6,204) imin
       write(6,205) symjmin
       write(6,206) jmin
       write(6,207) symimax
       write(6,208) imax
       write(6,209) symjmax
       write(6,210) jmax
       write(6,211)
 202    format (' IJ CYCLE SEGMENTED          : YES         ')
 203    format (' SYMI minimal                : ',i4)
 204    format ('    I minimal                : ',i4)
 205    format (' SYMJ minimal                : ',i4)
 206    format ('    J minimal                : ',i4)
 207    format (' SYMI maximal                : ',i4)
 208    format ('    I maximal                : ',i4)
 209    format (' SYMJ maximal                : ',i4)
 210    format ('    J maximal                : ',i4)
 211    format (' Be very careful in using IJSEgmentation technique')
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
       subroutine cct3_divfok (wrk,wrksize,
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
#include "t31.fh"
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
       call cct3_fokunpck5 (symp,wrk(possfoka),wrk(possfokb),
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
       call cct3_fokunpck1 (wrk(possfok),wrk(possdp),norb(symp))
       if (nva(symp).gt.0) then
       call cct3_fokunpck2 (wrk(possfok),wrk(possfaa),norb(symp),
     & nva(symp),noa(symp))
       end if
       if ((noa(symp)*nva(symp)).gt.0) then
       call cct3_fokunpck3 (wrk(possfok),wrk(possfai),norb(symp),
     & nva(symp),noa(symp))
       end if
       if (noa(symp).gt.0) then
       call cct3_fokunpck4 (wrk(possfok),wrk(possfii),norb(symp),
     &  noa(symp))
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
       call cct3_fokunpck1 (wrk(possfok),wrk(possdp),norb(symp))
       if (nvb(symp).gt.0) then
       call cct3_fokunpck2 (wrk(possfok),wrk(possfaa),norb(symp),
     & nvb(symp),nob(symp))
       end if
       if ((nob(symp)*nvb(symp)).gt.0) then
       call cct3_fokunpck3 (wrk(possfok),wrk(possfai),norb(symp),
     & nvb(symp),nob(symp))
       end if
       if (nob(symp).gt.0) then
       call cct3_fokunpck4 (wrk(possfok),wrk(possfii),norb(symp),
     &  nob(symp))
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
       subroutine cct3_fokunpck1 (fok,dp,dimfok)
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
       subroutine cct3_fokunpck2 (fok,faa,dimfok,dimfa,shift)
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
       do 201 a=1,dimfa
       faa(a,b)=fok(shift+a,shift+b)
 201    continue
 200    continue
c
       return
       end
c
c     -------------
c
       subroutine cct3_fokunpck3 (fok,fai,dimfok,dimfa,dimfi)
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
       do 301 a=1,dimfa
       fai(a,i)=fok(dimfi+a,i)
 301    continue
 300    continue
c
       return
       end
c
c     -------------
c
       subroutine cct3_fokunpck4 (fok,fii,dimfok,dimfi)
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
       do 401 i=1,dimfi
       fii(i,j)=fok(i,j)
 401    continue
 400    continue
c
       return
       end
c
c     -------------
c
       subroutine cct3_fokunpck5 (symp,foka,fokb,dpa,dpb,dimfok,rc)
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
#include "t31.fh"
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
       return
       end
c
c     -------------------------------------
c
       subroutine t3reaintsta (wrk,wrksize)
c
c     this routine read integral file INTSTA (reorg), which contains
c     following integrals: foka,fokb,
c     <kl||ij>aaaa,<kl||ij>bbbb,<kl||ij>abab - naplano
c     <ka||ij>aaaa,<ka||ij>bbbb,<ka||ij>abab,<ka||ij>baab
c     <ab||ij>aaaa,<ab||ij>bbbb,<ab||ij>abab
c
c     two electron integrals are readed to their fix files,
c     foka,fokb are readed to N,P help files
c
c     use and destroy : N,P
c
#include "t31.fh"
#include "t32.fh"
#include "wrk.fh"
c
c     help variables
c
       integer lunsta,rc
c
c*    open INTSTA file
       lunsta=1
       if (iokey.eq.1) then
c      Fortran IO
c       open (unit=lunsta,file='INTSTA',form='unformatted')
       call molcas_binaryopen_vanilla(lunsta,'INTSTA')
c
       else
c      MOLCAS IO
       call daname (lunsta,'INTSTA')
       daddr(lunsta)=0
       end if
c
c1    read foka to N
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possn0,mapdn,mapin,rc)
c
c2    read fokb to P
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possp0,mapdp,mapip,rc)
c
c
c3    read <kl||ij>aaaa to W23 - naplano
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw230,mapdw23,mapiw23,rc)
c
c4    read <kl||ij>bbbb to W23 - naplano
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw230,mapdw23,mapiw23,rc)
c
c5    read <kl||ij>abab to W23 - naplano
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw230,mapdw23,mapiw23,rc)
c
c
c6    read <ie||mn>aaaa to W11
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw110,mapdw11,mapiw11,rc)
c
c7    read <ie||mn>bbbb to W12
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw120,mapdw12,mapiw12,rc)
c
c8    read <ie||mn>abab to W13
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw130,mapdw13,mapiw13,rc)
c
c9    read <ie||mn>baab to W14
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw140,mapdw14,mapiw14,rc)
c
c
c10   read <ab||ij>aaaa to W21
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw210,mapdw21,mapiw21,rc)
c
c11   read <ab||ij>bbbb to W22
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw220,mapdw22,mapiw22,rc)
c
c12   read <ab||ij>abab to W23
       call cct3_getmediate (wrk,wrksize,
     & lunsta,possw230,mapdw23,mapiw23,rc)
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
c     -----------------------
c
       subroutine t3reaccsd (wrk,wrksize,
     & eccsd)
c
c     this routine read CCSD results, it T1 and T2 amplitudes
c     and CCSD energy from saverst file
c
c     eccsd - Converged CCSD energy (O)
c
#include "t31.fh"
#include "t32.fh"
#include "wrk.fh"
c
       real*8 eccsd,dum(1)
c
c     help parameters
c
       integer lunrst,rc1,iteration
c
c
c1    open file savename
       lunrst=1
c
       if (iokey.eq.1) then
c      Fortran IO
c       open (unit=lunrst,file=filerst,form='unformatted')
        call molcas_binaryopen_vanilla(lunrst,filerst)
c
       else
c      MOLCAS IO
       call daname (lunrst,filerst)
       daddr(lunrst)=0
       end if
c
c2    get T1aa
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst110,mapdt11,mapit11,rc1)
c
c3    get T1bb
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst120,mapdt12,mapit12,rc1)
c
c4    get T2aaaa
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst210,mapdt21,mapit21,rc1)
c
c5    get T2bbbb
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst220,mapdt22,mapit22,rc1)
c
c6    get T2abab
       call cct3_getmediate (wrk,wrksize,
     & lunrst,posst230,mapdt23,mapit23,rc1)
c
c7    get energy,niter
       if (iokey.eq.1) then
c      Fortran IO
       read (lunrst,end=1) eccsd,iteration
       else
c
c      MOLCAS IO
       call ddafile (lunrst,2,dum,1,daddr(lunrst))
       eccsd=dum(1)
       end if
c
       goto 999
c
 1      write(6,*) ' ENERGY AND NIT WAS NOT IN SAVE FILE, CHANGED TO 0'
       write(6,*) ' USE CCSD ENERGY FROM CCSD OUTPUT FILE'
       eccsd=0.0d0
       iteration=0
c
 999   if (iokey.eq.1) then
c      Fortran IO
       close (lunrst)
c
       else
c      MOLCAS IO
       call daclos (lunrst)
       end if
c
       return
       end
c
c     -------------------------
c
        subroutine t3wresult (symi,symj,i,j,eaaa,eaab,eabb,ebbb)
c
c       this routine write:
c       0) value od SymiMin,Imin,SymJmin,Jmin from which
c          accumilation started
c       1) value of SymI,SymJ,I,J
c       2) present stage of energies
c       into T3tEne file and overwrite previous values
c
#include "t31.fh"
c
        integer symi,symj,i,j
        real*8 eaaa,eaab,eabb,ebbb
c
c       help variable
c
        integer lun
c
        lun=1
        Call Molcas_Open(lun,'T3tEne')
*       open (unit=lun,file='T3tEne')
c
        write (lun,97) symimin,imin,symjmin,jmin
        write (lun,98) symi,symj
        write (lun,98) i,j
        write (lun,99) eaaa
        write (lun,99) eaab
        write (lun,99) eabb
        write (lun,99) ebbb
c
97      format (2x,4(i4,2x))
98      format (2x,2(i4,2x))
99      format (2x,f22.16)
c
        close (lun)
c
        return
        end
c
c     -------------------------
c
        subroutine t3gresult (symi,symj,i,j,eaaa,eaab,eabb,ebbb)
c
c       this routine get:
c       1) value of SymI,SymJ,I,J
c       2) present stage of energies
c       from T3tEne file
c
        integer symi,symj,i,j
        real*8 eaaa,eaab,eabb,ebbb
c
c       help variable
c
        integer lun,bullshit
c
        lun=1
        Call Molcas_Open(Lun,'T3tEne')
*       open (unit=lun,file='T3tEne')
c
c       read blank, since there is Symimin,imin,Symjmin,jmin
c       on first line
c
        read (lun,*) bullshit
c
        read (lun,*) symi,symj
        read (lun,*) i,j
        read (lun,*) eaaa
        read (lun,*) eaab
        read (lun,*) eabb
        read (lun,*) ebbb
c
c98      format (2x,2(i4,2x))
c99      format (2x,f22.16)
c
        close (lun)
c
        return
        end
