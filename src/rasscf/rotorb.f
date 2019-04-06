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
      SUBROUTINE rotorb (cmoo, cmon, c, x, x2, y, vec, a, b, thmax,FA)
c
c     RASSCF program: version IBM-3090: SX section
c
c     PURPOSE: Calculation of the rotation matrix X from the super-CI
c              coefficient matrix C, formation of exp(x), and rotation
c              of the orbitals with this matrix.
c              These orbitals are used in the next
c              RASSCF iteration as starting orbitals.
c              Called from SXCTL
c
c          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
c

      IMPLICIT REAL*8 (a-h,o-z)

#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
#include "gas.fh"
      Parameter (ROUTINE='ROTORB  ')
#include "rasscf.fh"
#include "WrkSpc.fh"

      DIMENSION cmoo(*), cmon(*), c(*), x(*), vec(*),  x2(*)
      DIMENSION y(*), a(*), b(*), FA(*)
      DIMENSION DAMPGAS(10),COREGAS(10)
      INTEGER IDAMPGAS(0:4,0:4),ICOREGAS(0:4,0:4)
      LOGICAL iFrzAct
c
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
         WRITE(LF,*)' Entering ',ROUTINE


        Write(LF,*)
        Write(LF,*)'FI+FA in RotOrb bf Unitary transform'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do

      END IF

      Acc=1.0D-13
      istbm=1
      istmo1=0
      ib=0
      rotmax=0.0D0
      cmax=0.0D0
      thmax=0.0d0
      iOff = 1
* A long loop over symmetry
      DO isym=1,nsym
        no=norb(isym)
        nb=nbas(isym)
        nf=nfro(isym)
        io=ib+nf
        nio=nish(isym)
        nao=nash(isym)
        noc=nio+nao
        neo=nssh(isym)
        nae=nao+neo
        IF (noc.eq.0.or.nae.eq.0 .or. nio+neo.eq.0) THEN
c
        IF(IPRLEV.ge.VERBOSE)
     & write(6,*) 'No rotations active for symmetry=',isym
c      No rotations active in this symmetry
c      Move orbitals from old to new set
c
          call dcopy_(nb*nb, cmoo(istmo1+1), 1, cmon(istmo1+1), 1)
          GO TO 100
        END IF
        istmo=istmo1+nf*nb
c
c      Form the quadratic matrix X from the coefficients C
c
        call dcopy_(no*no,[0.0D0],0,x,1)
        DO nr=1,noc
         DO np=max(nr+1,nio+1),no
          jpr=no*(nr-1)+np
          jrp=no*(np-1)+nr
          xx=c(istbm+nr+noc*(np-nio-1))
*
*         Any numerical information smaller than Acc is considered
*         numerical noise and is ignored.
*
          If (Abs(xx).lt.Acc) xx=0.0D0
*
          x(jpr)=xx
          x(jrp)=-xx
         END DO
        END DO
c
c      Set zero to matrix elements corresponding to orbitals
c      not allowed to rotate into each other.
c
        ij=0
        DO ni=1,no
         DO nj=1,no
          ij=ij+1
          IF (ixsym(ni+io).ne.ixsym(nj+io)) x(ij)=0.0D0
         END DO
        END DO
c
c  Freeze Active orbitals for difficult orbital optimization
c  This is activated by FRAC keyword in RASSCF module
c
        iFrzAct = .false.
        if(iFrzAct)then
          ij=0
          DO ni=1,no
           DO nj=1,no
            ij=ij+1
c io=ib+nf ! offeset counting all nbas of previous sym and nfro of current sym.
            IF( (ni.gt.nio) .and. (ni.lt.nio+nao) .or.
     &          (nj.gt.nio) .and. (nj.lt.nio+nao) ) x(ij)=0.0D0
           END DO
          END DO
        end if
c
c      For optimization of the core-hole we want to damp/eliminate
c      certain rotations to find the local minimum.
c
c      Damp/eliminate certain rotations inside RAS/GAS.
c
c      Only loop over the active spaces
c
c     DIMENSION DAMPGAS(10),COREGAS(10)
c     INTEGER IDAMPGAS(0:4,0:4),ICOREGAS(0:4,0:4)
      IDAMP = 0
      IDAMPGAS = 0
      DAMPGAS = 0.0D0
      ICORE = 0
      ICOREGAS = 0
      COREGAS = 0.1D0
C     COREGAS(1) = 0.0D0
C     COREGAS(2) = 0.0D0
      ICOREGAS(0,1) = 1
      ICOREGAS(1,0) = 2
C     ICOREGAS(0,2) = 1
C     ICOREGAS(2,0) = 2
      ICOREGAS(2,1) = 3
      ICOREGAS(1,2) = 4
      ICOREGAS(3,1) = 5
      ICOREGAS(1,3) = 6
      ICOREGAS(1,4) = 1
      ICOREGAS(4,1) = 2
C      write(6,*) 'NOC,NO',NOC,NO
       IF(ICORE.EQ.1.OR.IDAMP.EQ.1) THEN ! New keywords
c
         IJ = 0
         DO NACI = 1, NO ! Notice only loop over inactive and active
           IF(IDOGAS) THEN ! Find which GAS this index belong to
             ISPACE = 0
             IF(NACI.LE.NIO) THEN
               ISPACE = 0
             ELSE IF(NACI.GT.NOC) THEN
               ISPACE = NGAS + 1
             ELSE
               ISUM = 0
               DO IGAS = 1, NGAS
                 IF((NGSSH(IGAS,ISYM)+NIO+ISUM).GE.NACI) THEN
                   ISPACE = IGAS
                   EXIT
                 ELSE
                   ISUM = ISUM + NGSSH(IGAS,ISYM)
                 END IF
               END DO
             END IF
           ELSE ! Find which RAS this index belong to
             IF(NACI.LE.NIO) THEN
               ISPACE = 0
             ELSE IF(NACI.GT.NOC) THEN
               ISPACE = 4
             ELSE IF((NRS1(ISYM)+NIO).GE.NACI) THEN
               ISPACE = 1
             ELSE IF((NRS1(ISYM)+NRS2(ISYM)+NIO).GE.NACI) THEN
               ISPACE = 2
             ELSE
               ISPACE = 3
             END IF
           END IF
c
           DO NACJ = 1, NO ! Over all orbitals due to counting
C            write(6,*) 'NACI,NACJ',NACI,NACJ
             IF(IDOGAS) THEN ! Find which GAS this index belong to
               JSPACE = 0
               IF(NACJ.LE.NIO) THEN
                 JSPACE = 0
               ELSE IF(NACJ.GT.NOC) THEN
                 JSPACE = NGAS + 1
C                IJ = IJ + 1
C                CYCLE
               ELSE
                 ISUM = 0
                 DO IGAS = 1, NGAS
                   IF((NGSSH(IGAS,ISYM)+NIO+ISUM).GE.NACJ) THEN
                     JSPACE = IGAS
                     EXIT
                   ELSE
                     ISUM = ISUM + NGSSH(IGAS,ISYM)
                   END IF
                 END DO
               END IF
             ELSE ! Find which RAS this index belong to
               JSPACE = 0
               IF(NACJ.LE.NIO) THEN
                 JSPACE = 0
               ELSE IF(NACJ.GT.NOC) THEN
                 JSPACE = 4
C                IJ = IJ + 1
C                write(6,*) 'IJ',IJ
C                CYCLE
               ELSE IF((NRS1(ISYM)+NIO).GE.NACJ) THEN
                 JSPACE = 1
               ELSE IF((NRS1(ISYM)+NRS2(ISYM)+NIO).GE.NACJ) THEN
                 JSPACE = 2
               ELSE
                 JSPACE = 3
               END IF
             END IF
             IJ = IJ + 1
             write(6,*) 'NACI,NACJ',NACI,NACJ
             write(6,*) 'ISPACE,JSPACE',ISPACE,JSPACE
             write(6,*) 'IJ',IJ
             IF(IDAMP.EQ.1) THEN
               IF(IDAMPGAS(ISPACE,JSPACE).GT.0) THEN
                 X(IJ) = X(IJ)*DAMPGAS(IDAMPGAS(ISPACE,JSPACE))
               END IF
             ELSE IF(ICORE.EQ.1) THEN
               IF(ICOREGAS(ISPACE,JSPACE).GT.0) THEN
                 write(6,*) ' dampning'
             write(6,*)'ICOREGAS(ISPACE,JSPACE)',ICOREGAS(ISPACE,JSPACE)
             write(6,*) 'fac',COREGAS(ICOREGAS(ISPACE,JSPACE))
                 X(IJ) = X(IJ)*COREGAS(ICOREGAS(ISPACE,JSPACE))
               END IF
             END IF
           END DO
c
         END DO
c temp print
        ij=0
        DO ni=1,no
         DO nj=1,no
          ij=ij+1
          write(6,*) 'ij,x(ij)',ij,x(ij)
         END DO
        END DO
c
       END IF
c
c      Now form the unitary matrix exp(X)
c
        CALL expx (x, x2, y, vec, a, b, thm, no)
        thmax=max(thmax,thm)
c
c      Check for largest non diagonal element
c
        ij=0
        DO ni=1,no
          DO nj=1,no
            ij=ij+1
            IF (ni.eq.nj) GO TO 40
            IF (abs(x(ij)).gt.abs(rotmax)) rotmax=x(ij)
  40        CONTINUE
          END DO
        END DO
c
c      Check for large rotations and phase of new orbital
c
 1010   FORMAT (6x,'Molecular orbital',i4,
     &             ' of symmetry',i2,' MO space',i2,'  weight is',f12.6)
        ii=1
        ist=0
        DO ni=1,no
          IF (ni.le.nio) THEN
            motype=1
            xn=0.0D0
            DO nii=1,nio
              xn=xn+x(ist+nii)**2
            END DO
            IF (xn.lt.0.5D0) Then
             IF(IPRLEV.GE.TERSE) THEN
              Call WarningMessage(1,'Large orbital rotation.')
              Write(LF,1010) ni,isym,motype,xn
             END IF
            End If
          END IF
          IF (ni.gt.nio.and.ni.le.noc) THEN
            motype=2
            xn=0.0D0
            nio1=nio+1
            DO nii=nio1,noc
              xn=xn+x(ist+nii)**2
            END DO
            IF (xn.lt.0.5D0) Then
             IF(IPRLEV.GE.TERSE) THEN
              Call WarningMessage(1,'Large orbital rotation.')
              Write(LF,1010) ni,isym,motype,xn
             END IF
            End If
          END IF
          IF (ni.gt.noc) THEN
            motype=3
            noc1=noc+1
            xn=0.0D0
            DO nii=noc1,no
              xn=xn+x(ist+nii)**2
            END DO
            IF (xn.lt.0.5D0) Then
             IF(IPRLEV.GE.TERSE) THEN
              Call WarningMessage(1,'Large orbital rotation.')
              Write(LF,1010) ni,isym,motype,xn
             END IF
            End If
          END IF
          ii=ii+no+1
          ist=ist+no
        END DO
c
C Transformation of the Fock matrix according to the orbital rotation matrix.
C This step is not required other than for FCIQMC as orbital energies could be
c used for choosing the reference determinant.

        IF (iprlev.ge.debug) THEN
          CALL recprt ('X in RotOrb', ' ', x, no, no)
          write(6,*) 'FA for sym = ', iSym
          write(6,*) 'iOff is set to = ', iOff
          Call TriPrt(' ',' ',FA(iOff),no)
        END IF


        Call GetMem('Unit','Allo','Real',LUNIT,no*no)
        Call GetMem('SqFA','Allo','Real',iSqFA,no*no)
        call fzero(Work(Lunit),no*no)
        call fzero(Work(iSqFA),no*no)
        Call Square(FA(iOff),Work(iSqFA),1,no,no)
        IF (iprlev.ge.debug) THEN
          CALL recprt ('Square FA in RotOrb', ' ', Work(iSqFA), no, no)
        END IF
        CALL DGEMM_('N','N',
     &              no,no,no,
     &              1.0d0,x,no,
     &              Work(iSqFA),no,
     &              0.0d0,Work(LUNIT),no)

        CALL DGEMM_('N','T',
     &              no,no,no,
     &              1.0d0,Work(LUNIT),no,
     &              x,no,
     &              0.0d0,Work(iSqFA),no)

        call Fold_Mat(1,[no],Work(iSqFA),FA(iOff))

        iOff = iOff + (no*no+no)/2

        Call GetMem('Unit','Free','Real',LUNIT,no*no)
        Call GetMem('SqFA','Free','Real',iSqFA,no*no)

c
c     Print output in MO-basis
c
c      Transform new orbitals to AO-basis
c
        CALL DGEMM_('N','N',
     &              nb,no,no,
     &              1.0d0,cmoo(istmo+1),nb,
     &              x,no,
     &              0.0d0,cmon(istmo+1),nb)
c
c      Calculate max change in occupied molecular orbital coefficient
c
        nbo=noc*nb
        DO i=1,nbo
          term=cmon(istmo+i)-cmoo(istmo+i)
          IF (abs(term).gt.abs(cmax)) cmax=term
        END DO
c
 100    CONTINUE
c
c      Move frozen orbitals from old to new set of orbitals
c
        IF (nf.ne.0) THEN
          nfb=nf*nb
          call dcopy_(nfb, cmoo(istmo1+1), 1, cmon(istmo1+1), 1)
        END IF
c
c      MOVE DELETED ORBITALS.
c
        nd=ndel(isym)
        IF (nd.ne.0) THEN
          ist=istmo1+nb*(nf+no)+1
          ndb=nd*nb
          call dcopy_(ndb, cmoo(ist), 1, cmon(ist), 1)
        END IF
c
        istbm=istbm+noc*nae
        istmo1=istmo1+nb*nb
        ib=ib+nb
      END DO
*     ^ End of the long loop over symmetry

      IF(PURIFY(1:4).eq.'ATOM') CALL SPHPUR(CMON)
      IF(PURIFY(1:6).eq.'LINEAR') CALL LINPUR(CMON)
c
c     Orthogonalize new MO's and move them back to CMOO
c
      CALL supsch (x2, cmoo, cmon)
* PAM07: In ortho, cmoo is scratch and will be destroyed.
      CALL ortho_rasscf (x2, cmoo, cmon, y)
      call dcopy_(ntot2, cmon, 1, cmoo, 1)
c
      IF (iprlev.ge.debug) THEN
        Write(LF,*)
        Write(LF,*) ' >>> Exit RotOrb <<< '
        Write(LF,*)
      END IF
c
      RETURN
      END
