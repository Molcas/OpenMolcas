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
      SUBROUTINE rotorb (cmoo, cmon, c, x, x2, y, thmax,FA)
!
!     RASSCF program: version IBM-3090: SX section
!
!     PURPOSE: Calculation of the rotation matrix X from the super-CI
!              coefficient matrix C, formation of exp(x), and rotation
!              of the orbitals with this matrix.
!              These orbitals are used in the next
!              RASSCF iteration as starting orbitals.
!              Called from SXCTL
!
!          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
!
      use stdalloc, only: mma_allocate, mma_deallocate
      use gas_data, only: iDoGAS, NGAS, NGSSH
      use rasscf_global, only: PURIFY, CMAX, ROTMAX, iXSym
      use PrintLevel, only: DEBUG,VERBOSE,TERSE
      use output_ras, only: LF,IPRLOC
      use general_data, only: NSYM,NASH,NBAS,NDEL,NFRO,NISH,NORB,NRS1,  &
     &                        NRS2,NSSH,NTOT2


      IMPLICIT None

      Character(LEN=16), Parameter :: ROUTINE='ROTORB  '

      Real*8 cmoo(*), cmon(*), c(*), x(*), x2(*)
      Real*8 y(*), FA(*), THMAX

      Real*8 DAMPGAS(10),COREGAS(10)
      INTEGER IDAMPGAS(0:4,0:4),ICOREGAS(0:4,0:4)
      LOGICAL iFrzAct
      REAL*8, PARAMETER  :: Thrs=1.0D-14
      Real*8, Allocatable:: Unit(:), SqFA(:)
      INTEGER I, IB, ICORE, IDAMP, IGAS, II, IJ, IO, iOff, iOrb,        &
     &        iPrLev, iSpace, IST, ISTBM, ISTMO, ISTMO1, ISUM, ISYM,    &
     &        jPr, jSPace, MOType, NACI, NACJ, NAE, NAO, NB, NBO, ND,   &
     &        NDB, NEO, NF, NFB, NI, NII, NIO, NIO1, NJ, NO, NOC, NOC1, &
     &        NP, NR
      REAL*8 TERM, THM, Xn, XX
      REAL*8, Parameter :: ACC=1.0D-13
!
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
         WRITE(LF,*)' Entering ',ROUTINE


        Write(LF,*)
        Write(LF,*)'FI+FA in RotOrb by Unitary transform'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do

      END IF

      istbm=1
      istmo1=0
      ib=0
      rotmax=0.0D0
      cmax=0.0D0
      thmax=0.0d0
      iOff = 1
! A long loop over symmetry
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
!
        IF(IPRLEV.ge.VERBOSE)                                           &
     & write(6,*) 'No rotations active for symmetry=',isym
!      No rotations active in this symmetry
!      Move orbitals from old to new set
!
          call dcopy_(nb*nb, cmoo(istmo1+1), 1, cmon(istmo1+1), 1)
          GO TO 100
        END IF
        istmo=istmo1+nf*nb
!
!      Form the quadratic matrix X from the coefficients C
!
        call dcopy_(no*no,[0.0D0],0,x,1)
        DO nr=1,noc
         DO np=max(nr+1,nio+1),no
          jpr=no*(nr-1)+np
          !jrp=no*(np-1)+nr
          xx=c(istbm+nr+noc*(np-nio-1))
!
!         Any numerical information smaller than Acc is considered
!         numerical noise and is ignored.
!
          If (Abs(xx).lt.Acc) Cycle
!
          x(jpr)=xx
          !x(jrp)=-xx
         END DO
        END DO
!
!      Set zero to matrix elements corresponding to orbitals
!      not allowed to rotate into each other.
!
        ij=0
        DO ni=1,no
         DO nj=1,no
          ij=ij+1
          IF (ixsym(ni+io).ne.ixsym(nj+io)) x(ij)=0.0D0
         END DO
        END DO
!
!  Freeze Active orbitals for difficult orbital optimization
!  This is activated by FRAC keyword in RASSCF module
!
        iFrzAct = .false.
        if(iFrzAct)then
          ij=0
          DO ni=1,no
           DO nj=1,no
            ij=ij+1
! io=ib+nf ! offset counting all nbas of previous sym and nfro of current sym.
            IF( (ni.gt.nio) .and. (ni.lt.nio+nao) .or.                  &
     &          (nj.gt.nio) .and. (nj.lt.nio+nao) ) x(ij)=0.0D0
           END DO
          END DO
        end if
!
!      For optimization of the core-hole we want to damp/eliminate
!      certain rotations to find the local minimum.
!
!      Damp/eliminate certain rotations inside RAS/GAS.
!
!      Only loop over the active spaces
!
!     DIMENSION DAMPGAS(10),COREGAS(10)
!     INTEGER IDAMPGAS(0:4,0:4),ICOREGAS(0:4,0:4)
      IDAMP = 0
      IDAMPGAS = 0
      DAMPGAS = 0.0D0
      ICORE = 0
      ICOREGAS = 0
      COREGAS = 0.1D0
!     COREGAS(1) = 0.0D0
!     COREGAS(2) = 0.0D0
      ICOREGAS(0,1) = 1
      ICOREGAS(1,0) = 2
!     ICOREGAS(0,2) = 1
!     ICOREGAS(2,0) = 2
      ICOREGAS(2,1) = 3
      ICOREGAS(1,2) = 4
      ICOREGAS(3,1) = 5
      ICOREGAS(1,3) = 6
      ICOREGAS(1,4) = 1
      ICOREGAS(4,1) = 2
!      write(6,*) 'NOC,NO',NOC,NO
       IF(ICORE.EQ.1.OR.IDAMP.EQ.1) THEN ! New keywords
!
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
!
           DO NACJ = 1, NO ! Over all orbitals due to counting
!            write(6,*) 'NACI,NACJ',NACI,NACJ
             IF(IDOGAS) THEN ! Find which GAS this index belong to
               JSPACE = 0
               IF(NACJ.LE.NIO) THEN
                 JSPACE = 0
               ELSE IF(NACJ.GT.NOC) THEN
                 JSPACE = NGAS + 1
!                IJ = IJ + 1
!                CYCLE
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
!                IJ = IJ + 1
!                write(6,*) 'IJ',IJ
!                CYCLE
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
                 write(6,*) ' damping'
             write(6,*)'ICOREGAS(ISPACE,JSPACE)',ICOREGAS(ISPACE,JSPACE)
             write(6,*) 'fac',COREGAS(ICOREGAS(ISPACE,JSPACE))
                 X(IJ) = X(IJ)*COREGAS(ICOREGAS(ISPACE,JSPACE))
               END IF
             END IF
           END DO
!
         END DO
! temp print
        ij=0
        DO ni=1,no
         DO nj=1,no
          ij=ij+1
          write(6,*) 'ij,x(ij)',ij,x(ij)
         END DO
        END DO
!
       END IF
!
!      Now form the unitary matrix exp(X)
!
        CALL exp_eig(no,x,thm)
        thmax=max(thmax,thm)
!
!      Check for largest non diagonal element
!
        ij=0
        DO ni=1,no
          DO nj=1,no
            ij=ij+1
            IF (abs(x(ij)).lt.Thrs) THEN
              x(ij) = 0.0D0
              GO TO 40
            END IF
            IF (ni.eq.nj) GO TO 40
            IF (abs(x(ij)).gt.abs(rotmax)) rotmax=x(ij)
  40        CONTINUE
          END DO
        END DO
!
!      Check for large rotations and phase of new orbital
!
 1010   FORMAT (6x,'Molecular orbital',i4,                              &
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
!
! Transformation of the Fock matrix according to the orbital rotation matrix.
! This step is not required other than for FCIQMC as orbital energies could be
! used for choosing the reference determinant.

        IF (iprlev.ge.debug) THEN
          CALL recprt ('X in RotOrb', ' ', x, no, no)
          write(6,*) 'FA for sym = ', iSym
          write(6,*) 'iOff is set to = ', iOff
          Call TriPrt(' ',' ',FA(iOff),no)
        END IF


        Call mma_allocate(UNIT,no*no,Label='UNIT')
        Call mma_allocate(SqFA,no*no,Label='SqFA')
        Unit(:)=0.0D0
        SqFA(:)=0.0D0
        Call Square(FA(iOff),SqFA,1,no,no)
        IF (iprlev.ge.debug) THEN
          CALL recprt ('Square FA in RotOrb', ' ', SqFA, no, no)
        END IF
        CALL DGEMM_('N','N',                                            &
     &              no,no,no,                                           &
     &              1.0d0,x,no,                                         &
     &                    SqFA,no,                                      &
     &              0.0d0,UNIT,no)

        CALL DGEMM_('N','T',                                            &
     &              no,no,no,                                           &
     &              1.0d0,UNIT,no,                                      &
     &                    x,no,                                         &
     &              0.0d0,SqFA,no)

        call Fold_Mat(1,[no],SqFA,FA(iOff))

        iOff = iOff + (no*no+no)/2

        Call mma_deallocate(Unit)
        Call mma_deallocate(SqFA)
!
!     Print output in MO-basis
!
!      Transform new orbitals to AO-basis
!
        CALL DGEMM_('N','N',                                            &
     &              nb,no,no,                                           &
     &              1.0d0,cmoo(istmo+1),nb,                             &
     &              x,no,                                               &
     &              0.0d0,cmon(istmo+1),nb)
!
!      Calculate max change in occupied molecular orbital coefficient
!
        nbo=noc*nb
        DO i=1,nbo
          term=cmon(istmo+i)-cmoo(istmo+i)
          IF (abs(term).gt.abs(cmax)) cmax=term
        END DO
!
 100    CONTINUE
!
!      Move frozen orbitals from old to new set of orbitals
!
        IF (nf.ne.0) THEN
          nfb=nf*nb
          call dcopy_(nfb, cmoo(istmo1+1), 1, cmon(istmo1+1), 1)
        END IF
!
!      MOVE DELETED ORBITALS.
!
        nd=ndel(isym)
        IF (nd.ne.0) THEN
          ist=istmo1+nb*(nf+no)+1
          ndb=nd*nb
          call dcopy_(ndb, cmoo(ist), 1, cmon(ist), 1)
        END IF
!
        istbm=istbm+noc*nae
        istmo1=istmo1+nb*nb
        ib=ib+nb
      END DO
!     ^ End of the long loop over symmetry

      IF(PURIFY(1:4).eq.'ATOM') CALL SPHPUR(CMON)
      IF(PURIFY(1:6).eq.'LINEAR') CALL LINPUR(CMON)
!
!     Orthogonalize new MO's and move them back to CMOO
!
      CALL supsch (x2, cmoo, cmon)
! PAM07: In ortho, cmoo is scratch and will be destroyed.
      CALL ortho_rasscf (x2, cmoo, cmon, y)
      call dcopy_(ntot2, cmon, 1, cmoo, 1)
!
      IF (iprlev.ge.debug) THEN
        Write(LF,*)
        Write(LF,*) ' >>> Exit RotOrb <<< '
        Write(LF,*)
      END IF
!
      END SUBROUTINE rotorb
