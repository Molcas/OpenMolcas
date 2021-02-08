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
      Subroutine Vqr(LUQRP,MPLBL,ISIM,ZETA,NZ,OPMAT)
C
C...  Matrix of the atomic quasi-relativistic potential in the
C     1-center basis set. Mass-velocity and Darwin potentials,
*     which are read from the library in atomic units.
C
      Implicit Real*8 (A-H,O-Z)
      Character OrbLab*4, MPLbl*20, Keyw*40,
     &          RdName*40, DumCha*40, Line*80
      Real*8        ZETA(NZ),OPMAT(NZ*(NZ+1)/2)

*     Internal
      Parameter (MPoint=250)
      Real*8        R(Mpoint),V(Mpoint),FI(Mpoint),FIVFJ(Mpoint)

C...  auxiliar constant pool:    only DFAC is used here.
      PARAMETER (lp1=5,lp12=lp1*lp1,lp13=(lp1*lp1+lp1)/2)
      COMMON/CONST/RCA(lp1,lp13),DFAC(lp12),KOSUU(lp13),NYU(lp1,lp13)
      EXTERNAL      RDNAME,SIMPLM
      DATA          PI/3.141592653589793D0/
C...  Flags: MV bit 2, DW bit 3
*     Data iMVPot/2/,iDWPot/4/
C
*     LU6=6
      N=ISIM
C
C...  Reads the Mass Velocity and Darwin potentials for the
C     valence orbitals
C
*     WRITE (LU6,600) ISIM,LUQRP
C
C               locate the beginning of the data
      KEYW=MPLBL
        DUMCHA=RDNAME(LUQRP,KEYW)
        IF (DUMCHA.EQ.' ') THEN
          WRITE (6,'(1X,A20'//
     &       ','' MV & DW potentials not found in unit '',I3)')
     &       KEYW,LUQRP
          Call Quit_OnUserError()
        ENDIF
*     Write (*,*) 'Data found in libray. Start reading.'
      Read (LUQRP,*) Npoint
      If (NPOINT.GT.MPoint) Then
        Write (6,*) 'VQR: nPoint',nPoint
        Write (6,*) 'Abend: Increade mPoint'
        Call Quit_OnUserError()
      Endif
      Read (LUQRP,'(4d20.13)') (R(k),k=1,Npoint)
*
1     Read (LUQRP,'(a80)',End=999) Line
      IF (Line(1:8).eq.MPLbl(1:8)) THEN
C...                             there is no potential for this symmetry
*       WRITE (*,603)
        Do  IJ=1,NZ*(NZ+1)/2
          OPMAT(IJ)=0D0
        End Do
        RETURN
      ELSE
        Read (Line,'(a4)') OrbLab
        Read (LUQRP,'(4d20.13)') (V(k),k=1,Npoint)
C
C...    check the symmetry
        IF (    ((INDEX(ORBLAB,'S').NE.0).AND.(ISIM.NE.1))
     .      .OR.((INDEX(ORBLAB,'P').NE.0).AND.(ISIM.NE.2))
     .      .OR.((INDEX(ORBLAB,'D').NE.0).AND.(ISIM.NE.3))
     .      .OR.((INDEX(ORBLAB,'F').NE.0).AND.(ISIM.NE.4))
     .         ) THEN
          GO TO 1
        ENDIF
      ENDIF
*     WRITE (*,601) ORBLAB,NPOINT,R(1),R(NPOINT)
C
C...  Matrix of a numerical local operator V(r)
C     given in a logaritmic mesh
C
      PREN=2D0**(N+1)/(sqrt(DFAC(2*N  )*sqrt(2D0*PI)))
      N2P1=2*N+1
      IJ=0
C
      DO 20 I=1,NZ
        ZI=ZETA(I)
        RNORI=PREN*sqrt(sqrt(ZI**N2P1))
        DO 210 K=1,NPOINT
          ERRE=R(K)
          FI(K)=RNORI*ERRE**(N-1)*EXP(-ZI*ERRE*ERRE)
210     CONTINUE
C
        DO 21 J=1,I
          ZJ=ZETA(J)
          RNORJ=PREN*sqrt(sqrt(ZJ**N2P1))
          DO 220 K=1,NPOINT
            ERRE=R(K)
            ERRE2=ERRE*ERRE
            FIVFJ(K)=RNORJ*ERRE**(N-1)*EXP(-ZJ*ERRE2)*FI(K)*V(K)*ERRE2
220       CONTINUE
C
          RES=SIMPLM(NPOINT,FIVFJ,R)
C
C...      add up the contribution from 0 to R(1)
          IJ=IJ+1
          OPMAT(IJ)=RES+FIVFJ(1)*R(1)/3D0
21      CONTINUE
20    CONTINUE
C
      RETURN
999   CONTINUE
      Write (6,*)
     &   ' Troubles reading MV and DW potentials from QRP library'
      Call Quit_OnUserError()
C
*600   FORMAT (/1X,'Symmetry ',I1,' : '
*     .        ,'Mass Velocity and Darwin potentials are read'
*     .        ,'  from unit ',I2)
*601   FORMAT (/' Orbital ',A4,' : ',I5,' points in the logarithmic ',
*     &         'mesh from ',F10.6,' to ',F10.6)
*602   FORMAT (/1X,' The symmetry label of one of the relativistic'
*     .,' potentials is not S, P, D, or F')
*603   FORMAT (/1X,' There is no potential for this symmetry.')
C
      END
