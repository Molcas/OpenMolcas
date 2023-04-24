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
!
!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
      SUBROUTINE ALFas(NDP,YMIN,YH,NCN,V,SWF,VLIM,KVMAX,AFLAG,ZMU,EPS,  &
     &                                    GV,BFCT,INNODE,INNR,IWR)
!     USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
      USE LEVEL_COMMON
!***********************************************************************
!** The subroutine ALF (Automatic vibrational Level Finder) will
!   automatically generate the eigenvalues from the first vibrational
!   level (v=0) to a user specified level (v=KVMAX) or the highest
!   allowed vibrational level of a given smooth single (or double)
!   minimum potential (V). These energies are stored and returned to the
!   calling program in the molecular constants array GV(v=0-KVMAX).
!** For any errors that cannot be resolved within the subroutine, ALF
!   returns AFLAG with a value that defines which error had occured.
!** Uses the Schrodinger solver subroutine SCHRQas.
!
!** On entry:
!    NDP   is the number of datapoints used for the potential.
!    YMIN  is the innermost dimensionless radial distance
!    YH    is the dimensionless radial meshvalue
!    NCN   is the (integer) inverse power defining the linmiting attractive
!          long-range behaviour of the potential.  For a barrier, set NCN=99
!    V(i)  is the scaled input potential in 'AS' units
!    VLIM  is the potential asymptote (cm-1).
!    KVMAX is v for the highest vibrational level we wish to find.
!    AFLAG is rot.quantum J for the (centrifugally distorted) potential
!    ZMU   is the reduced mass of the diatom (amu).
!    EPS   is the energy convergence criterion (cm-1).
!    BFCT  it the internal unit scaling factor (2*mu/hbar^2)*RH^2.
!    INNODE specifies whether wave fx. initiation @ RMIN starts with a
!        note (normal case: INNODE > 0) or zero slope (when INNODE.le.0)
!    IWR    specifies the level of printing inside SCHRQ
!           <> 0 : print error & warning descriptions.
!           >= 1 : also print final eigenvalues & node count.
!           >= 2 : also show end-of-range wave function amplitudes.
!           >= 3 : print also intermediate trial eigenvalues, etc.
!
!** On exit:
!    KVMAX   is vib.quantum number for the highest vibrational level
!            found (may be less than the input value of KVMAX).
!    AFLAG   returns calculation outcome to calling program.
!            >=  0 : found all levels to v=KVMAX{input} & AFLAG= J
!             = -1 : KVMAX larger than number of levels found.
!    GV(v)   contains the vibrational energy levels found for v=0-KVMAX
!    INNR(v) labels each level as belonging to the inner (INNR = 1) or
!            outer (INNR = 0) well.
!
!** Flags: Modify only when debugging.
!    AWO   specifies the level of printing inside ALF
!          <> 0 : print error & warning descriptions.
!          >  0 : also print intermediate ALF messages.
!    INNER specifies wave function matching (& initiation) conditions.
!        .le.0 : Match inward & outward solutions at outermost well t.p.
!          > 0 : Match at innermost well inner turning point
!        For most normal cases set INNER = 0,  but ......
!            To find "inner-well-dominated" solutions of an asymmetric
!            double minimum potential, set  INNER > 0.
!    LPRWF specifies option of printing out generated wavefunction
!          > 0 : print wave function every LPRWF-th  point.
!          < 0 : compactly write to channel-7 every |LPRWF|-th wave
!                function value.
!          A lead "card" identifies the level, gives the position of
!          1-st point and radial mesh, & states No. of  points.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** The dimensioning parameters must be consistant with the sizes of the
!   arrays used in the calling program.
!
!    NVIBMX  is the maximum number of vibrational levels considered.
!            Note: NVIBMX should be larger than KVMAX.
!
      IMPLICIT NONE
      INTEGER NVIBMX
      PARAMETER (NVIBMX= 400)
!!!
      INTEGER NDIMR
!     PARAMETER (NDIMR= 200001)
! A limit set by the -fmax-stack-var-size in OpenMolcas is making arrays
! of the above size too large. If we can't get that increased, we could
! use an ALLOCATABLE array or use -frecursive.
!     PARAMETER (NDIMR= 131074)
!     REAL*8 PRV,ARV,RFN(NDIMR),YVB(NDIMR),DRDY2(NDIMR),FAS(NDIMR),
!    1                                         SDRDY(NDIMR),VBZ(NDIMR)
      REAL*8 PRV,ARV
!     REAL*8, ALLOCATABLE :: RFN(:),YVB(:),DRDY2(:),FAS(:),
!    1                                         SDRDY(:),VBZ(:)
      COMMON /BLKAS/PRV,ARV!,RFN,YVB,DRDY2,SDRDY,FAS,VBZ
!!!
!** NF counts levels found in automatic search option
!
      INTEGER NDP,KVMAX,KV,KVB,KVBB,AFLAG,NF,NBEG,NEND,NBEGG(0:NVIBMX), &
     &  NENDD(0:NVIBMX),INNR(0:NVIBMX),ICOR,IWR,IPMIN,IPMINN,           &
     &  I,LTRY,AWO,INNODE,INNER,LPRWF,JROT,NPMIN,NPMAX,NCN
!
      REAL*8 YMIN,YMAX,YH,V(NDP),SWF(NDP),VLIM,EO,ZMU,EPS,BFCT,GAMA,    &
     &  VMIN,VMAX,VME1,VME2,VME3,RE,PMAX, ESAV, ZPEHO, DGDV2, BMAX,     &
     &  GV(0:KVMAX),VPMIN(10),YPMIN(10),VPMAX(10),YPMAX(10)
      DATA AWO/1/,LPRWF/0/,KVB/-1/,KVBB/-2/
      NDIMR= 131074
!     CALL MMA_ALLOCATE(RFN,NDIMR,LABEL='RFN')
!     CALL MMA_ALLOCATE(YVB,NDIMR,LABEL='YVB')
!     CALL MMA_ALLOCATE(DRDY2,NDIMR,LABEL='DRDY2')
!     CALL MMA_ALLOCATE(FAS,NDIMR,LABEL='FAS')
!     CALL MMA_ALLOCATE(SDRDY,NDIMR,LABEL='SDRDY')
!     CALL MMA_ALLOCATE(VBZ,NDIMR,LABEL='VBZ')
      ipminn=0
! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
!     WRITE(6,*) ''
!     WRITE(6,*) 'NPP=',NDP
!     WRITE(6,*) 'YMIN=',YMIN
!     WRITE(6,*) 'YH=',YH
!     WRITE(6,*) 'NCN1=',NCN
!     DO I=1,3
!      WRITE(6,*) 'RVB=',RVB(I)
!      WRITE(6,*) 'VJ=',V(I)
!      WRITE(6,*) 'WF1=',SWF(I)
!      WRITE(6,*) 'GV=',GV(I)
!      WRITE(6,*) 'INNR=',INNR(I)
!     ENDDO
!     WRITE(6,*) 'VLIM1=',VLIM
!     WRITE(6,*) 'VMAX=',KVMAX
!     WRITE(6,*) 'AFLAG=',AFLAG
!!!!!! Don't comment the ZMU write statement, unless you want to remove
!!!!!! ZMU altogether:
      WRITE(6,*) 'ZMU=',ZMU
!!!!!!
!!!!!!
!     WRITE(6,*) 'EPS=',EPS
!     WRITE(6,*) 'BFCT=',BFCT
!     WRITE(6,*) 'INNOD1=',INNODE
!     WRITE(6,*) 'IWR=',IWR
!     WRITE(6,*) ''
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Check that the array dimensions are adequate.
      IF(KVMAX.GT.NVIBMX) THEN
          WRITE(6,602) KVMAX, NVIBMX
!         STOP
          CALL ABEND()
          ENDIF
!
!** Initialize remaining variables and flags. NF is label of level being sought
      NF= 0
      KVB= -1
      KV= 0
      INNER= 0
      LTRY= 0
!** Initialize level counters for each well.
      DO  I= 0,KVMAX
          INNR(I)= -1
          ENDDO
!** Store input rotational quantum number.
      JROT= AFLAG
      AFLAG= -1
!
!** YMAX is the outer radial distance over which potential is defined.
      YMAX= YMIN + DBLE(NDP-1)*YH
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Locate the potential minima.
      NPMIN= 0
      IPMIN= 2
      VMIN= 1.d99
      VME2= VBZ(2)
      VME3= VBZ(3)
      DO  I= 4,NDP-1
          VME1= VME2
          VME2= VME3
          VME3= VBZ(I)
          IF((VME2.LT.VME1).AND.(VME2.LT.VME3)) THEN
              NPMIN= NPMIN + 1
              YPMIN(NPMIN)= YVB(I)
              VPMIN(NPMIN)= VME2/BFCT
              IF(NPMIN.EQ.1) THEN
                  IPMIN= I
                  ENDIF
              IF(VPMIN(NPMIN).LT.VMIN) THEN
                  RE= YPMIN(NPMIN)
                  VMIN= VPMIN(NPMIN)
                  IPMINN= I
                  ENDIF
              IF(NPMIN.EQ.10) GOTO 10
              ENDIF
          END DO
   10 IF(NPMIN.EQ.0) THEN
          IF(V(2).LE.V(1)) THEN
!** If NO minimum & potential has negative slope, print a warning and stop.
              WRITE(6,608) JROT
              KVMAX= -1
              RETURN
              ENDIF
!...  but if potl. alway has positive slope, mesh point #1 is minimum
          NPMIN= 1
          IPMIN= 1
          YPMIN(NPMIN)= YVB(1)
          RE= YVB(1)
          VPMIN(NPMIN)= VBZ(1)
          VMIN= YPMIN(NPMIN)
          WRITE(6,618) VPMIN(1),YMIN
!         WRITE(6,*) 'Stuff about minima but error'
          ENDIF
!** Locate any potential maxima (if they exist).
      NPMAX= 0
      VMAX= -9.d99
      VME2= VBZ(IPMIN)
      VME3= VBZ(IPMIN+1)
      DO  I= IPMIN+2, NDP-1
          VME1= VME2
          VME2= VME3
          VME3= VBZ(I)
          IF((VME2.GT.VME1).AND.(VME2.GT.VME3)) THEN
              NPMAX= NPMAX + 1
              YPMAX(NPMAX)= YVB(I)
              VPMAX(NPMAX)= VME2/BFCT
              IF(VPMAX(NPMAX).GT.VMAX) VMAX= VPMAX(NPMAX)
              IF(NPMAX.EQ.10) GOTO 150
              ENDIF
          END DO
  150 IF((NPMAX.EQ.0).OR.                                               &
     &         ((NPMAX.GT.0).AND.(YPMAX(NPMAX).LT.YPMIN(NPMIN)))) THEN
!** If no maxima found or there is no barrier past outermost minimum,
!   set an energy maximum to be the value at the end of the radial range.
          NPMAX= NPMAX+ 1
          YPMAX(NPMAX)= YVB(NDP-1)
!?? should this limit be set at  VLIM ??
          VPMAX(NPMAX)= VBZ(NDP-1)/BFCT
          IF(VPMAX(NPMAX).GT.VMAX) VMAX= VPMAX(NPMAX)
          ENDIF
!
!** If innermost maximum lies inside innermost minimum, the potential
!   turns over in short range region OR have a minimim at mesh point #1:
!   PRINT a Warning
      IF(YPMAX(1).LT.YPMIN(1)) THEN
          WRITE(6,610) YPMAX(1)
          ENDIF
!
!** Otherwise, print out potential extrema count
      IF(NPMIN.GT.0) THEN
!         WRITE(6,*) 'Stuff about minima but gives error'
!         WRITE(6,614) NPMIN, (VPMIN(I),I= 1,NPMIN)
!         WRITE(6,616) (YPMIN(I), I= 1,NPMIN)
!         WRITE(6,*) 'Stuff about maximum but gives error'
!         WRITE(6,618) NPMAX, (VPMAX(I),I= 1,NPMAX)
!         WRITE(6,616) (YPMAX(I), I= 1,NPMAX)
          IF(NPMIN.GT.2) THEN
!** If potential has more than two minima - print warning & stop
              WRITE(6,620)
!cc           STOP
              ENDIF
          ENDIF
!** Set BMAX as barrier height of double-minimum potential
      BMAX= -9.d+09
      IF(NPMIN.GT.1) THEN
          DO  I= 1,NPMAX
              IF((YPMAX(I).GT.YPMIN(1)).AND.(YPMAX(I).LT.YPMIN(2)))     &
     &        BMAX= VPMAX(I)
              ENDDO
          ENDIF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*** Use harmonic approximation to estimate zero point energy.
      ZPEHO= DSQRT((VBZ(IPMINN+20)-VBZ(IPMINN))/400.d0)/BFCT
      EO= VMIN + ZPEHO
      WRITE(6,*) ''
      WRITE(6,634) BFCT
!     WRITE(6,*) 'IPMINN:                                 ',IPMINN
!     WRITE(6,*) 'VBZ(1):                                 ',VBZ(1)
!     WRITE(6,*) 'VBZ(IPMINN):                            ',VBZ(IPMINN)
!     WRITE(6,*) 'VBZ(IPMINN+20)',VBZ(IPMINN+20)
      WRITE(6,632) ZPEHO
      WRITE(6,*) 'Trial energy obtained from harmonic oscillator:  ',EO
!
!=========== Begin Actual Eigenvalue Calculation Loop Here =============
!** Compute eigenvalues ... etc. up to the KVMAX'th vibrational level.
!** When attempts to find the next eigenvalue fails, then perhaps the
!   next level is located in a second (inner) well. If so, then the
!   subroutine will set INNER = 1, and attempt to find that level.
!
      ICOR= 0
  100 KVBB= KVB
      KVB= KVBB ! Make sure it's "refernced"
      KVB= KV
      KV= NF
  110 ESAV= EO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Call subroutine SCHRQ to find eigenvalue EO and eigenfunction SWF(I).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
      WRITE(6,*) ''
      WRITE(6,*) 'Exiting alf.f'
      WRITE(6,*) 'Entering schrq.f'
!     WRITE(6,*) 'Entering schrq.f with the following parameters:'
!     WRITE(6,*) ''
!     WRITE(6,*) 'KV=',KV
!     WRITE(6,*) 'JROT=',JROT
!     WRITE(6,*) 'EO=',EO
!     WRITE(6,*) 'GAMA=',GAMA
!     WRITE(6,*) 'VMAX=',PMAX
!     WRITE(6,*) 'VLIM=',VLIM
!     DO I=1,3
!      WRITE(6,*) 'V=',V(I)
!      WRITE(6,*) 'WF=',SWF(I)
!     ENDDO
!     WRITE(6,*) 'BFCT=',BFCT
!     WRITE(6,*) 'EEPS=',EPS
!     WRITE(6,*) 'YMIN=',YMIN
!     WRITE(6,*) 'YH=',YH
!     WRITE(6,*) 'NPP=',NDP
!     WRITE(6,*) 'NBEG=',NBEG
!     WRITE(6,*) 'NEND=',NEND
!     WRITE(6,*) 'INNODE=',INNODE
!     WRITE(6,*) 'INNER=',INNER
!     WRITE(6,*) 'IWR=',IWR
!     WRITE(6,*) 'LPRWF=',LPRWF
      WRITE(6,*) ''
      CALL SCHRQas(KV,JROT,EO,GAMA,PMAX,VLIM,V,SWF,BFCT,EPS,YMIN,YH,NDP,&
     &                               NBEG,NEND,INNODE,INNER,IWR,LPRWF)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.LT.0) THEN
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** The SCHRQ error condition is KV < 0.  Allow for 3 cases:
!     EO > VMAX : energy from previous trial above potential maximum
!     NF = 0 : Looking for the first vibrational level (v = 0)
!     NF > 0 : Looking for the other vibrational levels (v > 0)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF(EO.GT.VMAX) THEN
!** For the case when the previous trial gave energy above the potential
!   maximum, make one last ditch attempt to find the highest bound level
!   (quasi or otherwise) in the potential.
              IF(LTRY.LT.1) THEN
                  LTRY= 1
                  KV= 999
                  EO= VMAX - 0.0001d0
                  GOTO 110
!... if that was unsuccessful, then print out a warning and exit.
                ELSE
                  WRITE(6,622) NF, EO, VMAX
                  KV= NF-1
                  GOTO 200
                ENDIF
              ENDIF
          WRITE(6,624) NF,JROT,ESAV
!.. eigenvalue of -9.9d9 signifies that eigenvalue search failed completely
          KVMAX= NF-1
          EO= -9.9d9
          RETURN
          ENDIF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** If calculated vibrational level is the desired level, NF, then ...
!   call SCECOR to calculate dG/dv and predict next higher level
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.EQ.NF) THEN
          NBEGG(KV)= NBEG
          NENDD(KV)= NEND
          GV(NF)= EO
          INNR(NF)= INNER
  120     NF= NF + 1
          IF(NF.LE.KVMAX) THEN
              IF(INNR(NF).GT.0) GOTO 120
!... if the next level was found earlier in overshoot ...
            ELSE
              IF(AWO.GT.0) WRITE(6,626) JROT,KVMAX
              AFLAG= JROT
              RETURN
            ENDIF
          ICOR= 0
          CALL SCECORas(KV,NF,JROT,INNER,ICOR,IWR,EO,YH,BFCT,NDP,       &
     &                                  NCN,VBZ,SDRDY,BMAX,VLIM,DGDV2)
          IF(EO.GT.VPMAX(NPMAX)) THEN
!... if estimated energy above highest barrier, set value below it
              EO=  VPMAX(NPMAX) - 0.05d0*DGDV2
              ICOR= 20
              ENDIF
          LTRY= 0
          KV= NF
          GOTO 100
          ENDIF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.NE.NF) THEN
!*** If last level found is not the desired one ...
          IF(INNR(KV).EQ.-1) THEN
!... Record vibrational level (if haven't already) for posterity.
              GV(KV)= EO
              INNR(KV)= INNER
              ENDIF
          ICOR= ICOR+1
          IF(ICOR.LE.20) THEN
!... Call subroutine using semiclassical methods to estimate correct energy
              CALL SCECORas(KV,NF,JROT,INNER,ICOR,IWR,EO,YH,BFCT,NDP,   &
     &                                  NCN,VBZ,SDRDY,BMAX,VLIM,DGDV2)
              IF(EO.GT.VPMAX(NPMAX)) THEN
!... if estimated energy above highest barrier, set value below it
                  KV= 999
                  EO=  VPMAX(NPMAX) - 0.05d0*DGDV2
              ENDIF
              GOTO 100
          ENDIF
!** If the calculated wavefunction is still for the wrong vibrational
!   level, then write out a warning return
          WRITE(6,628) NF,JROT
          KVMAX= NF-1
      ENDIF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  200 IF(AFLAG.LT.0) THEN
!** If unable to find all KVMAX+1 levels requested, then return KVMAX as
!  v for the highest vibrational level actually found, and print out the
!  the energy of that level.
          IF(AWO.NE.0) WRITE(6,630) KVMAX, GV(KVMAX)
! Make sure the following variables are "referenced":
          WRITE(6,*) NBEGG,NDIMR,NENDD,RE,YMAX
      ENDIF
!     CALL MMA_DEALLOCATE(RFN)
!     CALL MMA_DEALLOCATE(YVB)
!     CALL MMA_DEALLOCATE(DRDY2)
!     CALL MMA_DEALLOCATE(FAS)
!     CALL MMA_DEALLOCATE(SDRDY)
!     CALL MMA_DEALLOCATE(VBZ)
      RETURN
!-----------------------------------------------------------------------
  602 FORMAT(/'  *** ALF ERROR ***'/4X,'Number of vib levels requested='&
     & ,i4,' exceeds internal ALF array dimension  NVIBMX=',i4)
! 604 FORMAT(/' *** ALF ERROR ***   Find NO potential minima for   J=',
!    1  i4)
! 606 FORMAT(/'  ALF  finds onee potential minimum of',1PD15.7,
!    1  '  at  R(1)=',0Pf9.6)
  608 FORMAT(/'  *** ALF ERROR ***   Unable to find a potential minimum &
     & for   J=',i4)
  610 FORMAT(/'  *** ALF CAUTION ***'/ 4X,'The potential turns over in t&
     &he short range region at  y= ',G15.8)
! 614 FORMAT(' Find',F3.5,'  potential minima:   Vmin=',8F11.3)
! 616 FORMAT(19x,'located at   y =',8f11.5)
  618 FORMAT(' Find',I3,'  potential maxima:   Vmax=',8F11.3)
  620 FORMAT(' *** So  STOP !!!!')
  622 FORMAT(/' ALF search finds next estimated trial energy  E(v=',I3, &
     & ')=',G15.8/8X,'lies above potential maximum or asymptote at  VMAX&
     &=',G15.8)
  624 FORMAT(/' *** SCHRQ FAILS in ALF when searching for  v=',i3,      &
     &  ' J=',i3,'   with   EO=',f9.3/5x,'Check range and/or contact    &
     &. Nike Dattani [nike@hpqc.org,ndattani@uwaterloo.ca]')
  626 FORMAT(/' ALF successfully finds all (J=',i3,') vibrational levels&
     & up to   v= KVMAX=',I3)
  628 FORMAT(4x,'ALF fails to find level   v=',i3,', J=',i3)
  630 FORMAT(' Highest calculated level found by ALF is   E(v=',I3,     &
     &  ')=',1PD17.9 /)
  632 FORMAT(' Zero point energy (measured from VLIM) approximated using&
     & a harmonic osccilator:        ',8F11.3)
  634 FORMAT(' Mult. V(R) by this factor (BFCT) for solving the SE in   &
     & dimensionless units: ',E20.13)
      END
