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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine ALFas(NDP,YMIN,YH,NCN,V,SWF,VLIM,KVMAX,AFLAG,ZMU,EPS,GV,BFCT,INNODE,INNR,IWR)
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
!        note (normal case: INNODE > 0) or zero slope (when INNODE <= 0)
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
!         <= 0 : Match inward & outward solutions at outermost well t.p.
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

!use stdalloc, only: mma_allocate, mma_deallocate
use LEVEL_COMMON

implicit none
integer NVIBMX
parameter(NVIBMX=400)
integer NDIMR
!parameter (NDIMR=200001)
! A limit set by the -fmax-stack-var-size in OpenMolcas is making arrays
! of the above size too large. If we can't get that increased, we could
! use an ALLOCATABLE array or use -frecursive.
!parameter (NDIMR=131074)
!real*8 PRV, ARV, RFN(NDIMR), YVB(NDIMR), DRDY2(NDIMR), FAS(NDIMR), SDRDY(NDIMR), VBZ(NDIMR)
real*8 PRV, ARV
!real*8, allocatable :: RFN(:), YVB(:), DRDY2(:), FAS(:), SDRDY(:), VBZ(:)
common/BLKAS/PRV,ARV!,RFN,YVB,DRDY2,SDRDY,FAS,VBZ
!** NF counts levels found in automatic search option
integer NDP, KVMAX, KV, KVB, KVBB, AFLAG, NF, NBEG, NEND, NBEGG(0:NVIBMX), NENDD(0:NVIBMX), INNR(0:NVIBMX), ICOR, IWR, IPMIN, &
        IPMINN, I, LTRY, AWO, INNODE, INNER, LPRWF, JROT, NPMIN, NPMAX, NCN
real*8 YMIN, YMAX, YH, V(NDP), SWF(NDP), VLIM, EO, ZMU, EPS, BFCT, GAMA, VMIN, VMAX, VME1, VME2, VME3, RE, PMAX, ESAV, ZPEHO, &
       DGDV2, BMAX, GV(0:KVMAX), VPMIN(10), YPMIN(10), VPMAX(10), YPMAX(10)
data AWO/1/,LPRWF/0/,KVB/-1/,KVBB/-2/

NDIMR = 131074
!call mma_allocate(RFN,NDIMR,LABEL='RFN')
!call mma_allocate(YVB,NDIMR,LABEL='YVB')
!call mma_allocate(DRDY2,NDIMR,LABEL='DRDY2')
!call mma_allocate(FAS,NDIMR,LABEL='FAS')
!call mma_allocate(SDRDY,NDIMR,LABEL='SDRDY')
!call mma_allocate(VBZ,NDIMR,LABEL='VBZ')
ipminn = 0
! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
!write(6,*) ''
!write(6,*) 'NPP=',NDP
!write(6,*) 'YMIN=',YMIN
!write(6,*) 'YH=',YH
!write(6,*) 'NCN1=',NCN
!do I=1,3
!  write(6,*) 'RVB=',RVB(I)
!  write(6,*) 'VJ=',V(I)
!  write(6,*) 'WF1=',SWF(I)
!  write(6,*) 'GV=',GV(I)
!  write(6,*) 'INNR=',INNR(I)
!end do
!write(6,*) 'VLIM1=',VLIM
!write(6,*) 'VMAX=',KVMAX
!write(6,*) 'AFLAG=',AFLAG
! Don't comment the ZMU write statement, unless you want to remove
! ZMU altogether:
write(6,*) 'ZMU=',ZMU
!
!write(6,*) 'EPS=',EPS
!write(6,*) 'BFCT=',BFCT
!write(6,*) 'INNOD1=',INNODE
!write(6,*) 'IWR=',IWR
!write(6,*) ''
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Check that the array dimensions are adequate.
if (KVMAX > NVIBMX) then
  write(6,602) KVMAX,NVIBMX
  !stop
  call ABEND()
end if

! Initialize remaining variables and flags. NF is label of level being sought
NF = 0
KVB = -1
KV = 0
INNER = 0
LTRY = 0
! Initialize level counters for each well.
do I=0,KVMAX
  INNR(I) = -1
end do
! Store input rotational quantum number.
JROT = AFLAG
AFLAG = -1

! YMAX is the outer radial distance over which potential is defined.
YMAX = YMIN+dble(NDP-1)*YH
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Locate the potential minima.
NPMIN = 0
IPMIN = 2
VMIN = 1.d99
VME2 = VBZ(2)
VME3 = VBZ(3)
do I=4,NDP-1
  VME1 = VME2
  VME2 = VME3
  VME3 = VBZ(I)
  if ((VME2 < VME1) .and. (VME2 < VME3)) then
    NPMIN = NPMIN+1
    YPMIN(NPMIN) = YVB(I)
    VPMIN(NPMIN) = VME2/BFCT
    if (NPMIN == 1) then
      IPMIN = I
    end if
    if (VPMIN(NPMIN) < VMIN) then
      RE = YPMIN(NPMIN)
      VMIN = VPMIN(NPMIN)
      IPMINN = I
    end if
    if (NPMIN == 10) goto 10
  end if
end do
10 continue
if (NPMIN == 0) then
  if (V(2) <= V(1)) then
    ! If NO minimum & potential has negative slope, print a warning and stop.
    write(6,608) JROT
    KVMAX = -1
    return
  end if
  !...  but if potl. alway has positive slope, mesh point #1 is minimum
  NPMIN = 1
  IPMIN = 1
  YPMIN(NPMIN) = YVB(1)
  RE = YVB(1)
  VPMIN(NPMIN) = VBZ(1)
  VMIN = YPMIN(NPMIN)
  write(6,618) VPMIN(1),YMIN
  !write(6,*) 'Stuff about minima but error'
end if
! Locate any potential maxima (if they exist).
NPMAX = 0
VMAX = -9.d99
VME2 = VBZ(IPMIN)
VME3 = VBZ(IPMIN+1)
do I=IPMIN+2,NDP-1
  VME1 = VME2
  VME2 = VME3
  VME3 = VBZ(I)
  if ((VME2 > VME1) .and. (VME2 > VME3)) then
    NPMAX = NPMAX+1
    YPMAX(NPMAX) = YVB(I)
    VPMAX(NPMAX) = VME2/BFCT
    if (VPMAX(NPMAX) > VMAX) VMAX = VPMAX(NPMAX)
    if (NPMAX == 10) goto 150
  end if
end do
150 continue
if ((NPMAX == 0) .or. ((NPMAX > 0) .and. (YPMAX(NPMAX) < YPMIN(NPMIN)))) then
  ! If no maxima found or there is no barrier past outermost minimum,
  ! set an energy maximum to be the value at the end of the radial range.
  NPMAX = NPMAX+1
  YPMAX(NPMAX) = YVB(NDP-1)
  !?? should this limit be set at  VLIM ??
  VPMAX(NPMAX) = VBZ(NDP-1)/BFCT
  if (VPMAX(NPMAX) > VMAX) VMAX = VPMAX(NPMAX)
end if

! If innermost maximum lies inside innermost minimum, the potential
! turns over in short range region OR have a minimim at mesh point #1:
! PRINT a Warning
if (YPMAX(1) < YPMIN(1)) then
  write(6,610) YPMAX(1)
end if

! Otherwise, print out potential extrema count
if (NPMIN > 0) then
  !write(6,*) 'Stuff about minima but gives error'
  !write(6,614) NPMIN,(VPMIN(I),I= 1,NPMIN)
  !write(6,616) (YPMIN(I),I=1,NPMIN)
  !write(6,*) 'Stuff about maximum but gives error'
  !write(6,618) NPMAX,(VPMAX(I),I=1,NPMAX)
  !write(6,616) (YPMAX(I),I=1,NPMAX)
  if (NPMIN > 2) then
    ! If potential has more than two minima - print warning & stop
    write(6,620)
    !stop
  end if
end if
!** Set BMAX as barrier height of double-minimum potential
BMAX = -9.d+09
if (NPMIN > 1) then
  do I=1,NPMAX
    if ((YPMAX(I) > YPMIN(1)) .and. (YPMAX(I) < YPMIN(2))) BMAX = VPMAX(I)
  end do
end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*** Use harmonic approximation to estimate zero point energy.
ZPEHO = dsqrt((VBZ(IPMINN+20)-VBZ(IPMINN))/400.d0)/BFCT
EO = VMIN+ZPEHO
write(6,*) ''
write(6,634) BFCT
!write(6,*) 'IPMINN:                                 ',IPMINN
!write(6,*) 'VBZ(1):                                 ',VBZ(1)
!write(6,*) 'VBZ(IPMINN):                            ',VBZ(IPMINN)
!write(6,*) 'VBZ(IPMINN+20)',VBZ(IPMINN+20)
write(6,632) ZPEHO
write(6,*) 'Trial energy obtained from harmonic oscillator:  ',EO

!=========== Begin Actual Eigenvalue Calculation Loop Here =============
! Compute eigenvalues ... etc. up to the KVMAX'th vibrational level.
! When attempts to find the next eigenvalue fails, then perhaps the
! next level is located in a second (inner) well. If so, then the
! subroutine will set INNER = 1, and attempt to find that level.

ICOR = 0
100 continue
KVBB = KVB
KVB = KVBB ! Make sure it's "referenced"
KVB = KV
KV = NF
110 continue
ESAV = EO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Call subroutine SCHRQ to find eigenvalue EO and eigenfunction SWF(I).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
write(6,*) ''
write(6,*) 'Exiting alf.f'
write(6,*) 'Entering schrq.f'
!write(6,*) 'Entering schrq.f with the following parameters:'
!write(6,*) ''
!write(6,*) 'KV=',KV
!write(6,*) 'JROT=',JROT
!write(6,*) 'EO=',EO
!write(6,*) 'GAMA=',GAMA
!write(6,*) 'VMAX=',PMAX
!write(6,*) 'VLIM=',VLIM
!do I=1,3
!  write(6,*) 'V=',V(I)
!  write(6,*) 'WF=',SWF(I)
!end do
!write(6,*) 'BFCT=',BFCT
!write(6,*) 'EEPS=',EPS
!write(6,*) 'YMIN=',YMIN
!write(6,*) 'YH=',YH
!write(6,*) 'NPP=',NDP
!write(6,*) 'NBEG=',NBEG
!write(6,*) 'NEND=',NEND
!write(6,*) 'INNODE=',INNODE
!write(6,*) 'INNER=',INNER
!write(6,*) 'IWR=',IWR
!write(6,*) 'LPRWF=',LPRWF
write(6,*) ''
call SCHRQas(KV,JROT,EO,GAMA,PMAX,VLIM,V,SWF,BFCT,EPS,YMIN,YH,NDP,NBEG,NEND,INNODE,INNER,IWR,LPRWF)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (KV < 0) then
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! The SCHRQ error condition is KV < 0.  Allow for 3 cases:
  !   EO > VMAX : energy from previous trial above potential maximum
  !   NF = 0 : Looking for the first vibrational level (v = 0)
  !   NF > 0 : Looking for the other vibrational levels (v > 0)
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (EO > VMAX) then
    ! For the case when the previous trial gave energy above the potential
    ! maximum, make one last ditch attempt to find the highest bound level
    ! (quasi or otherwise) in the potential.
    if (LTRY < 1) then
      LTRY = 1
      KV = 999
      EO = VMAX-0.0001d0
      goto 110
    else
      !... if that was unsuccessful, then print out a warning and exit.
      write(6,622) NF,EO,VMAX
      KV = NF-1
      goto 200
    end if
  end if
  write(6,624) NF,JROT,ESAV
  !.. eigenvalue of -9.9d9 signifies that eigenvalue search failed completely
  KVMAX = NF-1
  EO = -9.9d9
  return
end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! If calculated vibrational level is the desired level, NF, then ...
! call SCECOR to calculate dG/dv and predict next higher level
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (KV == NF) then
  NBEGG(KV) = NBEG
  NENDD(KV) = NEND
  GV(NF) = EO
  INNR(NF) = INNER
  120 continue
  NF = NF+1
  if (NF <= KVMAX) then
    if (INNR(NF) > 0) goto 120
  else
    !... if the next level was found earlier in overshoot ...
    if (AWO > 0) write(6,626) JROT,KVMAX
    AFLAG = JROT
    return
  end if
  ICOR = 0
  call SCECORas(KV,NF,JROT,INNER,ICOR,IWR,EO,YH,BFCT,NDP,NCN,VBZ,SDRDY,BMAX,VLIM,DGDV2)
  if (EO > VPMAX(NPMAX)) then
    !... if estimated energy above highest barrier, set value below it
    EO = VPMAX(NPMAX)-0.05d0*DGDV2
    ICOR = 20
  end if
  LTRY = 0
  KV = NF
  goto 100
end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (KV /= NF) then
  ! If last level found is not the desired one ...
  if (INNR(KV) == -1) then
    !... Record vibrational level (if haven't already) for posterity.
    GV(KV) = EO
    INNR(KV) = INNER
  end if
  ICOR = ICOR+1
  if (ICOR <= 20) then
    !... Call subroutine using semiclassical methods to estimate correct energy
    call SCECORas(KV,NF,JROT,INNER,ICOR,IWR,EO,YH,BFCT,NDP,NCN,VBZ,SDRDY,BMAX,VLIM,DGDV2)
    if (EO > VPMAX(NPMAX)) then
      !... if estimated energy above highest barrier, set value below it
      KV = 999
      EO = VPMAX(NPMAX)-0.05d0*DGDV2
    end if
    goto 100
  end if
  ! If the calculated wavefunction is still for the wrong vibrational
  ! level, then write out a warning return
  write(6,628) NF,JROT
  KVMAX = NF-1
end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
200 continue
if (AFLAG < 0) then
  ! If unable to find all KVMAX+1 levels requested, then return KVMAX as
  ! v for the highest vibrational level actually found, and print out the
  ! the energy of that level.
  if (AWO /= 0) write(6,630) KVMAX,GV(KVMAX)
  ! Make sure the following variables are "referenced":
  write(6,*) NBEGG,NDIMR,NENDD,RE,YMAX
end if
!call mma_deallocate(RFN)
!call mma_deallocate(YVB)
!call mma_deallocate(DRDY2)
!call mma_deallocate(FAS)
!call mma_deallocate(SDRDY)
!call mma_deallocate(VBZ)
return
!-----------------------------------------------------------------------
602 format(/'  *** ALF ERROR ***'/4X,'Number of vib levels requested=',i4,' exceeds internal ALF array dimension  NVIBMX=',i4)
!604 format(/' *** ALF ERROR ***   Find NO potential minima for   J=',i4)
!606 format(/'  ALF  finds onee potential minimum of',1PD15.7,'  at  R(1)=',0Pf9.6)
608 format(/'  *** ALF ERROR ***   Unable to find a potential minimum for   J=',i4)
610 format(/'  *** ALF CAUTION ***'/4X,'The potential turns over in the short range region at  y= ',G15.8)
!614 format(' Find',F3.5,'  potential minima:   Vmin=',8F11.3)
!616 format(19x,'located at   y =',8f11.5)
618 format(' Find',I3,'  potential maxima:   Vmax=',8F11.3)
620 format(' *** So  STOP !!!!')
622 format(/' ALF search finds next estimated trial energy  E(v=',I3,')=',G15.8/8X, &
           'lies above potential maximum or asymptote at  VMAX=',G15.8)
624 format(/' *** SCHRQ FAILS in ALF when searching for  v=',i3,' J=',i3,'   with   EO=',f9.3/5x, &
           'Check range and/or contact Nike Dattani [nike@hpqc.org,ndattani@uwaterloo.ca]')
626 format(/' ALF successfully finds all (J=',i3,') vibrational levels up to   v= KVMAX=',I3)
628 format(4x,'ALF fails to find level   v=',i3,', J=',i3)
630 format(' Highest calculated level found by ALF is   E(v=',I3,')=',1PD17.9/)
632 format(' Zero point energy (measured from VLIM) approximated using a harmonic osccilator:        ',8F11.3)
634 format(' Mult. V(R) by this factor (BFCT) for solving the SE in dimensionless units: ',E20.13)

end subroutine ALFas
