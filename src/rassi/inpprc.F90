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

subroutine INPPRC()

use Index_Functions, only: nTri_Elem
use rasdef, only: NRS1, NRS1T, NRS2, NRS2T, NRS3, NRS3T
use rassi_global_arrays, only: ESHFT, HAM, HDIAG, JBNUM, LROOT
use rassi_aux, only: AO_Mode, CMO1, CMO2, DMAB, ipglob, jDisk_TDM, JOB_INDEX, mTRA
use rassi_data, only: NASH, NASHT, NBASF, NBMX, NBSQ, NBSQPR, NBST, NBTRI, NCMO, NISH, NISHT, NOSH, NSSH, NSSHT, NTDMAB, NTDMZZ, &
                      NTRA
use kVectors, only: nk_Vector
use Lebedev_quadrature, only: order_table
use OneDat, only: sOpSiz, sRdFst, sRdNxt
use Cntrl, only: BINA, Coor, Do_SK, DO_TMOM, DQVD, FORCE_NON_AO_TDM, HAVE_DIAG, HAVE_HEFF, HEff, IBINA, ICOMP, IFDCPL, IFEJOB, &
                 IFGCAL, IFHAM, IFHCOM, IFHDIA, IFHEFF, IFHEXT, IFJ2, IFJZ, IFMCAL, IFSHFT, IFSO, IFTDM, IFTRD1, IFXCAL, IPUSED, &
                 IRREP, ISOCMP, L_Eff, LuTDM, MLTPLT, MXPROP, NATO, nAtoms, NBINA, NJOB, NOHAM, NPROP, NQUAD, NrNATO, NSOPR, &
                 nSOThr_PRT, nState, ONLY_OVERLAPS, PNAME, PRCI, PRMEE, PRMER, PRMES, PRORB, PRSXY, PRTRA, PRXVE, PRXVR, PRXVS, &
                 PTYPE, RefEne, RFPert, SAVEDENS, SONATNSTATE, SONTOSTATES, SOPRNM, SOPRTP, SOThr_PRT, ToFile, TRACK
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: I, I1, I2, IADD, IBYTE, ICMP, ICMPLST(MXPROP), IDUM(1), IERR, II, III, IMISS, IOPT, IPROP, IPRP, IRC, ISOPR, &
                     ISYLAB, J, JOB1, JOB2, MISSAMX, MISSAMY, MISSAMZ, MPROP, MSOPR, N, NATOM, NPRPLST
real(kind=wp) :: XAXIS, ZAXIS
logical(kind=iwp) :: IsAvail(MXPROP), IsAvailSO(MXPROP), JOBMATCH
character(len=8) :: LABEL, LABEL2, PRPLST(MXPROP)
character(len=3) :: lIrrep(8)
integer(kind=iwp), external :: IsFreeUnit

! Analysing and post-processing the input that was read in readin_rassi.

call mma_allocate(jDisk_TDM,2,nTri_Elem(nState),Label='jDisk_TDM')
jDisk_TDM(1,:) = -1
jDisk_TDM(2,:) = 0
! PAM07: The printing of spin-orbit Hamiltonian matrix elements:
! If no value for SOTHR_PRT was given in the input, it has a
! negative value that was set in init_rassi:
if (SOTHR_PRT < Zero) then
  ! Assign default settings. SOTHR_PRT is in cm-1:
  if (IPGLOB >= 4) then
    NSOTHR_PRT = 10000
    SOTHR_PRT = 1.0e-4_wp
  else if (IPGLOB >= 3) then
    NSOTHR_PRT = 100
    SOTHR_PRT = 0.01_wp
  else if (IPGLOB >= 2) then
    NSOTHR_PRT = 20
    SOTHR_PRT = One
  else
    NSOTHR_PRT = 0
  end if
end if

! Some sizes:
NBMX = maxval(NBASF(1:nIrrep))
NBSQ = 0
IPRP = 0
do I=1,nIrrep
  NBSQPR(I) = NBSQ
  NBSQ = NBSQ+NBASF(I)**2
end do
NBTRI = (NBSQ+NBST)/2
! Sizes of some data sets:
! ACTUAL SIZES OF TDMAB AND TDMZZ DEPENDS ON BOTH LSYM1 AND LSYM2.
! HOWEVER, MAX POSSIBLE SIZE IS WHEN LSYM1=LSYM2.
NCMO = sum(NOSH(1:nIrrep)*NBASF(1:nIrrep))
NTRA = sum(NOSH(1:nIrrep)**2)
NTDMZZ = sum(NBASF(1:nIrrep)**2)
NTDMAB = NTRA
SaveDens = (IFTRD1 .or. IFTDM) .or. (SONATNSTATE > 0) .or. (SONTOSTATES > 0) .or. NATO .or. Do_TMOM
if (SaveDens) then
  write(u6,*)
  write(u6,*) ' Info: creating TDMFILE'
  LUTDM = IsFreeUnit(21)
  call DANAME_MF(LUTDM,'TDMFILE')
  AO_Mode = .true.
  iByte = 8*3*nTri_Elem(nstate-1)*nTDMZZ

  ! For the time we will move over to compact mode if the required
  ! estimate of disk space if more than 1 Gb.

  if (iByte > 1024**3) AO_Mode = .false.
  ! Force for debugging purpose.
  if (Force_NON_AO_TDM) AO_Mode = .false.
  write(u6,*) '       estimated file size ',iByte/1024,'kB'

  ! For small basis set with symmetry we might not benefit from
  ! doing this.

  if (NASHT**2+1 > nTDMAB) AO_Mode = .true.

  if (.not. AO_Mode) then
    write(u6,*) '       TDMs in reduced format'
    call mma_allocate(JOB_INDEX,nState,Label='JOB_INDEX')
    JOB_INDEX(:) = JBNUM(1:nState)
    !write(u6,*) 'Job_Index=',Job_Index
    call mma_allocate(CMO1,nCMO,Label='CMO1')
    call mma_allocate(CMO2,nCMO,Label='CMO2')
    call mma_allocate(DMAB,nTDMAB,Label='DMAB')
    mTRA = nTRA
  else
    write(u6,*) '       TDMs in AO format'
  end if
  write(u6,*)
end if

! Upcase property names in lists of requests:
do IPROP=1,NPROP
  call UPCASE(PNAME(IPROP))
end do
do ISOPR=1,NSOPR
  call UPCASE(SOPRNM(ISOPR))
end do

! Which properties are available in the ONEINT file?
! (IPUSED will be set later, set it to zero now.)
!write(u6,*) 'Which properties are available in the ONEINT file?'
IPRP = 0
IRC = -1
IOPT = ibset(ibset(0,sOpSiz),sRdFst)
LABEL = 'UNDEF'
call iRDONE(IRC,IOPT,LABEL,ICMP,IDUM,ISYLAB)
if (IRC == 0) then
  IPRP = 1
  call UPCASE(LABEL)
  PRPLST(1) = LABEL
  ICMPLST(1) = ICMP
  IPUSED(1) = 0

  do I=1,MXPROP
    if (IPRP >= MXPROP) exit
    IRC = -1
    IOPT = ibset(ibset(0,sOpSiz),sRdNxt)
    call iRDONE(IRC,IOPT,LABEL,ICMP,IDUM,ISYLAB)
    if (IRC /= 0) exit
    IPRP = IPRP+1
    call UPCASE(LABEL)
    PRPLST(IPRP) = LABEL
    ICMPLST(IPRP) = ICMP
    IPUSED(IPRP) = 0

    ! Copy the EF2 integral label for non-rel hyperfine calculations
    if (LABEL(1:3) == 'EF2') then
      if (IPRP >= MXPROP) exit
      IPRP = IPRP+1
      LABEL2 = LABEL
      LABEL2(1:4) = 'ASDO'
      PRPLST(IPRP) = LABEL2
      ICMPLST(IPRP) = ICMP
      IPUSED(IPRP) = 0
    end if

    ! Now the ASD are calculated from X2C magnetic integrals
    if ((LABEL(1:5) == 'MAGXP') .and. (ICMP <= 6)) then
      if (IPRP >= MXPROP) exit
      IPRP = IPRP+1
      LABEL2 = LABEL
      LABEL2(1:5) = 'ASD  '
      PRPLST(IPRP) = LABEL2
      ICMPLST(IPRP) = ICMP
      IPUSED(IPRP) = 0
    end if

    if (LABEL(1:4) == 'PSOI') then
      if (IPRP >= MXPROP) exit
      IPRP = IPRP+1
      LABEL2 = LABEL
      LABEL2(1:4) = 'PSOP'
      PRPLST(IPRP) = LABEL2
      ICMPLST(IPRP) = ICMP
      IPUSED(IPRP) = 0
    end if
    if (LABEL(1:6) == 'DMS  1') then
      if (IPRP >= MXPROP) exit
      IPRP = IPRP+1
      LABEL2 = LABEL
      LABEL2(1:6) = 'DMP   '
      PRPLST(IPRP) = LABEL2
      ICMPLST(IPRP) = ICMP
      IPUSED(IPRP) = 0
    end if
  end do
end if
NPRPLST = IPRP

! Add empty slots for on-the-fly TM integrals.

! If the RASSI code is run several instances on the same job some
! of these labels will already be available on the file and need
! not to be added to the list.

if (Do_TMOM .and. (PRPLST(IPRP)(1:4) /= 'TMOM')) then
  PRPLST(IPRP+1) = 'TMOM0  R'
  ICMPLST(IPRP+1) = 1
  IPUSED(IPRP+1) = 0
  PRPLST(IPRP+2) = 'TMOM0  I'
  ICMPLST(IPRP+2) = 1
  IPUSED(IPRP+2) = 0
  IPRP = IPRP+2

  PRPLST(IPRP+1) = 'TMOM  RS'
  ICMPLST(IPRP+1) = 1
  IPUSED(IPRP+1) = 0
  PRPLST(IPRP+2) = 'TMOM  RS'
  ICMPLST(IPRP+2) = 2
  IPUSED(IPRP+2) = 0
  PRPLST(IPRP+3) = 'TMOM  RS'
  ICMPLST(IPRP+3) = 3
  IPUSED(IPRP+3) = 0
  PRPLST(IPRP+4) = 'TMOM  RA'
  ICMPLST(IPRP+4) = 1
  IPUSED(IPRP+4) = 0
  PRPLST(IPRP+5) = 'TMOM  RA'
  ICMPLST(IPRP+5) = 2
  IPUSED(IPRP+5) = 0
  PRPLST(IPRP+6) = 'TMOM  RA'
  ICMPLST(IPRP+6) = 3
  IPUSED(IPRP+6) = 0
  PRPLST(IPRP+7) = 'TMOM  IS'
  ICMPLST(IPRP+7) = 1
  IPUSED(IPRP+7) = 0
  PRPLST(IPRP+8) = 'TMOM  IS'
  ICMPLST(IPRP+8) = 2
  IPUSED(IPRP+8) = 0
  PRPLST(IPRP+9) = 'TMOM  IS'
  ICMPLST(IPRP+9) = 3
  IPUSED(IPRP+9) = 0
  PRPLST(IPRP+10) = 'TMOM  IA'
  ICMPLST(IPRP+10) = 1
  IPUSED(IPRP+10) = 0
  PRPLST(IPRP+11) = 'TMOM  IA'
  ICMPLST(IPRP+11) = 2
  IPUSED(IPRP+11) = 0
  PRPLST(IPRP+12) = 'TMOM  IA'
  ICMPLST(IPRP+12) = 3
  IPUSED(IPRP+12) = 0
  IPRP = IPRP+12

  ! Not currently in use.
  !PRPLST(IPRP+1) = 'TMOM2  R'
  !ICMPLST(IPRP+1) = 1
  !IPUSED(IPRP+1) = 0
  !PRPLST(IPRP+ 2) = 'TMOM2  R'
  !ICMPLST(IPRP+2) = 2
  !IPUSED(IPRP+2) = 0
  !PRPLST(IPRP+3) = 'TMOM2  R'
  !ICMPLST(IPRP+3) = 3
  !IPUSED(IPRP+3) = 0
  !PRPLST(IPRP+4) = 'TMOM2  I'
  !ICMPLST(IPRP+4) = 1
  !IPUSED(IPRP+4) = 0
  !PRPLST(IPRP+5) = 'TMOM2  I'
  !ICMPLST(IPRP+5) = 2
  !IPUSED(IPRP+5) = 0
  !PRPLST(IPRP+6) = 'TMOM2  I'
  !ICMPLST(IPRP+6) = 3
  !IPUSED(IPRP+6) = 0
  !IPRP = IPRP+6
!else if (Do_TMOM) then
!  PRPLST(IPRP+1) = 'TMOM0  R'
!  ICMPLST(IPRP+1) = 1
!  IPUSED(IPRP+1) = 0
!  PRPLST(IPRP+2) = 'TMOM0  I'
!  ICMPLST(IPRP+2) = 1
!  IPUSED(IPRP+2) = 0
!  IPRP = IPRP+2
!
!  ! Not currently in use.
!  PRPLST(IPRP+1) = 'TMOM2  R'
!  ICMPLST(IPRP+1) = 2
!  IPUSED(IPRP+1) = 0
!  PRPLST(IPRP+2) = 'TMOM2  R'
!  ICMPLST(IPRP+2) = 3
!  IPUSED(IPRP+2) = 0
!  PRPLST(IPRP+3) = 'TMOM2  I'
!  ICMPLST(IPRP+3) = 2
!  IPUSED(IPRP+3) = 0
!  PRPLST(IPRP+4) = 'TMOM2  I'
!  ICMPLST(IPRP+4) = 3
!  IPUSED(IPRP+4) = 0
!  IPRP = IPRP+4
end if

NPRPLST = IPRP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
! Add some property names by defaults, if no input:                    C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

if (NPROP == 0) then

  if (NSOPR == 0) then
    ! If no input at all, use this selection:
    do IPRP=1,NPRPLST

      if ((PRPLST(IPRP) == 'MLTPL  0') .or. (PRPLST(IPRP) == 'MLTPL  1') .or. (PRPLST(IPRP) == 'MLTPL  2') .or. &
          (PRPLST(IPRP) == 'MLTPL  3') .or. (PRPLST(IPRP) == 'OMQ') .or. (PRPLST(IPRP) == 'ANGMOM') .or. &
          (PRPLST(IPRP)(1:4) == 'TMOM') .or. (PRPLST(IPRP) == 'VELOCITY') .or. (PRPLST(IPRP) == 'MLTPV  2') .or. &
          (PRPLST(IPRP)(1:4) == 'EMFR')) then

        NSOPR = NSOPR+1
        SOPRNM(NSOPR) = PRPLST(IPRP)

        ! SET the proper default type label. That is, which WE-reduced density
        ! should the integral be contracted with. If not set here it will
        ! default to a Hermitian singlet reduced density.
        if (PRPLST(IPRP) == 'ANGMOM') then
          SOPRTP(NSOPR) = 'ANTISING'
        else if (PRPLST(IPRP) == 'VELOCITY') then
          SOPRTP(NSOPR) = 'ANTISING'
        else if (PRPLST(IPRP) == 'OMQ') then
          SOPRTP(NSOPR) = 'ANTISING'
        else if (PRPLST(IPRP) == 'MLTPV  2') then
          SOPRTP(NSOPR) = 'ANTISING'
        else if (PRPLST(IPRP) == 'TMOM  RS') then
          SOPRTP(NSOPR) = 'HERMSING'
        else if (PRPLST(IPRP) == 'TMOM  RA') then
          SOPRTP(NSOPR) = 'ANTISING'
        else if (PRPLST(IPRP) == 'TMOM  IS') then
          SOPRTP(NSOPR) = 'HERMSING'
        else if (PRPLST(IPRP) == 'TMOM  IA') then
          SOPRTP(NSOPR) = 'ANTISING'
        else if (PRPLST(IPRP) == 'TMOM0  R') then
          SOPRTP(NSOPR) = 'HERMTRIP'
        else if (PRPLST(IPRP) == 'TMOM0  I') then
          SOPRTP(NSOPR) = 'HERMTRIP'
        else
          SOPRTP(NSOPR) = 'HERMSING'
        end if

        ! Set the number of elements of the property
        ISOCMP(NSOPR) = ICMPLST(IPRP)

        ! Add properties for the explicit spin part of the transition moments
        ! in the case of spin-orbit coupled wave functrions.

        if (IFSO) then
          if (PRPLST(IPRP) == 'MLTPL  0') then
            NSOPR = NSOPR+1
            SOPRNM(NSOPR) = PRPLST(IPRP)
            ISOCMP(NSOPR) = ICMPLST(IPRP)
            SOPRTP(NSOPR) = 'ANTITRIP'
          else if (PRPLST(IPRP) == 'MLTPL  1') then
            NSOPR = NSOPR+1
            SOPRNM(NSOPR) = PRPLST(IPRP)
            ISOCMP(NSOPR) = ICMPLST(IPRP)
            SOPRTP(NSOPR) = 'ANTITRIP'
          ! Uncomment to activate more options!!!
          !else if (PRPLST(IPRP) == 'MLTPL  2') then
          !  NSOPR = NSOPR+1
          !  SOPRNM(NSOPR) = PRPLST(IPRP)
          !  ISOCMP(NSOPR) = ICMPLST(IPRP)
          !  SOPRTP(NSOPR) = 'HERMTRIP'
          !else if (PRPLST(IPRP) == 'ANGMOM') then
          !  NSOPR = NSOPR+1
          !  SOPRNM(NSOPR) = PRPLST(IPRP)
          !  ISOCMP(NSOPR) = ICMPLST(IPRP)
          !  SOPRTP(NSOPR) = 'ANTITRIP'
          end if
        end if

      end if

      ! Add some properties if DQVD is requested
      if (DQVD) then
        ! 'MLTPL  2' already there by default
        if (PRPLST(IPRP) == 'EF0    1') then
          NSOPR = NSOPR+1
          SOPRNM(NSOPR) = PRPLST(IPRP)
          ISOCMP(NSOPR) = ICMPLST(IPRP)
        end if
      end if

    end do
  end if
  ! If no PROP input, copy the SOPR selection:
  NPROP = NSOPR
  PNAME(1:NPROP) = SOPRNM(1:NPROP)
  PTYPE(1:NPROP) = SOPRTP(1:NPROP)
  ICOMP(1:NPROP) = ISOCMP(1:NPROP)

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !                                                                    C
else
  !                                                                    C
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! If no SOPR input, copy the PROP selection:
  if (NSOPR == 0) then
    NSOPR = NPROP
    SOPRNM(1:NSOPR) = PNAME(1:NSOPR)
    SOPRTP(1:NSOPR) = PTYPE(1:NSOPR)
    ISOCMP(1:NSOPR) = ICOMP(1:NSOPR)
  end if

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !                                                                    C
end if
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Lists (above) now contain either a default choice, or
! a selection by the user. Check that integrals are
! available on the oneint file.

! Check if we need to activate IFJ2/IFJZ automatically
if (natoms == 1) ifj2 = 1
xaxis = sum(abs(coor(1,1:natoms))+abs(coor(2,1:natoms)))
zaxis = sum(abs(coor(3,1:natoms)))
if ((xaxis < 1.0e-10_wp) .and. (zaxis > 0.1_wp)) ifjz = 1

! Check if angular momentum integrals have been computed
MISSAMX = 1
MISSAMY = 1
MISSAMZ = 1
do IPRP=1,NPRPLST
  if (PRPLST(IPRP) == 'ANGMOM') then
    if (ICMPLST(IPRP) == 1) MISSAMX = 0
    if (ICMPLST(IPRP) == 2) MISSAMY = 0
    if (ICMPLST(IPRP) == 3) MISSAMZ = 0
  end if
end do
if (MISSAMX+MISSAMY+MISSAMZ > 0) then
  if (IFJ2 == 1) then
    write(u6,*) ' J2 values cannot be computed.'
    write(u6,*) ' Reason: Angular momentum integrals are missing.'
    IFJ2 = 0
  end if
  if (IFGCAL .or. IFXCAL) then
    write(u6,*) ' Neither GCAL or XCAL can be computed.'
    write(u6,*) ' Reason: Angular momentum integrals are missing.'
    IFJ2 = 0
  end if
end if
if (MISSAMZ > 0) then
  if (IFJZ == 1) then
    write(u6,*) ' Omega values will not be computed.'
    write(u6,*) ' Reason: Angular momentum integrals are missing.'
    IFJZ = 0
  end if
end if

! Make sure we don't check SO properties if no SO requested
if (.not. IFSO) NSOPR = 0
! Is everything available that we may need?
IMiss = 0
IsAvail(:) = .false.
do IPROP=1,NPROP
  do IPRP=1,NPRPLST
    if ((PNAME(IPROP) == PRPLST(IPRP)) .and. (ICOMP(IPROP) == ICMPLST(IPRP))) then
      IsAvail(IPROP) = .true.
      IPUSED(IPRP) = 1
      exit
    end if
  end do

  if (IPRP > NPRPLST) then
    IMiss = IMiss+1
    if (IMiss == 1) then
      write(u6,*)
      call WarningMessage(1,'Requested integrals are missing.')
      write(u6,*) ' Property name, and component:',PNAME(IPROP),ICOMP(IPROP)
      write(u6,*) ' This record cannot be found. Some of the requested'
      write(u6,*) ' properties cannot be computed. Suggested fix: Try'
      write(u6,*) ' recomputing one-electron integrals with keyword'
      write(u6,*) " 'OneOnly', and additional keywords for the"
      write(u6,*) ' properties needed.'
    else
      write(u6,*) ' Also missing:',PNAME(IPROP),ICOMP(IPROP)
    end if
  end if
end do
IMiss = 0
IsAvailSO(:) = .false.
do ISOPR=1,NSOPR
  do IPRP=1,NPRPLST
    if ((SOPRNM(ISOPR) == PRPLST(IPRP)) .and. (ISOCMP(ISOPR) == ICMPLST(IPRP))) then
      IPUSED(IPRP) = 1
      IsAvailSO(ISOPR) = .true.
      exit
    end if
  end do

  if (IPRP > NPRPLST) then
    IMiss = IMiss+1
    if (IMiss == 1) then
      write(u6,*)
      call WarningMessage(1,'Requested integrals are missing.')
      write(u6,*) ' SO-Property name, and component:',SOPRNM(ISOPR),ISOCMP(ISOPR)
      write(u6,*) ' This record cannot be found. Some of the requested'
      write(u6,*) ' properties cannot be computed. Suggested fix: Try'
      write(u6,*) ' recomputing one-electron integrals with keyword'
      write(u6,*) " 'OneOnly', and additional keywords for the"
      write(u6,*) ' properties needed.'
    else
      write(u6,*) ' Also missing:',SOPRNM(ISOPR),ISOCMP(ISOPR)
    end if
  end if
end do
!nf
if (IfDCpl) then
  call Get_iScalar('Unique atoms',natom)
  IErr = 3*natom
  do IPrp=1,NPrpLst
    if (PrpLst(IPrp)(1:3) == 'EF1') then
      IErr = IErr-1
      IPUsed(IPrp) = 1
    end if
  end do
  if (IErr /= 0) then
    write(u6,*)
    write(u6,'(A,i5,A)') '  Approx derivative couplings require',3*natom,' field integrals.'
    write(u6,*) ' Add the keywords EFLD=0 and ONEONLY to the SEWARD input and recompute ONEINT.'
    write(u6,*) ' Program will continue, but DerCpl calculations are disabled.'
    IfDCpl = .false.
  end if
end if
!nf
! Temporarily use IPUSED to mark operators that may be needed:
do IPRP=1,NPRPLST
  do IPROP=1,NPROP
    if ((PNAME(IPROP) == PRPLST(IPRP)) .and. (ICOMP(IPROP) == ICMPLST(IPRP))) then
      IPUSED(IPRP) = 1
      exit
    end if
  end do
end do

if (.not. IFHAM) then
  if (IFSHFT) then
    call WarningMessage(1,'Ignored user request.')
    write(u6,*) ' INPCTL Warning: Both keywords ONEL and SHIFT.'
    write(u6,*) ' User-supplied diagonal energy shifts are meaningless since no Hamiltonian will be'
    write(u6,*) ' used/computed anyway. Ignored.'
    IFSHFT = .false.
  end if
  if (IFHDIA) then
    call WarningMessage(1,'Ignored user request.')
    write(u6,*) ' INPCTL Warning: Both keywords ONEL and HDIAG.'
    write(u6,*) ' User-supplied diagonal H elements are meaningless since no Hamiltonian will be'
    write(u6,*) ' used/computed anyway. Ignored.'
    IFHDIA = .false.
  end if
end if
! In any Spin-Orbit calculation, we need SO integrals as well.
! Right now, only AMFI is available. Check that first.
if (IFSO .or. IFGCAL .or. IFXCAL .or. IFMCAL) then
  IERR = 3
  do IPRP=1,NPRPLST
    if (PRPLST(IPRP) == 'AMFI') then
      IERR = IERR-1
      IPUSED(IPRP) = 1

    end if
    !if (PRPLST(IPRP) == 'PSOI') then
    !  IERR = IERR-1
    !  IPUSED(IPRP) = 1
    !end if

  end do
  if (IERR /= 0) then
    call WarningMessage(1,'Incomplete integrals.')
    write(u6,*) ' Spin-Orbit interaction calculation was requested'
    write(u6,*) ' but this requires three components of AMFI'
    write(u6,*) ' integrals. Add the keywords AMFI and ONEONLY to'
    write(u6,*) ' the SEWARD input and recompute ONEINT.'
    write(u6,*) ' Program will continue, but SO/EPRG/MAGN'
    write(u6,*) ' calculations are disabled.'
    IFSO = .false.
    NSOPR = 0
    IFGCAL = .false.
    IFXCAL = .false.
    IFMCAL = .false.
  end if
end if
! SVC2009 Check for the presence of the ANGMOM integrals needed for G factor
! or Magnetic Moment calculations.
if (IFGCAL .or. IFXCAL .or. IFMCAL) then
  IERR = 3
  do IPRP=1,NPRPLST
    if (PRPLST(IPRP) == 'ANGMOM') then
      IERR = IERR-1
      IPUSED(IPRP) = 1
    end if

    !if (PRPLST(IPRP) == 'PSOI') then
    !  write(u6,*) '5*****rassi/inpprc ANGMOM'
    !  IERR = IERR-1
    !  IPUSED(IPRP) = 1
    !end if

  end do
  if (IERR /= 0) then
    call WarningMessage(1,'Incomplete integrals.')
    write(u6,*) ' EPRG or MAGN keyword was requested'
    write(u6,*) ' but this requires three components of ANGMOM'
    write(u6,*) ' integrals. Add the keywords ANGMOM and ONEONLY to'
    write(u6,*) ' the SEWARD input and recompute ONEINT.'
    write(u6,*) ' Program will continue, but EPRG/MAGN calculations'
    write(u6,*) ' are disabled.'
    IFGCAL = .false.
    IFXCAL = .false.
    IFMCAL = .false.
  end if
end if
! Similarly, check integrals for which we want matrix elements over
! SO eigenstates.
do IPRP=1,NPRPLST
  do ISOPR=1,NSOPR
    if ((SOPRNM(ISOPR) == PRPLST(IPRP)) .and. (ISOCMP(ISOPR) == ICMPLST(IPRP))) then
      IPUSED(IPRP) = 1
      exit
    end if
  end do
end do

! Reassemble the PNAME, ICOMP arrays.

do IPRP=1,NPRPLST
  if (IPUSED(IPRP) == 0) then
    do IPROP=1,NPROP
      if (PNAME(IPROP) == PRPLST(IPRP)) PNAME(IPROP) = 'REMOVE'
    end do
  else
    IADD = 1
    do IPROP=1,NPROP
      if ((PNAME(IPROP) == PRPLST(IPRP)) .and. (ICOMP(IPROP) == ICMPLST(IPRP))) IADD = 0
    end do
    if (IADD == 1) then
      NPROP = NPROP+1
      PNAME(NPROP) = PRPLST(IPRP)
      PTYPE(NPROP) = 'UNDEF.'
      ICOMP(NPROP) = ICMPLST(IPRP)
      IsAvail(NPROP) = .true.
    end if
  end if
end do
do IPROP=1,NPROP
  if (.not. IsAvail(IPROP)) PNAME(IPROP) = 'REMOVE'
end do
MPROP = NPROP
do IPROP=1,NPROP
  if (PNAME(IPROP) == 'DONE') exit
  if (PNAME(IPROP) == 'REMOVE') then
    do N=1,MPROP-IPROP
      if (PNAME(IPROP+N) == 'DONE') exit
      if (PNAME(IPROP+N) /= 'REMOVE') exit
    end do
    PNAME(IPROP:MPROP-N) = PNAME(IPROP+N:MPROP)
    PTYPE(IPROP:MPROP-N) = PTYPE(IPROP+N:MPROP)
    ICOMP(IPROP:MPROP-N) = ICOMP(IPROP+N:MPROP)
    PNAME(MPROP-N+1:MPROP) = 'DONE'
    MPROP = MPROP-N
  end if
end do
NPROP = MPROP

! Reassemble the SOPRNM, ISOCMP arrays:

do IPRP=1,NPRPLST
  if (IPUSED(IPRP) == 0) then
    do ISOPR=1,NSOPR
      if (SOPRNM(ISOPR) == PRPLST(IPRP)) SOPRNM(ISOPR) = 'REMOVE'
    end do
  else
    IADD = 1
    do ISOPR=1,NSOPR
      if ((SOPRNM(ISOPR) == PRPLST(IPRP)) .and. (ISOCMP(ISOPR) == ICMPLST(IPRP))) IADD = 0
    end do
    if (IADD == 1) then
      NSOPR = NSOPR+1
      SOPRNM(NSOPR) = PRPLST(IPRP)
      SOPRTP(NSOPR) = 'UNDEF.'
      ISOCMP(NSOPR) = ICMPLST(IPRP)
      IsAvailSO(NSOPR) = .true.
    end if
  end if
end do
do ISOPR=1,NSOPR
  if (.not. IsAvailSO(ISOPR)) SOPRNM(ISOPR) = 'REMOVE'
end do
MSOPR = NSOPR
do ISOPR=1,NSOPR
  if (SOPRNM(ISOPR) == 'DONE') exit
  if (SOPRNM(ISOPR) == 'REMOVE') then
    do N=1,MSOPR-ISOPR
      if (SOPRNM(ISOPR+N) == 'DONE') exit
      if (SOPRNM(ISOPR+N) /= 'REMOVE') exit
    end do
    SOPRNM(ISOPR:MSOPR-N) = SOPRNM(ISOPR+N:MSOPR)
    SOPRTP(ISOPR:MSOPR-N) = SOPRTP(ISOPR+N:MSOPR)
    ISOCMP(ISOPR:MSOPR-N) = ISOCMP(ISOPR+N:MSOPR)
    SOPRNM(MSOPR-N+1:MSOPR) = 'DONE'
    MSOPR = MSOPR-N
  end if
end do
NSOPR = MSOPR
if (.not. IFSO) NSOPR = 0

! IPUSED is used later for other purposes, and should be initialized
! to zero.
IPUSED(1:NPRPLST) = 0

! PTYPE and SOPRTP is set here if not already set above. Note that this
! is a fallback procedure that should not be used actively. This
! fallback typically assigns the type of WE-reduced density to be
! used for user-specified lists of properties.

do IPROP=1,NPROP
  if (PTYPE(IPROP) /= 'UNDEF.') cycle
  select case (PNAME(IPROP))
    case ('ANGMOM','MLTPV  2','OMQ','EMFR  RA','EMFR  IA','TMOM  RA','TMOM  IA','TMOM2  I')
      PTYPE(IPROP) = 'ANTISING'
    case ('AMFI','EMFR0  I')
      PTYPE(IPROP) = 'ANTITRIP'
    case ('TMOM0  R','TMOM0  I')
      PTYPE(IPROP) = 'HERMTRIP'
    case default
      ! check partial labels
      if (PNAME(IPROP)(1:4) == 'PSOP') then
        PTYPE(IPROP) = 'ANTISING'
      else if (PNAME(IPROP)(1:3) == 'ASD') then
        PTYPE(IPROP) = 'HERMTRIP'
      else
        ! default value
        PTYPE(IPROP) = 'HERMSING'
      end if
  end select
end do

do ISOPR=1,NSOPR
  if (SOPRTP(ISOPR) /= 'UNDEF.') cycle
  select case (SOPRNM(ISOPR))
    case ('VELOCITY','ANGMOM','MLTPV  2','OMQ','EMFR  RA','EMFR  IA','TMOM  RA','TMOM  IA','TMOM2  I')
      SOPRTP(ISOPR) = 'ANTISING'
    case ('AMFI','EMFR0  I')
      SOPRTP(ISOPR) = 'ANTITRIP'
    case ('TMOM0  R','TMOM0  I')
      SOPRTP(ISOPR) = 'HERMTRIP'
    case default
      ! check partial labels
      if (SOPRNM(ISOPR)(1:4) == 'PSOP') then
        SOPRTP(ISOPR) = 'ANTISING'
      else if (SOPRNM(ISOPR)(1:3) == 'ASD') then
        SOPRTP(ISOPR) = 'HERMTRIP'
      else
        ! default value
        SOPRTP(ISOPR) = 'HERMSING'
      end if
  end select
end do

! Write out various input data:

call Get_cArray('Irreps',lIrrep,24)
lIrrep(1:nIrrep) = adjustr(lIrrep(1:nIrrep))

! determine if there are any matching wavefunctions
JOBMATCH = .false.
do job1=1,njob
  do job2=1,job1-1
    if ((MLTPLT(JOB1) == MLTPLT(JOB2)) .and. (IRREP(JOB1) == IRREP(JOB2))) JOBMATCH = .true.
  end do
end do
! make decision regarding the use of input hamiltonian/diagonal values
if (ifheff) then
  if (have_heff) then
    do J=1,NSTATE
      HAM(j,j) = HEFF(j,j)
      do I=1,J-1
        HAM(i,j) = Half*(HEFF(i,j)+HEFF(j,i))
        HAM(j,i) = HAM(i,j)
      end do
    end do
    if (jobmatch) call WarningMessage(1,'HEFF used for a situation where possible extra interaction between states is ignored!')
  else
    call WarningMessage(2,'HEFF used but none is available!')
    call Quit_OnUserError()
  end if
else if (ifejob) then
  if (have_heff) call WarningMessage(1,'EJOB used when HEFF is available, possible extra interaction between states is ignored!')
  if (have_diag) then
    do I=1,NSTATE
      HAM(i,i) = REFENE(i)
    end do
  else if (have_heff) then
    do I=1,NSTATE
      HAM(i,i) = HEFF(i,i)
    end do
  else
    call WarningMessage(2,'EJOB used but no energies available!')
    call Quit_OnUserError()
  end if
  if (jobmatch) call WarningMessage(1,'EJOB used for a situation where possible extra interaction between states is ignored!')
else if (.not. (ifhext .or. ifhdia .or. ifshft .or. ifhcom)) then
  ! the user has selected no procedure...
  if (have_heff .and. (.not. jobmatch)) then
    ifheff = .true.
    do J=1,NSTATE
      HAM(j,j) = HEff(j,j)
      do I=1,J-1
        HAM(i,j) = Half*(HEff(i,j)+HEff(j,i))
        HAM(j,i) = HAM(i,j)
      end do
    end do
  else if (have_diag) then
    ifhdia = .true.
    HDIAG(1:NSTATE) = REFENE(1:NSTATE)
  end if
end if

! Enable Hamiltonian if available/requested and not explicitly disabled
if ((.not. NOHAM) .and. (IFHEXT .or. IFHEFF .or. IFHCOM .or. IFEJOB)) IFHAM = .true.

if (IPGLOB >= 2) then
  write(u6,*)
  write(u6,*) '  The following data are common to all the states:'
  write(u6,*) '  ------------------------------------------------'
  write(u6,*) '  (note: frozen counts as inactive, deleted as secondary)'
  write(u6,*)
  write(u6,'(6X,A,I2)') 'Nr of irreps:',nIrrep
  write(u6,*)
  write(u6,'(6X,A)') '           Total     No./Irrep'
  write(u6,'(6X,A,8X,8I4)') 'Irrep       ',(I,I=1,nIrrep)
  write(u6,'(6X,A,8X,8(1X,A))') '            ',(lIrrep(I),I=1,nIrrep)
  write(u6,*)
  write(u6,'(6X,A,I4,4X,8I4)') 'INACTIVE    ',NISHT,(NISH(I),I=1,nIrrep)
  write(u6,'(6X,A,I4,4X,8I4)') 'ACTIVE      ',NASHT,(NASH(I),I=1,nIrrep)
  write(u6,'(6X,A,I4,4X,8I4)') 'SECONDARY   ',NSSHT,(NSSH(I),I=1,nIrrep)
  write(u6,'(6X,A,I4,4X,8I4)') 'BASIS       ',NBST,(NBASF(I),I=1,nIrrep)
  write(u6,*)
  write(u6,'(6X,A,I4,4X,8I4)') 'RAS1        ',NRS1T,(NRS1(I),I=1,nIrrep)
  write(u6,'(6X,A,I4,4X,8I4)') 'RAS2        ',NRS2T,(NRS2(I),I=1,nIrrep)
  write(u6,'(6X,A,I4,4X,8I4)') 'RAS3        ',NRS3T,(NRS3(I),I=1,nIrrep)
  write(u6,*)
  if (.not. (TRACK .or. ONLY_OVERLAPS)) then
    write(u6,*) '       MATRIX ELEMENTS WILL BE COMPUTED FOR THE FOLLOWING ONE-ELECTRON OPERATORS, UNLESS ZERO BY SYMMETRY.'
    write(u6,*) '  (Herm=Hermitian, Anti=Antihermitian, Sing=Singlet operator, Trip=Triplet operator)'
    do i1=1,nProp,3
      i2 = min(i1+2,nProp)
      write(u6,'(3(5X,A8,1X,I3,1X,A1,A8,A1))') (PNAME(i),ICOMP(i),'(',PTYPE(i),')',i=i1,i2)
    end do
  end if
  if (IFHAM) then
    write(u6,*)
    write(u6,*) '      EIGENSTATES OF A SPIN-FREE HAMILTONIAN WILL BE COMPUTED BASED ON:'
    write(u6,*)
    ! which kind of base hamiltonian is taken?
    if (IFHEXT) then
      write(u6,*) ' a Hamiltonian matrix that was supplied in the input.'
    else if (IFHEFF) then
      write(u6,*) ' a (effective) Hamiltonian matrix that was read from the wavefunction file(s).'
    else if (IFEJOB) then
      write(u6,*) ' a Hamiltonian matrix assumed to be diagonal with energies read from the wavefunction file(s).'
    else
      write(u6,*) ' A Hamiltonian matrix computed by RASSI.'
    end if
    ! which kind of corrections are applied?
    if (IFHDIA) &
      write(u6,*) ' In addition, the diagonal energies of the hamiltonian matrix will be replaced by either the user (HDIAG '// &
                  'keyword) or read from the wavefunction file(s).'
    if (IFSHFT) &
      write(u6,*) ' In addition, the diagonal energies of the hamiltonian matrix will be shifted by the user (SHIFT keyword).'
    if (NSOPR > 0) then
      write(u6,*)
      write(u6,*) ' SO coupling elements will be added.'
    end if
  end if
end if
if (NSOPR > 0) then
  if (.not. (TRACK .or. ONLY_OVERLAPS)) then
    write(u6,*)
    write(u6,*) '       MATRIX ELEMENTS OVER SPIN EIGENSTATES FOR:'
    do i1=1,NSOPR,3
      i2 = min(i1+2,NSOPR)
      write(u6,'(3(5X,A8,1X,I3,1X,A1,A8,A1))') (SOPRNM(i),ISOCMP(i),'(',SOPRTP(i),')',i=i1,i2)
    end do
  end if
  if (IFHAM) then
    write(u6,*)
    write(u6,*) '      EIGENSTATES OF SPIN-ORBIT HAMILTONIAN WILL BE COMPUTED'
  end if
end if
if (IPGLOB >= 4) then
  write(u6,*) 'Initial default flags are:'
  write(u6,*) '     PRSXY :',PRSXY
  write(u6,*) '     PRORB :',PRORB
  write(u6,*) '     PRTRA :',PRTRA
  write(u6,*) '     PRCI  :',PRCI
  write(u6,*) '     IFHAM :',IFHAM
  write(u6,*) '     IFHEXT:',IFHEXT
  write(u6,*) '     IFHEFF:',IFHEFF
  write(u6,*) '     IFEJOB:',IFEJOB
  write(u6,*) '     IFHCOM:',IFHCOM
  write(u6,*) '     IFSHFT:',IFSHFT
  write(u6,*) '     IFHDIA:',IFHDIA
  write(u6,*) '     IFSO  :',IFSO
  write(u6,*) '     NATO  :',NATO
  write(u6,*) '     RFPERT:',RFPERT
  write(u6,*) '     TOFILE:',ToFile
  write(u6,*) '     PRXVR :',PRXVR
  write(u6,*) '     PRXVE :',PRXVE
  write(u6,*) '     PRXVS :',PRXVS
  write(u6,*) '     PRMER :',PRMER
  write(u6,*) '     PRMEE :',PRMEE
  write(u6,*) '     PRMES :',PRMES
end if
if (IPGLOB >= 2) then
  if (NATO .and. (NRNATO > 0)) write(u6,*) ' Natural orbitals will be computed for the lowest eigenstates. NRNATO=',NRNATO
  if (BINA) then
    write(u6,*) ' Bi-natural orbitals will be computed for the following pairs of states:'
    write(u6,'(5X,8(2X,A1,I2,A1,I2,A1))') ('(',IBINA(1,I),',',IBINA(2,I),')',I=1,NBINA)
  end if
  write(u6,*)
  write(u6,*) ' Nr of states:',NSTATE
  do II=1,NSTATE,20
    III = min(II+19,NSTATE)
    write(u6,*)
    write(u6,'(1X,A8,5x,20I4)') '  State:',(I,I=II,III)
    write(u6,'(1X,A8,5x,20I4)') ' JobIph:',(JBNUM(I),I=II,III)
    write(u6,'(1X,A8,5x,20I4)') 'Root nr:',(LROOT(I),I=II,III)
  end do
  if (IFSHFT) then
    write(u6,*)
    write(u6,*) 'Each input state will be shifted with an individual'
    write(u6,*) 'amount of energy. These energy shifts are (a.u.):'
    write(u6,'(1X,5F16.8)') (ESHFT(I),I=1,NSTATE)
  end if
end if

! Added by Ungur Liviu on 04.11.2009
! Addition of NSTATE, JBNUM, and LROOT to RunFile.

call Put_iscalar('NSTATE_SINGLE',NSTATE)
call Put_iArray('JBNUM_SINGLE',JBNUM,NSTATE)
call Put_iArray('LROOT_SINGLE',LROOT,NSTATE)

! Generate the quadrature points for isotropic integration of the exponential operator

if (Do_TMOM) then
  if (Do_SK) then
    nQuad = 1
  else
    nk_Vector = 1
    nQuad = order_table(4,(L_Eff-1)/2)
  end if
else
  nk_Vector = 0
  nQuad = 0
end if

call XFLUSH(u6)

end subroutine INPPRC
