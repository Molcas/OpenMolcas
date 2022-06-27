!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
! 2021: Remove GOTOs

subroutine CONFIG(NREF,IOCR,nIOCR,L0,L1,L2,L3,JSY,INTNUM,LSYM,JJS,LV,IFCORE,ICOR,NONE_,JONE,JREFX,NFREF)

use guga_global, only: IB, ICASE, IFIRST, ILIM, IRC, ISPIN, IV0, IWAY, J2, JNDX, JRC, LN, MXVERT, N, NIORB, NSM, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NREF, nIOCR, IOCR(nIOCR), L0(*), L1(*), L2(*), L3(*), INTNUM, LV, IFCORE, ICOR(*), NONE_, JONE(*)
integer(kind=iwp), intent(out) :: JSY(3000), LSYM, JJS(18), JREFX(9000), NFREF
integer(kind=iwp) :: I, I1, IBS, IDIF, IEL, IFEXC, IFREF, IIJ, IJJ, IND, INHOLE, IOC(55), IPART, IRC1, IRC2, IREF, IRR, ISP(55), &
                     ISTA, ITEMP, ITU, IX1, IX2, IX3, IX4, J, JHOLE, JJ1, JND, JPART, JRC1, JRC2, JRC21, JRX, JSYL, K, KM, KM1, L, &
                     LMN, LMN0, LNS, M, MND, NCORR, NSJ, NSJ1
logical(kind=iwp) :: first
integer(kind=iwp), allocatable :: ISO(:,:), JSYM(:), TEMP(:)

JSYL = 30000
JRX = size(JREFX)
IBS = 0
IEL = 2
LSYM = 1
NFREF = 0
LNS = NIORB+LV+1
! Compute wave function symmetry
IRR = 0
do I=LNS,LN
  IRR = IRR+1
  if (IOCR(IRR) == 1) LSYM = Mul(LSYM,NSM(I))
end do
write(u6,'(6X,A,I3)') 'WAVE-FUNCTION SYMMETRY LABEL:',LSYM

! Initialize arrays JJS, JNDX, JSY
JJS(:) = 0
! CONSTRUCT JNDX
! ILIM=2 or ILIM=4 was set by input: Normally 4, but 2 if keyword FIRST
! has been given.
JNDX(:) = 0
JSY(:) = 0
JND = 0
LMN = 0

call mma_allocate(JSYM,JSYL,label='JSYM')
call mma_allocate(ISO,LN,JRX,label='ISO')
call mma_allocate(TEMP,LN,label='TEMP')
do IIJ=1,ILIM
  ! Special cases:
  JRC(IIJ) = LMN
  if ((N == 0) .and. (IIJ > 1)) cycle
  if ((N == 1) .and. (IIJ > 2)) cycle
  if ((N == 2) .and. (ISPIN == 1) .and. (IIJ == 3)) cycle
  if ((N == 2) .and. (ISPIN == 3) .and. (IIJ == 4)) cycle
  ! ISTA is actually an offset, not a start point.
  ISTA = (IIJ-1)*MXVERT
  ! Element in row IJJ=IV0+1-IIJ of the tables is the top vertex, for
  ! IIJ=1 (Valence), 2 (Singles), 3 (Doubles T), 4 (Doubles S).
  ! It is found as element ISTA+IJJ of arrays L0, L1 etc.
  IJJ = IV0+1-IIJ
  KM = 1
  J2(KM) = IJJ
  first = .true.
  do
    if (first) then
      KM = KM+1
      IWAY(KM) = 0
      first = .false.
    end if
    KM1 = KM-1
    ! At KM, trying for a way down.
    if ((L0(ISTA+J2(KM1)) /= 0) .and. (IWAY(KM) < 1)) then
      ! IWAY is 0, and the next vertex L0(ISTA+J2(KM1)) is actually there:
      J2(KM) = L0(ISTA+J2(KM1))
      IWAY(KM) = 1
      IOC(KM1) = 0
      ISP(KM1) = 0
    else if ((L1(ISTA+J2(KM1)) /= 0) .and. (IWAY(KM) < 2)) then
      ! IWAY is 1, and the next vertex is actually there:
      J2(KM) = L1(ISTA+J2(KM1))
      IWAY(KM) = 2
      IOC(KM1) = 1
      ISP(KM1) = 1
    else if ((L2(ISTA+J2(KM1)) /= 0) .and. (IWAY(KM) < 3)) then
      ! IWAY is 2, and the next vertex is actually there:
      J2(KM) = L2(ISTA+J2(KM1))
      IWAY(KM) = 3
      IOC(KM1) = 1
      ISP(KM1) = 2
    else if ((L3(ISTA+J2(KM1)) /= 0) .and. (IWAY(KM) < 4)) then
      ! IWAY is 3, and the next vertex is actually there:
      J2(KM) = L3(ISTA+J2(KM1))
      IWAY(KM) = 4
      IOC(KM1) = 2
      ISP(KM1) = 3
    else
      KM = KM-1
      ! No more ways to try from this vertex. Up to higher vertex:
      if (KM == 1) exit
      ! If KM was 1, then finish this IIJ value and do next one.
      ! Else, take new route from the higher vertex.
      cycle
    end if

    ! While trying to find a new wav through the graph, we found
    ! a feasible edge from a vertex J2(KM) to a vertex at level KM1.
    if (KM1 == NIORB+LV) IBS = IB(J2(KM))
    if (KM /= LN+1) then
      first = .true.
      cycle
    end if
    ! KM has reached top of the graph.

    JND = JND+1
    ! A formally legal way was found. Should it be accepted?
    NSJ = 1
    INHOLE = 0
    do I=1,LN
      if (IOC(I) == 1) NSJ = Mul(NSJ,NSM(I))
      if ((I <= NIORB+LV) .and. (I > LV)) INHOLE = INHOLE+2-IOC(I)
    end do
    ! STRIKE OUT INTERNAL CONFIGURATIONS
    IPART = 0
    if (JND > IRC(1)) IPART = IPART+1
    if (JND > IRC(2)) IPART = IPART+1
    if ((IPART == 0) .and. (NSJ /= LSYM)) cycle
    ! TEST IF TOO HIGHLY EXCITED
    ! TEST ALSO IF REFERENCE STATE
    IFEXC = 0
    IFREF = 0
    JJ1 = 0
    do IREF=1,NREF
      JHOLE = 0
      JPART = IPART
      do I=1,LN
        if (I <= LV) then
          IDIF = IOC(I)
        else if (I <= NIORB+LV) then
          IDIF = IOC(I)-2
        else
          JJ1 = JJ1+1
          if (IOC(I) == IOCR(JJ1)) cycle
          IDIF = IOC(I)-IOCR(JJ1)
        end if
        if (IDIF <= 0) then
          JHOLE = JHOLE-IDIF
        else
          JPART = JPART+IDIF
        end if
      end do
      if (JPART /= JHOLE) then
        write(u6,*) 'Config: JPART /= JHOLE'
        write(u6,*) 'JPART,JHOLE=',JPART,JHOLE
        call Abend()
      end if
      if (JPART <= IEL) IFEXC = 1
      if (JPART == 0) IFREF = 1
    end do
    if (IFEXC == 0) cycle
    if ((IPART == 2) .and. (INTNUM /= 0)) then
      ! INTERACTING SPACE
      if ((INHOLE == 2) .and. (IBS /= 0)) cycle
    end if
    ! NO CORE-CORE CORRELATION
    if (IFCORE /= 0) then
      NCORR = 0
      do I=1,LN
        if (ICOR(I) == 0) cycle
        NCORR = NCORR+2-IOC(I)
      end do
      if (NCORR > 1) cycle
    end if
    ! SINGLY OCCUPIED ORBITALS
    do I=1,NONE_
      if (IOC(JONE(I)) /= 1) cycle
    end do
    LMN = LMN+1
    L = JND
    IND = LMN
    JNDX(L) = IND
    if (IIJ == 1) then
      ! CONSTRUCT INDEX LIST FOR REFERENCE STATES
      if (LMN > JRX) then
        write(u6,*) 'Config: LMN > JRX'
        write(u6,*) 'LMN,JRX=',LMN,JRX
        write(u6,*) 'This error is almost certainly that the problem'
        write(u6,*) ' at hand requires larger arrays than GUGA is'
        write(u6,*) ' compiled for. Please check your input against'
        write(u6,*) ' the manual. If you are certain that this'
        write(u6,*) ' calculation should be possible, please report'
        write(u6,*) ' this as a bug.'
        call Abend()
      end if
      JREFX(LMN) = 0
      if (IFREF /= 0) then
        NFREF = NFREF+1
        JREFX(LMN) = NFREF
      end if
    end if

    JRC(IIJ) = LMN
    if (LMN <= JSYL) then
      JSYM(LMN) = NSJ
      if (IIJ > 2) then
        NSJ1 = NSJ+1
        if (IIJ == 3) JJS(NSJ1) = JJS(NSJ1)+1
        if (IIJ == 4) JJS(NSJ1+9) = JJS(NSJ1+9)+1
      end if
    end if

    if (LMN > JRX) then
      write(u6,*) 'Config: LMN > JRX'
      write(u6,*) 'LMN,JRX=',LMN,JRX
      write(u6,*) 'This error may be due to a bug.'
      write(u6,*) ' Please check your input against the manual.'
      write(u6,*) ' If there is no input errors, please report'
      write(u6,*) ' this as a bug.'
      call Abend()
    end if
    ISO(:,LMN) = ISP(1:LN)
  end do
end do

JRC(ILIM) = LMN
write(u6,214)
IX1 = JRC(1)
IX2 = JRC(2)-JRC(1)
if (IFIRST == 0) then

  !if (N == 2) JRC(3) = JRC(2)

  IX3 = JRC(3)-JRC(2)
  IX4 = JRC(4)-JRC(3)
  write(u6,215) IX1,IX2,IX3,IX4
  if (LMN > JSYL) then
    write(u6,*) 'Config: LMN > JSYL'
    write(u6,*) 'LMN,JSYL=',LMN,JSYL
    call Abend()
  end if
  if ((IX1 >= 8192) .or. (IX2 >= 8192) .or. (IX3 >= 8192) .or. (IX4 >= 8192)) then
    write(u6,*) 'Config: IX? >= 8192'
    write(u6,*) 'IX1,IX2,IX3,IX4=',IX1,IX2,IX3,IX4
    write(u6,*) 'This error is almost certainly that the problem'
    write(u6,*) ' at hand requires larger arrays than GUGA is'
    write(u6,*) ' compiled for. Please check your input against'
    write(u6,*) ' the manual. If you are certain that this'
    write(u6,*) ' calculation should be possible, please report'
    write(u6,*) ' this as a bug.'
    call Abend()
  end if
  ! SORT BY SYMMETRY
  if (NSYM /= 1) then
    ITU = 2
    do
      IRC1 = IRC(ITU)+1
      IRC2 = IRC(ITU+1)
      JRC1 = JRC(ITU)+1
      JRC2 = JRC(ITU+1)
      JRC21 = JRC2-1
      do I=JRC1,JRC21
        I1 = I+1
        do J=I1,JRC2
          if (JSYM(J) >= JSYM(I)) cycle
          ITEMP = JSYM(I)
          JSYM(I) = JSYM(J)
          JSYM(J) = ITEMP
          TEMP(:) = ISO(:,I)
          ISO(:,I) = ISO(:,J)
          ISO(:,J) = TEMP
          do K=IRC1,IRC2
            if (JNDX(K) == I) then
              JNDX(K) = J
            else if (JNDX(K) == J) then
              JNDX(K) = I
            end if
          end do
        end do
      end do
      if (ITU == 3) exit
      ITU = 3
    end do

    ! JJS(2..9): JJS(I+1)=Nr of internal triplet states per symmetry I.
    ! JJS(11..18): JJS(I+10)=Nr of internal singlet states per symmetry I.
    write(u6,407) (JJS(I+1),I=1,NSYM)
    write(u6,408) (JJS(I+10),I=1,NSYM)
    do I=2,NSYM
      I1 = I+1
      JJS(I1) = JJS(I)+JJS(I1)
      JJS(I1+9) = JJS(I+9)+JJS(I1+9)
    end do
    ! Now, JJS is changed to contain instead the corresponding cumulative sum
    ! summed over the symmetries.
  end if
else
  write(u6,216) IX1,IX2
  if ((IX1 >= 8192) .or. (IX2 >= 8192)) then
    write(u6,*) 'Config: IX? >= 8192'
    write(u6,*) 'IX1,IX2=',IX1,IX2
    call Abend()
  end if
end if

! PACK OCCUPATION AND SYMMETRY VECTORS FOR CI
LMN0 = JRC(ILIM)*LN
!PAM97 M = (LMN0+29)/30
M = (LMN0+14)/15
call mma_allocate(ICASE,M,label='ICASE')
ICASE(:) = 0
M = 0
do L=1,LMN
  do K=1,LN
    M = M+1
    MND = ISO(K,L)
    !IOCC((M+29)/30) = OR(IOCC((M+29)/30),1SHIFT(MND,2*((M+29)/30*30-M)))
    !PAM97 QOCC((M+29)/30) = PACK(QOCC((M+29)/30),MND,2*M-(2*M-1)/60*60,2)
    call ICPCK(ICASE,M,MND)
  end do
  !PAM97 NSJ = JSYM(L)-1
  !PAM97 JSY((L+9)/10) = ior(JSY((L+9)/10),ishft(NSJ,29-3*mod(L-1,10)))
  call JSPCK(JSY,L,JSYM(L))
end do
call mma_deallocate(JSYM)
call mma_deallocate(ISO)
call mma_deallocate(TEMP)

return

214 format(//,6X,'INTERNAL CONFIGURATIONS (REAL)')
215 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7,/,6X, &
           'NUMBER OF TRIPLET COUPLED DOUBLES',I7,/,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
216 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
407 format(/6X,'INTERNAL TRIPLET STATES PER SYMMETRY:',6X,8I5)
408 format(6X,'INTERNAL SINGLET STATES PER SYMMETRY:',6X,8I5)

end subroutine CONFIG
