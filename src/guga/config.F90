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
!***********************************************************************

subroutine CONFIG(NREF,IOCR,nIOCR,L0,L1,L2,L3,JSYM,JSY,INTNUM,LSYM,JJS,ISO,LV,IFCORE,ICOR,NONE_,JONE,JREFX,NFREF)

use Definitions, only: iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NREF, nIOCR, IOCR(nIOCR), L0(*), L1(*), L2(*), L3(*), INTNUM, LV, IFCORE, ICOR(*), NONE_, JONE(*)
integer(kind=iwp), intent(_OUT_) :: JSYM(*), JSY(*), JJS(*), ISO(*), JREFX(*)
integer(kind=iwp), intent(out) :: LSYM, NFREF
#include "real_guga.fh"
#include "integ.fh"
#include "d.fh"
integer(kind=iwp) :: I, I1, IBS, IDIF, IEL, IFEXC, IFREF, IIJ, IJJ, IND, INHOLE, INTOT, IOC(55), IPART, IRC1, IRC2, IREF, IRR, &
                     ISP(55), ISTA, ITEMP, ITU, IX1, IX2, IX3, IX4, J, JHOLE, JJ1, JND, JPART, JRC1, JRC2, JRC21, JRX, JSYL, &
                     JSYLL, K, KM, KM1, L, LMN, LMN0, LNS, M, M1, M2, MND, NCORR, NSJ, NSJ1

JSYL = 30000
JSYLL = 3000
JRX = 9000
IBS = 0
IEL = 2
LSYM = 1
NFREF = 0
LNS = NIORB+LV+1
! Compute wave function symmetry
IRR = 0
do I=LNS,LN
  IRR = IRR+1
  if (IOCR(IRR) == 1) LSYM = MUL(LSYM,NSM(I))
end do
write(IW,'(6X,A,I3)') 'WAVE-FUNCTION SYMMETRY LABEL:',LSYM

! Initialize arrays JJS, JNDX, ICASE, JSY
do I=1,18
  JJS(I) = 0
end do
! CONSTRUCT JNDX
! ILIM=2 or ILIM=4 was set by input: Normally 4, but 2 if keyword FIRST
! has been given. ILIM is in integ.fh
INTOT = IRC(ILIM)
do I=1,INTOT
  JNDX(I) = 0
end do
do J=1,MXCASE
  ICASE(J) = 0
end do
do I=1,JSYLL
  JSY(I) = 0
end do
JND = 0
LMN = 0

do IIJ=1,ILIM
  ! Special cases:
  JRC(IIJ) = LMN
  if ((N == 0) .and. (IIJ > 1)) goto 10
  if ((N == 1) .and. (IIJ > 2)) goto 10
  if ((N == 2) .and. (ISPIN == 1) .and. (IIJ == 3)) goto 10
  if ((N == 2) .and. (ISPIN == 3) .and. (IIJ == 4)) goto 10
  ! ISTA is actually an offset, not a start point.
  ISTA = (IIJ-1)*MXVERT
  ! Element in row IJJ=IV0+1-IIJ of the tables is the top vertex, for
  ! IIJ=1 (Valence), 2 (Singles), 3 (Doubles T), 4 (Doubles S).
  ! It is found as element ISTA+IJJ of arrays L0, L1 etc.
  IJJ = IV0+1-IIJ
  KM = 1
  J2(KM) = IJJ
11 KM = KM+1
  IWAY(KM) = 0
12 KM1 = KM-1
  ! At KM, trying for a way down.
  if ((L0(ISTA+J2(KM1)) == 0) .or. (IWAY(KM) >= 1)) GO TO 14
  ! IWAY is 0, and the next vertex L0(ISTA+J2(KM1)) is actually there:
  J2(KM) = L0(ISTA+J2(KM1))
  IWAY(KM) = 1
  IOC(KM1) = 0
  ISP(KM1) = 0
  GO TO 20

14 if ((L1(ISTA+J2(KM1)) == 0) .or. (IWAY(KM) >= 2)) GO TO 15
  ! IWAY is 1, and the next vertex is actually there:
  J2(KM) = L1(ISTA+J2(KM1))
  IWAY(KM) = 2
  IOC(KM1) = 1
  ISP(KM1) = 1
  GO TO 20

15 if ((L2(ISTA+J2(KM1)) == 0) .or. (IWAY(KM) >= 3)) GO TO 16
  ! IWAY is 2, and the next vertex is actually there:
  J2(KM) = L2(ISTA+J2(KM1))
  IWAY(KM) = 3
  IOC(KM1) = 1
  ISP(KM1) = 2
  GO TO 20

16 if ((L3(ISTA+J2(KM1)) == 0) .or. (IWAY(KM) >= 4)) GO TO 17
  ! IWAY is 3, and the next vertex is actually there:
  J2(KM) = L3(ISTA+J2(KM1))
  IWAY(KM) = 4
  IOC(KM1) = 2
  ISP(KM1) = 3
  GO TO 20

17 KM = KM-1
  ! No more ways to try from this vertex. Up to higher vertex:
  if (KM == 1) GO TO 10
  ! If KM was 1, then finish this IIJ value and do next one.
  ! Else, goto 12, to take new route from the higher vertex.
  GO TO 12

20 continue
  ! While trying to find a new wav through the graph, we found
  ! a feasible edge from a vertex J2(KM) to a vertex at level KM1.
  if (KM1 == NIORB+LV) IBS = IB(J2(KM))
  if (KM /= LN+1) GO TO 11
  ! KM has reached top of the graph.

  JND = JND+1
  ! A formally legal way was found. Should it be accepted?
  NSJ = 1
  INHOLE = 0
  do I=1,LN
    if (IOC(I) == 1) NSJ = MUL(NSJ,NSM(I))
    if ((I <= NIORB+LV) .and. (I > LV)) INHOLE = INHOLE+2-IOC(I)
  end do
  ! STRIKE OUT INTERNAL CONFIGURATIONS
  IPART = 0
  if (JND > IRC(1)) IPART = IPART+1
  if (JND > IRC(2)) IPART = IPART+1
  if ((IPART == 0) .and. (NSJ /= LSYM)) GO TO 12
  ! TEST IF TOO HIGHLY EXCITED
  ! TEST ALSO IF REFERENCE STATE
  IFEXC = 0
  IFREF = 0
  JJ1 = 0
  do IREF=1,NREF
    JHOLE = 0
    JPART = IPART
    do I=1,LN
      if (I > LV) GO TO 250
      IDIF = IOC(I)
      GO TO 251
250   if (I > NIORB+LV) GO TO 252
      IDIF = IOC(I)-2
      GO TO 251
252   JJ1 = JJ1+1
      if (IOC(I) == IOCR(JJ1)) GO TO 112
      IDIF = IOC(I)-IOCR(JJ1)
251   if (IDIF > 0) GO TO 114
      JHOLE = JHOLE-IDIF
      GO TO 112
114   JPART = JPART+IDIF
112   continue
    end do
    if (JPART /= JHOLE) then
      write(u6,*) 'Config: JPART /= JHOLE'
      write(u6,*) 'JPART,JHOLE=',JPART,JHOLE
      call Abend()
    end if
    if (JPART <= IEL) IFEXC = 1
    if (JPART == 0) IFREF = 1
  end do
  if (IFEXC == 0) GO TO 12
  if ((IPART /= 2) .or. (INTNUM == 0)) GO TO 115
  ! INTERACTING SPACE
  if ((INHOLE == 2) .and. (IBS /= 0)) GO TO 12
  ! NO CORE-CORE CORRELATION
115 if (IFCORE == 0) GO TO 116
  NCORR = 0
  do I=1,LN
    if (ICOR(I) == 0) GO TO 117
    NCORR = NCORR+2-IOC(I)
117 continue
  end do
  if (NCORR > 1) GO TO 12
  ! SINGLY OCCUPIED ORBITALS
116 if (NONE_ == 0) GO TO 118
  do I=1,NONE_
    if (IOC(JONE(I)) /= 1) GO TO 12
  end do
118 LMN = LMN+1
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
    if (IFREF == 0) GO TO 210
    NFREF = NFREF+1
    JREFX(LMN) = NFREF
  end if

210 continue
  JRC(IIJ) = LMN
  if (LMN > JSYL) GO TO 985
  JSYM(LMN) = NSJ
  if (IIJ <= 2) GO TO 985
  NSJ1 = NSJ+1
  if (IIJ == 3) JJS(NSJ1) = JJS(NSJ1)+1
  if (IIJ == 4) JJS(NSJ1+9) = JJS(NSJ1+9)+1

985 continue
  M = (LMN-1)*LN
  if (M+LN > ISPA) then
    write(u6,*) 'Config: M+LN > ISPA'
    write(u6,*) 'M,LN,ISPA=',M,LN,ISPA
    write(u6,*) 'This error may be due to a bug.'
    write(u6,*) ' Please check your input against the manual.'
    write(u6,*) ' If there is no input errors, please report'
    write(u6,*) ' this as a bug.'
    call Abend()
  end if
  do K=1,LN
    ISO(M+K) = ISP(K)
  end do
  GO TO 12
10 continue
end do

JRC(ILIM) = LMN
write(IW,214)
214 format(//,6X,'INTERNAL CONFIGURATIONS (REAL)')
IX1 = JRC(1)
IX2 = JRC(2)-JRC(1)
if (IFIRST /= 0) GO TO 205

!if (N == 2) JRC(3) = JRC(2)

IX3 = JRC(3)-JRC(2)
IX4 = JRC(4)-JRC(3)
write(IW,215) IX1,IX2,IX3,IX4
215 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7,/,6X, &
           'NUMBER OF TRIPLET COUPLED DOUBLES',I7,/,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
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
if (NSYM == 1) GO TO 410
ITU = 2
400 IRC1 = IRC(ITU)+1
IRC2 = IRC(ITU+1)
JRC1 = JRC(ITU)+1
JRC2 = JRC(ITU+1)
JRC21 = JRC2-1
if (JRC21 < JRC1) GO TO 401
do I=JRC1,JRC21
  I1 = I+1
  do J=I1,JRC2
    if (JSYM(J) >= JSYM(I)) GO TO 403
    ITEMP = JSYM(I)
    JSYM(I) = JSYM(J)
    JSYM(J) = ITEMP
    M1 = (I-1)*LN
    M2 = (J-1)*LN
    do K=1,LN
      ITEMP = ISO(M1+K)
      ISO(M1+K) = ISO(M2+K)
      ISO(M2+K) = ITEMP
    end do
    do K=IRC1,IRC2
      if (JNDX(K) /= I) GO TO 421
      JNDX(K) = J
      GO TO 420
421   if (JNDX(K) /= J) GO TO 420
      JNDX(K) = I
420   continue
    end do
403 continue
  end do
end do
401 if (ITU == 3) GO TO 406
ITU = 3
GO TO 400

406 continue
! JJS(2..9): JJS(I+1)=Nr of internal triplet states per symmetry I.
! JJS(11..18): JJS(I+10)=Nr of internal singlet states per symmetry I.
write(IW,407) (JJS(I+1),I=1,NSYM)
407 format(/6X,'INTERNAL TRIPLET STATES PER SYMMETRY:',6X,8I5)
write(IW,408) (JJS(I+10),I=1,NSYM)
408 format(6X,'INTERNAL SINGLET STATES PER SYMMETRY:',6X,8I5)
do I=2,NSYM
  I1 = I+1
  JJS(I1) = JJS(I)+JJS(I1)
  JJS(I1+9) = JJS(I+9)+JJS(I1+9)
end do
! Now, JJS is changed to contain instead the corresponding cumulative sum
! summed over the symmetries.
GO TO 410

205 write(IW,216) IX1,IX2
216 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
if ((IX1 >= 8192) .or. (IX2 >= 8192)) then
  write(u6,*) 'Config: IX? >= 8192'
  write(u6,*) 'IX1,IX2=',IX1,IX2
  call Abend()
end if

! PACK OCCUPATION AND SYMMETRY VECTORS FOR CI
410 LMN0 = JRC(ILIM)*LN
!PAM97 M1 = (LMN0+29)/30
M1 = (LMN0+14)/15
if (M1 > MXCASE) then
  write(u6,*) 'Config: M1 > MXCASE'
  write(u6,*) 'M1,MXCASE=',M1,MXCASE
  call Abend()
end if
M = 0
do L=1,LMN
  do K=1,LN
    M = M+1
    MND = ISO(M)
    !IOCC((M+29)/30) = OR(IOCC((M+29)/30),1SHIFT(MND,2*((M+29)/30*30-M)))
    !PAM97 QOCC((M+29)/30) = PACK(QOCC((M+29)/30),MND,2*M-(2*M-1)/60*60,2)
    call ICPCK(ICASE,M,MND)
  end do
  !PAM97 NSJ = JSYM(L)-1
  !PAM97 JSY((L+9)/10) = ior(JSY((L+9)/10),ishft(NSJ,29-3*mod(L-1,10)))
  call JSPCK(JSY,L,JSYM(L))
end do

return

end subroutine CONFIG
