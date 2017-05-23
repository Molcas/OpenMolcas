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
      Subroutine ChkVec(OrbFileName,iVer,NSYM_L,NBAS_L,NORB_L,
     &                    InfoLbl,iRc)
      Character*(*) OrbFileName
      Character*8 InfoLbl
      Character*80 line
      Dimension NBAS_L(8)
      Dimension NORB_L(8)
#include "warnings.fh"
#include "inporbfmt.fh"
      Logical lExists

* Purpose: To check if OrbFileName is a valid orbital file, and which
* information it contains.
      Call f_Inquire(OrbFileName,lExists)
      If (.not.lExists) Go To 921
      LU=99
      LU=IsFreeUnit(LU)
      Call Molcas_Open(LU, OrbFileName)
* Check version!
      READ(LU,'(A80)',ERR=920,END=920) Line
      iVer=0
      Do jVer=1,mxVer
        if(Magic(jVer).eq.Line(1:len(Magic(jVer)))) iVer=jVer
      End Do
      iRc=_RC_IO_ERROR_READ_
      IF(iVer.eq.0) Return

* Find and read information section:
  10  CONTINUE
      READ(LU,'(A80)',ERR=920) LINE
      IF(LINE.NE.'#INFO') GOTO 10
  11  CONTINUE
      READ(LU,'(A80)',ERR=920) LINE
      IF(LINE(1:1).EQ.'*') GOTO 11
      READ(LINE,*) IUHF,NSYM_L
  12  CONTINUE
      READ(LU,'(A80)',ERR=920) LINE
      IF(LINE(1:1).EQ.'*') GOTO 12
      READ(LINE,*) (NBAS_L(I),I=1,NSYM_L)
  13  CONTINUE
      READ(LU,'(A80)',ERR=920) LINE
      IF(LINE(1:1).EQ.'*') GOTO 13
      READ(LINE,*) (NORB_L(I),I=1,NSYM_L)

* Create a InfoLbl telling what data sets are available:
      iORB=0
      iOCC=0
      iONE=0
      iIND=0
* Find section InfoLbls:
      If (IUHF.eq.0) Then
  21    CONTINUE
        READ(LU,'(A80)',END=900,ERR=900) LINE
        IF(LINE(1:4).EQ.'#ORB') iORB=1
        IF(LINE(1:4).EQ.'#OCC') iOCC=1
        IF(LINE(1:4).EQ.'#ONE') iONE=1
        IF(LINE(1:4).EQ.'#IND') iIND=1
        GOTO 21
      Else
  22    CONTINUE
        READ(LU,'(A80)',END=900,ERR=900) LINE
        IF(LINE(1:5).EQ.'#UORB') iORB=1
        IF(LINE(1:5).EQ.'#UOCC') iOCC=1
        IF(LINE(1:5).EQ.'#UONE') iONE=1
        IF(LINE(1:4).EQ.'#IND') iIND=1
        GOTO 22
      End If

* Intentionally reached end of file:
 900  CONTINUE
      InfoLbl=' '
      ip=0
      ip=ip+iorb
      if(iorb.eq.1) InfoLbl(ip:ip)='C'
      ip=ip+iocc
      if(iocc.eq.1) InfoLbl(ip:ip)='O'
      ip=ip+ione
      if(ione.eq.1) InfoLbl(ip:ip)='E'
      ip=ip+iind
      if(iind.eq.1) InfoLbl(ip:ip)='I'

* -----------------------------------------------------------------
      CLOSE(LU)
      iRC=_RC_ALL_IS_WELL_
      RETURN
* -----------------------------------------------------------------
 920  CONTINUE
      CLOSE(LU)
 921  CONTINUE
      iRC=_RC_IO_ERROR_READ_
      RETURN
* -----------------------------------------------------------------
      End
