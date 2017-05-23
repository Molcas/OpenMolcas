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
      Subroutine SingP(nCalls,iQ_Atoms,ipStoreCoo,nPart2)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "WrkSpc.fh"

*
*-- If this is first call, issue a warning.
*
      If(nCalls.eq.0) then
        Write(6,*)
        Write(6,*)
        Write(6,*)'---->>>  WARNING  <<<----'
        Write(6,*)
        Write(6,*)'You have specified that a set of single-point '
     &//'calculations are to be preformed.'
        Write(6,*)'This means that the input will be given to some '
     &//'extent a new meaning.'
        Write(6,*)

*
*-- Put coordinates in a new vector if first call.
*
        kaunter=0
        nPart2=nPart
        Call GetMem('Store','Allo','Real',ipStoreCoo,nPart2*nCent*3)
        Do 11, iPart=1,nPart2
          Do 12, iCent=1,nCent
            kaunter=kaunter+1
            Work(ipStoreCoo+3*(kaunter-1))=Cordst(kaunter,1)
            Work(ipStoreCoo+3*(kaunter-1)+1)=Cordst(kaunter,2)
            Work(ipStoreCoo+3*(kaunter-1)+2)=Cordst(kaunter,3)
12        Continue
11      Continue

*
*-- Put dummies that will be substituted for the qm-region.
*
        nAllQm=(((iQ_Atoms-1)/nAtom)+1)*nCent
        Do 16, i=1,nAllQm
          Do 17, j=1,3
            Cordst(i,j)=0.0d0
17        Continue
16      Continue

*
*-- Put the coordinates of first iteration.
*
        Do 18, iCent=1,nCent
          Cordst(nAllQm+iCent,1)=Work(ipStoreCoo+3*(iCent-1))
          Cordst(nAllQm+iCent,2)=Work(ipStoreCoo+3*(iCent-1)+1)
          Cordst(nAllQm+iCent,3)=Work(ipStoreCoo+3*(iCent-1)+2)
18      Continue

*
*-- Set new value on some variables.
*
        nMicro=1
        nMacro=1
        DelX=0
        DelFi=0
        DelR=0
        QmEq=.true.
        nPart=(nAllQm/nCent)+1
        Write(6,*)
        Write(6,*)'Resetting for FIT:'
        Write(6,*)'Number of macrosteps:',nMacro
        Write(6,*)'Number of microsteps:',nMicro
        Write(6,*)'No translation, rotation or radie modification.'
        Write(6,*)'Take the QmEq path.'

*
*-- If not first call, then collect relevant coordinates.
*
      Else
        Initial1=(((iQ_Atoms-1)/nAtom)+1)*nCent
        Initial2=3*nCent*nCalls-1
        Do 21, iCent=1,nCent
        Cordst(Initial1+iCent,1)=Work(ipStoreCoo+Initial2+(iCent-1)*3+1)
        Cordst(Initial1+iCent,2)=Work(ipStoreCoo+Initial2+(iCent-1)*3+2)
        Cordst(Initial1+iCent,3)=Work(ipStoreCoo+Initial2+(iCent-1)*3+3)
21      Continue
      Endif

*
*-- Up-date nCalls.
*
      nCalls=nCalls+1


      Return
      End
