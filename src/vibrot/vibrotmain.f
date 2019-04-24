************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1986, Per-Olof Widmark                                 *
*               1986, Bjorn O. Roos                                    *
************************************************************************
************************************************************************
*                                                                      *
* WRITTEN IN 1986 BY                                                   *
* PER-OLOF WIDMARK AND BJOERN O. ROOS                                  *
* DEPARTMENT OF THEORETICAL CHEMISTRY                                  *
* UNIVERSITY OF LUND                                                   *
* SWEDEN                                                               *
*                                                                      *
************************************************************************
      subroutine vibrotmain(ireturn)
      Implicit real*8 (A-H,O-Z)
#include "dimensions.fh"
#include "intinp.fh"
#include "observ.fh"
      Dimension R(npoint+4),PotR(npoint+4)
C
      Call qEnter('Main')
C
C  Logical units
C  Vibwvs can be saved for later Transition moment calculations
C  where it will be redefined as Vibwvs1 or Vibwvs2:
      Vibwvs=12
      Vibwvs1=Vibwvs+1
      Vibwvs2=Vibwvs+2
      Call Daname(Vibwvs,'VIBWVS')
      Call Daname(Vibwvs1,'VIBWVS1')
      Call Daname(Vibwvs2,'VIBWVS2')
      ncase=0
      Call Vibinp(ncase,ngrid,nvib,Umin,Umax,R,PotR,E0,dE0,Redm,
     *            Teas,Req,scale,temp)
      If ( ncase.eq.1 ) then
        Call Vibrot(ngrid,nvib,Umin,Umax,R,PotR,E0,dE0,
     *              Redm,Req,scale,temp)
        If(IfPrWf.gt.0) then
          Call PrWf_VibRot(ngrid,R)
        End If
      Else If ( ncase.eq.2 ) then
        Write(6,*) iallrot
        Do 10 i=1,iobs
          Call Vibtrm(ngrid,Umin,Umax,Teas,R,EoutO(1,i),Titobs(i))
10      Continue
      Else
        Write(6,1000)
1000    Format(/1x,'No ROVIbrational or TRANsition keywords'//
     & 'specified in input.')
      End If
      Call Daclos(Vibwvs)
      Call Daclos(Vibwvs1)
      Call Daclos(Vibwvs2)
      Call qExit('Main')
      ireturn=0
      return
      End
