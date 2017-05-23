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
      Subroutine GF_Print(EVal,EVec,dDipM,iel,nX,nDim,ictl,IRInt,Lu_10,
     &                    iOff)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "constants2.fh"
#include "WrkSpc.fh"
      Real*8 EVal(nDim), EVec(2,nX,nDim),dDipM(ndim,iel),IRInt(nDim)
      Parameter(Inc=6)
      Real*8 Tmp(Inc)
      Character*80 Format, Line*120
      Character*(LENIN6) ChDisp(3*MxAtom),Label
*
      LUt=6
*
      Call Get_iScalar('nChDisp',nChDisp)
      Call Get_cArray('ChDisp',ChDisp,(LENIN6)*nChDisp)
*
      iIRInt=0
      Do iHarm = 1, nDim, Inc
         Jnc=Min(Inc,nDim-iHarm+1)
         Label=' '
         Write(Format,'(A,I3,A)') '(5X,A,1x,',Jnc,'(I7,3X))'
         Write (LUt,Format) Label,(i,i=iHarm,iHarm+Jnc-1)
         Write (LUt,*)
*
         Label='Frequency:'
         Write(Format,'(A,I3,A)') '(5X,A,1x,',Jnc,'F10.2)'
         Line=' '
         Write (Line,Format) Label,(EVal(i),i=iHarm,iHarm+Jnc-1)
*
*------- Replace minus signs with sign for imaginary unit.
*
         Do i = 1, 120
            If (Line(i:i).eq.'-') Line(i:i)='i'
         End Do
         Write (LUt,'(A)') Line
         Write (LUt,*)
*
         If (ictl.ne.0) Then
            Label='Intensity:'
            Write(Format,'(A,I3,A)') '(5X,A,1x,',Jnc,'E10.3)'
            call dcopy_(Jnc,0.0d0,0,Tmp,1)
            Do k=1,Jnc
              Do l=1,iel
               Tmp(k)=tmp(k)+dDipM(k+iHarm-1,l)**2
              End Do
            End Do
            Write (6,Format) Label,(RF*tmp(i),i=1,Jnc)
            write(6,*)
            Do i=1,Jnc
               iIRInt=iIRInt+1
               IRInt(iIRInt)=RF*tmp(i)
            enddo
         Else
            Do i=1,Jnc
               iIRInt=iIRInt+1
               IRInt(iIRInt)=Zero
            enddo
         End if
*
         Write(Format,'(A,I3,A)') '(5X,A,1x,',Jnc,'F10.5)'
         Do iInt = 1, nX
            Write (LUt,Format) ChDisp(iInt+iOff)(1:LENIN6),
     &            (EVec(1,iInt,i),
     &            i=iHarm,iHarm+Jnc-1)
         End Do
         Write (LUt,*)
         Write (LUt,*)
      End Do
*
      Call GetMem('Temp','ALLO','REAL',ipT,nDim*nX)
      ij=-1
      Do i=1,nDim
         Do j=1,nX
            ij=ij+1
            Work(ipT+ij)=Evec(1,j,i)
         End Do
      End Do
      Call WRH(Lu_10,1,nX,nDim,Work(ipT),EVAL,1,'*FREQUENCIES')
      Call GetMem('Temp','FREE','REAL',ipT,nDim*nX)
*
      If (ictl.ne.0) Then
*        Write(Lu_10,*) '*BEGIN PROJECTED DIPOLE TRANSITIONS'
         Do j=1,iel
            Call WRH(Lu_10,1,ndim,ndim,rdum,dDipM(1,j),2,
     &               '*DIPOLE TRANSITIONS')
         End Do
*        Write(Lu_10,*) '*END PROJECTED DIPOLE TRANSITIONS'
      End If
*
      Return
      End
