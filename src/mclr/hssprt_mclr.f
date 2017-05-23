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
      SubRoutine HssPrt_MCLR(ideg,Hess,ldisp)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

#include "Input.fh"
      Integer  kDisp(8),ldisp(nsym)
      Character Title*39
      Real*8     Hess(*)
      integer ideg(*)
*
      Ind(idisp,jdisp)=idisp*(idisp-1)/2+jdisp
*
      iDisp=0
      Do iIrrep=1,nIrrep
           kDisp(iIrrep)=iDisp
           iDisp=iDisp+lDisp(iIrrep)
           Write (6,*) lDisp(iIrrep)
      End Do
*
      Call GetMem('Temp','ALLO','REAL',ipT,iDisp**2)
      iaa=0
      Do iIrrep=1,nIrrep
        If (ldisp(iirrep).ne.0) Then
        Write(title,'(A,A)') 'Hessian in Irrep ',
     &                chirr(iIrrep)

        Do i=1,lDisp(iirrep)
         Do j=1,i
          ii=ind(i,j)-1
          jj=iaa+ind(i,j)
          Work(ipT+ii)=Hess(jj)*
     &        sqrt(DBLE(ideg(i+kdisp(iirrep))*
     &                    ideg(j+kdisp(iirrep)) ))
         End Do
        End Do
        Call TriPrt(title,' ',Work(ipT),ldisp(iirrep))
        iaa=iaa+ind(ldisp(iirrep),ldisp(iirrep))
        End If
       End Do
       Call GetMem('Temp','FREE','REAL',ipT,idum)

      Return
      End
