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
      Subroutine RdDens(d1,nd1,d2,nd2)
      Implicit Real*8(a-h,o-z)

#include "Input.fh"

#include "SysDef.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "Files_mclr.fh"
      Real*8 D1(nd1),d2(nd2)
      Dimension rdum(1)
*                                                                      *
************************************************************************
*                                                                      *
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      call dcopy_(nd1,[0.0d0],0,d1,1)
      Call Getmem('TDEND','ALLO','REAL',ipG2tt,nd2)
      Call Getmem('TDEND','ALLO','REAL',ipD2,nd2)
      Call Getmem('TDEND','ALLO','REAL',ipD1,nd1)
      call dcopy_(nd2,[0.0d0],0,Work(ipg2tt),1)
      jDisk = ITOC(3)
      Do i=1,lroots
         W=0.0d0
         Do j=1,nroots
            If (iroot(j).eq.i) W=Weight(j)
         End Do
         Call dDaFile(LUJOB ,2,Work(ipD1),nd1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,nd1,jDisk)
         Call dDaFile(LUJOB ,2,Work(ipD2),ND2,jDisk)
         Call dDaFile(LUJOB ,0,rdum,ND2,jDisk)
         If (W.ne.0.0d0) Then
            call daxpy_(nd2,w,Work(ipD2),1,Work(ipG2tt),1)
            call daxpy_(nd1,w,Work(ipD1),1,D1,1)
         End If
!         Write(*,*) i,w,"LUJOB",LUJOB
!        Call Triprt('D',' ',D1,ntash)
      End Do
      Call Put_D2AV(Work(ipG2tt),nd2)
      Call Put_D1AV(D1,nd1)
*
      Do iB=1,ntash
       Do jB=1,iB
        iDij=iTri(ib,jB)
        Do kB=1,ib
         Do lB=1,kB
          iDkl=iTri(kB,lB)
          fact=1.0d00
          if(iDij.ge.iDkl .and. kB.eq.lB) fact=2.0d00
          if(iDij.lt.iDkl .and. iB.eq.jB) fact=2.0d00
          iijkl=itri(iDij,iDkl)
          D2(iijkl)=Fact*Work(ipG2tt+iijkl-1)
         End Do
        End Do
       End Do
      End Do
c
      Call Getmem('TDEND','FREE','REAL',ipG2tt,nd2)
      Call Getmem('TDEND','FREE','REAL',ipD2,nd2)
      Call Getmem('TDEND','FREE','REAL',ipD1,nd1)
      Return
      End
