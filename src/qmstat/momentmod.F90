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
      Subroutine MomentMod(ipRe,ipNRe,iCmo,nBRe,nBNRe,LindMOs,iS1,iS2   &
     &                    ,First,DiffMax)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "files_qmstat.fh"
#include "qminp.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension DipRe(3),DipNRe(3)

      Logical LindMOs(MxBas),First

      If(First.and.iPrint.ge.5) then
        Write(6,*)
        Write(6,*)'     Modifications of dipoles by renormalization and'&
     &//' basis reduction.'
        Write(6,*)
        Write(6,*)'     State pair    |  Difference '
        Write(6,*)'     --------------|---------------------'
        First=.false.
      Endif

      nSize1=nBNRe*(nBNRe+1)/2
      nSize2=nBRe*(nBRe+1)/2
      Call GetMem('DipX','Allo','Real',ipDx,nSize1+4)
      Call GetMem('DipY','Allo','Real',ipDy,nSize1+4)
      Call GetMem('DipZ','Allo','Real',ipDz,nSize1+4)
      Call GetMem('DipXre','Allo','Real',ipDxRe,nSize2)
      Call GetMem('DipYre','Allo','Real',ipDyRe,nSize2)
      Call GetMem('DipZre','Allo','Real',ipDzRe,nSize2)
      Call GetMem('DipXsq','Allo','Real',ipDxsq,nBNRe**2)
      Call GetMem('DipYsq','Allo','Real',ipDysq,nBNRe**2)
      Call GetMem('DipZsq','Allo','Real',ipDzsq,nBNRe**2)
      Call GetMem('DipXm','Allo','Real',ipDxM,nBNRe**2)
      Call GetMem('DipYm','Allo','Real',ipDyM,nBNRe**2)
      Call GetMem('DipZm','Allo','Real',ipDzM,nBNRe**2)
      Call GetMem('TEMP','Allo','Real',ipTEMP,nBNRe**2)
      irc=-1
      iopt=0
      iSmLbl=0
!--- X
      icomp=1
      Call RdOne(irc,iopt,'Mltpl  1',icomp,Work(ipDx),iSmLbl)
      Call Square(Work(ipDx),Work(ipDxsq),iONE,nBNRe,nBNRe)
      Call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,ONE,Work(iCmo)              &
     &          ,nBNRe,Work(ipDxsq),nBNRe,ZERO,Work(ipTEMP)             &
     &          ,nBNRe)
      Call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,ONE,Work(ipTEMP)            &
     &          ,nBNRe,Work(iCmo),nBNRe,ZERO,Work(ipDxM)                &
     &          ,nBNRe)
!--- Y
      icomp=2
      Call RdOne(irc,iopt,'Mltpl  1',icomp,Work(ipDy),iSmLbl)
      Call Square(Work(ipDy),Work(ipDysq),iONE,nBNRe,nBNRe)
      Call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,ONE,Work(iCmo)              &
     &          ,nBNRe,Work(ipDysq),nBNRe,ZERO,Work(ipTEMP)             &
     &          ,nBNRe)
      Call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,ONE,Work(ipTEMP)            &
     &          ,nBNRe,Work(iCmo),nBNRe,ZERO,Work(ipDyM)                &
     &          ,nBNRe)
!--- Z
      icomp=3
      Call RdOne(irc,iopt,'Mltpl  1',icomp,Work(ipDz),iSmLbl)
      Call Square(Work(ipDz),Work(ipDzsq),iONE,nBNRe,nBNRe)
      Call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,ONE,Work(iCmo)              &
     &          ,nBNRe,Work(ipDzsq),nBNRe,ZERO,Work(ipTEMP)             &
     &          ,nBNRe)
      Call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,ONE,Work(ipTEMP)            &
     &          ,nBNRe,Work(iCmo),nBNRe,ZERO,Work(ipDzM)                &
     &          ,nBNRe)
!--- Triangualize and reduce.
      kaunt1=0
      kaunt2=0
      Do 2001, i=1,nBNRe
        Do 2002, j=1,nBNRe
          If(j.le.i) then
            If(LindMOs(i).and.LindMOs(j)) then
              Work(ipDxRe+kaunt1)=Work(ipDxM+kaunt2)
              Work(ipDyRe+kaunt1)=Work(ipDyM+kaunt2)
              Work(ipDzRe+kaunt1)=Work(ipDzM+kaunt2)
              kaunt1=kaunt1+1
            Endif
          Endif
          kaunt2=kaunt2+1
2002    Continue
2001  Continue
!--- Density
      DipNRe(1)=Ddot_(nBNRe**2,Work(ipDxM),iONE,Work(ipNRe),iONE)
      DipNRe(2)=Ddot_(nBNRe**2,Work(ipDyM),iONE,Work(ipNRe),iONE)
      DipNRe(3)=Ddot_(nBNRe**2,Work(ipDzM),iONE,Work(ipNRe),iONE)
      DipRe(1)=Ddot_(nSize2,Work(ipDxRe),iONE,Work(ipRe),iONE)
      DipRe(2)=Ddot_(nSize2,Work(ipDyRe),iONE,Work(ipRe),iONE)
      DipRe(3)=Ddot_(nSize2,Work(ipDzRe),iONE,Work(ipRe),iONE)
      Diffx=abs(DipRe(1)-DipNRe(1))
      Diffy=abs(DipRe(2)-DipNRe(2))
      Diffz=abs(DipRe(3)-DipNRe(3))
      If(iPrint.ge.5) then
        Write(6,99)iS1,iS2,'(',Diffx,',',Diffy,',',Diffz,')'
      Endif
99    Format('     ',2I3,'          ',3(A,F10.7),A)
!--- Return number
      DiffMax=Diffy
      If(Diffx.ge.Diffy) DiffMax=Diffx
      If(Diffz.ge.Diffx.and.Diffz.ge.Diffy) DiffMax=Diffz
!--- Deallocate en masse.
      Call GetMem('DipX','Free','Real',ipDx,nSize1+4)
      Call GetMem('DipY','Free','Real',ipDy,nSize1+4)
      Call GetMem('DipZ','Free','Real',ipDz,nSize1+4)
      Call GetMem('DipXre','Free','Real',ipDxRe,nSize2)
      Call GetMem('DipYre','Free','Real',ipDyRe,nSize2)
      Call GetMem('DipZre','Free','Real',ipDzRe,nSize2)
      Call GetMem('DipXsq','Free','Real',ipDxsq,nBNRe**2)
      Call GetMem('DipYsq','Free','Real',ipDysq,nBNRe**2)
      Call GetMem('DipZsq','Free','Real',ipDzsq,nBNRe**2)
      Call GetMem('DipXm','Free','Real',ipDxM,nBNRe**2)
      Call GetMem('DipYm','Free','Real',ipDyM,nBNRe**2)
      Call GetMem('DipZm','Free','Real',ipDzM,nBNRe**2)
      Call GetMem('TEMP','Free','Real',ipTEMP,nBNRe**2)

      Return
      End
