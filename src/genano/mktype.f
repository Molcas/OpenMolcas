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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
* This routine creates a list of basis function labels.                *
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine MkType
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
      Parameter (MxLst=MxLqn*(MxLqn*(MxLqn+6)+11)/6+1)
#include "common.fh"
#include "symlab.fh"
      Character*4 CrtLst(0:MxLst)
      Character*4 SphLst(0:MxLst)
*----------------------------------------------------------------------*
* Get labels from utility routine                                      *
*----------------------------------------------------------------------*
      Call Make_Labels(CrtLst,SphLst,MxLst,MxLqn)
*     Write(*,*) 'Cartesian label'
*     ind=0
*     Do l=0,MxLqn
*        Write(*,'(15(1x,a4))') (CrtLst(ind+m),m=0,l*(l+3)/2)
*        ind=ind+l*(l+3)/2+1
*     End Do
*     Write(*,*) 'Spherical label'
*     ind=0
*     Do i=0,MxLqn
*        Do l=i,0,-2
*           Write(*,'(15(1x,a4))') (SphLst(ind+m),m=0,2*l)
*           ind=ind+2*l+1
*        End Do
*     End Do
*----------------------------------------------------------------------*
* Copy labels                                                          *
*----------------------------------------------------------------------*
      iFrom=0
      iTo=1
      Do l=0,Mxlqn
         if(l.lt.2) Then
            Do m=0,2*l
               Type(iTo+m)=CrtLst(iFrom+m)
            End Do
         Else
            Do m=0,2*l
               Type(iTo+m)=SphLst(iFrom+m)
            End Do
         End If
         iFrom=iFrom+l*(l+3)/2+1
         iTo=iTo+2*l+1
      End Do
*     Write(*,*) 'Local labels'
*     ind=1
*     Do l=0,MxLqn
*        Write(*,'(15(1x,a4))') (Type(ind+m),m=0,2*l)
*        ind=ind+2*l+1
*     End Do
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
