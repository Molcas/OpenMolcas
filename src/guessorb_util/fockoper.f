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
* Copyright (C) 2004, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine creates a model fock operator based on atomic orbital   *
* energies.                                                            *
*                                                                      *
* This is a very preliminary routine for testing purposes. It relies   *
* on the basis set being of ANO type.                                  *
*                                                                      *
* Absolutely NOT to be used for production!!!!                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: Oct 2004                                                    *
*                                                                      *
************************************************************************
      Subroutine FockOper(RC,Fock)
      Implicit None
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "commgo.fh"
*----------------------------------------------------------------------*
* Parameters                                                           *
*----------------------------------------------------------------------*
      Integer MxComp
      Parameter ( MxComp = 4 )
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Real*8  Fock(*)
      Integer RC
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Logical Debug
      Logical Trace
      Logical Found
      Integer iSym
      Integer iBas
      Integer iOff
      Integer iNuc
      Integer nBasTot
      Integer iUse(MxAtom,MxComp)
      Integer nData
      Integer i,k
      Real*8  energy
*----------------------------------------------------------------------*
* Some setup                                                           *
*----------------------------------------------------------------------*
      Debug=.false.
      Trace=.false.
      If(Trace) Write(6,*) '>>> Entering fockoper'
      RC=0
*----------------------------------------------------------------------*
* Setup various counters.                                              *
*----------------------------------------------------------------------*
      nBasTot=0
      Do iSym=1,nSym
         nBasTot=nBasTot+nBas(iSym)
      End Do
*----------------------------------------------------------------------*
* Is Fock operator on disk?                                            *
*----------------------------------------------------------------------*
      Call Qpg_dArray('Eorb',Found,nData)
      If(Found) Then
         Call Get_dArray('Eorb',Fock,nData)
         If(Debug) Then
            Write(6,*)
            Write(6,*) 'Found Eorb'
            Write(6,'(10F12.6)') (Fock(i),i=1,nData)
            Write(6,*)
         End If
         If(.true.) Return
      End If
      If(.true.) Then
         RC=1
         return
      End If
      Write(6,*) '***'
      Write(6,*) '*** Warning: using built in fock operator'
      Write(6,*) '***'
*----------------------------------------------------------------------*
* Create model Fock operator.                                          *
*----------------------------------------------------------------------*
      iOff=0
      Do iSym=1,nSym
         If(Debug) Then
            Write(6,*) '***'
            Write(6,*) '*** Symmetry',iSym
            Write(6,*) '***'
         End If
         Do i=1,MxAtom
            Do k=1,MxComp
               iUse(i,k)=0
            End Do
         End Do
         Do iBas=1,nBas(iSym)
            iNuc=0
            Do i=1,nNuc
               If(Name(i).eq.Label(iBas+iOff)(1:LENIN)) iNuc=i
            End Do
            If(Debug) Then
               Write(6,'(2(a,i3),3a,i3,f6.2)')
     &            'iSym:',iSym,' iBas:',iBas,
     &            ' = ',Label(iBas+iOff)(1:LENIN),
     &            Label(iBas+iOff)(LENIN1:LENIN8),
     &            iNuc,xCharge(iNuc)
            End If
            If(iNuc.eq.0) Then
               Call SysAbendMsg('fockoper','Fatal','001')
            End If
            energy=0.0d0
            If(Abs(xCharge(iNuc)-1.0d0).lt.1.0d-3) Then
               If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'1s      ') Then
                  iUse(iNuc,1)=iUse(iNuc,1)+1
                  If(iUse(iNuc,1).eq.1) energy=-0.50000d0
               End If
            Else If(Abs(xCharge(iNuc)-3.0d0).lt.1.0d-3) Then
               If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'1s      ') Then
                  iUse(iNuc,1)=iUse(iNuc,1)+1
                  If(iUse(iNuc,1).eq.1) energy=-2.47773d0
                  If(iUse(iNuc,1).eq.2) energy=-0.19632d0
               End If
            Else If(Abs(xCharge(iNuc)-6.0d0).lt.1.0d-3) Then
               If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'1s      ') Then
                  iUse(iNuc,1)=iUse(iNuc,1)+1
                  If(iUse(iNuc,1).eq.1) energy=-11.32554d0
                  If(iUse(iNuc,1).eq.2) energy=-0.70563d0
               Else If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'2px     ')
     &         Then
                  iUse(iNuc,2)=iUse(iNuc,2)+1
                  If(iUse(iNuc,2).eq.1) energy=-0.43335d0
               Else If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'2py     ')
     &         Then
                  iUse(iNuc,3)=iUse(iNuc,3)+1
                  If(iUse(iNuc,3).eq.1) energy=-0.43335d0
               Else If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'2pz     ')
     &         Then
                  iUse(iNuc,4)=iUse(iNuc,4)+1
                  If(iUse(iNuc,4).eq.1) energy=-0.43335d0
               End If
            Else If(Abs(xCharge(iNuc)-7.0d0).lt.1.0d-3) Then
               If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'1s  ') Then
                  iUse(iNuc,1)=iUse(iNuc,1)+1
                  If(iUse(iNuc,1).eq.1) energy=-15.62909d0
                  If(iUse(iNuc,1).eq.2) energy=-0.94531d0
               Else If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'2px     ')
     &         Then
                  iUse(iNuc,2)=iUse(iNuc,2)+1
                  If(iUse(iNuc,2).eq.1) energy=-0.56758d0
               Else If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'2py     ')
     &         Then
                  iUse(iNuc,3)=iUse(iNuc,3)+1
                  If(iUse(iNuc,3).eq.1) energy=-0.56758d0
               Else If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'2pz     ')
     &         Then
                  iUse(iNuc,4)=iUse(iNuc,4)+1
                  If(iUse(iNuc,4).eq.1) energy=-0.56758d0
               End If
            Else If(Abs(xCharge(iNuc)-8.0d0).lt.1.0d-3) Then
               If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'1s  ') Then
                  iUse(iNuc,1)=iUse(iNuc,1)+1
                  If(iUse(iNuc,1).eq.1) energy=-20.66866d0
                  If(iUse(iNuc,1).eq.2) energy=-1.24433d0
               Else If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'2px     ')
     &         Then
                  iUse(iNuc,2)=iUse(iNuc,2)+1
                  If(iUse(iNuc,2).eq.1) energy=-0.63192d0
               Else If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'2py     ')
     &         Then
                  iUse(iNuc,3)=iUse(iNuc,3)+1
                  If(iUse(iNuc,3).eq.1) energy=-0.63192d0
               Else If(Label(iBas+iOff)(LENIN1:LENIN8).eq.'2pz     ')
     &         Then
                  iUse(iNuc,4)=iUse(iNuc,4)+1
                  If(iUse(iNuc,4).eq.1) energy=-0.63192d0
               End If
            Else
               Call SysAbendMsg('fockoper','Fatal','002')
            End If
            Fock(iOff+iBas)=energy
         End Do
         iOff=iOff+nBas(iSym)
      End Do
*----------------------------------------------------------------------*
* Done, deallocate the rest.                                           *
*----------------------------------------------------------------------*
      If(trace) Write(6,*) '<<< Exiting fockoper'
      Return
      End
