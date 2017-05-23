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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_Set2CL(iAtomPair,n2CF_Added)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: Set list of two-center functions for atom pair iAtomPair
C     based on diagonal integrals. On exit, n2CF_Added is the number of
C     two-center funtions added (not the total number; some might be
C     present already).
C
      Implicit None
      Integer iAtomPair
      Integer n2CF_Added
#include "WrkSpc.fh"
#include "localdf.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom

      Character*8 Label

      Integer iAtom, jAtom
      Integer ip_D, l_D
      Integer nBas_iAtom, nBas_jAtom
      Integer uv, nuv
      Integer n2CF_Old
      Integer ip_OldList, l_OldList
      Integer ip_NewList, l_NewList
      Integer ipi, ipj
      Integer iS, jS
      Integer iShell, jShell
      Integer u, v
      Integer i2CFun

      Real*8  Thr

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_2CFunctions
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Init counter
      n2CF_Added=0

      ! Get requested accuracy
      Thr=Thr_Accuracy

      ! Get atoms of atom pair
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)

      ! Get dimensions
      nBas_iAtom=LDF_nBas_Atom(iAtom)
      nBas_jAtom=LDF_nBas_Atom(jAtom)

      ! Get pointer to diagonals
      If (iAtom.eq.jAtom) Then
         l_D=nBas_iAtom*(nBas_iAtom+1)/2
         Call GetMem('DiaLT','Allo','Real',ip_D,l_D)
         Call LDF_Q2LT(iAtom,Work(iWork(ip_AP_Diag-1+iAtomPair)),
     &                       Work(ip_D))
      Else
         l_D=0
         ip_D=iWork(ip_AP_Diag-1+iAtomPair)
      End If
      ip_D=ip_D-1

      ! Count diagonals greater than threshold
      If (iAtom.eq.jAtom) Then
         nuv=nBas_iAtom*(nBas_iAtom+1)/2
      Else If (iAtom.gt.jAtom) Then
         nuv=nBas_iAtom*nBas_jAtom
      Else
         Call WarningMessage(2,'LDF_Set2CL: iAtom<jAtom')
         Call LDF_Quit(1)
         nuv=0
      End If
      Do uv=1,nuv
         If (abs(Work(ip_D+uv)).gt.Thr) Then
            n2CF_Added=n2CF_Added+1
         End If
      End Do

      ! Set new list
      If (n2CF_Added.gt.0) Then
         ! If two-center functions are already present, get pointer to
         ! old list
         n2CF_Old=AP_2CFunctions(1,iAtomPair)
         If (n2CF_Old.gt.0) Then
            ip_OldList=AP_2CFunctions(2,iAtomPair)
            l_OldList=4*n2CF_Old
         Else
            ip_OldList=0
            l_OldList=0
         End If
         ! Store total number of two-center functions
         iWork(ip_AP_2CFunctions+2*(iAtomPair-1))=n2CF_Old+n2CF_Added
         ! Allocate new list
         l_NewList=4*AP_2CFunctions(1,iAtomPair)
         Write(Label,'(A,I5.5)') '2CF',iAtomPair-1
         Call GetMem(Label,'Allo','Inte',ip_NewList,l_NewList)
         ! Store pointer to new list
         iWork(ip_AP_2CFunctions+2*(iAtomPair-1)+1)=ip_NewList
         ! Copy contents of old list to new list
         ! Deallocate old list
         If (n2CF_Old.gt.0) Then
            Call iCopy(l_OldList,iWork(ip_OldList),1,
     &                           iWork(ip_NewList),1)
            Call GetMem(Label,'Free','Inte',ip_OldList,l_OldList)
         End If
         ! Store data for newly added two-center functions
         If (iAtom.eq.jAtom) Then
            ipi=LDF_lShell_Atom(iAtom)-1
            uv=ip_D
            i2CFun=n2CF_Old
            Do iS=1,LDF_nShell_Atom(iAtom)
               iShell=iWork(ipi+iS)
               Do jS=1,iS-1
                  jShell=iWork(ipi+jS)
                  Do v=1,nBasSh(jShell)
                     Do u=1,nBasSh(iShell)
                        uv=uv+1
                        If (abs(Work(uv)).gt.Thr) Then
                           iWork(ip_NewList+4*i2CFun)=iS
                           iWork(ip_NewList+4*i2CFun+1)=u
                           iWork(ip_NewList+4*i2CFun+2)=jS
                           iWork(ip_NewList+4*i2CFun+3)=v
                           i2CFun=i2CFun+1
                        End If
                     End Do
                  End Do
               End Do
               Do v=1,nBasSh(iShell)
                  Do u=1,v
                     uv=uv+1
                     If (abs(Work(uv)).gt.Thr) Then
                        iWork(ip_NewList+4*i2CFun)=iS
                        iWork(ip_NewList+4*i2CFun+1)=u
                        iWork(ip_NewList+4*i2CFun+2)=iS
                        iWork(ip_NewList+4*i2CFun+3)=v
                        i2CFun=i2CFun+1
                     End If
                  End Do
               End Do
            End Do
         Else
            ipi=LDF_lShell_Atom(iAtom)-1
            ipj=LDF_lShell_Atom(jAtom)-1
            uv=ip_D
            i2CFun=n2CF_Old
            Do jS=1,LDF_nShell_Atom(jAtom)
               jShell=iWork(ipj+jS)
               Do iS=1,LDF_nShell_Atom(iAtom)
                  iShell=iWork(ipi+iS)
                  Do v=1,nBasSh(jShell)
                     Do u=1,nBasSh(iShell)
                        uv=uv+1
                        If (abs(Work(uv)).gt.Thr) Then
                           iWork(ip_NewList+4*i2CFun)=iS
                           iWork(ip_NewList+4*i2CFun+1)=u
                           iWork(ip_NewList+4*i2CFun+2)=jS
                           iWork(ip_NewList+4*i2CFun+3)=v
                           i2CFun=i2CFun+1
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End If
#if defined (_DEBUG_)
         If (i2CFun.ne.AP_2CFunctions(1,iAtomPair)) Then
            Call WarningMessage(2,'LDF_Set2CL: count error!')
            Write(6,'(A,3I10)') 'iAtomPair,iAtom,jAtom:',iAtomPair,
     &                          iAtom,jAtom
            Write(6,'(A,2I10)') 'i2CFun,AP_2CFunctions:',i2CFun,
     &                          AP_2CFunctions(1,iAtomPair)
            Write(6,'(A,2I10)') 'n2CF_Old,n2CF_Added:',n2CF_Old,
     &                          n2CF_Added
            Call LDF_Quit(1)
         End If
#endif
      End If

      If (l_D.gt.0) Then
         ip_D=ip_D+1
         Call GetMem('DiaLT','Free','Real',ip_D,l_D)
      End If

      End
