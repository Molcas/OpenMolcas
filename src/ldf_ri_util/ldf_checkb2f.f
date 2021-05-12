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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_CheckB2F(n,Full)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: check extraction blocked to full and back.
C
      Implicit None
      Integer n
      Real*8  Full(n,n)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Character*12 SecNam
      Parameter (SecNam='LDF_CheckB2F')

      Real*8 Tol
      Parameter (Tol=1.0d-14)

      Logical  BlockMatricesAreIdentical

      Logical  isSymmetric
      External isSymmetric

      Integer ip_Packed, l_Packed
      Integer ip_Packed2, l_Packed2
      Integer ip_Full2, l_Full2
      Integer u, v, ip0, uv, iCount
      Integer ip_BlocksFromFull, ip_BlocksFromPacked

      Integer i, j
      Integer iTri
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      ! Check dimension
      If (n.ne.nBas_Valence) Then
         Call WarningMessage(2,SecNam//': n != nBas_Valence')
         Call LDF_Quit(1)
      End If

      ! Check that Full is symmetric
      If (.not.isSymmetric(Full,n,1.0d-14)) Then
         Call WarningMessage(2,SecNam//': Full is not symmetric')
         Call LDF_Quit(1)
      End If

      ! Allocate extra Packed and Full matrices
      l_Packed=n*(n+1)/2
      Call GetMem('Packed','Allo','Real',ip_Packed,l_Packed)
      l_Packed2=l_Packed
      Call GetMem('Packed2','Allo','Real',ip_Packed2,l_Packed2)
      l_Full2=n**2
      Call GetMem('Full2','Allo','Real',ip_Full2,l_Full2)

      ! Allocate block matrices
      Call LDF_AllocateBlockMatrix('BfF',ip_BlocksFromFull)
      Call LDF_AllocateBlockMatrix('BfP',ip_BlocksFromPacked)

      ! Get packed matrix
      ip0=ip_Packed-1
      Do v=1,n
         Do u=v,n
            Work(ip0+iTri(u,v))=Full(u,v)
         End Do
      End Do

      ! Get extra Packed and Full matrices
      Call dCopy_(l_Packed,Work(ip_Packed),1,Work(ip_Packed2),1)
      Call dCopy_(l_Full2,Full,1,Work(ip_Full2),1)

      ! Get block matrix from Full and from Packed
      Call LDF_Full2Blocked(Full,.False.,ip_BlocksFromFull)
      Call LDF_Full2Blocked(Work(ip_Packed),.True.,ip_BlocksFromPacked)

      ! Check that block matrices are identical
      If (.not.BlockMatricesAreIdentical(ip_BlocksFromFull,
     &                                   ip_BlocksFromPacked,
     &                                   1.0d-14)) Then
         Call WarningMessage(2,SecNam//': Error in extraction')
         Call LDF_Quit(1)
      End If

      ! Extract Packed and Full matrices from blocked ones
      Call LDF_Blocked2Full(ip_BlocksFromFull,.True.,Work(ip_Packed2))
      Call LDF_Blocked2Full(ip_BlocksFromFull,.False.,Work(ip_Full2))

      ! Check that they are identical to input matrix
      iCount=0
      uv=0
      Do v=1,n
         Do u=1,n
            uv=uv+1
            If (abs(Work(ip_Full2-1+uv)-Full(u,v)).gt.Tol) Then
               iCount=iCount+1
            End If
         End Do
      End Do
      If (iCount.ne.0) Then
         Call WarningMessage(2,SecNam//': Error in Full extraction')
         Call LDF_Quit(1)
      End If
      iCount=0
      Do uv=1,l_Packed
         If (abs(Work(ip_Packed2-1+uv)-Work(ip_Packed-1+uv)).gt.Tol)
     &   Then
            iCount=iCount+1
         End If
      End Do
      If (iCount.ne.0) Then
         Call WarningMessage(2,SecNam//': Error in Packed extraction')
         Call LDF_Quit(1)
      End If

      ! Deallocate block matrices
      Call LDF_DeallocateBlockMatrix('BfF',ip_BlocksFromFull)
      Call LDF_DeallocateBlockMatrix('BfP',ip_BlocksFromPacked)

      ! Deallocate extra Packed and Full matrices
      Call GetMem('Full2','Free','Real',ip_Full2,l_Full2)
      Call GetMem('Packed2','Free','Real',ip_Packed2,l_Packed2)
      Call GetMem('Packed','Free','Real',ip_Packed,l_Packed)

      Write(6,'(/,A,A,/)')
     & SecNam,': Full2Blocked and Blocked2Full checked - all seems OK!'
      Call xFlush(6)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Logical Function BlockMatricesAreIdentical(ip_Block1,ip_Block2,
     &                                           Tol)
      Implicit None
      Integer ip_Block1
      Integer ip_Block2
      Real*8  Tol
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer iCount
      Integer AB
      Integer ip1, ip2
      Integer k, n

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      iCount=0
      Do AB=1,NumberOfAtomPairs
         ip1=iWork(ip_Block1-1+AB)
         ip2=iWork(ip_Block2-1+AB)
         n=LDF_nBas_Atom(AP_Atoms(1,AB))*LDF_nBas_Atom(AP_Atoms(2,AB))
         Do k=0,n-1
            If (abs(Work(ip1+k)-Work(ip2+k)).gt.Tol) Then
               iCount=iCount+1
            End If
         End Do
      End Do
      BlockMatricesAreIdentical=iCount.eq.0

      End
