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
      Subroutine LDF_ComputeZVec(iAtomPair,ip_V,l_V,ip_G,l_G,ip_Z,l_Z,
     &                           irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: Compute Z vectors by CD of the G matrix:
C
C        G[JK] = (J|K) = sum(I) Z(J,I)*Z(K,I)
C
C     The Z vectors are ordered according to auxiliary shells (plus
C     shell pairs if 2-center functions are included), as specified
C     through the atom pair list.
C
C     Linear dependence and 2-center functions may be redefined by this
C     routine (if linear dependence is detected). In that case, M may be
C     changed on exit from this routine.
C
C     The array V may be used as work space (it serves no other purpose)
C
C     irc is a return code (0 if succesful).
C
      Implicit None
      Integer iAtomPair
      Integer ip_V, l_V
      Integer ip_G, l_G
      Integer ip_Z, l_Z
      Integer irc
#include "WrkSpc.fh"

      Character*15 SecNam
      Parameter (SecNam='LDF_ComputeZVec')

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair
#if defined (_DEBUGPRINT_)
      Logical  isSymmetric, hasNonnegativeDiagonal, obeysCauchySchwarz
      External isSymmetric, hasNonnegativeDiagonal, obeysCauchySchwarz
#endif

      Real*8  Thr

      Integer M
      Integer ip_ID, l_ID
      Integer ip_Z_, l_Z_
      Integer ipG0, ipZ0
      Integer nZ

      Integer i, j
      Integer iTri
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      ! Init
      irc=0
      ip_Z=0
      l_Z=0

      ! Get number of auxiliary functions
      M=LDF_nBasAux_Pair(iAtomPair)
      If (M.lt.1) Return

#if defined (_DEBUGPRINT_)
      If (l_G.ne.M*M) Then
         Write(6,'(A,A,2I8)') SecNam,': l_G,M:',l_G,M
         irc=-1
         Return
      End If
      ! Check G matrix for symmetry and non-negative diagonals
      ! and Cauchy-Schwarz inequality
      irc=0
      If (.not.isSymmetric(Work(ip_G),M,1.0d-15)) Then
         Write(6,'(A,A)')
     &   SecNam,': G matrix is not symmetric!'
         irc=irc+1
      End If
      If (.not.hasNonnegativeDiagonal(Work(ip_G),M)) Then
         Write(6,'(A,A)')
     &   SecNam,': G matrix has negative diagonals!'
         irc=irc+1
      End If
      If (.not.obeysCauchySchwarz(Work(ip_G),M,1.0d-14)) Then
         Write(6,'(A,A)')
     &   SecNam,
     &   ': G matrix does not obey the Cauchy-Schwarz inequality!'
         irc=irc+1
      End If
      If (irc.ne.0) Then
         Write(6,'(A,A)')
     &   SecNam,': => G matrix not PSD!!'
         Call LDF_Quit(1)
      End If
#endif

      ! Incomplete CD of G matrix
      Thr=1.0d-14
      l_ID=M
      Call GetMem('CD_ID','Allo','Inte',ip_ID,l_ID)
      l_Z_=M*M
      If (l_V.ge.l_Z_) Then
         ip_Z_=ip_V
      Else
         Call GetMem('Z_','Allo','Real',ip_Z_,l_Z_)
      End If
      nZ=0
      Call CD_InCore_P(Work(ip_G),M,Work(ip_Z_),M,
     &                 iWork(ip_ID),nZ,Thr,irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': CD_InCore_P returned code',irc
         irc=1
         Call GetMem('CD_ID','Free','Inte',ip_ID,l_ID)
         If (ip_Z_.ne.ip_V) Then
            Call GetMem('Z_','Free','Real',ip_Z_,l_Z_)
         End If
         Return
      End If

      ! Handle linear dependence
      Call LDF_RemoveLinDep(iAtomPair,Work(ip_Z_),iWork(ip_ID),M,nZ)

      ! Compute linearly independent G matrix
      Call dGeMM_('N','T',nZ,nZ,nZ,
     &            1.0d0,Work(ip_Z_),M,Work(ip_Z_),M,
     &            0.0d0,Work(ip_G),nZ)

      ! Deallocate ID and Z_ arrays
      Call GetMem('CD_ID','Free','Inte',ip_ID,l_ID)
      If (ip_Z_.ne.ip_V) Then
         Call GetMem('Z_','Free','Real',ip_Z_,l_Z_)
      End If

      ! Compute Z vectors without pivoting by complete CD of G (with
      ! linear dependence removed!)
      Call CCD_InCore(Work(ip_G),nZ,irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': CCD_InCore returned code',irc
         irc=1
         Return
      End If

      ! Allocate memory for Z vectors
      l_Z=nZ*(nZ+1)/2
      Call GetMem('ZVec','Allo','Real',ip_Z,l_Z)

      ! Copy Z vectors from G
      ipZ0=ip_Z-1
      Do J=1,nZ
         ipG0=ip_G-1+nZ*(J-1)
         Do I=J,nZ
            Work(ipZ0+iTri(I,J))=Work(ipG0+I)
         End Do
      End Do

#ifndef _DEBUGPRINT_
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(l_G)
#endif
      End
