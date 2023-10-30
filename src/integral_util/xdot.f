      Real*8 Function XDot(A,n,m,k,l)
      use Constants, only: Zero, One
      Implicit None
      Integer n, m, k, l
      Real*8 A(n*m+1,k,l)

      Integer ik, il
      Real*8, external :: DDot_
*
      XDot=Zero
      Do ik= 1, k
         Do il = 1, l
            XDot=XDot+DDot_(n*m,One,0,A(1,ik,il),1)
         End Do
      End Do
*
      End Function XDot

