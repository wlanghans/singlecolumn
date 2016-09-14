! ==================================================================
      SUBROUTINE tridiag(a,b,c,d,nsize)

      use grid
      implicit none

!! to solve system of linear eqs on tridiagonal matrix n times n
!! after Peaceman and Rachford, 1955
!! a,b,c,d - are vectors of order n 
!! a,b,c - are coefficients on the LHS
!! d - is initially RHS on the output becomes a solution vector

!-------------------------------------------------------------------

        INTEGER, INTENT(in) :: nsize
        REAL, DIMENSION(nsize), INTENT(in)    :: a,b
        REAL, DIMENSION(nsize), INTENT(inout) :: c,d

        INTEGER :: i
        REAL :: p
        REAL, DIMENSION(nsize) :: q

        c(nsize)=0.
        q(1)=-c(1)/b(1)
        d(1)=d(1)/b(1)

        DO i=2,nsize
           p=1./(b(i)+a(i)*q(i-1))
           q(i)=-c(i)*p
           d(i)=(d(i)-a(i)*d(i-1))*p
        ENDDO

        DO i=nsize-1,1,-1
           d(i)=d(i)+q(i)*d(i+1)
        ENDDO

      END SUBROUTINE tridiag

