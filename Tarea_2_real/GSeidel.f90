subroutine Gseid(a, b, n, x, imax, es, lambda)
	integer :: i, j, n
	n = 10
	double precision :: a(n, n), b(n), x(n,n), dummy, sum
	Do i = 1, n
		Do j = 1, n
			dummy = a(i,i)/dummy
		End Do
		b(i) = b(i)/dummy
	End Do
	Do i = 1, n
		sum = b(i)
		Do j = 1, n
			if (i /= j)
				sum = sum - a(i,j)*x(i,j)
