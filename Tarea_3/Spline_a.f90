program Gancho

implicit none

integer :: n, i, j
integer, parameter :: n = 10
double precision :: x(n), a(n), b(n), c(n), d(n), h(n), alfa(n), mu(n), l(n), z(n)

Open(unit = 2, file = 'matrix1.txt')
Do i = 1,n
	read(10,*) x(i)
End Do

Do i = 2,n-1
	h(i) = x(i+1) - x(i)
	alfa(i) = (3/h(i))*(a(i+1) - a(i)) - (3/h(i-1))*(a(i) - a(i-1))
End Do

l(1) = 1
mu(1) = 0
z(1) = 0

Do i = 2,n-1
	l(i) = 2*(x(i+1) - x(i-1)) - h(i-1)*mu(i-1)
	mu(i) = h(i)/l(i)
	z(i) = (alfa(i) - h(i-1)*z(i-1))/l(i)
End Do

l(n) = 1
z(n) = 0
 c(n) = 0

Do j = n-1,1,-1
	c(j) = z(j) - u(j)*c(j+1)
	b(j) = (a(j+1) - a(j))/h(j) - h(j)*(c(j+1) + 2*c(j))/3
	d(j) = (c(j+1) - c(j))/(3*h(j))
End Do



end programm
