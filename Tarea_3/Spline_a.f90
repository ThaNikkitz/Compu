program Gancho

implicit none

integer :: i, j
integer, parameter :: n = 22
double precision :: blah(n,2), x(n), a(n), b(n), c(n), d(n), h(n), alfa(n), mu(n), l(n), z(n)

Open(unit = 2, file = 'matrix1.txt')
Do i = 1,n
	read(2,*) (blah(i,j), j = 1,2)
End Do

Do j = 1,n
	x(j) = blah(j,1)
	a(j) = blah(j,2)
End Do

 Close(unit = 2)

Do i = 1,n-1
	h(i) = x(i+1) - x(i)
End Do
Do i = 2, n-1
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
	c(j) = z(j) - mu(j)*c(j+1)
	b(j) = (a(j+1) - a(j))/h(j) - h(j)*(c(j+1) + 2*c(j))/3
	d(j) = (c(j+1) - c(j))/(3*h(j))
End Do

Open(unit = 10, file = 'coeficientes.txt')
	Do i = 1,n
	write(10,*) a(i), b(i), c(i), d(i)
	End Do
 Close(unit = 10)

call Spline(a, b, c, d, e, x, xx, y)

write(*,*) y

end program

!----------------------------------------------------------------------------------

Subroutine Spline(a, b, c, d, e, x, xx, y)

implicit none

integer :: i, j
integer, parameter :: n = 22
double precision :: e(n,4), a, b, c, d, x(n), xx, y

write(*,*) 'Ingrese el valor al que quiere aproximarle una solución'
read(*,*) xx

Open(unit = 11, file = 'coeficientes.txt')
Do i = 1,n	
	read(10,*) (e(i,j), j = 1,4)
End Do 

Open(unit = 3, file = 'matrix1.txt')
Do j = 1, n-1
	read(3,*) x(j)
	If(xx <= x(j+1) .and. xx >= x(j)) then
		a = e(j,1)
		b = e(j,2)
		c = e(j,3)
		d = e(j,j)
		y = a + b*(xx-x(j)) + c*(xx - x(j))**2 + d*(xx - x(j))**3
	End If
End Do
 Close(Unit = 3)
 Close(unit = 11)

End Subroutine
