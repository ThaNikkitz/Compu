program Gancho

implicit none

integer :: i, j
integer, parameter :: n = 22
double precision :: blah(n,2), x(n), a(n), b(n), c(n), d(n), e(n,4), h(n), alfa(n), mu(n), l(n), z(n), y, xx, aa, bb, cc, dd

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

Open(Unit = 1, file = 'resultados.txt')
Do i = 1, 5*n
	xx = (i - 1.0d0)*212.8d0/(5.0d0*n-1.0d0)
	call Spline(aa, bb, cc, dd, e, x, xx, y)
	write(1,*) y
End Do
 Close(Unit = 1)

end program

!----------------------------------------------------------------------------------

Subroutine Spline(aa, bb, cc, dd, e, x, xx, y)

implicit none

integer :: i, j
integer, parameter :: n = 22
double precision :: e(n,4), aa, bb, cc, dd, x(n), xx, y

!write(*,*) 'Ingrese el valor al que quiere aproximarle una solución'
!read(*,*) xx

Open(Unit = 11, file = 'coeficientes.txt')
Do i = 1,n	
	read(11,*) (e(i,j), j = 1,4)
End Do 

Open(Unit = 3, file = 'matrix1.txt')
Do j = 1, n-1
	read(3,*) x(j)
	If(xx < x(j+1) .and. xx >= x(j)) then
		aa = e(j,1)
		bb = e(j,2)
		cc = e(j,3)
		dd = e(j,4)
		y = aa + bb*(xx - x(j)) + cc*(xx - x(j))**2 + dd*(xx - x(j))**3
	Else If (xx == x(n)) then
		aa = e(n,1)
		bb = e(n,2)
		cc = e(n,3)
		dd = e(n,4)
		y = aa + bb*(xx - x(j)) + cc*(xx - x(j))**2 + dd*(xx - x(j))**3
	End If
	write(*,*) xx, x(n)
End Do

 Close(Unit = 3)
 Close(Unit = 11)

End Subroutine
