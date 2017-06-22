program Tarea4
implicit none

integer :: i, j, k, n
double precision :: y, m, t, dt, alfa, beta
double precision, allocatable :: xp(:), yi(:)

x = xi
m = 1
xp(m) = x

Do i = 1, n
	yp(i,m) = yi
	y(i) = yi(i)
End Do

Do
	xend = x + xout
	If (xend > xf) Then
		xend = xf
	End If
	h ) dx
	Call Integrator(x, y, n, h, xend)
	m = m + 1
	xp(m) = x
	Do i = 1,n
		yp(i, m) = yi
	End Do
	If (x >= xf) Then
		Exit
	End If
End Do

end program

!-------------------------------------------------------------

Subroutine Integrator(x, y, n, h, xend)

integer :: n
double precision :: x, y, h, xend

Do
	If (xend - x < h) Then
		h = xend - x
	End If
	Call RK4(x, y, n, h)
	If (x >= xend) Then
		Exit
	End If
End Do
End Subroutine

!-------------------------------------------------------------

Subroutine RK4(x, y, n, h)

integer :: n
double precision :: x, y, h

Call Derivs(x, y, k1)
Do i = 1, n
	ym(i) = y(i) + k1(i)*h/2
End Do
Call Derivs(x + h/2, ym, k2)
Do i = 1, n
	ym(i) = y(i) + k2(i)*h/2
End Do
Call Derivs(x + h/2, ym, k3)
Do i = 1, n
	ye(i) = y(i) + k3(i)*h
End Do
Call Derivs(x + h, ye, k4)
Do i = 1, n
	slope(i) = (k1(i) + 2*(k2(i) + k3(i)) + k4(i))/6
	y(i) = y(i) + slope(i)*h
End Do
x = x + h
End Subroutine

!-------------------------------------------------------------

Subroutine Derivs(x, y, dy)






