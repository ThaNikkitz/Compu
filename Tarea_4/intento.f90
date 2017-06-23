program Tarea4
implicit none

integer :: i, j, k, n
double precision :: y(2), t, dt, alfa(3), yi(2), beta(7), a, b, h, tf, ti, tp

Open(unit = 10, file = 'input.txt')

read(10, *) h

read(10, *) (alfa(j), j = 1, 3)

read(10, *) (beta(i), i = 1, 7)

 Close(unit = 10)

Open(unit = 11, file = 'oli.txt')

!write(11, '(3(1x,a15))') 'Tau', 'Tau_gorrito', 'y', 'dy'

Do j = 1, 3
Do k = 1, 7
b = beta(k)
a = alfa(j)

yi(1) = 0.
yi(2) = 1.
ti = 0.
tf = 8.

t = ti
n = (tf - ti)/h + 1
tp = t

y = yi

!write(11,*) 'a, b'
!write(11,*) a, b
!write(11,*) 

Do i = 1, n
	write(11, '(4(1x,1pe15.7))') t, t-log(b), y(1), y(2)
	Call RK4(t, y, n, h, a, b)
End Do
End Do
End Do

 Close(unit = 11)

end program

!-------------------------------------------------------------

Subroutine Integrator(t, y, n, h, tend, a, b)

double precision :: t, y(2), h, tend, a, b

Do
	If (tend - t < h) Then
		h = tend - t
	End If
	Call RK4(t, y, n, h, a, b)
	If (t >= tend) Then
		Exit
	End If
End Do

End Subroutine

!-------------------------------------------------------------

Subroutine RK4(t, y, n, h, a, b)

double precision :: t, y(2), h, ym(2), ye(2), slope(2), k1(2), k2(2), k3(2), k4(2), a, b

Call Derivs(t, y, k1, a, b)
	ym = y + k1*h/2

Call Derivs(t + h/2, ym, k2, a, b)
	ym = y + k2*h/2

Call Derivs(t + h/2, ym, k3, a, b)
	ye = y + k3*h

Call Derivs(t + h, ye, k4, a, b)
	slope = (k1 + 2*(k2 + k3) + k4)/6
	y = y + slope*h

t = t + h

End Subroutine
!-------------------------------------------------------------

Subroutine Derivs(t, y, dy, a, b)

double precision :: y(2), dy(2), t, b, a

dy(1) = y(2)
dy(2) = -3.*y(2) - y(1)*(2 - (a/(b*exp(t)))*(1 - 1/((b**2)*exp(2*t))))

End Subroutine


