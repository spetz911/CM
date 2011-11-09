#! /usr/bin/env python3
# -*- coding: utf-8 -*-






start_f = lambda x: sin(x)
border1 = lambda x: 0
border2 = lambda x: 0

x0 = 0
x1 = 1

h = 0.1
tau = 0.5

T = 7.0

u0 = [start_f(x) for x in frange(x0, x1, h)]
U1 = [u0]
U2 = [u0]






s0 = "15*u_t = 0.1*u_xx + 0.2*u_x + -3.0*u + x*t"




