#!/usr/bin/env python
# -*- coding: utf-8 -*-
import PyQt4
import PyQt4.Qt
import PyQt4.QtGui
import PyQt4.QtCore
from PyQt4 import *
from PyQt4.Qt import *
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *
from PyQt4.QtCore import *

class Timer(QTime):
	def __init__(self):
		QTime.__init__(self)
		self.n = 0
		self.p = 0
		self.start()
	def tick(self):#милисекунды
		self.p = self.n
		self.n = self.elapsed()
		self.d = abs(self.p-self.n)
	def dt(self):
		return self.d*0.0001
