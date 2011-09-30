#!/usr/bin/env python
# -*- coding: utf-8 -*-
from timer import Timer

import pdb

import sys
import PyQt4
import PyQt4.Qt
import PyQt4.QtGui
import PyQt4.QtCore
from PyQt4 import *
from PyQt4.Qt import *
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *
from PyQt4.QtCore import *
from OpenGL.GL import *
from OpenGL.GLU import *
import Image

VERTEX_ATTRIBUTE = 0
TEXCOORD_ATTRIBUTE = 1

def loadImage(imageName):
	glEnable(GL_TEXTURE_2D)
	im = Image.open(imageName)
	try:
		ix, iy, image = im.size[0], im.size[1], im.tostring("raw", "RGBA", 0, -1)
	except SystemError:
		ix, iy, image = im.size[0], im.size[1], im.tostring("raw", "RGBX", 0, -1)
	ID=0
	ID = glGenTextures(1)
	glBindTexture(GL_TEXTURE_2D, ID)
	glPixelStorei(GL_UNPACK_ALIGNMENT,1)
	glTexImage2D(GL_TEXTURE_2D, 0, 3, ix, iy, 0, GL_RGBA, GL_UNSIGNED_BYTE, image)
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE)
	return ID


class Scene(QGLWidget):
	clicked = QtCore.pyqtSignal()
	def __init__(self, parent = None):
		QGLWidget.__init__(self)

		#~ self.setGeometry(0, 0, 640, 480)

		self.xRot = 0
		self.yRot = 0
		self.zRot = 0
		
		self.timer = QTimer(self)
		self.connect(self.timer, SIGNAL('timeout()'), self, SLOT('updateGL()'))
		self.timer.start(0)
		self.shader = None
		
		self.coords = (
		(( +1, -1, -1 ), ( -1, -1, -1 ), ( -1, +1, -1 ), ( +1, +1, -1 )),
		(( +1, +1, -1 ), ( -1, +1, -1 ), ( -1, +1, +1 ), ( +1, +1, +1 )),
		(( +1, -1, +1 ), ( +1, -1, -1 ), ( +1, +1, -1 ), ( +1, +1, +1 )),
		(( -1, -1, -1 ), ( -1, -1, +1 ), ( -1, +1, +1 ), ( -1, +1, -1 )),
		(( +1, -1, +1 ), ( -1, -1, +1 ), ( -1, -1, -1 ), ( +1, -1, -1 )),
		(( -1, -1, +1 ), ( +1, -1, +1 ), ( +1, +1, +1 ), ( -1, +1, +1 ))
		)
		
		
		
		self.vertices2 = [(-1.0, 0.5, 0.0, 0.5), (1.0, 1.0, 0.0, 0.5),
						 (1.0, -1.0, 0.0, 0.5), (1.0, 0.0, 0.0, 0.0),
						 (0.0,  0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0)]

		#~ glVertex2f(-1.,0.5)
		#~ glVertex2f(1.,1.)
		#~ glVertex2f(1.,-1.)


		self.texCoords2 = [(0.0, 0.1, 0.0, 0.0), (1.0, 0.1, 0.0, 0.0),
						  (0.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0)]

		

		
	def makeObject(self):
		self.texCoords = []
		self.vertices = []

		for i in range(6):
			#~ self.textures.append(
					#~ self.bindTexture(
							#~ QtGui.QPixmap(':/images/side%d.png' % (i + 1))))
			for j in range(4):
				self.texCoords.append(((j == 0 or j == 3), (j == 0 or j == 1)))

				x, y, z = self.coords[i][j]
				self.vertices.append((0.2 * x, 0.2 * y, 0.2 * z))
			
		
		
	def initializeGL(self):
		print("INIT")
		self.rotateBy(+42 * 16, +42 * 16, -21 * 16)
		self.makeObject()
		self.qglClearColor(Qt.red)
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);
		
		if not(QGLShaderProgram.hasOpenGLShaderPrograms()):
			print("shaders not support")
			exit(1)
		
		vshader = QGLShader(QGLShader.Vertex, self);
		
		vsrc = """
			 attribute highp vec4 vertex;
			 attribute mediump vec4 texCoord;
			 varying mediump vec4 texc;
			 uniform mediump mat4 matrix;
			 void main(void)
			 {
				 gl_Position = matrix * vertex;
				 texc = texCoord;
			 }
		"""
		vshader.compileSourceCode(vsrc);

		fshader = QGLShader(QGLShader.Fragment, self);
		fsrc = """
		uniform sampler2D texture;
		uniform mediump vec4 color;
		varying mediump vec4 texc;
		void main(void)
		{
			float x = texc.x;
			float y = texc.y;
			float z = texc.z;
			float cl = x*x+y*y;
			if (cl > 10.0) cl = 0.5;
			if (cl < 0.0001) cl = 1.0;
			
			gl_FragColor = vec4(0.5); //texture2D(texture, texc.st);
		}
		"""
		fshader.compileSourceCode(fsrc);

		program = QGLShaderProgram(self);
		program.addShader(vshader);
		program.addShader(fshader);
		#~ program.bindAttributeLocation("vertex", VERTEX_ATTRIBUTE);
		#~ program.bindAttributeLocation("texCoord", TEXCOORD_ATTRIBUTE);
		program.link();
		
		vertexLocation = program.attributeLocation("vertex");
		texCoordLocation = program.attributeLocation("texCoord");
		#~ textureLocation = program.uniformLocation("texture");
		#~ self.matrixLocation = program.uniformLocation("matrix");
		#~ self.colorLocation = program.uniformLocation("color");
		
		program.bind();
		
		program.setUniformValue("texture", 0);
		color = QColor(0, 255, 0, 0)
		program.setUniformValue("color", color);
		
		
		program.enableAttributeArray(vertexLocation);
		program.enableAttributeArray(texCoordLocation);
		program.setAttributeArray(vertexLocation, self.vertices);
		program.setAttributeArray(texCoordLocation, self.texCoords);
		
		
		
		self.shader = program
		
		texture = loadImage("cubelogo.png")
		glBindTexture(GL_TEXTURE_2D, texture);
		
	#	self.texture = QGLContext.bindTexture(self.context,
	#										  QImage("cubelogo.png"));

		
	def paintGL(self):
		program = self.shader
		if( program == None):
			print("paintGL error")
			exit(1)

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

		m = QMatrix4x4()
		m.ortho(-0.5, 0.5, 0.5, -0.5, 4.0, 15.0)
		m.translate(0.0, 0.0, -10.0)
		m.rotate(self.xRot / 16.0, 1.0, 0.0, 0.0)
		m.rotate(self.yRot / 16.0, 0.0, 1.0, 0.0)
		m.rotate(self.zRot / 16.0, 0.0, 0.0, 1.0)

		program.setUniformValue('matrix', m)


		for i in range(6):
			#~ glBindTexture(GL_TEXTURE_2D, self.textures[i])
			glDrawArrays(GL_TRIANGLE_FAN, i * 4, 4)

		#~ glDrawArrays(GL_TRIANGLES, 0, 3);
		
		#~ program.disableAttributeArray(self.vertexLocation);
		#~ program.disableAttributeArray(self.texCoordLocation);

		return
		
		
		#~ glDrawArrays(GL_QUADS, 0, 4);
		
		#~ glBindTexture(GL_TEXTURE_2D, self.texture);
		
		#~ for i in range(6):
			#~ glDrawArrays(GL_TRIANGLE_FAN, i * 4, 4);
		
	

		#~ glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		

		#~ glColor(1.0,0.5,1.0)
		#~ 
		#~ glBegin(GL_QUADS)
#~ 
		#~ glTexCoord2f(0.,1.)
		#~ glVertex2f(-1.,0.5)
#~ 
		#~ glTexCoord2f(1.,1.)
		#~ glVertex2f(1.,1.)
		#~ 
		#~ glTexCoord2f(1.,0.)
		#~ glVertex2f(1.,-1.)
#~ 
		#~ glTexCoord2f(0.,0.)
		#~ glVertex2f(-1.,-1.)
		#~ glEnd()
		
	def minimumSizeHint(self):
		return QtCore.QSize(50, 50)

	def sizeHint(self):
		return QtCore.QSize(200, 200)
		
	def rotateBy(self, xAngle, yAngle, zAngle):
		self.xRot += xAngle
		self.yRot += yAngle
		self.zRot += zAngle
		self.updateGL()
				
	def resizeGL(self, width, height):
		side = min(width, height)
		glViewport((width - side) / 2, (height - side) / 2, side, side)

	def mousePressEvent(self, event):
		self.lastPos = event.pos()

	def mouseMoveEvent(self, event):
		dx = event.x() - self.lastPos.x()
		dy = event.y() - self.lastPos.y()

		if event.buttons() & QtCore.Qt.LeftButton:
			self.rotateBy(8 * dy, 8 * dx, 0)
		elif event.buttons() & QtCore.Qt.RightButton:
			self.rotateBy(8 * dy, 0, 8 * dx)

		self.lastPos = event.pos()

	def mouseReleaseEvent(self, event):
		self.clicked.emit()

		


	
class Interface(QWidget):
	def __init__(self):
		self.app = QApplication(sys.argv)
		QWidget.__init__(self)
		self.setWindowTitle('Fractal')
		self.scene = Scene(self)

		
	def run(self):
		self.scene.show()
		sys.exit(self.app.exec_())


def main():
	app = QApplication(sys.argv)
	scene = Scene()
	scene.resize(640, 480);
	scene.setWindowTitle('Fractal')
	scene.show()
	sys.exit(app.exec_())
	return 0


if (__name__ == "__main__"):
	main()
