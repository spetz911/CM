#!/usr/bin/env python



import sys
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
from PyQt4 import QtCore, QtGui, QtOpenGL

try:
	from OpenGL.GL import *
except ImportError:
	app = QtGui.QApplication(sys.argv)
	QtGui.QMessageBox.critical(None, "OpenGL textures",
			"PyOpenGL must be installed to run this example.")
	sys.exit(1)

import textures_rc


class GLWidget(QtOpenGL.QGLWidget):

	clicked = QtCore.pyqtSignal()

	VERTEX_ATTRIBUTE, TEXCOORD_ATTRIBUTE = range(2)

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
	#read from file "fs.glsl"
	fsrc = """
"""

	coords = (
		(( +1, -1, -1 ), ( -1, -1, -1 ), ( -1, +1, -1 ), ( +1, +1, -1 )),
		(( +1, +1, -1 ), ( -1, +1, -1 ), ( -1, +1, +1 ), ( +1, +1, +1 )),
		(( +1, -1, +1 ), ( +1, -1, -1 ), ( +1, +1, -1 ), ( +1, +1, +1 )),
		(( -1, -1, -1 ), ( -1, -1, +1 ), ( -1, +1, +1 ), ( -1, +1, -1 )),
		(( +1, -1, +1 ), ( -1, -1, +1 ), ( -1, -1, -1 ), ( +1, -1, -1 )),
		(( -1, -1, +1 ), ( +1, -1, +1 ), ( +1, +1, +1 ), ( -1, +1, +1 ))
	)

	def __init__(self, parent=None, shareWidget=None):
		super(GLWidget, self).__init__(parent, shareWidget)

		self.clearColor = QtCore.Qt.black
		self.xRot = 0
		self.yRot = 0
		self.zRot = 0

		self.clearColor = QtGui.QColor()
		self.lastPos = QtCore.QPoint()

		self.program = None

	def minimumSizeHint(self):
		return QtCore.QSize(50, 50)

	def sizeHint(self):
		return QtCore.QSize(200, 200)

	def rotateBy(self, xAngle, yAngle, zAngle):
		self.xRot += xAngle
		self.yRot += yAngle
		self.zRot += zAngle
		self.updateGL()

	def setClearColor(self, color):
		self.clearColor = color
		self.updateGL()

	def initializeGL(self):
		self.makeObject()

	    fbo = QGLFramebufferObject(width(), height(),
	            QGLFramebufferObject::Depth);
    	fbo2 = QGLFramebufferObject(width(), height() );


		glEnable(GL_DEPTH_TEST)
		glEnable(GL_CULL_FACE)

		vshader = QtOpenGL.QGLShader(QtOpenGL.QGLShader.Vertex, self)
		vshader.compileSourceCode(self.vsrc)

		fshader = QtOpenGL.QGLShader(QtOpenGL.QGLShader.Fragment, self)
		fshader.compileSourceCode(self.fsrc)

		program = QtOpenGL.QGLShaderProgram(self)
		program.addShader(vshader)
		program.addShader(fshader)
		program.link()
		VERTEX_ATTRIBUTE = program.attributeLocation("vertex");
		TEXCOORD_ATTRIBUTE = program.attributeLocation("texCoord");
		

		program.bind()
		program.setUniformValue('texture', 0)

		program.enableAttributeArray(VERTEX_ATTRIBUTE)
		program.enableAttributeArray(TEXCOORD_ATTRIBUTE)
		program.setAttributeArray(VERTEX_ATTRIBUTE, self.vertices)
		program.setAttributeArray(TEXCOORD_ATTRIBUTE, self.texCoords)
		
		self.program = program

	def paintGL(self):
		self.qglClearColor(self.clearColor)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

		m = QtGui.QMatrix4x4()
		m.ortho(-0.5, 0.5, 0.5, -0.5, 4.0, 15.0)
		m.translate(0.0, 0.0, -5.0)
		m.scale(0.5)
		self.program.setUniformValue('matrix', m)
		
		glDrawArrays(GL_TRIANGLE_FAN, 0, 4)


     float aspect = width() / height();

     QMatrix4x4 projMatrix;
     projMatrix.setToIdentity();

     QMatrix4x4 viewMatrix;
     viewMatrix.setToIdentity();

     QMatrix4x4 modelMatrix;
     modelMatrix.setToIdentity();

     pmvMatrix.setToIdentity();

//     float near = 1.0;
//     float far = 100.0;

     projMatrix.perspective(45.0, aspect, 1.0, 100.0);


     viewMatrix.translate(posX, posY, posZ);
   //  viewMatrix.rotate(90, 1.0, 0.0, 0.0);

     viewMatrix.translate(0, +0.5, -2);

     pmvMatrix *= projMatrix;
     pmvMatrix *= viewMatrix;

     # render to the framebuffer object
     fbo->bind();
     program->bind();

     glViewport(0, 0, fbo->size().width(), fbo->size().height());


     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

     glBindTexture(GL_TEXTURE_2D, cubeTexture);
  //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);


         program->setUniformValue(colorLocation, color1);
         glDrawArrays(GL_TRIANGLES, 0, 6);


     fbo->release();



     glBindTexture(GL_TEXTURE_2D, fbo->texture());
 

		


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

	def makeObject(self):
		self.textures = []
		self.texCoords = []
		self.vertices = []

		for i in range(6):
			self.textures.append(
					self.bindTexture(
							QtGui.QPixmap(':/images/side%d.png' % (i + 1))))

			for j in range(4):
				self.texCoords.append(( (j == 0 or j == 3) * 1, (j == 0 or j == 1) * 1))

				x, y, z = self.coords[i][j]
				self.vertices.append(( x, y, z))



if __name__ == '__main__':
#~ def main():
	app = QApplication(sys.argv)
	scene = GLWidget()
	scene.resize(640, 480);
	scene.setWindowTitle('Fractal')
	scene.rotateBy(+42 * 16, +42 * 16, -21 * 16)

	scene.show()
	sys.exit(app.exec_())

