/****************************************************************************
**
** Copyright (C) 2010 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#include "glwidget.h"
#include <QtGui/QImage>

#include <math.h>

#include "global.h"


GLWidget::GLWidget(QWidget *parent)
  : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    // create the framebuffer object - make sure to have a current
    // context before creating it
    makeCurrent();

    timerId = startTimer(20);
    setWindowTitle(tr("OpenGL lololo"));
    xxxx = 3.0;

    posX = 0;
    posY = 0;
    posZ = 0;

    rotX = -10;
    rotY = 120;

}

GLWidget::~GLWidget()
{
    delete fbo;
    delete fbo2;
}

void GLWidget::setupShaders()
{
    program = new QGLShaderProgram (this);
    program2 = new QGLShaderProgram (this);
    DOFprogram = new QGLShaderProgram (this);
    DOF2program = new QGLShaderProgram (this);

    program->addShaderFromSourceCode(QGLShader::Vertex,
         "attribute vec4 vertex;\n"
         "uniform mat4 matrix;\n"

         "attribute vec4 texCoord;\n"
         "varying vec4 texc;\n"
         "out vec4 vertex1;\n"

         "void main(void)\n"
         "{\n"
         "   vertex.w =  1.0;\n"

         "   gl_Position =  matrix * vertex;\n"
         "   gl_Position.w *= -1.0; \n"

         "   texc = texCoord;\n"
         "   gl_Position /=  gl_Position.w ;\n"

         "   vertex1 = gl_Position;\n"
         "}");
    program2->addShaderFromSourceCode(QGLShader::Vertex,
         "attribute vec4 vertex;\n"
         "uniform mat4 matrix;\n"
         "out vec4 vertex1;\n"

         "void main(void)\n"
         "{\n"
         "   vertex.w =  1.0;\n"
         "   gl_Position =  matrix * vertex;\n"
         "   gl_Position.w *= -1.0; \n"
         "   gl_Position /=  gl_Position.w ;\n"

         "   vertex1 = gl_Position;\n"
         "}");

    DOFprogram->addShaderFromSourceCode(QGLShader::Vertex,
         "attribute vec4 vertex;\n"
         "varying vec4 texc;\n"
         "void main(void)\n"
         "{\n"
         "   vertex.w =  1.0;\n"
         "   mat4 I =  mat4(1.0);\n"
         "   gl_Position =  2.0 * vertex - 1.0;\n"
         "   gl_Position.z = 0.0; \n"
         "   gl_Position.w = 1.0; \n"
         "   texc = (gl_Position + 1.0)*0.5; \n"
         "}");

    DOF2program->addShaderFromSourceCode(QGLShader::Vertex,
         "attribute vec4 vertex;\n"
         "varying vec4 texc;\n"
         "void main(void)\n"
         "{\n"
         "   vertex.w =  1.0;\n"
         "   mat4 I =  mat4(1.0);\n"
         "   gl_Position =  2.0 * vertex - 1.0;\n"
         "   gl_Position.z = 0.0; \n"
         "   gl_Position.w = 1.0; \n"
         "   texc = (gl_Position + 1.0)*0.5; \n"
         "}");

    program->addShaderFromSourceCode(QGLShader::Fragment,
        "uniform vec4 color;\n"
        "uniform sampler2D texture;\n"
        "varying vec4 texc;\n"
        "in vec4 vertex1;\n"
        "void main(void)\n"
        "{\n"
        "   vec4 cl = texture2D(texture, texc.st); \n"
        "   gl_FragColor = cl;\n"
        "   gl_FragColor.a = vertex1.z*0.5 + 0.5;\n"

        "}");
    program2->addShaderFromSourceCode(QGLShader::Fragment,
        "uniform vec4 color;\n"
        "in vec4 vertex1;\n"
        "void main(void)\n"
        "{\n"
        "   gl_FragColor = color;\n"
        "   gl_FragColor.a = vertex1.z*0.5 + 0.5;\n"
        "}");
    DOFprogram->addShaderFromSourceCode(QGLShader::Fragment,
        "uniform sampler2D texture;  \n"
        "uniform sampler2D texture1; \n"
//!        "uniform sampler2D diafragm; \n"
        "uniform float width;\n"
        "uniform float height;\n"
        "uniform float focus;\n"
        "uniform float zero;\n"
        "varying vec4 texc;\n"

        "void main(void)\n"
        "{\n"
        "   vec2  dx = vec2(2.0/width, 0.0);      \n"
        "   vec2  dy = vec2(0.0, 2.0/height);     \n"
        "   vec4  Idxdy = vec4(0.0);            \n"
        "   vec4  Sdxdy = vec4(0.0);            \n"
        "   float rmin = length(dx + dy)/2.0;   \n"
        "   float rx = rmin*7;                  \n"
        "   float PI = radians(180);            \n"
        "   float z = 0.0;            \n"

        "   gl_FragColor = texture2D(texture, texc.st);\n"
        "   float depth = gl_FragColor.a*2.0 - 1.0; "


        "     gl_FragColor = vec4(0.0);\n"
        "                                  \n"
        "     for(int i = -15; i<15; i++)     \n"
        "     {                            \n"
        "       for(int j = -15; j<15; j++)     \n"
        "       {                                 \n"
        "           if (i*i+j*j > 225) {continue;} \n"
        "           vec2 delta = dx*i + dy*j;     \n"
        "           Idxdy = texture2D(texture,  texc.st + delta); \n"
        "           Sdxdy = texture2D(texture1, texc.st + delta); \n"
        "           z  = Idxdy.a*2.0 - 1.0;         \n"
        "           rx = abs(0.05*(z - focus)/(z - zero));  \n"
      //!  "           if (rx < rmin) { rx = rmin;}                  \n"

        "  float    ss = Sdxdy.r + Sdxdy.g/256.0;    \n"

        "           if (length(delta)< rx)        \n"
        "           {                             \n"
        "              float xx = length(delta)/rx; \n"
        "              if (depth < z + 2*rx+11.01) {        \n"
        "              gl_FragColor += Idxdy*vec4(ss)*(1-0.3*(xx*xx));        \n"
        "              gl_FragColor.a = 0.0; }        \n"
        "           }                            \n"


        "       }                              \n"
        "     }                              \n"

        "}");

    DOF2program->addShaderFromSourceCode(QGLShader::Fragment,
        "uniform sampler2D texture;\n"
        "uniform sampler2D diafragm;\n"
        "uniform float width;\n"
        "uniform float height;\n"
        "uniform float focus;\n"
        "uniform float zero;\n"
        "varying vec4 texc;\n"

        "void main(void)\n"
        "{\n"
        "   vec2  dx = vec2(2.0/width, 0.0);    \n"
        "   vec2  dy = vec2(0.0, 2.0/height);   \n"
        "   vec4  Idxdy = vec4(0.0);            \n"
        "   float Sdxdy = 0.0;                  \n"
        "   float rmin = length(dx + dy)/2.0;   \n"
        "   float PI = radians(180);            \n"
        "   float z = 0.0;                      \n"

        "   gl_FragColor = texture2D(texture, texc.st); \n"
        "   float depth = gl_FragColor.a*2.0 - 1.0;     \n"
        "   float rx = abs(0.05*(depth - focus)/(depth - zero));  \n"
     //!   "   if (rx < rmin) { rx = rmin;}                  \n"

        "     gl_FragColor = vec4(0.0);            \n"
        "                                          \n"
        "     for(int i = -15; i<15; i++)          \n"
        "     {                                    \n"
        "       for(int j = -15; j<15; j++)        \n"
        "       {                                  \n"
        "           if (i*i+j*j > 225) {continue;} \n"
        "           vec2 delta = dx*i + dy*j;      \n"
        "           if (length(delta)< rx)         \n"
        "           {                              \n"
        "              float xx = length(delta)/rx; \n"
        //! It is functoin(delta) or stencil texture
        "              Sdxdy += 1-0.3*xx*xx;       \n"
        "           }                              \n"
        "       }                                  \n"
        "     }                                    \n"
        "     if (Sdxdy < 1.0) { Sdxdy = 1.0;}     \n"
        "     float ff = 1.0/Sdxdy;                \n"
        "     ff = 1.0/Sdxdy;                        \n"
        "     gl_FragColor.r = floor(256.0*ff)/256.0; \n"
        "     gl_FragColor.g = floor (256.0*256.0*(ff - gl_FragColor.r))/256.0;   \n"


        "     gl_FragColor.b = 256.0*256.0/Sdxdy;          \n"

        //!        "     gl_FragColor.r = rx;          \n"

        "}");


    DOFprogram->link();
    DOF2program->link();
    program->link();
    program2->link();

}


void GLWidget::initializeGL()
{

    glViewport(0,0,width(), height());

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_TEXTURE_2D);
  //  glDisable(GL_BLEND);
  //  glDisable(GL_ALPHA_TEST);

    cubeTexture = bindTexture(QImage(":res/cubelogo.png"));
//    diafTexture = bindTexture(QImage(":res/diafragm.png"));

    fbo = new QGLFramebufferObject(width(), height() , QGLFramebufferObject::Depth);
    fbo2 = new QGLFramebufferObject(width(), height() );

    setWindowTitle(tr("failolo"));

     setupShaders();

//    setWindowTitle(QString( program->log() ));

     vertexLocation = program->attributeLocation("vertex");
     vertex2Location = program2->attributeLocation("vertex");
     texCoordLocation = program->attributeLocation("texCoord");

     DOFvertexLocation = DOFprogram->attributeLocation("vertex");
     DOFwidthLocation = DOFprogram->uniformLocation("width");
     DOFheightLocation = DOFprogram->uniformLocation("height");
     DOFfocusLocation = DOFprogram->uniformLocation("focus");
     DOFzeroLocation = DOFprogram->uniformLocation("zero");
     DOF2vertexLocation = DOF2program->attributeLocation("vertex");
     DOF2widthLocation = DOF2program->uniformLocation("width");
     DOF2heightLocation = DOF2program->uniformLocation("height");
     DOF2focusLocation = DOF2program->uniformLocation("focus");
     DOF2zeroLocation = DOF2program->uniformLocation("zero");


     matrixLocation = program->uniformLocation("matrix");
     matrix2Location = program2->uniformLocation("matrix");
     colorLocation = program->uniformLocation("color");
     color2Location = program2->uniformLocation("color");
     program->setUniformValue("texture", 0);
     DOFprogram->setUniformValue("texture", 0);
     DOFprogram->setUniformValue("texture1", 1);
//     DOFprogram->setUniformValue("diafragm", 2);
     DOF2program->setUniformValue("texture",  0);
//     DOF2program->setUniformValue("diafragm", 1);

     if ( QGLShaderProgram::hasOpenGLShaderPrograms() ) {
         setWindowTitle(QString( "NE FAIL" ));

     }

}


void GLWidget::resizeGL(int w, int h)
{
   glViewport(0, 0, w, h);
   float aspect = (float)h/w;
   pmvMatrix.perspective(45.0, aspect, 1.0, 100.0);
   fbo  = new QGLFramebufferObject(w, h , QGLFramebufferObject::Depth);
   fbo2 = new QGLFramebufferObject(w, h);

}

void GLWidget::paintGL()
{


    QColor clearColor(0,0,0);
    qglClearColor(clearColor);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


     QColor color1 (255, 0, 0, 0);

     glViewport(0, 0, fbo->size().width(), fbo->size().height());


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



     viewMatrix.translate(0, 0, -1);
     viewMatrix.rotate(rotX, 1.0, 0.0, 0.0);
     viewMatrix.rotate(rotY, 0.0, 1.0, 0.0);
     viewMatrix.translate(0, 0, 1);

     viewMatrix.translate(posX, posY, posZ);
   //  viewMatrix.rotate(90, 1.0, 0.0, 0.0);

     viewMatrix.translate(0, +0.5, -2);

     pmvMatrix *= projMatrix;
     pmvMatrix *= viewMatrix;

     // render to the framebuffer object
     //============================================================================
     fbo->bind();
     program->bind();

     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

     glBindTexture(GL_TEXTURE_2D, cubeTexture);
  //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  //   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

     QMatrix4x4 tempe = pmvMatrix;
     int i,j;
     for (i=0; i<7; i++)
         for (j=0; j<7; j++)
         {
         //if (i==j) continue;
         program->enableAttributeArray(vertexLocation);
         program->setAttributeArray(vertexLocation, triangleVertices, 3);
         program->enableAttributeArray(texCoordLocation);
         program->setAttributeArray(texCoordLocation, texCoordArray, 2);

         pmvMatrix = tempe;
         pmvMatrix.translate(i, 0, j);
         //     pmvMatrix.translate(1.0,0.0,1.0);
         program->setUniformValue(matrixLocation, pmvMatrix);

         program->setUniformValue(colorLocation, color1);
         glDrawArrays(GL_TRIANGLES, 0, 6);

         program->disableAttributeArray(vertexLocation);

         program->disableAttributeArray(texCoordLocation);
     }

     for (i=0; i<7; i++)
     {
         program->enableAttributeArray(vertexLocation);
         program->setAttributeArray(vertexLocation, triangleVertices, 3);
         program->enableAttributeArray(texCoordLocation);
         program->setAttributeArray(texCoordLocation, texCoordArray, 2);

         pmvMatrix = tempe;
         pmvMatrix.translate(i,0.4,i);
         glPointSize(GLfloat(7.0));
         program->setUniformValue(matrixLocation, pmvMatrix);

         program->setUniformValue(colorLocation, color1);
         glDrawArrays(GL_POINTS, 0, 6);

         program->disableAttributeArray(vertexLocation);

         program->disableAttributeArray(texCoordLocation);
     }

//     program2->bind();
//     program2->enableAttributeArray(vertex2Location);
//     program2->setAttributeArray(vertex2Location, triangleVertices, 3);

//     pmvMatrix = tempe;
//     program2->setUniformValue(matrix2Location, pmvMatrix);

//     program2->setUniformValue(color2Location, color1);
//     glDrawArrays(GL_POINTS, 0, 6);

//     program2->disableAttributeArray(vertex2Location);



     fbo->release();

     //============================================================================
     fbo2->bind();
     glViewport(0, 0, width(), height());

//     glActiveTexture( 0 );
     glBindTexture(GL_TEXTURE_2D, fbo->texture());

//     glActiveTexture(GL_TEXTURE1);
//     glBindTexture(GL_TEXTURE_2D, diafTexture);

     glClear(GL_COLOR_BUFFER_BIT);

     DOF2program->bind();
     DOF2program->enableAttributeArray(DOF2vertexLocation);
     DOF2program->setAttributeArray(DOF2vertexLocation, quadVertices, 6);

     DOF2program->setUniformValue(DOF2widthLocation, GLfloat(width())   );
     DOF2program->setUniformValue(DOF2heightLocation, GLfloat(height()) );
     DOF2program->setUniformValue(DOF2focusLocation, GLfloat(0.0) );
//     float xxx = -(101)/(99) ;
     DOF2program->setUniformValue(DOF2zeroLocation, GLfloat( -101.0/99.0 ) );

     glDrawArrays(GL_TRIANGLES, 0, 6);

     DOF2program->disableAttributeArray(DOFvertexLocation);
     fbo2->release();


     //============================================================================
//     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

     glViewport(0, 0, width(), height());


     glBindTexture(GL_TEXTURE_2D, fbo->texture());
     glBindTexture(GL_TEXTURE_2D, fbo2->texture());
//     glBindTexture(GL_TEXTURE_2D, diafTexture);

     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

     DOFprogram->bind();
     DOFprogram->enableAttributeArray(DOFvertexLocation);
     DOFprogram->setAttributeArray(DOFvertexLocation, quadVertices, 6);

     DOFprogram->setUniformValue(DOFwidthLocation, GLfloat(width())   );
     DOFprogram->setUniformValue(DOFheightLocation, GLfloat(height()) );
     DOFprogram->setUniformValue(DOFfocusLocation, GLfloat(0.0) );
//     float xxx = -(101)/(99) ;
     DOFprogram->setUniformValue(DOFzeroLocation, GLfloat( -101.0/99.0 ) );


     glDrawArrays(GL_TRIANGLES, 0, 6);

     DOFprogram->disableAttributeArray(DOFvertexLocation);



}


void  GLWidget::mousePressEvent(QMouseEvent *event) {
  lastPos = event->pos();
  killTimer(timerId);
}
void  GLWidget::mouseMoveEvent(QMouseEvent *event) {
  GLfloat dx = (GLfloat)(event->x() - lastPos.x()) / width();
  GLfloat dy = (GLfloat)(event->y() - lastPos.y()) / height();

  GLfloat dd = 0.1;

  if (event->buttons()==Qt::LeftButton) {
    rotX -= 180 * dy;
    rotY -= 180 * dx;
    updateGL();
  } else if (event->buttons()==Qt::RightButton) {
      posX += dd * sin(rotY);
      posZ += dd * cos(rotY);
    updateGL();
  }
  lastPos = event->pos();
}

void GLWidget::wheelEvent(QWheelEvent *event) {
    int numDegrees = event->delta() / 8;
    int numSteps = numDegrees / 15;
    GLfloat dd = numSteps*0.1;

    posX += dd * sin(-rotY/180*3.14);
    posZ += dd * cos(-rotY/180*3.14);

    updateGL();
//    event->accept();
}
