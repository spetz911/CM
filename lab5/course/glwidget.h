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

#include <QtOpenGL>
#include <QMouseEvent>
#include <QWheelEvent>


class GLWidget : public QGLWidget
{
public:
    GLWidget(QWidget *parent);
    ~GLWidget();
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void timerEvent(QTimerEvent *) { update(); }

    void mousePressEvent(QMouseEvent *event); //{ killTimer(timerId); }
    void mouseReleaseEvent(QMouseEvent *event) { timerId = startTimer(20); }
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

    void drawCube(int i, GLfloat z, GLfloat ri, GLfloat jmp, GLfloat amp);
    void setupShaders();



private:
    GLfloat rot[3], xOffs[3], yOffs[3], xInc[3];
    GLuint pbufferList;
    GLuint cubeTexture;
    GLuint diafTexture;
    int timerId;

    QGLFramebufferObject *fbo;
    QGLFramebufferObject *fbo2;
    QGLShaderProgram *program;
    QGLShaderProgram *program2;
    QGLShaderProgram *DOFprogram;
    QGLShaderProgram *DOF2program;

    QMatrix4x4 pmvMatrix;
    QMatrix4x4 DOFMatrix;

    QPoint lastPos;

    float xxxx;

    int vertexLocation ;
    int vertex2Location ;
    int DOFvertexLocation ;
    int DOF2vertexLocation ;
    int texCoordLocation ;
    int matrixLocation ;
    int matrix2Location ;
    int colorLocation ;
    int color2Location ;
    int DOFwidthLocation;
    int DOFheightLocation;
    int DOF2widthLocation;
    int DOF2heightLocation;
    int DOFfocusLocation;
    int DOFzeroLocation;
    int DOF2focusLocation;
    int DOF2zeroLocation;

    //Camera:
       GLfloat rotX;
       GLfloat rotY;

       GLfloat posX;
       GLfloat posY;
       GLfloat posZ;


};

