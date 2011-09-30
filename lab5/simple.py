program.addShaderFromSourceCode(QGLShader.Vertex,
	"attribute highp vec4 vertex;\n"
	"attribute mediump mat4 matrix;\n"
	"void main(void)\n"
	"{\n"
	"   gl_Position = matrix * vertex;\n"
	"}");
program.addShaderFromSourceCode(QGLShader.Fragment,
	"uniform mediump vec4 color;\n"
	"void main(void)\n"
	"{\n"
	"   gl_FragColor = color;\n"
	"}");
program.link();
program.bind();


With the above shader program active, we can draw a green triangle as follows:

static GLfloat const triangleVertices[] = {
	60.0f,  10.0f,  0.0f,
	110.0f, 110.0f, 0.0f,
	10.0f,  110.0f, 0.0f
};

QColor color(0, 255, 0, 255);

QMatrix4x4 pmvMatrix;
pmvMatrix.ortho(rect());

program.enableAttributeArray(vertexLocation);
program.setAttributeArray(vertexLocation, triangleVertices, 3);
program.setUniformValue(matrixLocation, pmvMatrix);
program.setUniformValue(colorLocation, color);

glDrawArrays(GL_TRIANGLES, 0, 3);

program.disableAttributeArray(vertexLocation);
