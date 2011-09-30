
class MyShader :

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

	fsrc = """
uniform sampler2D texture;
varying mediump vec4 texc;
void main(void)
{
	gl_FragColor = texture2D(texture, texc.st);
}
"""
	
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

	fsrc = """
uniform sampler2D texture;
varying mediump vec4 texc;
void main(void)
{
	gl_FragColor = texture2D(texture, texc.st);
}
"""
	
	
	
	
