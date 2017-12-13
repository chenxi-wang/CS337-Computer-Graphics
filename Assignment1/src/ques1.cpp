/*
	question 1: draw a red cube
	name: 王晨曦
	student id: 515030910467

	the main function is at the bottom of this file
*/

#include <Windows.h>
#include <GL/glut.h>

// an array containing cordinates of 8 vertices of a cube
static const float vertex_list[][3] =
{
	-0.5f, -0.5f, -0.5f,
	0.5f, -0.5f, -0.5f,
	-0.5f, 0.5f, -0.5f,
	0.5f, 0.5f, -0.5f,
	-0.5f, -0.5f, 0.5f,
	0.5f, -0.5f, 0.5f,
	-0.5f, 0.5f, 0.5f,
	0.5f, 0.5f, 0.5f,
};

// an array containing the vertex indices in the same plain
static const GLint index_list[][4] = {
	{0, 1, 3, 2},
	{0, 1, 5, 4},
	{2, 3, 7, 6},
	{1, 3, 7, 5},
	{0, 2, 6, 4},
	{4, 5, 7, 6},
};

// draw a red cube and rotate it for viewing
void drawCube(void)
{
	glClearColor(0.5, 0.5, 0.5, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glPushMatrix();

	// rotate the cube for viewing
	GLfloat rotate1 = -9;
	GLfloat rotate2 = 15;
	glRotatef(rotate1, 0, 1, 0);
	glRotatef(rotate2, 1, 0, 0);

	glColor3f(1, 0, 0);

	// draw 6 plains of a cube
	int i, j;
	for (i = 0; i < 6; ++i) {
		glBegin(GL_QUADS);
		for (j = 0; j < 4; ++j)
			glVertex3fv(vertex_list[index_list[i][j]]);
		glEnd();
	}

	glPopMatrix();

	// add lighting to the environment
	GLfloat light_position[] = { -0.8,0.4,0.5,1 };
	GLfloat light_ambient[] = { 0.3,0,0,1 };
	GLfloat light_diffuse[] = { 1,0,0,1 };
	GLfloat light_specular[] = { 1,0.1,0,1 };
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	
	glutSwapBuffers();
}

void init() {
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(500, 500);
	glutCreateWindow("question 1");
	glClearColor(0.5, 0.5, 0.5, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void main(int argc, char **argv)
{
	
	// initial settings
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	
	// draw a red cube
	init();
	glutDisplayFunc(drawCube);
	glutIdleFunc(drawCube);
	glutMainLoop();
}