/*
	question 2: draw a blue triangle
	name: Íõ³¿êØ
	student id: 515030910467

	the main function is at the bottom of this file
*/

#include <Windows.h>
#include <GL/glut.h>

// initial function of the window
void init() {
	glutInitWindowSize(480, 480);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("question 2");
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
}

// draw a triangle
void drawTriangle() {
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_TRIANGLES);
	glColor3f(0, 0, 1);
	glVertex3f(-0.5, -0.4, 0);
	glVertex3f(0.5, -0.4, 0);
	glVertex3f(0, 0.4, 0);
	glEnd();
	glFlush();
}

int main(int argc, char** argv) {

	// initial settings
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);

	// draw a triangle
	init();
	glutDisplayFunc(drawTriangle);
	glutMainLoop();

	return 0;
}