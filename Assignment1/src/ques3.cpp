/*
	question 3: draw and move a car
	name: Íõ³¿êØ
	student id: 515030910467

	the main function is at the bottom of this file
	the class definition is at the top
*/

#include <GL/glut.h>
#include <cmath>
#include <Windows.h>

const double PI = 3.1415926;

// implementation of class Car
// use function draw() to draw the shape
// use function move() to updata the position
// to construct a car, just give an initial position (x0,y0)
class Car {
private:
	GLfloat x0, y0;		// a base position for drawing
	bool toRight;		// a flag to judge whether move to right or left

	void drawTriangle(GLfloat x, GLfloat y);
	void drawRectangle(GLfloat x, GLfloat y);
	void drawWheel(GLfloat x, GLfloat y, GLfloat* color);

public:
	Car(GLfloat* pos);
	~Car() {}
	void draw();
	void move();
};

Car::Car(GLfloat* pos) {
	x0 = pos[0];	y0 = pos[1];
	toRight = true;
}

// draw every components of the car
// including a triangle, a rectangle and two wheels
void Car::draw() {
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	drawTriangle(x0, y0);
	drawRectangle(x0, y0 - 0.2);
	GLfloat frontwheelColor[3] = { 0.0, 1.0, 0.0 };
	GLfloat rearwheelColor[3] = { 0.9, 0.9, 0.0 };
	drawWheel(x0 - 0.1, y0 - 0.3, frontwheelColor);
	drawWheel(x0 + 0.1, y0 - 0.3, rearwheelColor);
	glFlush();
}

// update the position of a car
void Car::move() {
	if (x0 >= 0.81 || x0 <= -0.81)
		toRight = !toRight;
	if (toRight)
		x0 += 0.01;
	else
		x0 -= 0.01;
	glutPostRedisplay();
}

// draw a triangle for the car
void Car::drawTriangle(GLfloat x, GLfloat y) {
	glBegin(GL_TRIANGLES);
	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(x - 0.1, y - 0.1, 0);
	glVertex3f(x + 0.1, y - 0.1, 0);
	glVertex3d(x, y + 0.15, 0);
	glEnd();
}

// draw a rectangle for the car
void Car::drawRectangle(GLfloat x, GLfloat y) {
	glBegin(GL_QUADS);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(x + 0.2, y + 0.1, 0);
	glVertex3f(x + 0.2, y - 0.1, 0);
	glVertex3f(x - 0.2, y - 0.1, 0);
	glVertex3f(x - 0.2, y + 0.1, 0);
	glEnd();
}

// draw a wheel for the car
// the color of a wheel needs to be specified
// use a polygon with many edges to represent a wheel
void Car::drawWheel(GLfloat x, GLfloat y, GLfloat* color) {
	GLfloat r = 0.1;		// radius of a circle
	double n = 200;			// number of polygon edges 
	glBegin(GL_POLYGON);
	glColor3fv(color);
	for (int i = 0; i < n; ++i) {
		glVertex3f(x + r / 2 * cos(2 * PI*i / n), y + r*sin(2 * PI*i / n), 0.0);
	}
	glEnd();
}

// construct a global variable car
GLfloat pos[2] = { -0.8,0.1 };
Car car(pos);

// initial fucntion 
void init() {
	glutInitWindowSize(1000, 500);
	glutInitWindowPosition(50, 50);
	glutCreateWindow("question 3");
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

// a function to call the car to draw
void display() {
	car.draw();
}

// a function to call a car to move
void timeControler(GLint v) {
	car.move();
	glutTimerFunc(20, timeControler, 1);
}


// the main function
int main(int argc, char** argv) {

	// initial settings
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);

	// draw and move a car
	init();
	glutDisplayFunc(display);
	glutTimerFunc(33, timeControler, 1);
	glutMainLoop();

	return 0;
}