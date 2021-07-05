#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <initializer_list>
#include "bitmap_image.hpp"
#include "1605089_Object.cpp"



// #include <GL/glext.h>

#define pi (2*acos(0.0))
using namespace std;

double cameraMove;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double uppAngle;
double cylinderAngle;
double rotateAngle;
bool intersect;
double x3, y3,z3;

int levelOfRecursion;
int numberOfObjects;
int numberOfpixels;
int numberOflights;

vector <Object*> objects;
vector <Light> lights;


point u;//upDirection
point r;//rightDirection
point l;//lookDirection
point pos;//position





void drawSS()
{

	for(int i =0; i< objects.size(); i++){
		objects[i]->draw();
	}
	
}

void capture(){

	cout <<"capturing start"<<endl;

	bitmap_image image(numberOfpixels, numberOfpixels);

	double viewAngle = (80 * pi) / 180;
	double planeDistance = (500/2.0) /tan(viewAngle/2.0); //windoeheight and windowwidth 500
	struct point topleft;
	topleft.x = pos.x + (l.x* planeDistance - r.x * (500 / 2.0) + u.x * (500 / 2.0));
	topleft.y = pos.y + (l.y* planeDistance - r.y * (500 / 2.0) + u.y * (500 / 2.0));
	topleft.z = pos.x + (l.z* planeDistance - r.z * (500 / 2.0) + u.z * (500 / 2.0));

	// double du = 500.0/numberOfpixels;
	// double dv = 500.0/numberOfpixels;
	// int nearest = -1;
	double du = 0.74; // trial and error
	double dv = 0.74;
	

	// Choose middle of the grid cell
	topleft.x = topleft.x + r.x*(0.5*du) - u.x*(0.5*dv);
	topleft.y = topleft.y + r.y*(0.5*du) - u.y*(0.5*dv);
	topleft.z = topleft.z + r.z*(0.5*du) - u.z*(0.5*dv);

	for(int i = 0; i < numberOfpixels; i++){
		for(int j = 0; j< numberOfpixels; j++){
			
			double t;
			double tMin = 999999;
			int nearest = -1;

			struct point curPixel;
			curPixel.x = topleft.x + r.x*i*du - u.x*j*dv;
			curPixel.y = topleft.y + r.y*i*du - u.y*j*dv;
			curPixel.z = topleft.z + r.z*i*du - u.z*j*dv;

			struct point dir;

			dir.x = curPixel.x - pos.x;
			dir.y = curPixel.y - pos.y;
			dir.z = curPixel.z - pos.z;

			dir = normalize(dir);

			Ray ray = Ray(pos, dir);

			double *color = new double[3];

			for(int k = 0; k < objects.size(); k++){

				t = objects[k]->intersect(ray, color, 0);
				// if( k == 7){
				// 	cout << " for "<< k <<"th object t is " << t << endl;
				// }
				
				if(t>0 && t < tMin){
					nearest = k;
					tMin = t;
					
				}


			}
			if( nearest!= -1){

				double temp;
				temp = objects[nearest]->intersect(ray, color, 1);
				
				
				

			}
			image.set_pixel(i,j,color[0]*255,color[1]*255,color[2]*255);
			

			// free (color);

			// for(int l = 0; l < 3; l++){
			// 	delete []color;
			// }



		}

		
	}



	image.save_image("output.bmp");
	image.clear();



}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '0':
			
			capture();
			cout <<"capturing done"<<endl;
			break;

		case '1':
			//drawgrid=1-drawgrid;
			l = rotation(u, l, cameraAngle); //rotate/look left
			r = crossProduct(l, u);
			break;
		case '2':
			l = rotation(u, l, -cameraAngle); //rotate/look right
			r = crossProduct(l, u);
			break;
		case '3':
			l = rotation(r, l, cameraAngle); //look up
			u = crossProduct(r, l);
			break;
		case '4':
			l = rotation(r, l, -cameraAngle); //look down
			u = crossProduct(r, l);
			break;
		case '5':
			r = rotation(l, r, -cameraAngle); //tilt clockwise
			u = crossProduct(r, l);
			break;
		case '6':
			r = rotation(l, r, cameraAngle); //tilt counter-clockwise
			u = crossProduct(r, l);
			break;
		case 'p':
		
			
			if(angle<=50){
				angle = angle + 2;

			}
			break;

		case 'w':
		
			
			if(angle>=-50){
				angle = angle - 2;

			}
			break;

		case 'e':
		
			
			if(uppAngle<=60){
				uppAngle = uppAngle + 2;

			}
			break;
		case 'r':
		
			
			if(uppAngle>=-60){
				uppAngle = uppAngle - 2;

			}
			break;

		case 'a':
		
			
			if(cylinderAngle<=20){
				cylinderAngle = cylinderAngle + 2;

			}
			break;
		case 's':
		
			
			if(cylinderAngle>=-20){
				cylinderAngle = cylinderAngle - 2;

			}
			break;

		case 'd':
		
			rotateAngle = rotateAngle + 2;

			
			break;
		case 'f':
		
		    rotateAngle = rotateAngle - 2;
			break;
			

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			//cameraMove -= 3.0;
			pos.x = pos.x - cameraMove*l.x; //move backward
			pos.y = pos.y - cameraMove*l.y;
			pos.z = pos.z - cameraMove*l.z;
			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;
			pos.x = pos.x + cameraMove*l.x; //move forward
			pos.y = pos.y + cameraMove*l.y;
			pos.z = pos.z + cameraMove*l.z;
			break;

		case GLUT_KEY_RIGHT:
			//cameraAngle += 0.03;
			pos.x = pos.x + cameraMove*r.x; //move right
			pos.y = pos.y + cameraMove*r.y;
			pos.z = pos.z + cameraMove*r.z;
			break;
		case GLUT_KEY_LEFT:
			//cameraAngle -= 0.03;

			pos.x = pos.x - cameraMove*r.x; //move left
			pos.y = pos.y - cameraMove*r.y;
			pos.z = pos.z - cameraMove*r.z;
			break;

		case GLUT_KEY_PAGE_UP:

			pos.x = pos.x + cameraMove*u.x; //move up
			pos.y = pos.y + cameraMove*u.y;
			pos.z = pos.z + cameraMove*u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:

			pos.x = pos.x - cameraMove*u.x; //move down
			pos.y = pos.y - cameraMove*u.y;
			pos.z = pos.z - cameraMove*u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			struct point A;
			struct point B;
			A.x = cos(90-uppAngle)*40;
			A.y = cos(angle)*40;
			A.z = sqrt((40*40) - (A.x*A.x) - (A.y*A.y));

			B.x = cos(90-(uppAngle+cylinderAngle))*220;
			B.y = cos(angle)*220;
			B.z = sqrt((220*220) - (B.x*B.x) - (B.y*B.y));

			double x2,y2,z2;
			x2 = B.x - A.x;
			y2 = B.y - A.y;
			z2 = B.z - A.z;

			double t;

			t = (400 - A.y)/y2; // not sure
			
			x3 = A.x + t*x2;
			y3 = A.y + t*y2;
			z3 = A.z + t*z2;

			printf("%f %f %f",&x3,&y3,&z3);

			if(x3>=-200.0 && x3<=200.0){
				printf("\n%d",intersect);
				if(z3>=-200.0 && z3<=200.0){
					intersect = true;
				}
			}

			printf("\n%d",intersect);

			
			
			break;

		case GLUT_RIGHT_BUTTON:

			//........
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	// resets display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

    //glColor3f(1,0,0);
    //drawSquare(10);

    drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	//angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraMove=50.0;
	cameraAngle=5.0;
	angle=0;
	uppAngle = 0;
	cylinderAngle = 0;
	rotateAngle = 0;
	intersect = false;

	pos.x = 100;
	pos.y = 100;
	pos.z = 100;

	u.x = 0;
	u.y = 0;
	u.z = 1;
    
	r.x = -1/sqrt(2);
	r.y = 1/sqrt(2);
	r.z = 0;

	l.x = -1/sqrt(2);
	l.y = -1/sqrt(2);
	l.z = 0;



	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

void loadData(){

	// ifstream newfile("/home/eliza/Documents/Offline 3/INPUT.txt");
	fstream newfile;
	newfile.open("INPUT.txt", ios::in);

	if(newfile.is_open()){
		string lines;
		newfile >> levelOfRecursion;
		cout << "level of recursion  " << levelOfRecursion<<endl;
		newfile >> numberOfpixels;
		cout<<"pixels "<<numberOfpixels<<endl;
		newfile >> numberOfObjects;
		cout<<"objectss "<<numberOfObjects<<endl;

		for (int i =0; i<numberOfObjects;i++){

			newfile >> lines;
					
			if(lines == string("triangle")){
				cout<<"triangle got"<<endl;

				struct point A,B,C;
				double color[3];
				double ambient, diffuse, specular, recursive_reflection_coefficient;
				int shine;



				newfile >> A.x >> A.y >> A.z;
				newfile >> B.x >> B.y >> B.z;
				newfile >> C.x >> C.y >> C.z;
				newfile >> color[0] >> color[1] >> color[2];
				newfile >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
				newfile >> shine;

				Object *triangle = new Triangle(A,B,C);

				triangle->setColor(color[0],color[1],color[2]);
				triangle->setCoEfficients(ambient,diffuse,specular,recursive_reflection_coefficient);
				triangle->setShine(shine);

				objects.push_back(triangle);


			}

			if(lines == string("sphere")){
				cout<<"sphere got" <<endl;

				struct point center;
				double radius;
				double color[3];
				double ambient, diffuse, specular, recursive_reflection_coefficient;
				int shine;



				newfile >> center.x >> center.y >> center.z;
				newfile >> radius;
				newfile >> color[0] >> color[1] >> color[2];
				newfile >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
				newfile >> shine;

				Object *sphere = new Sphere(center, radius);

				sphere->setColor(color[0],color[1],color[2]);
				sphere->setCoEfficients(ambient,diffuse,specular,recursive_reflection_coefficient);
				sphere->setShine(shine);

				objects.push_back(sphere);


			}

			if(lines == string("general")){
				cout<<"general got"<<endl;

				double A, B, C, D, E,F,G,H,I,J;
				struct point reference;
				double length, width, height;
				double color[3];
				double ambient, diffuse, specular, recursive_reflection_coefficient;
				int shine;

				newfile >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
				newfile >> reference.x >> reference.y >> reference.z >> length >> width >> height;
				newfile >> color[0] >> color[1] >> color[2];
				newfile >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
				newfile >> shine;

				Object *general = new General(A,B,C,D,E,F,G,H,I,J,reference,length,width,height);

				general->setColor(color[0],color[1],color[2]);
				general->setCoEfficients(ambient,diffuse,specular,recursive_reflection_coefficient);
				general->setShine(shine);

				general->printA();
				
				objects.push_back(general);



			}


		}

		newfile >> numberOflights;
		cout <<"number of light sources " << numberOflights<<endl;

		for (int i = 0; i < numberOflights; i++){

			struct point position;
			double color[3];

			newfile >> position.x >> position.y >> position.z;
			newfile >> color[0] >> color[1] >> color[2];

			Light light = Light(position, color[0], color[1], color[2]);

			lights.push_back(light);
			cout << lights[i].position.y<<endl;

			
		}
		cout << lights.size() << endl;
		

		newfile.close();
	}

}



int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	// GLenum err = glewInit();
	// if (GLEW_OK != err)
	// {
	// 	fprintf(stderr, "Error %s\n", glewGetErrorString(err));
	// 	exit(1);
	// }
	// fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));

	// if (GLEW_ARB_vertex_program)
	// 	fprintf(stdout, "Status: ARB vertex programs available.\n");

	// if (glewGetExtension("GL_ARB_fragment_program"))
	// 	fprintf(stdout, "Status: ARB fragment programs available.\n");

	// if (glewIsSupported("GL_VERSION_1_4  GL_ARB_point_sprite"))
	// 	fprintf(stdout, "Status: ARB point sprites available.\n");

	init();

	Object *floor = new Floor(1000, 20);

	// floor->setColor(1,1,1);
	floor->setCoEfficients(0.4 , 0.2 , 0.1 ,0.3);
	floor->setShine(5);

	objects.push_back(floor);

	cout << "hi"<<endl;
	loadData();
	cout << "hi2"<<endl;

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	// // loadData();

	glutMainLoop();		//The main loop of OpenGL


	// memory management

	for (int i=0; i< objects.size(); i++){
		delete objects[i];
	}



	return 0;
}
